!!****if* source/process_fetch_list
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!   process_fetch_list
!!
!! MODIFICATIONS
!!
!!  2022-05-20 K. Weide  Added pattern argument
!!***

#include "paramesh_preprocessor.fh"

Subroutine process_fetch_list(pattern, fetch_list,                     &
                                    istack,                            &
                                    mype,                              &
                                    nprocs,                            &
                                    n_to_left,                         &
                                    tag_offset)

   use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t
   use gr_pmCommPatternData, ONLY: gr_pmCommPatternPtr, &
        gr_pmPrintCommPattern
   Use tree, ONLY: lnblocks, last_buffer
   Use paramesh_Dimensions, ONLY: maxblocks_alloc
   Use paramesh_mpi_interfaces, only : compress_fetch_list
   Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

#include "Flashx_mpi_implicitNone.fh"
#include "FortranLangFeatures.fh"

   TYPE(gr_pmCommPattern_t),POINTER_INTENT_IN :: pattern
   Integer, Intent(inout), Dimension(:,:) :: fetch_list
   Integer, Intent(in)    :: istack, mype, nprocs
   Integer, Intent(in)    :: n_to_left(0:nprocs-1)
   Integer, Intent(inout) :: tag_offset

   Integer :: no_of_remote_neighs, num_sending_pes
   integer :: strtBuffer
   Integer,Dimension (:),  allocatable :: recvrequest
   Integer,Dimension (:,:),allocatable :: recvstatus
   Integer :: i, j, k, ierror, ierrorcode, iprocs, ii, jj, jm1, jp
   Integer :: ll, kk, isrc, idest, itag, isize

!-----Compress the list of possible off processor blocks by eliminating 
!-----redundant entries (this routine also sorts the list)

!!$   associate(commatrix_send => pattern % commatrix_send, &
!!$                commatrix_recv => pattern % commatrix_recv, &
!!$                to_be_sent => pattern % to_be_sent, &
!!$                to_be_received => pattern % to_be_received &
!!$                )
   associate(p => pattern)

      no_of_remote_neighs = 0
      If (istack > 0) Then
       Call compress_fetch_list(fetch_list,                            &
                                istack,                                &
                                no_of_remote_neighs,                   &
                                mype,                                  &
                                nprocs,                                &
                                n_to_left)
      End If

      Do i = 1, no_of_remote_neighs
        fetch_list(2,i) = fetch_list(2,i) + 1
      End Do

!-----Construct commatrix_recv (the number of blocks to receive from each
!-----processor)

      p % commatrix_send(:) = 0
      p % commatrix_recv(:) = 0
      do i = 1,no_of_remote_neighs
         p % commatrix_recv(fetch_list(2,i)) =                             &
            p % commatrix_recv(fetch_list(2,i)) + 1
      End Do
      
!-----Constrcut commatrix_send (the of blocks to send to each other processor) 
!-----by providing the complete commatrix_recv to all processors


      Call MPI_ALLTOALL (p % commatrix_recv,1,MPI_INTEGER,                 &
                         p % commatrix_send,1,MPI_INTEGER,                 &
                         amr_mpi_meshComm,ierror)

!-----Compute the maximum no. of bloacks which any processor
!-----is going to receive.

       iprocs = 0
       do i = 1,nprocs
          iprocs = iprocs + min(1,p % commatrix_recv(i))
       End Do
       num_sending_pes = iprocs ! max(1,iprocs)

!------Compute the maximum no. of bloacks which any processor
!------is going to receive.

       iprocs = 0
       do i = 1,nprocs
          iprocs = iprocs + min(1,p % commatrix_send(i))
       End Do
       p % num_recipient_pes = iprocs ! max(1,iprocs)

!-----Evaluate smallest guard block starting index over all pe
!-----store this into variable strt_buffer which is Used in amr_1blk_guardcell

      last_buffer = maxblocks_alloc

      k = maxblocks_alloc + 1
      Do i = 0,nprocs-1
         k = k - p % commatrix_recv(i+1)
      End Do
      strtBuffer = min(k,last_buffer)
      p % strt_buffer = strtBuffer
#ifdef DEBUG_XTRA
      write(*,*) 'pe ',mype,' begin process_fetch_list: ' &
           ,'p % strt_buffer',p % strt_buffer &
           ,' p % commatrix_recv', p % commatrix_recv
#endif

      If (strtBuffer <= lnblocks) Then
        Write(*,*)  & 
        'ERROR in process_fetch_list : guard block starting index',    & 
        strtBuffer,' not larger than lnblocks',lnblocks,               &
        ' processor no. ',mype,' maxblocks_alloc ',                    & 
        maxblocks_alloc
        Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierror)
      End If

!-----Dynamically allocate memory to store the lists of blocks to be
!-----sent and received.

      If (Allocated(p % to_be_sent))     Deallocate(p % to_be_sent)
      If (Allocated(p % to_be_received)) Deallocate(p % to_be_received)

      Allocate ( p % to_be_sent(3,                                         &
                            max(1,maxval(p % commatrix_send)),             &
                            max(1,p % num_recipient_pes) ) )
      Allocate ( p % to_be_received(3,                                     &
                                max(1,maxval(p % commatrix_recv)),         &
                                max(1,num_sending_pes) ) )

!-----Construct arrays to_be_sent and to_be_received which contain
!-----the lists of blocks to be packaged.

      p % to_be_sent(:,:,:) = -1
      p % to_be_received(:,:,:) = -1
#ifdef DEBUG_XTRA
      write(*,*) 'pe ',mype,' in process_fetch_list: ' &
           ,'p%laddres BEFORE:',p%laddress(:,max(1,size(p%laddress,2)-9):)
#endif
      p % laddress(:,:) = 0

!-----Set up the array to_be_received on each processor
      If (no_of_remote_neighs > 0) Then

         jp = 1
         ii = 0
         jm1 = 1
         Do jj = 1, no_of_remote_neighs
            If (jj == 1 .or.                                          &
            fetch_list(2,jj) == fetch_list(2,jm1)) Then
               ii = ii + 1
            Else
               jp = jp + 1
               ii = 1
            End If
            p % to_be_received(:,ii,jp) = fetch_list(:,jj)
            jm1 = jj
         End Do

         jj = no_of_remote_neighs
         p % laddress(1,strtBuffer:strtBuffer+jj-1) = fetch_list(1,1:jj)
         p % laddress(2,strtBuffer:strtBuffer+jj-1) = fetch_list(2,1:jj)-1
         
      End If

#ifdef DEBUG_XTRA
      write(*,*) 'pe ',mype,' in process_fetch_list: ' &
           ,'no_of_remote_neighs',no_of_remote_neighs  &
           , 'p%laddress AFTER:',p%laddress(:,max(1,size(p%laddress,2)-9):)
#endif
      If (Allocated(recvrequest)) Deallocate( recvrequest )
      Allocate ( recvrequest(nprocs) )
      If (Allocated(recvstatus)) Deallocate( recvstatus )
      Allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )

! Post receives 
      kk = 0
      Do i = 1,nprocs
         
         isrc = i-1
         idest= mype
#ifdef PM_UNIQUE_MPI_TAGS
         itag = isrc*nprocs + idest+1 + tag_offset
#else
         itag = tag_offset !OK because 1 msg per process.
#endif
                                ! receive to pe=j
         If (p % commatrix_send(i).gt.0) Then
            kk = kk+1
            isize = 3*p % commatrix_send(i)
            call MPI_IRECV(p % to_be_sent(1,1,kk),isize,                 &
                 MPI_INTEGER,isrc ,itag,amr_mpi_meshComm,                & 
                 recvrequest(kk),ierror)
         End If
      End Do

! Post sends

      ll = 0
      Do j = 1,nprocs
         
         isrc = mype
         idest= j-1
#ifdef PM_UNIQUE_MPI_TAGS
         itag = isrc*nprocs + idest+1 + tag_offset
#else
         itag = tag_offset !OK because 1 msg per process.
#endif
                                ! send from mype=i
         If (p % commatrix_recv(j).gt.0) Then
            ll = ll+1
            isize = 3*p % commatrix_recv(j)
            Call MPI_SSEND(p % to_be_received(1,1,ll),isize,MPI_INTEGER, &
                           idest,itag,amr_mpi_meshComm,ierror)
         End If
      End Do

#ifdef PM_UNIQUE_MPI_TAGS
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset
#endif

      If (kk.gt.0) Then
         Call MPI_Waitall(kk,recvrequest,recvstatus,ierror)
      End If
     
      If (Allocated(recvrequest)) Deallocate( recvrequest )
      If (Allocated(recvstatus)) Deallocate( recvstatus )

      end associate
#ifdef DEBUG_XTRA
      write(*,*) 'pe ',mype,' end process_fetch_list: ' &
           ,'pattern % strt_buffer',pattern % strt_buffer &
           ,' pattern % commatrix_recv', pattern % commatrix_recv
#endif
#ifdef DEBUG
      call gr_pmPrintCommPattern(pattern,'p_f_l:pattern',mype)
#endif
End Subroutine process_fetch_list

