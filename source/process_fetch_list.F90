!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      Subroutine process_fetch_list(fetch_list,                        &
                                    istack,                            &
                                    mype,                              &
                                    nprocs,                            &
                                    n_to_left,                         &
                                    tag_offset)


      Use tree
      Use paramesh_Dimensions
      Use mpi_morton
      Use paramesh_mpi_interfaces, only : compress_fetch_list 
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

      Include 'mpif.h'

      Integer, Intent(inout), Dimension(:,:) :: fetch_list
      Integer, Intent(in)    :: istack, mype, nprocs
      Integer, Intent(in)    :: n_to_left(0:nprocs-1)
      Integer, Intent(inout) :: tag_offset

      Integer :: no_of_remote_neighs, max_no_to_be_received
      Integer,Dimension (:),  allocatable :: recvrequest
      Integer,Dimension (:,:),allocatable :: recvstatus
      Integer :: i, j, k, ierror, ierrorcode, iprocs, ii, jj, jm1, jp
      Integer :: ll, kk, isrc, idest, itag, isize

!-----Compress the list of possible off processor blocks by eliminating 
!-----redundant entries (this routine also sorts the list)

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

      commatrix_send(:) = 0
      commatrix_recv(:) = 0
      do i = 1,no_of_remote_neighs
         commatrix_recv(fetch_list(2,i)) =                             &
            commatrix_recv(fetch_list(2,i)) + 1
      End Do
      
!-----Constrcut commatrix_send (the of blocks to send to each other processor) 
!-----by providing the complete commatrix_recv to all processors


      Call MPI_ALLTOALL (commatrix_recv,1,MPI_INTEGER,                 & 
                         commatrix_send,1,MPI_INTEGER,                 & 
                         amr_mpi_meshComm,ierror)

!-----Compute the maximum no. of bloacks which any processor
!-----is going to receive.

       iprocs = 0
       do i = 1,nprocs
          iprocs = iprocs + min(1,commatrix_recv(i))
       End Do
       max_no_to_be_received = max(1,iprocs)

!------Compute the maximum no. of bloacks which any processor
!------is going to receive.

       iprocs = 0
       do i = 1,nprocs
          iprocs = iprocs + min(1,commatrix_send(i))
       End Do
       max_no_to_send = max(1,iprocs)

!-----Evaluate smallest guard block starting index over all pe
!-----store this into variable strt_buffer which is Used in amr_1blk_guardcell

      last_buffer = maxblocks_alloc

      k = maxblocks_alloc + 1
      Do i = 0,nprocs-1
         k = k - commatrix_recv(i+1)
      End Do
      strt_buffer = min(k,last_buffer)

      If (strt_buffer <= lnblocks) Then
        Write(*,*)  & 
        'ERROR in process_fetch_list : guard block starting index',    & 
        strt_buffer,' not larger than lnblocks',lnblocks,              & 
        ' processor no. ',mype,' maxblocks_alloc ',                    & 
        maxblocks_alloc
        Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierror)
      End If

!-----Dynamically allocate memory to store the lists of blocks to be
!-----sent and received.

      If (Allocated(to_be_sent))     Deallocate(to_be_sent)
      If (Allocated(to_be_received)) Deallocate(to_be_received)

      Allocate ( to_be_sent(3,                                         & 
                            max(1,maxval(commatrix_send)),             & 
                            max(1,max_no_to_send) ) )
      Allocate ( to_be_received(3,                                     & 
                                max(1,maxval(commatrix_recv)),         & 
                                max(1,max_no_to_be_received) ) )

!-----Construct arrays to_be_sent and to_be_received which contain
!-----the lists of blocks to be packaged.

      to_be_sent = -1
      to_be_received = -1
      laddress = 0

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
            to_be_received(:,ii,jp) = fetch_list(:,jj)
            jm1 = jj
         End Do

         jj = no_of_remote_neighs
         laddress(1,strt_buffer:strt_buffer+jj-1) = fetch_list(1,1:jj)
         laddress(2,strt_buffer:strt_buffer+jj-1) = fetch_list(2,1:jj)-1
         
      End If

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
         If (commatrix_send(i).gt.0) Then
            kk = kk+1
            isize = 3*commatrix_send(i)
            call MPI_IRECV(to_be_sent(1,1,kk),isize,                   & 
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
         If (commatrix_recv(j).gt.0) Then
            ll = ll+1
            isize = 3*commatrix_recv(j)
            Call MPI_SSEND(to_be_received(1,1,ll),isize,MPI_INTEGER,   & 
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

      End Subroutine process_fetch_list

