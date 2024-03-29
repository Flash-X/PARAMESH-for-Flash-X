!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/mpi_morton_bnd
!! NAME
!!
!!   mpi_morton_bnd
!!
!! SYNOPSIS
!!
!!   Call mpi_morton_bnd(mype, nprocs, tag_offset)
!!   Call mpi_morton_bnd(integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in)    :: mype       Local processor id.
!!   Integer, Intent(in)    :: nprocs     Number of processors.
!!   Integer, Intent(inout) :: tag_offset A unique id used in marking messages.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!    paramesh_dimensions
!!    physicaldata
!!    tree
!!    timings
!!    mpi_morton
!!    constants
!!
!! CALLS
!!
!!    mpi_amr_write_guard_comm
!!    process_fetch_list
!!
!! RETURNS
!!
!!    Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine does a communications analysis for guardcell filling and
!!   constructs and stores lists of off-processor blocks which are to be 
!!   communicated during guardcell filling.  Also stored are which sections of 
!!   the blocks to be fetched.
!!
!! AUTHORS
!!
!!    Written by Peter MacNeice  and Michael Gehmeyr, February 2000.
!!    Major simplification and rewrite by Kevin Olson August 2007.
!!
!!***

#include "paramesh_preprocessor.fh"

      Subroutine mpi_morton_bnd(mype,nprocs,tag_offset)

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use mpi_morton
      Use constants

      Use paramesh_mpi_interfaces, only : mpi_amr_write_guard_comm,    & 
                                          process_fetch_list
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Include Statements.
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in)    ::  mype,nprocs
      Integer, intent(inout) ::  tag_offset

!-----Local variables and arrays.
      Real    :: eps,accuracy
      Real    :: pbsize(3),pcoord(3),pbndbox(2,3)

      Integer :: lb,i,j,k,j00
      integer :: ia,ib,ja,jb,ka,kb
      Integer :: ierror
      Integer :: max_no_of_blocks
      Integer :: istack, ioff, joff, koff, itemp
      Integer :: isize, isrc, idest, itag, kk
      Integer :: nguarda, iproc
      Integer :: npts_neigh1,npts_neigh2
      Integer,Dimension(:),        Allocatable :: n_to_left
      Integer,Dimension(:,:,:,:,:),Allocatable :: psurr_blks
      Integer,Dimension (:),       Allocatable :: recvrequest
      Integer,Dimension (:,:),     Allocatable :: recvstatus
      Integer,Dimension (:,:),     Allocatable :: fetch_list
      Integer,Dimension (:,:),     Allocatable :: tfetch_list

!-----Begin Executable Code

      accuracy = 100./10.**precision(accuracy)
      eps = accuracy                                                                                                     
      nguarda = max(nguard,nguard_work)

      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      allocate(fetch_list(3,npts_neigh2))
      allocate(n_to_left(0:nprocs-1))

!-----store the max no of blocks on any one processor
      Call MPI_ALLREDUCE(lnblocks,                                     & 
                         max_no_of_blocks,                             & 
                         1,                                            & 
                         MPI_INTEGER,                                  & 
                         MPI_MAX,                                      & 
                         amr_mpi_meshComm,                               & 
                         ierror)

!-----COMPUTE the number of blocks to the 'left' (ie. stored on processors with
!-----smaller process ids) of every other processor
      Call MPI_ALLGATHER(lnblocks,                                     &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         n_to_left,                                    &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         amr_mpi_meshComm,                               &
                         ierror)
                        
      Do iproc = nprocs-1, 1, -1
         n_to_left(iproc) = n_to_left(iproc-1)
      End Do
      n_to_left(iproc) = 0
      Do iproc = 2, nprocs-1
         n_to_left(iproc) = n_to_left(iproc) + n_to_left(iproc-1)
      End Do

!-----FETCH surr_blks lists of parents
      i = size(surr_blks,dim=2)
      j = size(surr_blks,dim=3)
      k = size(surr_blks,dim=4)
      Allocate(psurr_blks(3,i,j,k,lnblocks))
      psurr_blks(:,:,:,:,:) = -1

      if(allocated(recvrequest)) deallocate( recvrequest )
      allocate ( recvrequest(maxblocks) )
      if(allocated(recvstatus)) deallocate( recvstatus )
      allocate ( recvstatus(MPI_STATUS_SIZE,maxblocks) )

      isize = 3*size(surr_blks,dim=2)*                                 &
                size(surr_blks,dim=3)*                                 &
                size(surr_blks,dim=4)

!-----Post receives (children receive from parent block)
      kk = 0
      Do lb = 1,lnblocks
         If (parent(1,lb) > 0) Then
         If (parent(2,lb) .ne. mype) then
          isrc  = parent(2,lb)
#ifdef PM_UNIQUE_MPI_TAGS
          itag  = mype*max_no_of_blocks + lb
#else
          itag  = lb !Locally unique, not globally unique.
#endif
          kk = kk+1
          Call MPI_IRECV(psurr_blks(1,1,1,1,lb),                       &
                         isize,                                        & 
                         MPI_INTEGER,                                  &
                         isrc,                                         &
                         itag,                                         &
                         amr_mpi_meshComm,                               & 
                         recvrequest(kk),                              &
                         ierror)
          ! Commented out next two lines (Else ...), not needed - KW 2018-11-05
!!$         Else
!!$            psurr_blks(:,:,:,:,lb) = surr_blks(:,:,:,:,parent(1,lb))
         End If  ! End If (parent(2,lb) .ne. mype) 
         End If  ! End If (parent(1,lb) > 0)
      End Do  ! End Do lb = 1, lnblocks

!-----Post sends from parents to their children
      Do lb = 1, lnblocks
        Do j = 1,nchild
          If (child(1,j,lb) > 0) Then
          If (child(2,j,lb) .ne. mype) Then
           idest = child(2,j,lb)
#ifdef PM_UNIQUE_MPI_TAGS
           itag  = child(2,j,lb)*max_no_of_blocks + child(1,j,lb)
#else
           itag  = child(1,j,lb) !Locally unique, not globally unique.
#endif
           Call MPI_SSEND(surr_blks(1,1,1,1,lb),                       &
                          isize,                                       &
                          MPI_INTEGER,                                 &
                          idest,                                       &
                          itag,                                        &
                          amr_mpi_meshComm,                              &
                          ierror)
          End If  ! End If (child(2,j,lb) .ne. mype)
          End If  ! End If (child(1,j,lb) > 0)
        End Do  ! End Do j = 1,nchild
      End Do  ! End Do lb = 1, lnblocks

      If (kk.gt.0)                                                     & 
        Call MPI_WAITALL(kk,recvrequest,recvstatus,                    & 
                         ierror)

!-----Construct a list of potential neighbors of all blocks on this
!-----processor, and potential neighbors of their parents.
!-----Exclude any which are on this processor.

      istack = 0

      Do lb = 1, lnblocks

      If (nodetype(lb) <= 2 .or. advance_all_levels) Then

!-------Compute geometry information for parent if spherical_pm is defined
       If (parent(1,lb) > 0 .and. spherical_pm) Then

        pbsize(:) = bsize(:,lb)*2.             ! size of parent block
        ioff = mod(which_child(lb)-1,2)        ! coord for parent block
        joff = mod((which_child(lb)-1)/2,2)
        koff = mod((which_child(lb)-1)/4,2)
        If (ioff == 0) Then
          pcoord(1) = bnd_box(2,1,lb)
        Else
          pcoord(1) = bnd_box(1,1,lb)
        End If
        If (joff == 0) Then
          pcoord(2) = bnd_box(2,2,lb)
        Else
          pcoord(2) = bnd_box(1,2,lb)
        End If
        If (ndim < 2) pcoord(2) = coord(2,lb)
        If(koff == 0) Then
          pcoord(3) = bnd_box(2,3,lb)
        Else
          pcoord(3) = bnd_box(1,3,lb)
        End If
        If (ndim < 3) pcoord(3) = coord(3,lb)
        pbndbox(1,:) = pcoord(:) - bsize(:,lb)
        pbndbox(2,:) = pcoord(:) + bsize(:,lb)
        If (ioff == 0) Then
          pbndbox(2,1) = pcoord(1) + pbsize(1)
        Elseif(ioff == 1) Then
          pbndbox(1,1) = pcoord(1) - pbsize(1)
        End If
        If (joff == 0) then
          pbndbox(2,2) = pcoord(2) + pbsize(2)
        Elseif(joff == 1) Then
          pbndbox(1,2) = pcoord(2) - pbsize(2)
        End If
        If (koff == 0) Then
          pbndbox(2,3) = pcoord(3) + pbsize(3)
        Elseif(koff == 1) Then
          pbndbox(1,3) = pcoord(3) - pbsize(3)
        End If

       End If  ! End If (parent(1,lb) > .0 .and. spherical_pm)

!------ADD OFF PROCESSOR NEIGHBORS OF BLOCK 'lb' TO FETCH LIST
       Do k = 1,1+2*k3d
        Do j = 1,1+2*k2d
        Do i = 1,3
        If (.NOT. diagonals) then
           if (abs(i-2)+abs(j-(1+k2d))+abs(k-(1+k3d)) .NE. 1) cycle
        end If
        If (surr_blks(1,i,j,k,lb) > 0 .and.                            & 
            surr_blks(2,i,j,k,lb) .ne. mype) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = surr_blks(1,i,j,k,lb)
            fetch_list(2,istack) = surr_blks(2,i,j,k,lb)
            j00 = j

!-----------if this block is a polar block then change the way j is applied
!-----------in the formula
            
            j00 = j00 + 1 - k2d
            If(spherical_pm) Then
               If(lsingular_line) Then
                  If(abs(pbndbox(1,2)) < eps.and.j+1-k2d == 1) Then
                     j00 = 3
                  Else If(abs(pbndbox(2,2)-pi) < eps .and.             &
                          j+1-k2d == 3) Then
                     j00 = 1
                  End If
               End If
            End If

!-----------compute message type - note this index is computed to reflect the part
!-----------of the remote block to be acquired, not the part of the local blocks
!-----------guardcells which will be filled.
            kk = k + 1 - k3d
            fetch_list(3,istack) = (4-i)+((4-j00)-1)*3+((4-kk)-1)*9
            if(nguarda.gt.nmax_lays) fetch_list(3,istack) = 14

         End If  ! End If surr_blks(1,i,j,k,lb) > 0 .and. ...
         End Do  ! End Do i = 1,3
         End Do  ! End Do j = 1,1+2*k2d
        End Do  ! End Do k = 1,1+2*k3d

!-------ADD PARENT TO FETCH LIST (if off processor)
        If (parent(1,lb) > 0 .and. parent(2,lb) .ne. mype) Then
            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = parent(1,lb)
            fetch_list(2,istack) = parent(2,lb)
            fetch_list(3,istack) = 14
        End If

!-------ADD PARENT'S surrounding blocks to fetch list if the 
!-------block 'lb' is a leaf and it is at a refinement jump.
        If (parent(1,lb) > 0 .and. nodetype(lb) == 1 .and. & 
            parent(2,lb) .ne. mype                   .and. &
            any(surr_blks(1,1:3,1:1+2*k2d,1:1+2*k3d,lb) > -20 .and.    & 
                surr_blks(1,1:3,1:1+2*k2d,1:1+2*k3d,lb) < 0)) Then

        ka = 1; kb = 1+2*k3d
        ja = 1; jb = 1+2*k2d
        ia = 1; ib = 3
#ifdef PM_OPTIMIZE_MORTONBND_FETCHLIST
        ! Do this optimization only if parent does not touch a domain boundary
        ! anywhere (otherwise boundary condition can be called with invalid
        ! input data):
        if (ALL(psurr_blks(1,1:3,1:1+2*k2d,1:1+2*k3d,lb) > -20)) then
           if (nguarda .LE. nzb/2) then
              koff = mod((which_child(lb)-1)/4,2)
              ka = 1+koff; kb = 1+k3d+koff
           end if
           if (nguarda .LE. nyb/2) then
              joff = mod((which_child(lb)-1)/2,2)
              ja = 1+joff; jb = 1+k2d+joff
           end if
           if (nguarda .LE. nxb/2) then
              ioff = mod(which_child(lb)-1,2)
              ia = 1+ioff; ib = 2+ioff
           end if
        end if
#endif

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
        If (psurr_blks(1,i,j,k,lb) > 0 .and.                           & 
            psurr_blks(2,i,j,k,lb) .ne. mype) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = psurr_blks(1,i,j,k,lb)
            fetch_list(2,istack) = psurr_blks(2,i,j,k,lb)
            j00 = j

!-----------if this block is a polar block then change the way j is applied
!-----------in the formula
            j00 = j00 + 1 - k2d
            If(spherical_pm) Then
               If(lsingular_line) Then
                  If(abs(pbndbox(1,2)) < eps.and.j+1-k2d == 1) Then
                     j00 = 3
                  Else If(abs(pbndbox(2,2)-pi) < eps .and.             &
                          j+1-k2d == 3) Then
                     j00 = 1
                  End If
               End If
            End If
        
!-----------compute message type - note this index is computed to reflect the part
!-----------of the remote block to be acquired, not the part of the local blocks
!-----------guardcells which will be filled.
            kk = k + 1 - k3d
            fetch_list(3,istack) = (4-i)+((4-j00)-1)*3+((4-kk)-1)*9
            If (nguarda > nmax_lays) fetch_list(3,istack) = 14

         End If  ! End If psurr_blks(1,i,j,k,lb) > 0 .and. ...
         End Do  ! End Do i = 1+ioff, ...
         End Do  ! End Do j = 1+joff, ...
         End Do  ! End Do k = 1+koff, ...

         End If  ! End If (nodetype(lb) == 1 ...

      End If  ! End If (nodetype(lb) <= 2 .or. advance_all_levels)

      End Do  ! End Do lb = 1, lnblocks

!-----Compress 'fetch_list' and eliminate any redundances
!-----Also, if a block appears in the list more than once this routine
!-----will adjust which section of the block is requested (i.e. corner, 
!-----face, or entire block)
      Call process_fetch_list(fetch_list,                              &
                              istack,                                  &
                              mype,                                    &
                              nprocs,                                  &
                              n_to_left,                               &
                              tag_offset)

!------Store communication info for future Use
       Call mpi_amr_write_guard_comm(nprocs)

!-----Deallocate any memory which was dynamically allocated for local 
!-----Use in this routine.
      If (Allocated(fetch_list)) Deallocate(fetch_list)
      If (Allocated(n_to_left)) Deallocate(n_to_left)
      If (Allocated(psurr_blks)) Deallocate(psurr_blks)
      If (Allocated(recvrequest)) Deallocate( recvrequest )
      If (Allocated(recvstatus)) Deallocate( recvstatus )

      Return

      Contains
        Subroutine expand_fetch_list

               If(allocated(tfetch_list)) deallocate(tfetch_list)
               allocate(tfetch_list(3,npts_neigh2))
               tfetch_list(:,:istack-1) = fetch_list(:,:istack-1)
               npts_neigh1 = npts_neigh1 + 3000
               npts_neigh2 = npts_neigh2 + 3000
               deallocate(fetch_list)
               allocate(fetch_list(3,npts_neigh2))
               fetch_list(:,:istack-1) = tfetch_list(:,:istack-1)
               deallocate(tfetch_list)

        End Subroutine expand_fetch_list

      End Subroutine mpi_morton_bnd



!--------------------------------------------------

!--------------------------------------------------

!--------------------------------------------------

