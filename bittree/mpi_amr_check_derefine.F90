!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* bittree/amr_check_derefine
!! NAME
!!
!!   amr_check_derefine
!! 
!! SYNOPSIS
!!
!!   call amr_check_derefine (mype)
!!
!!   call amr_check_derefine (integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: mype     
!!     The calling processor number.
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   bittree
!!
!! RETURNS
!!
!!   This routine does not return anything.  Upon return some derefine flags
!!   may be set to .false. to ensure a valid mesh.
!!
!! DESCRIPTION
!!
!!   This routine performs checks on the list of derefine flags to make sure that
!!   a valid mesh results when the new, child blocks are added.  It makes sure that
!!   a jump in refinement of no more than one level is ensured.
!!  
!!   The routine is called only my amr_refine_derefine and a user's application should
!!   not need to call this routine directly.
!!
!!   The Bittree version of this routine also updates the refine_delta array in
!!   Bittree (globally) for changed derefine flags.
!!
!! AUTHORS
!!
!!   Tom Klosterman (2019) 
!!
!!***


#include "paramesh_preprocessor.fh" 

      subroutine amr_check_derefine(mype)

      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_comm_data, ONLY : amr_mpi_meshComm

      use bittree
      
      implicit none

      include 'mpif.h'
      
      Integer, intent(in) :: mype

!-----local variables
      integer :: nprocs
      integer :: i,ierr
      integer :: loopcount
      integer :: li,lj,lk
      integer, dimension(mdim) :: gCell, gr
      integer :: lcoord(3),neighCoord(3)
      logical :: ref_mark,deref_mark,is_parent

      logical :: repeat,repeat_t
      logical :: deref_test(maxblocks_tr)
      logical :: first
      logical :: lreduce_datain(1),lreduce_dataout(1)
 
!----------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' entered amr_check_derefine'
#endif
      
      Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)

      repeat = .TRUE.
      loopcount=0
     
!-----Do not derefine non-leaf blocks
      Do i=1,lnblocks
        if (derefine(i).AND.nodetype(i).ne.1) derefine(i) = .FALSE.
      end do !i=1,lnblocks

!-----If any parents are marked for refine, make sure they are NOT
!-----marked for nodetype change on Bittree. Then turn off refine.
!-----NOTE: This should only be used for the behavior of marking parents
!-----for refine to signal children to not derefine.
      do i=1,lnblocks
        if (refine(i).AND.nodetype(i).gt.1) then
          call gr_getIntCoords(i,lcoord)
          call gr_btDerefineMark(lrefine(i)+1,lcoord*2,.FALSE.)
          refine(i) = .FALSE.
        end if
      end do !i=1,lnblocks

!----------------------------------------------------------------------

      do while (repeat)

#ifdef DEBUG
      write(*,*) 'amr_check_derefine : proc ',mype, &
                 ' derefine(1:lnblocks) ',derefine(1:lnblocks)
#endif

!-----Turn OFF deref_test if block can't be derefined
      deref_test = derefine 
      loopcount = loopcount+1

!-----Check surr_blks - if any neighbor is either a parent 
!-----or marked for refinement, do not derefine.
      do i=1,lnblocks
        if (derefine(i)) then

!---------If block doesn't have a parent, can't derefine
          if (parent(1,i).le.0) deref_test(i) = .FALSE.
          
!---------Loop over surr_blks
          kAxis: do lk = -k3d, k3d
          jAxis: do lj = -k2d, k2d
          iAxis: do li = -k1d, k1d
!-----------gr in range 1,3 while gCell in range -1,1
            gCell = (/li,lj,lk/)
            gr = gCell + (/k1d+1,k2d+1,k3d+1/)

!-----------Use surr_blks to check if a neighbor exists
            if (surr_blks(1,gr(1),gr(2),gr(3),i).ge.1) then

!-----------Calculate lcoord of neighbor and check parentage
            call gr_getNeighIntCoords(i,gCell,neighCoord)
            call gr_btIsParent(lrefine(i),neighCoord, is_parent)
            call gr_btGetRefine(lrefine(i), neighCoord,ref_mark)

!-----------Dont derefine if neighbor is marked for refine
            if (ref_mark) then
              deref_test(i) = .FALSE.
              exit kAxis
            end if !marked

!-----------Don't derefine if neighbor is a parent, unless it's 
!-----------an external neighbor with children marked for derefinement
            if (is_parent) then
              call gr_getIntCoords(i,lcoord)
              if (any((lcoord/2).ne.(neighCoord/2))) then
                call gr_btGetDerefine(lrefine(i)+1,  &
                                              neighCoord*2,  &
                                              deref_mark)
                if(.NOT.deref_mark) then
                  deref_test(i) = .FALSE.
                  exit kAxis
                end if
              else
                deref_test(i) = .FALSE.
                exit kAxis
              end if
            end if !is_parent

            end if !surr_blks>=1
          end do iaxis
          end do jaxis
          end do kaxis
        end if !derefine(i)
      end do !i=1,lnblocks

#ifdef DEBUG
      write(*,*) 'amr_check_derefine : proc ',mype,  &
                 ' deref_test(1:lnblocks) ',deref_test(1:lnblocks)
#endif

!-------------------------------------------------------------------------------
!-----Now: Set derefine flags based on deref_test and repeat process if necessary

      repeat = .FALSE.

      do i = 1,lnblocks

!-------If a block changed derefine status, signal loop to repeat 
        if (.NOT.deref_test(i).AND.nodetype(i).eq.1.AND.derefine(i) ) then
          repeat = .TRUE.
          derefine(i) = .FALSE.
        end if

!-------If a leaf is unmarked by deref_test, unmark (its parent) on bittree
        if (.NOT.deref_test(i).AND.nodetype(i).eq.1 &
             .AND. parent(1,i).ge.1 ) then
          call gr_getIntCoords(i,lcoord)
          call gr_btDerefineMark(lrefine(i), lcoord, .FALSE.)
        end if
      end do !i=1,lnblocks

!-----Reduce bittree by ANDing, unmarking blocks globally
      call bittree_refine_reduce_and(amr_mpi_meshComm)

!-----If any blocks are still marked for derefinement, check to make
!-----sure their parents are still marked on bittree. This ensures blocks
!-----only derefine if ALL siblings are marked for derefine.
      do i = 1,lnblocks
        if (derefine(i)) then
          call gr_getIntCoords(i,lcoord)
          call gr_btGetDerefine(lrefine(i),lcoord,deref_mark)
 
!---------If (parent) not marked, change status and request a loop repeat
          if (.NOT.deref_mark) then
            repeat = .TRUE.
            derefine(i) = .FALSE.
          end if
        end if 
      end do !i=1,lnblocks


!-----Check all processesors to see if a repeat is necessary 
      lreduce_datain(1) = repeat
      call mpi_logical_allreduce(lreduce_datain,lreduce_dataout,1, &
                        MPI_LOGICAL, MPI_LOR,amr_mpi_meshComm,ierr)
      repeat_t = lreduce_dataout(1)
      repeat = repeat_t

      end do !while(repeat)

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' exiting amr_check_derefine'
#endif

      end subroutine
