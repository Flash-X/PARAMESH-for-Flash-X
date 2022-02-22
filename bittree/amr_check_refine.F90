!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
! Copyright (C) 2019 The University of Chicago
! Copyright (C) 2022 UChicago Argonne, LLC and contributors
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* bittree/amr_check_refine
!! NAME
!!
!!   amr_check_refine
!! 
!! SYNOPSIS
!!
!!   call amr_check_refine (nprocs, mype, icontinue)
!!
!!   call amr_check_refine (integer, integer, integer)
!!
!! ARGUMENTS      
!!
!!   integer, intent(in) :: nprocs
!!     The number of processors being used.
!!
!!   integer, intent(in) :: mype     
!!     The calling processor number.
!!
!!   integer, intent(in) :: icontinue     
!!     A variable controling algorithm behavior. (unused in Bittree version)
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
!!   This routine does not return anything.  Upon return additional refine flags
!!   may be set to .true. to ensure a valid mesh.
!!
!! DESCRIPTION
!!
!!   This routine performs checks on the list of refine flags to make sure that
!!   a valid mesh results when the new, child blocks are added.  It makes sure that
!!   a jump in refinement of no more than one level is ensured.
!!  
!!   The routine is called only my amr_refine_derefine and a user's application should
!!   not need to call this routine directly.
!!
!!   The Bittree version of this routine also updates the refine_delta array in
!!   Bittree (globally) for new refine flags.
!!
!! AUTHORS
!!
!!   Tom Klosterman (2019) 
!!
!!***

#include "paramesh_preprocessor.fh"

      subroutine amr_check_refine(nprocs,mype,icontinue)
      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_comm_data, ONLY : amr_mpi_meshComm
      
      use bittree

      implicit none

      integer, intent(in) :: nprocs, mype
      integer, intent(in) :: icontinue   !unused

      include "mpif.h"

!-----local variables
      integer :: i,j
      integer :: ix,iy,iz
      integer :: li,lj,lk
      integer :: ierr
      integer, dimension(mdim) :: gCell,gr

      logical :: repeat,repeat_t
      logical :: ref_test(maxblocks_tr)
      logical :: lreduce_datain(1),lreduce_dataout(1)

      integer :: lcoord(3),neighCoord(3),childCoord(3)
      logical :: marked,is_parent,neigh_par

!---------------------------------------------------------------------------------
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' entered amr_check_refine'
#endif
      
      repeat = .TRUE.

!-----Repeat is left true if another round is needed
      do while (repeat)


#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype, &
                 ' refine(1:lnblocks) ',refine(1:lnblocks)
#endif

      ref_test(1:lnblocks) = .FALSE.

!-----Check adjacent children of neighbors of local leaf blocks to see if they
!-----are marked for refinement.

      do i=1,lnblocks
        if (nodetype(i).eq.1.AND..NOT.refine(i)) then

!---------Loop over surr_blks
          kAxis1: do lk = -k3d, k3d 
          jAxis1: do lj = -k2d, k2d 
          iAxis1: do li = -k1d, k1d 
!-----------gr in range 1,3 while gCell in range -1,1
            gCell = (/li,lj,lk/)
            gr = gCell + (/k1d+1,k2d+1,k3d+1/)

!-----------Use surr_blks to check if a neighbor exists
            if (surr_blks(1,gr(1),gr(2),gr(3),i).ge.1) then

!-----------Calculate lcoord of neighbor and check parentage
            call gr_getNeighIntCoords(i,gCell,neighCoord)
            call gr_btIsParent(lrefine(i),neighCoord,neigh_par)

!-----------Then if that neighbor is a parent, check its children 
            if (neigh_par) then
              
              do iz = 1,1+k3d
                do iy = 1,1+k2d
                  do ix = 1,1+k1d

!-----------------Check child if it's adjacent to block i, face j
                    if ( ((((1-li)/2)+1).eq.ix.OR.li.eq.0) .AND. &
                         ((((1-lj)/2)+1).eq.iy.OR.lj.eq.0) .AND. &
                         ((((1-lk)/2)+1).eq.iz.OR.lk.eq.0) ) then

!---------------------Calculate child's coordinates
                      childCoord = neighCoord*2 + (/ix-1,iy-1,iz-1/)

                      call gr_btGetRefine(lrefine(i)+1,   &
                                                  childCoord,marked)

!---------------------Set ref_test to TRUE if any of the children are marked
                      if(marked) then
                        ref_test(i) = .TRUE.
                        exit kAxis1
                      end if

                    end if !if adjacent

                  end do !ix
                end do !iy
              end do !iz
            
            end if !neigh_is_parent
          end if !(surr_blks(1,gr(1),gr(2),gr(3),i).ge.1)
         end do iAxis1
         end do jAxis1
         end do kAxis1
        end if ! nodetype(i).eq.1
      end do  ! i=1,lnblocks

#ifdef DEBUG
      write(*,*) 'amr_check_refine : proc ',mype,  &
                 ' ref_test(1:lnblocks) ',ref_test(1:lnblocks)
#endif

!-------------------------------------------------------------------------------
!-----Now: Set refine flags based on ref_test and repeat process if necessary

      repeat = .FALSE.
      
      do i = 1,lnblocks

!-------If a leaf is marked by ref_test but not by refine, change its 
!-------refine flag and mark it for refinement on bittree
        if (ref_test(i).AND.nodetype(i).eq.1 .AND. .NOT.refine(i) ) then

          repeat = .TRUE.
          refine(i) = .TRUE.
          derefine(i) = .FALSE.

          call gr_getIntCoords(i,lcoord)
          call gr_btRefineMark(lrefine(i),lcoord)
        end if !if bittree needs to be updated
      end do !i=1,lnblocks

!-----Check all processesors to see if a repeat is necessary, and update 
!-----Bittree's delta tree if necessary.
      lreduce_datain(1) = repeat
      call mpi_logical_allreduce(lreduce_datain,lreduce_dataout,1, &
                        MPI_LOGICAL, MPI_LOR,amr_mpi_meshComm,ierr)
      repeat_t = lreduce_dataout(1)
      repeat = repeat_t
      if (repeat) call bittree_refine_reduce(amr_mpi_meshComm)


      end do !while (repeat) 

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' exiting amr_check_refine'
#endif
      end subroutine
