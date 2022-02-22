!!****if* source/Grid/GridMain/paramesh/bittree/source/amr_morton_order_bittree.F90
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
!!  amr_morton_order_bittree
!!
!! SYNOPSIS
!!
!!  call amr_morton_order_bittree(nprocs, mype)
!!
!! DESCRIPTION
!!
!!  Replaces funciton of amr_refine_blocks, amr_derefine_blocks, and 
!!  amr_morton_order. Computes the tree data for new list of local blocks,
!!  Calculates new_loc to control where physical data of old blocks will go,
!!  and finally calls redist_blk to actually move physical data.
!!
!! ARGUMENTS
!!
!!  integer,intent(in)  :: nprocs
!!  integer,intent(in)  :: mype
!!
!!***

#include "paramesh_preprocessor.fh"

      subroutine amr_morton_order_bittree(nprocs,mype,ref_count)

!-----Use statements
      use paramesh_dimensions
      use physicaldata
      use tree
      use io
      use bittree, only: gr_btSortMortonBittree, old_localMortUB, &
                         localMortUB, old_localMortUB, gr_getIntCoords,&
                         amr_calculate_tree_data, gr_btIdentify, &
                         amr_exchange_work_bflags

      implicit none

!-----Include statements
      include 'mpif.h'

!-----Input/output statements
      integer, intent(in) :: nprocs, mype
      integer, intent(in) :: ref_count

!-----Local variables and arrays
      integer :: new_loc(2,maxblocks_tr)
      integer :: old_loc(2,maxblocks_tr)
      integer :: new_child(2,mchild,maxblocks_tr)
      integer :: lcoord(3),childCoord(3),lev,proc,locblk
      integer :: i, j, lnblocks_old
      logical :: marked
      logical, save :: first = .TRUE.


!-----Begin executable code
      lnblocks_old = lnblocks

!-----Store old version of localMortUB for use in amr_identify
      allocate(old_localMortUB(nprocs))
      old_localMortUB = localMortUB

!-----Adjust localMortUB (this sorts by work if applicable)
      if(ref_count.ne.0.OR..NOT.first) then
        call gr_btSortMortonBittree(nprocs, mype)
      end if
      first = .FALSE.

!-----Set new_lnblocks
      if(mype.eq.0) then
        new_lnblocks = localMortUB(mype+1)
      else
        new_lnblocks = localMortUB(mype+1) - localMortUB(mype)
      end if

!-----Get new locations of curent local blocks
      new_loc(:,:) = -1
      new_child(:,:,:) = -1
      do i=1,lnblocks_old
        if(.NOT.derefine(i)) then
          lev = lrefine(i)
          call gr_getIntCoords(i,lcoord)
          call gr_btIdentify(nprocs,lev,lcoord,   &
                                  proc,locblk,updated=.TRUE.)
          new_loc(1,i) = locblk
          new_loc(2,i) = proc

          if(refine(i)) then
            do j=1,nchild
              lev = lrefine(i)+1
              childCoord = lcoord*2 + (/ mod((j-1),2),       &
                                       mod((j-1)/2,2),       &
                                       mod((j-1)/4,2) /)
              call gr_btIdentify(nprocs,lev,childCoord,proc,&
                                      locblk,updated = .TRUE.)
              new_child(1,j,i) = locblk
              new_child(2,j,i) = proc
            end do
          end if

        end if
      end do

!---------------------------------------------------------------------
!-----Use Bittree to calculate new tree data for each new local block
      call amr_calculate_tree_data(nprocs,mype,lnblocks_old)

!-----Fill old_loc (use parent's location if newchild)
      old_loc(:,:) = -1
      do i=1,new_lnblocks
        if(.NOT.newchild(i)) then
          lev = lrefine(i)
          call gr_getIntCoords(i,lcoord)
        else
          lev = lrefine(i)-1
          call gr_getIntCoords(i,lcoord)
          lcoord = lcoord/2
        end if
        call gr_btIdentify(nprocs,lev,lcoord, &
                                proc,locblk,updated=.FALSE.)
        old_loc(1,i) = locblk
        old_loc(2,i) = proc
      end do

!---------------------------------------------------------------------
!-----Move physical data
      call amr_exchange_work_bflags(nprocs,mype,lnblocks_old, &
                                new_loc,old_loc,new_child)

      call amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)

!---------------------------------------------------------------------
!-----Update lnblocks and deallocate old_localMortUB
      lnblocks = new_lnblocks
      deallocate(old_localMortUB)

      return
      end subroutine
