!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_btSortMortonBittree.F90
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
!!  gr_btSortMortonBittree
!!
!! SYNOPSIS
!!
!!  call gr_btSortMortonBittree(nprocs, mype)
!!
!! DESCRIPTION
!!
!!  Sorts blocks by determining which section of the morton index
!!  will be distributed to each proc. Stores results in localMortUB.
!!  Must be called before blocks can be identified.
!!
!! ARGUMENTS
!!
!!  integer,intent(in)  :: nprocs
!!  integer,intent(in)  :: mype
!!  logical,intent(in),optional :: sort_by_work (overrides gr_btSortByWork)
!!
!!***

#include "paramesh_preprocessor.fh"

      subroutine gr_btSortMortonBittree(nprocs,mype,sort_by_work)

!-----Use statements
      use tree
      use bittree, only: bittree_block_count,localMortUB, &
                         gr_btDistributedSort
      use iso_c_binding, only: c_int,c_bool

      implicit none

!-----Input/output statements
      integer,intent(in) :: nprocs, mype
      logical,intent(in),optional :: sort_by_work

!-----Local variables and arrays
      integer :: j, over, over_blks
      integer(c_int) :: totblocks
      logical(c_bool) :: ctrue = .TRUE.
      logical :: lworksort

!-----Begin executable code
      if(present(sort_by_work)) then
        lworksort = sort_by_work
      else
        lworksort = gr_btSortByWork
      end if

      if(.NOT.lworksort) then
!-------If .NOT.lworksort, blocks distributed equally.
!-------Procs with j<=mod(totblocks,nprocs) have one extra block.
        call bittree_block_count(ctrue,totblocks)
        over = mod(totblocks,nprocs)
        over_blks = totblocks/nprocs + 1
        do j=1,nprocs
          if (j.lt.over) then
            localMortUB(j)=over_blks*j
          else
            localmortUB(j)=over_blks*over + totblocks/nprocs*(j-over)
          end if
        end do !1,nprocs

      else
!-------Sort by work, determined by runtime parameters
        if(gr_btCustomWork) then
          if(gr_btDistributedSort) then
            call bittree_sort_custom_distributed(nprocs,mype)
          else
            call bittree_sort_custom_consolidated(nprocs,mype)
          end if
        else
          call bittree_sort_nodetype(nprocs,mype)
        end if

      end if !.not.lworksort

      return
      end subroutine gr_btSortMortonBittree

!!****
!!
!! NAME
!!
!!  bittree_sort_nodetype
!!
!! SYNOPSIS
!!
!!  call bittree_sort_nodetype(nprocs, mype)
!!
!! DESCRIPTION
!!
!!  Sorts blocks by work. Each block is assigned a default work
!!  based on nodetype, which is determined by runtime parameter
!!
!! ARGUMENTS
!!
!!  integer,intent(in)  :: nprocs
!!  integer,intent(in)  :: mype
!!
!!***
      subroutine bittree_sort_nodetype(nprocs,mype)

!-----Use statements
      use gr_sortByWorkTools, only: gr_sortByWorkConsolidated
      use tree
      use bittree, only: bittree_block_count,localMortUB, &
                         bittree_get_bitid_list,&
                         bittree_is_parent
      use iso_c_binding, only: c_int,c_bool

      implicit none

!-----Input/output statements
      integer,intent(in) :: nprocs, mype

!-----Local variables and arrays
      integer :: i,totblocks
      logical(c_bool) :: ctrue = .TRUE.
      real,allocatable :: work(:)
      integer(c_int) :: totblocksc
      integer(c_int),allocatable :: bitids(:)
      logical(c_bool) :: is_par

!-----Get list of bitids corresponding to whole morton index
      call bittree_block_count(ctrue,totblocksc)
      totblocks = totblocksc
      allocate(bitids(totblocksc))
      call bittree_get_bitid_list(ctrue,0,totblocksc,bitids)

!-----Create work array for global block list, by nodetype
      allocate(work(totblocks))
      work(:) = 0.
      do i=1,totblocks
         call bittree_is_parent(ctrue,bitids(i),is_par)
         if(is_par) then
           work(i) = gr_btWorkDefaultPar
         else
           work(i) = gr_btWorkDefaultLeaf
         end if
      end do
      deallocate(bitids)

!-----Call consolidated sort-by-work routine
      call gr_sortByWorkConsolidated(nprocs,mype,work,&
                                         totblocks,localMortUB)

      deallocate(work)
      end subroutine

!!****
!!
!! NAME
!!
!!  bittree_sort_custom_consolidated
!!
!! SYNOPSIS
!!
!!  call bittree_sort_custom_consolidated(nprocs, mype)
!!
!! DESCRIPTION
!!
!!  Sorts blocks by work. Each block is assigned work
!!  based on the real-valued array work_block, which is a
!!  tree module member. New children are assigned work based
!!  on their parents. By default the scaling_factor is 1., but
!!  can be changed by the user.
!!
!! ARGUMENTS
!!
!!  integer,intent(in)  :: nprocs
!!  integer,intent(in)  :: mype
!!
!!***
      subroutine bittree_sort_custom_consolidated(nprocs,mype)

!-----Use statements
      use gr_sortByWorkTools, only: gr_sortByWorkConsolidated
      use tree
      use bittree, only: bittree_block_count,localMortUB, &
                         bittree_identify, gr_getIntCoords, &
                         bittree_is_parent, &
                         gr_btGetRefine,gr_btIsParent
      use paramesh_comm_data, ONLY : amr_mpi_meshComm, amr_mpi_real
      use iso_c_binding, only: c_int,c_bool

      implicit none

!-----Include statements
      include "mpif.h"

!-----Input/output statements
      integer,intent(in) :: nprocs, mype

!-----Local variables and arrays
      integer :: i,j,totblocks,ierr
      real,allocatable :: work(:)
      logical(c_bool) :: ctrue = .TRUE.
      integer(c_int):: lev,mort,bitid
      integer :: lcoord(3),lcoordc(3),childCoord(3),childCoordc(3)
      integer(c_int) :: totblocksc
      logical :: new_par
      logical(c_bool) :: is_par

      call bittree_block_count(ctrue,totblocksc)
      totblocks = totblocksc

!-----Put local blocks and their children into their (global)
!-----position in work array
      allocate(work(totblocks))
      work(:) = 0.
      do i=1,lnblocks
        call gr_getIntCoords(i,lcoord)
        lcoordc = int(lcoord,c_int)
        lev = int(lrefine(i)-1,c_int)
!-------Identify new loctation of block in global morton index
        call bittree_identify(ctrue,lev,lcoordc,mort,bitid)

        if(lev.eq.int((lrefine(i)-1),c_int)) then
!---------Check updated nodetype
          call bittree_is_parent(ctrue,bitid,is_par)

!---------Set work based on work_block, or to default if nodetype change
          if(is_par .eqv. (nodetype(i).gt.1)) then
            !nodetype unchanged
            work(mort+1) = work_block(i)
          else if (is_par .AND. (nodetype(i).eq.1)) then
            !new parent
            work(mort+1) = gr_btWorkDefaultPar
          else if (.NOT.is_par .AND. (nodetype(i).gt.1)) then
            !new leaf
            work(mort+1) = gr_btWorkDefaultLeaf
          end if

!---------If a new parent, set work for new children
          if(is_par .AND. (nodetype(i).eq.1)) then
            lev = lev + 1
            do j=1,nchild
              childCoord = lcoord*2 + (/ mod((j-1),2),         &
                                         mod((j-1)/2,2),       &
                                         mod((j-1)/4,2) /)
              childCoordc = int(childCoord,c_int)
              call bittree_identify(ctrue,lev,childCoordc,mort,bitid)
              work(mort+1) = work_block(i) * gr_btWorkChildScaling
            end do !j=1,nchild
          end if !new_par

        end if !lev.eq.(lrefine(i)-1)
      end do !i=1,lnblocks

!-----Sum work() across all processors in place. This creates a list with
!-----the work of all global blocks, ordered by morton index.
      call MPI_ALLREDUCE(MPI_IN_PLACE, work, totblocks, amr_mpi_real, &
                         MPI_SUM, amr_mpi_meshComm, ierr)


!-----Call consolidated sort-by-work routine
      call gr_sortByWorkConsolidated(nprocs,mype,work,&
                                         totblocks,localMortUB)
      deallocate(work)
      end subroutine

!!****
!!
!! NAME
!!
!!  bittree_sort_custom_distributed
!!
!! SYNOPSIS
!!
!!  call sort_by_work_optimize(nprocs, mype)
!!
!! DESCRIPTION
!!
!!  Sorts blocks by work. Each block is assigned work
!!  based on the real-valued array work_block, which is a
!!  tree module member. New children are assigned work based
!!  on their parents. By default the scaling_factor is 1., but
!!  can be changed by the user.
!!
!!  Unlike the above routine, this routine optimizes for memory.
!!  Instead of storing the work of every single block in existence,
!!  this routine only allocates enough memory to store the work
!!  of blocks on each proc in sequence.
!!
!! ARGUMENTS
!!
!!  integer,intent(in)  :: nprocs
!!  integer,intent(in)  :: mype
!!
!!***

      subroutine bittree_sort_custom_distributed(nprocs,mype)

!-----Use statements
      use gr_sortByWorkTools, only: gr_sortByWorkDistributed
      use tree
      use bittree, only: localMortUB, bittree_block_count, &
                         gr_btIsParent
      use paramesh_comm_data, ONLY : amr_mpi_meshComm, amr_mpi_real
      use iso_c_binding, only: c_int,c_bool

      implicit none

!-----Include statements
      include "mpif.h"

!-----Input/output statements
      integer,intent(in) :: nprocs, mype

!-----Local variables and arrays
      integer :: i,ii,j,ierr
      real    :: worktemp(maxblocks_tr)
      integer :: lnblockst, totblocks
      integer(c_int) :: totblocksc
      integer :: childrendone
      logical(c_bool), save :: ctrue = .TRUE.
      integer :: lcoord(3)
      logical :: is_par

      call bittree_block_count(ctrue,totblocksc)
      totblocks = totblocksc

!-----lnblockst = length of intermediate local block list
      lnblockst = lnblocks - count(derefine(1:lnblocks)) + &
        count(refine(1:lnblocks))*nchild

!-----Fill worktemp with work of local blocks.
!-----This list is "refined in place" (derefined blocks omitted,
!-----new children inserted in to preserve morton order).
      worktemp(:) = 0.
      ii = 1
      do i=1,lnblocks
        if(.NOT.derefine(i)) then

          childrendone=0
#ifdef ALT_MORTON_ORDER
          if(refine(i)) then
            do j=1,nchild/2
              worktemp(ii) = work_block(i) * gr_btWorkChildScaling
              ii = ii+1
              childrendone = childrendone + 1
            end do
          end if !refine(i)
#endif

          call gr_getIntCoords(i,lcoord)
          call gr_btIsParent(lrefine(i),lcoord, &
                                     is_par,updated=.TRUE.)
          if(is_par) then
            worktemp(ii) = gr_btWorkDefaultPar
          else
            if(nodetype(i).eq.1) worktemp(ii) = work_block(i)
            if(nodetype(i).gt.1) worktemp(ii) = gr_btWorkDefaultLeaf
          end if !is_par
          ii = ii+1

          if(refine(i)) then
            do j=childrendone+1,nchild
              worktemp(ii) = work_block(i) * gr_btWorkChildScaling
              ii = ii+1
            end do
          end if !refine(i)

        end if !.NOT.derefine(i)
      end do !1,lnblocks


!-----Call distributed sort-by-work routine
      call gr_sortByWorkDistributed(nprocs,mype,worktemp,&
                                    lnblockst,localMortUB,maxblocks_tr)

      end subroutine
