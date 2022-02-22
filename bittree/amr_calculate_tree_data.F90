!!****if* source/Grid/GridMain/AMR/Paramesh4/bittree/amr_calculate_tree_data
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
!!   amr_calculate_tree_data
!!
!! SYNOPSIS
!!
!!   call amr_calculate_tree_data (nprocs,mype,lnblocks_old)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: nprocs, mype
!!     number of procs. and proc. id.
!!   integer, intent(in) :: lnblocks_old
!!     old number of local blocks
!!
!! DESCRIPTION
!!
!!   The main workhorse of the Bittree implementation. This routine
!!   is called after the refinement pattern is settled, so Bittree
!!   call properly calculate tree data of new blocks. Unlike in the
!!   regular (non-bittree) implementation, no tree data (except bflags?)
!!   is actually communicated between processors. This routine just
!!   gets the new list of local blocks, and one-by-one computes
!!   their updated tree data. 
!!
!! NOTE
!!   In theory this routine could lead to slight rounding errors in
!!   bnd_box, bsize, coord, but this has not yet been observed.
!!   In general, bittree has not be well-tested on long runs with
!!   high levels of refinement.
!!
!!
!! HISTORY
!!
!!  2021 - 2022  Tom Klosterman
!!***

#include "paramesh_preprocessor.fh"

      subroutine amr_calculate_tree_data(nprocs,mype,lnblocks_old)

!-----Use statements
      use paramesh_dimensions
      use gr_specificData, ONLY : gr_nBlockX, gr_nBlockY, gr_nBlockZ
      use physicaldata
      use tree
      use io
      use bittree, only:gr_btIdentify, gr_btLocate, &
                 gr_btGetLocalBitids, gr_btGetRefine, &
                 bittree_is_parent, gr_btIsParent, &
                 gr_getIntCoords, gr_getNeighIntCoords
      use iso_c_binding, only: c_int,c_bool

      implicit none

!-----Include statements
      include "mpif.h"

!-----Input/output variables
      integer,intent(in)    :: nprocs,mype,lnblocks_old

!-----Local variables, arrays
      integer :: i,j,k
      integer :: li,lj,lk

      integer :: bitid_list(new_lnblocks)
      integer, dimension(mdim) :: gCell, gr
      integer :: lcoord(3),parCoord(3),childCoord(3),neighCoord(3)
      integer :: lev, cid, nid
      integer :: proc, locblk
      integer :: faceAxis, faceSide, irf
      logical(c_bool) :: is_par,c_is_par,n_is_par
      logical :: marked, logpar
      logical(c_bool) :: ctrue,cfalse
      logical :: forced, allChildRefine
      real :: grid_xyzmin(3),default_max(3)

      ctrue = logical(.TRUE.,c_bool)
      cfalse = logical(.FALSE.,c_bool)
     
!-----Get new list of local blocks
      call gr_btGetLocalBitids(mype,new_lnblocks, &
                                        bitid_list,updated=.TRUE.)
     
!-----Loop over new blocks
      do i=1,new_lnblocks
        
        call gr_btLocate(bitid_list(i), lev, lcoord, updated=.TRUE.)
        
        refine(i) = .FALSE.
        derefine(i) = .FALSE.

!--------------------------------------------------------------------
!-------Set lrefine
        lrefine(i) = lev

!--------------------------------------------------------------------
!-------Set bsize
        bsize(1,i) = (grid_xmax-grid_xmin)/real(ishft(gr_nBlockX,lev-1))
        bsize(2,i) = (grid_ymax-grid_ymin)/real(ishft(gr_nBlockY,lev-1))
        bsize(3,i) = (grid_zmax-grid_zmin)/real(ishft(gr_nBlockZ,lev-1))

!--------------------------------------------------------------------
!-------Set bnd_box
        if (ndim.eq.3) then
          default_max = (/grid_xmin,grid_ymin,grid_zmin/)
        elseif (ndim.eq.2) then
          default_max = (/grid_xmin,grid_ymin,grid_zmax/)
        else 
          default_max = (/grid_xmin,grid_ymax,grid_zmax/)
        end if

        grid_xyzmin = (/grid_xmin,grid_ymin,grid_zmin/)
        bnd_box(1,:,i) = real(lcoord)*bsize(:,i) + grid_xyzmin
        bnd_box(2,:,i) = max( real(lcoord+1)*bsize(:,i) + grid_xyzmin, &
                              default_max )

!--------------------------------------------------------------------
!-------Set coord
        coord(:,i) = (bnd_box(2,:,i) + bnd_box(1,:,i) )/2.0

!--------------------------------------------------------------------
!-------Set parent.
        if (lrefine(i).gt.1) then
          lev = lrefine(i)-1
          parCoord = lcoord/2
          call gr_btIdentify(nprocs,lev,parCoord,proc, &
                                  locblk,updated=.TRUE.)
          parent(1,i) = locblk
          parent(2,i) = proc
        else
          parent(:,i) = -1
        end if

!--------------------------------------------------------------------
!-------Set newchild by checking if parent WAS marked for refine.
        if(parent(1,i).gt.0) then
          call gr_btGetRefine(lrefine(i)-1,lcoord/2,marked)
          newchild(i) = marked
        else
          newchild(i) = .FALSE.
        end if

!--------------------------------------------------------------------
!-------Set which_child
        if(parent(1,i).gt.0) then
          which_child(i) = sum(mod(lcoord,2)*(/1,2,4/)) + 1
        else
          which_child(i) = -1
        end if

!--------------------------------------------------------------------
!-------Loop over children and set child.
        call bittree_is_parent(ctrue,int(bitid_list(i),c_int),is_par)
        allChildRefine = .TRUE.
        do j=1,nchild
          if (is_par) then
            lev = lrefine(i) + 1
            childCoord = lcoord*2 + (/ mod((j-1),2),         &
                                       mod((j-1)/2,2),       &
                                       mod((j-1)/4,2) /)
            call gr_btIdentify(nprocs,lev,childCoord,proc,&
                                      locblk,updated = .TRUE.,bitid=cid)
            child(1,j,i) = locblk
            child(2,j,i) = proc
            
            if(allChildRefine) then
              call bittree_is_parent(ctrue, int(cid,c_int), c_is_par)
              if (.NOT.c_is_par) allChildRefine=.FALSE.
            end if
          else
            child(1:2,j,i) = -1 
            allChildRefine = .FALSE.
          end if
        end do

!--------------------------------------------------------------------
!-------Set nodetype based on whether parent and whether all children
!-------are refined.
        if(is_par) then
          if (allChildRefine) then
            nodetype(i) = 3
          else
            nodetype(i) = 2
          end if
        else
          nodetype(i) = 1
        end if

!--------------------------------------------------------------------
!-------Set empty
        empty(i) = 0

      end do !i=1,new_lnblocks
          

!--------------------------------------------------------------------
!-----Set surr_blks. New loop so previously calculated nodetypes
!-----can be used.
      do i=1,new_lnblocks
        do lk = -k3d,k3d
          do lj = -k2d,k2d
            do li = -k1d,k1d
              gCell = (/li,lj,lk/)
              gr = gCell+(/k1d+1,k2d+1,k3d+1/)
              
              if ((0.ne.sum(abs(gCell)))) then
                
                lev = lrefine(i)
                call gr_getNeighIntCoords(i,gCell,neighCoord,.TRUE.)

!---------------If neighCoords are negative, fill surr_blks with BC
                if (any(neighCoord.lt.0)) then
                  surr_blks(3,gr(1),gr(2),gr(3),i) = -1
                  surr_blks(1:2,gr(1),gr(2),gr(3),i)=minval(neighCoord)

!---------------Else try to identify neighbor
                else
                  call gr_btIdentify(nprocs,lev,neighCoord,proc,&
                                        locblk,updated=.TRUE.,bitid=nid)
                
!-----------------If valid neighbor, set surr_blks accordingly
                  if (lev.eq.lrefine(i) ) then
                
                    surr_blks(1,gr(1),gr(2),gr(3),i) = locblk
                    surr_blks(2,gr(1),gr(2),gr(3),i) = proc
                
!-------------------If neigh is on processor, grab its nodetype.
!-------------------Otherwise, recalculate it.
                    if (proc.eq.mype) then
                      surr_blks(3,gr(1),gr(2),gr(3),i) = &
                        nodetype(locblk)
                    else
                      
                      call bittree_is_parent(ctrue,int(nid,c_int), &
                                             n_is_par)
                      if(n_is_par) then
                        allChildRefine = .TRUE.
                        do j=1,nchild
                          childCoord = neighCoord*2   +        &
                                       (/ mod((j-1),2),        &
                                       mod((j-1)/2,2),         &
                                       mod((j-1)/4,2) /)
                          if(allChildRefine) then
                            call gr_btIsParent(lrefine(i)+1, &
                                                       childCoord, &
                                                       logpar,.TRUE.)
                            if (.NOT.logpar) allChildRefine=.FALSE.
                          end if !allChildRefine
                        end do !j=1,nchild

                        if (allChildRefine) then
                          surr_blks(3,gr(1),gr(2),gr(3),i) = 3
                        else
                          surr_blks(3,gr(1),gr(2),gr(3),i) = 2
                        end if !allChildRefine
                      else
                        surr_blks(3,gr(1),gr(2),gr(3),i) = 1
                      end if !n_is_par
                      
                    end if !proc.eq.mype

                  else
                    surr_blks(:,gr(1),gr(2),gr(3),i) = -1
                  
                  end if !lev=lrefine(i)

                end if !any(neighCoord<0)



!---------------Copy into neigh
                if ((1.eq.sum(abs(gCell)))) then
                  faceAxis = sum((/1,2,3/) * abs(gCell)) !x:1, y:2, z:3
                  faceSide = sum(gCell)   !low side: -1, high side: +1
                  irf = faceAxis*2 - dim(0,faceSide)
                  neigh(1:2,irf,i) = surr_blks(1:2,gr(1),gr(2),gr(3),i)
                end if

              else !self
                surr_blks(1,gr(1),gr(2),gr(3),i) = i
                surr_blks(2,gr(1),gr(2),gr(3),i) = mype
                surr_blks(3,gr(1),gr(2),gr(3),i) = nodetype(i)
              
              end if !sum(abs(gCell))==0
            
            end do !iAxis
          end do !jAxis
        end do !kAxis
      end do !i=1,new_lnblocks

!--------------------------------------------------------------------
!--------------------------------------------------------------------

!-----Clear out old values
      if (lnblocks_old.gt.new_lnblocks) then
      do i = new_lnblocks+1,lnblocks_old
        derefine(i) = .FALSE.
        refine(i) = .FALSE.
        lrefine(i) = -1
        bsize(:,i) = -1.
        bnd_box(1:2,1:ndim,i) = -1.
        coord(1:ndim,i) = -1.
        parent(1:2,i) = -1
        newchild(i) = .FALSE.
        which_child(i) = -1
        nodetype(i) = -1
        child(1:2,1:nchild,i) = -1
        surr_blks(1:3,:,:,:,i) = -1
        neigh(1:2,:,i) = -1
      end do !i = new_lnblocks+1,lnblocks_old
      end if

      return
      end subroutine
