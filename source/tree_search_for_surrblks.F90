!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------


  subroutine tree_search_for_surrblks ()

      use local_tree_common
      use physicaldata
      use tree
      use paramesh_dimensions
      use paramesh_interfaces
      use paramesh_comm_data
      use mpi_morton, only : lperiodicx, lperiodicy, lperiodicz
      use constants

      implicit none

      include 'mpif.h'

      real :: neigh_coord(3), neigh_coord2(3)
      integer :: ierr, nprocs, mype
      integer :: neigh_lb,neigh_proc, neigh_nodetype
      integer :: i,j,k,ii,jj,kk,iboun,lb

      logical :: found

      real :: time_exe, time_max
      real :: eps,accuracy

      accuracy = 100./10.**precision(accuracy)
      if (accuracy > 1.0e-10) then
         eps = 1.e-10
      else
         eps = accuracy
      end if

      time_exe = MPI_WTIME()

      call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
      call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

! FIND SURROUNDING BLOCKS
      do lb = 1, lnblocks

         kk = -k3d
         do k = 1,1+2*k3d

         neigh_coord(3) = coord(3,lb) + kk*bsize(3,lb)
         if (lperiodicz.and.neigh_coord(3).lt.grid_zmin)               &
          neigh_coord(3) = neigh_coord(3) + (grid_zmax-grid_zmin)
         if (lperiodicz.and.neigh_coord(3).gt.grid_zmax)               &
          neigh_coord(3) = neigh_coord(3) - (grid_zmax-grid_zmin)

         jj = -k2d
         do j = 1,1+2*k2d

         neigh_coord(2) = coord(2,lb) + jj*bsize(2,lb)
         if (lperiodicy.and.neigh_coord(2).lt.grid_ymin)               &
          neigh_coord(2) = neigh_coord(2) + (grid_ymax-grid_ymin)
         if (lperiodicy.and.neigh_coord(2).gt.grid_ymax)               &
          neigh_coord(2) = neigh_coord(2) - (grid_ymax-grid_ymin)

         ii = -1
         do i = 1,3

         neigh_coord(1) = coord(1,lb) + ii*bsize(1,lb)
         if (lperiodicx.and.neigh_coord(1).lt.grid_xmin)               &
          neigh_coord(1) = neigh_coord(1) + (grid_xmax-grid_xmin)
         if (lperiodicx.and.neigh_coord(1).gt.grid_xmax)               &
          neigh_coord(1) = neigh_coord(1) - (grid_xmax-grid_xmin)

         neigh_coord2(:) = neigh_coord(:)

!--------Reset coordinates of neighbor in spherical coordinates
         If (spherical_pm) Then
            If ((((jj == -1).and.(abs(bnd_box(1,2,lb)) < eps)) .or.    &  
                 ((jj ==  1).and.(abs(bnd_box(2,2,lb)-pi) < eps)))     & 
                 .and. lsingular_line ) Then
               neigh_coord2(2) = coord(2,lb)
               If (neigh_coord2(3) < pi) Then
                  neigh_coord2(3) = neigh_coord2(3) + pi
               Elseif(neigh_coord2(3) > pi)  Then
                  neigh_coord2(3) = neigh_coord2(3) - pi
               End If
            End If
         End If

         neigh_lb = -1
         neigh_proc = -1
         neigh_nodetype = -1
         found = .false.

!--------First check boundaries

         do iboun = 1, nboundaries
            if (ndim == 1) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            elseif (ndim == 2) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
                neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
                neigh_coord2(2) < boundary_box(2,1+k2d,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            elseif (ndim == 3) then
            if (neigh_coord2(1) > boundary_box(1,1,iboun) .and.        &
                neigh_coord2(1) < boundary_box(2,1,iboun) .and.        &
                neigh_coord2(2) > boundary_box(1,1+k2d,iboun) .and.    &
                neigh_coord2(2) < boundary_box(2,1+k2d,iboun) .and.    &
                neigh_coord2(3) > boundary_box(1,1+2*k3d,iboun) .and.  &
                neigh_coord2(3) < boundary_box(2,1+2*k3d,iboun)) then
               found = .true.
               neigh_lb = boundary_index(iboun)
               neigh_proc = boundary_index(iboun)
            end if
            end if
         end do

         call search_sub_tree(local_tree,neigh_coord2,lrefine(lb),     &
                              neigh_lb,neigh_proc,neigh_nodetype,found)
         surr_blks(1,i,j,k,lb) = neigh_lb
         surr_blks(2,i,j,k,lb) = neigh_proc
         surr_blks(3,i,j,k,lb) = neigh_nodetype

         ii = ii + 1
      end do
      jj = jj + k2d
      end do
      kk = kk + k3d
      end do

      end do

      end subroutine tree_search_for_surrblks
