!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/gr_initParameshDomainBboxes
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
!!  gr_initParameshDomainBboxes
!!
!! SYNOPSIS
!!
!!  call gr_initParameshDomainBboxes(integer(IN) :: xlboundary,
!!                                   integer(IN) :: xrboundary,
!!                                   integer(IN) :: ylboundary,
!!                                   integer(IN) :: yrboundary,
!!                                   integer(IN) :: zlboundary,
!!                                   integer(IN) :: zrboundary)
!!
!! DESCRIPTION
!!
!!  Initialize 'boundary_box' and 'boundary_info' arrays for the
!!  domain boundaries.
!!
!!  The 'boundary_box' and 'boundary_info' arrays are used by Paramesh4
!!  to keep track of boundary conditions and where they apply.
!!  The first 2*NDIM slots in these arrays are reserved for information
!!  on external boundaries. (If there are obstacle blocks within the
!!  domain, then further slots may be used for those, see
!!  Simulation_defineDomain.)
!!
!!  The 'boundary_box' and 'boundary_info' arrays are used within
!!  the Paramesh4 routine find_surrblks.
!!
!! ARGUMENTS
!!
!!   xlboundary - boundary condition type of outer domain boundary in lower X direction.
!!   xrboundary - boundary condition type of outer domain boundary in upper X direction.
!!   ylboundary - boundary condition type of outer domain boundary in lower Y direction.
!!   yrboundary - boundary condition type of outer domain boundary in upper Y direction.
!!   zlboundary - boundary condition type of outer domain boundary in lower Z direction.
!!   zrboundary - boundary condition type of outer domain boundary in upper Z direction.
!!
!! SEE ALSO
!!
!!  Simulation_defineDomain
!!  Grid_initDomain
!!  find_surrblks
!!
!! HISTORY
!!
!!  2003 - 2022 Adapted for FLASH and Flash-X    U of Chicago
!!***
subroutine gr_initParameshDomainBboxes( xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use paramesh_dimensions, ONLY: ndim
!   use physicaldata
!   use workspace
   use tree, ONLY: boundary_index, boundary_box, &
                   grid_xmin, grid_xmax, grid_ymin, grid_ymax, &
                   grid_zmin, grid_zmax
   use paramesh_mpi_interfaces, ONLY : mpi_amr_global_domain_limits
   use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs, &
                         gr_imin, gr_imax, gr_jmin, gr_jmax, &
                         gr_kmin, gr_kmax
   use gr_specificData, ONLY : gr_nblockX, gr_nblockY, gr_nblockZ

   implicit none
#include "constants.h"

   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary
   integer,parameter :: five=5,six=6

   ! x boundaries
   boundary_box(1,2:3,1:2) = -1.e30
   boundary_box(2,2:3,1:2) =  1.e30
   boundary_box(1,1,1) = -1.e30
   boundary_box(2,1,1) = grid_xmin
   boundary_box(1,1,2) = grid_xmax
   boundary_box(2,1,2) = 1.e30
   boundary_index(1) = xlboundary
   boundary_index(2) = xrboundary
   if (grid_xmin - gr_imin > 0.5*gr_nblockX*(gr_imax-gr_imin)) then ! boundary blocks all along left domain face
#ifdef DEBUG_GRID
      print*,'****** LOW I *****',grid_xmin,gr_imin,0.5*gr_nblockX,(gr_imax-gr_imin)
#endif
      boundary_box(1:2,1:3,1) = 0.
      boundary_index(1) = 0
   end if
   if (gr_imax - grid_xmax > 0.5*gr_nblockX*(gr_imax-gr_imin)) then ! boundary blocks all along right domain face
#ifdef DEBUG_GRID
      print*,'****** HIGH I *****',grid_xmin,gr_imin,0.5*gr_nblockX,(gr_imax-gr_imin)
#endif
      boundary_box(1:2,1:3,2) = 0.
      boundary_index(2) = 0
   end if
   if (xlboundary.eq.PERIODIC) then ! periodic
      boundary_box(1,2:3,1:2) = 0.
      boundary_box(2,2:3,1:2) =  0.
      boundary_box(1,1,1) = 0.
      boundary_box(2,1,1) = 0.
      boundary_box(1,1,2) = 0.
      boundary_box(2,1,2) = 0.
      boundary_index(1) = 0
    end if
    if (xrboundary.eq.PERIODIC) then ! periodic
      boundary_box(1,2:3,1:2) = 0.
      boundary_box(2,2:3,1:2) =  0.
      boundary_box(1,1,1) = 0.
      boundary_box(2,1,1) = 0.
      boundary_box(1,1,2) = 0.
      boundary_box(2,1,2) = 0.
      boundary_index(2) = 0
    end if
! y boundaries
    if(ndim.ge.2) then

      boundary_box(1,1,3:4) = -1.e30
      boundary_box(2,1,3:4) =  1.e30
      boundary_box(1,3,3:4) = -1.e30
      boundary_box(2,3,3:4) =  1.e30
      boundary_box(1,2,3) = -1.e30
      boundary_box(2,2,3) = grid_ymin
      boundary_box(1,2,4) = grid_ymax
      boundary_box(2,2,4) = 1.e30
      boundary_index(3) = ylboundary
      boundary_index(4) = yrboundary
      if (grid_ymin - gr_jmin > 0.5*gr_nblockY*(gr_jmax-gr_jmin)) then ! boundary blocks all along south domain face
         boundary_box(1:2,1:3,3) = 0.
         boundary_index(3) = 0
      end if
      if (gr_jmax - grid_ymax > 0.5*gr_nblockY*(gr_jmax-gr_jmin)) then ! boundary blocks all along north domain face
         boundary_box(1:2,1:3,4) = 0.
         boundary_index(4) = 0
      end if
      if (ylboundary.eq.PERIODIC) then ! periodic
         boundary_box(1,1,3:4) = 0.
         boundary_box(2,1,3:4) = 0.
         boundary_box(1,3,3:4) = 0.
         boundary_box(2,3,3:4) = 0.
         boundary_box(1,2,3) = 0.
         boundary_box(2,2,3) = 0.
         boundary_box(1,2,4) = 0.
         boundary_box(2,2,4) = 0.
         boundary_index(3) = 0
      end if
      if (yrboundary.eq.PERIODIC) then ! periodic
         boundary_box(1,1,3:4) = 0.
         boundary_box(2,1,3:4) = 0.
         boundary_box(1,3,3:4) = 0.
         boundary_box(2,3,3:4) = 0.
         boundary_box(1,2,3) = 0.
         boundary_box(2,2,3) = 0.
         boundary_box(1,2,4) = 0.
         boundary_box(2,2,4) = 0.
         boundary_index(4) = 0
       end if
    endif
! z boundaries
    if(ndim.ge.3) then
       boundary_box(1,1:2,five:six) = -1.e30
       boundary_box(2,1:2,five:six) =  1.e30
       boundary_box(1,3,five) = -1.e30
       boundary_box(2,3,five) = grid_zmin
       boundary_box(1,3,six) = grid_zmax
       boundary_box(2,3,six) = 1.e30
       boundary_index(five) = zlboundary
       boundary_index(six) = zrboundary
       if (grid_zmin - gr_kmin > 0.5*gr_nblockZ*(gr_kmax-gr_kmin)) then ! boundary blocks all along bottom domain face
          boundary_box(1:2,1:3,5) = 0.
          boundary_index(5) = 0
       end if
       if (gr_kmax - grid_zmax > 0.5*gr_nblockZ*(gr_kmax-gr_kmin)) then ! boundary blocks all along top domain face
          boundary_box(1:2,1:3,6) = 0.
          boundary_index(6) = 0
       end if
       if (zlboundary.eq.PERIODIC) then ! periodic
          boundary_box(1,1:2,five:six) = 0.
          boundary_box(2,1:2,five:six) = 0.
          boundary_box(1,3,five) = 0.
          boundary_box(2,3,five) = 0.
          boundary_box(1,3,six) = 0.
          boundary_box(2,3,six) = 0.
          boundary_index(five) = 0
       end if
       if (zrboundary.eq.PERIODIC) then ! periodic
          boundary_box(1,1:2,five:six) = 0.
          boundary_box(2,1:2,five:six) =  0.
          boundary_box(1,3,five) = 0.
          boundary_box(2,3,five) = 0.
          boundary_box(1,3,six) = 0.
          boundary_box(2,3,six) = 0.
          boundary_index(six) = 0.
       end if
    end if
  
  return
end subroutine gr_initParameshDomainBboxes

