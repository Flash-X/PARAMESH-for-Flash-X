!!****if* source/Grid/GridMain/AMR/Paramesh4/bittree/gr_initParameshArrays
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
!!  gr_initParameshArrays
!!
!! SYNOPSIS
!!
!!  call gr_initParameshArrays(logical(IN) :: restart,
!!                             integer(IN) :: xlboundary,
!!                             integer(IN) :: xrboundary,
!!                             integer(IN) :: ylboundary,
!!                             integer(IN) :: yrboundary,
!!                             integer(IN) :: zlboundary,
!!                             integer(IN) :: zrboundary
!!                             )
!!
!! DESCRIPTION
!!
!!  Perform early initialization of some Grid data structures.
!!
!!  This routine prepares the Grid for being filled with
!!  meaningful data.
!!
!!  This Bittree version also initializes Bittree.
!!
!! ARGUMENTS
!!
!!   restart -   Is the grid being prepared for initialization with
!!               data from a checkpoint file?
!!   xlboundary - boundary condition type of outer domain boundary in lower X direction.
!!   xrboundary - boundary condition type of outer domain boundary in upper X direction.
!!   ylboundary - boundary condition type of outer domain boundary in lower Y direction.
!!   yrboundary - boundary condition type of outer domain boundary in upper Y direction.
!!   zlboundary - boundary condition type of outer domain boundary in lower Z direction.
!!   zrboundary - boundary condition type of outer domain boundary in upper Z direction.
!!
!! HISTORY
!!
!!  2003 - 2022 Adapted for FLASH and Flash-X and Bittree   U of Chicago
!!
!! SEE ALSO
!!
!!  gr_initParameshDomainBboxes
!!  gr_pmIoTreeMetadataIsValid
!!
!!***

subroutine gr_initParameshArrays(restart,&
                                     &  xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use paramesh_dimensions, ONLY: ndim, nvar
   use physicaldata, ONLY : surr_blks_valid, interp_mask_unk,&
                            lsingular_line, spherical_pm, polar_pm,&
                            use_flash_surr_blks_fill
   use workspace,    ONLY : interp_mask_work
   use tree,         ONLY : newchild, grid_analysed_mpi
   use paramesh_mpi_interfaces, ONLY : mpi_amr_global_domain_limits,&
                                       amr_morton_process,          &
                                       mpi_amr_boundary_block_info
   use paramesh_interfaces, ONLY : amr_refine_derefine
   use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
   use gr_interface, ONLY : gr_pmIoTreeMetadataIsValid
   use Logfile_interface, ONLY : Logfile_stampMessage
   use Driver_interface, only: Driver_abort
   use Simulation_interface, ONLY : Simulation_mapIntToStr

   implicit none
#include "constants.h"

   logical,intent(IN) :: restart
   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary

   integer :: i
   character(len=4) :: vname

   call mpi_amr_global_domain_limits()

   if(.NOT.use_flash_surr_blks_fill) &
     call Driver_abort("Error in initializing Bittree. Bittree is &
                       &only appropriate if use_flash_surr_blks_fill=True")
   if(lsingular_line .and. ndim > 1 .and. (spherical_pm .or. polar_pm)) &
     call Driver_abort("Error in initializing Bittree. Bittree is &
                       &not appropriate for spherical_pm, polar_pm")
   newchild(:) = .FALSE.
   call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                 &  ylboundary, yrboundary, &
                                 &  zlboundary, zrboundary)
   call amr_build_bittree()
   if (restart) then

      if (.NOT. gr_pmIoTreeMetadataIsValid()) then
         ! In this case, amr_refine_derefine must be called
         ! (with force_rebalance true) to properly initialize
         ! PARAMESH metainformation arrays.
         ! This will happen in particular if we are restarting from a
         ! checkpoint that was written by a run that used the Amrex
         ! Grid implementation.
         call amr_refine_derefine(force_rebalance=.TRUE.)
      else
         call amr_morton_process()
      endif

   else
     !call amr_reorder_grid() !needed when bittree used different order?
     call amr_refine_derefine()
   end if


  !CD: Inform PARAMESH that we no longer need orrery.
  surr_blks_valid = .true.


  if(restart) then
     grid_analysed_mpi=1
#ifdef CALL_BOUNDARY_BLOCK_INFO
     call mpi_amr_boundary_block_info(gr_meshMe, gr_meshNumProcs)
#endif
  end if
  ! reset for quadratic interpolation
  
  interp_mask_work(:) = 1

  ! AH: do not reset interp_mask_unk, but warn if not one
  !interp_mask_unk(:) = 1
  do i = 1, nvar
     if (interp_mask_unk(i) .NE. 1) then
        call Simulation_mapIntToStr(i,vname,MAPBLOCK_UNK)
        call Logfile_stampMessage("WARNING: interp_mask_unk  is not 1 for variable "//vname)
     end if
  end do
  
  return

end subroutine gr_initParameshArrays
