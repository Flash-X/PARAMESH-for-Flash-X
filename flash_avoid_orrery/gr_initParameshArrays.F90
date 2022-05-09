!!****if* source/Grid/GridMain/paramesh/flash_avoid_orrery/gr_initParameshArrays
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
!! SEE ALSO
!!
!!  gr_initParameshDomainBboxes
!!
!! HISTORY
!!
!!  2003 - 2022 Adapted for FLASH and Flash-X    U of Chicago
!!***
subroutine gr_initParameshArrays(restart,&
                                     &  xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use paramesh_dimensions
   use physicaldata
   use workspace
   use tree
   use paramesh_mpi_interfaces, ONLY : mpi_amr_global_domain_limits,&
                                       amr_morton_process,          &
                                       mpi_amr_boundary_block_info
   use paramesh_interfaces, ONLY : amr_refine_derefine
   use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs
   use gr_specificData, ONLY : gr_gidIsValid
   use gr_interface, ONLY : gr_pmIoTreeMetadataIsValid
   use Logfile_interface, ONLY : Logfile_stampMessage
   use Driver_interface, only: Driver_abort
   use Simulation_interface, ONLY : Simulation_mapIntToStr

   implicit none
#include "constants.h"
#include "Simulation.h"
   logical,intent(IN) :: restart
   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary

   integer :: i
   character(len=4) :: vname

   call mpi_amr_global_domain_limits()

! Make sure that blocks produced by divide_domain are in strict
! morton order
   if(restart) then
      if (.NOT. (gr_pmIoTreeMetadataIsValid() .OR. gr_gidIsValid)) then
         ! This will happen in particular if we are restarting from a
         ! checkpoint that was written by a run that used the Amrex
         ! Grid implementation.
         ! If only the surr_blks data was bad, then we can recover
         ! by a amr_morton_process call below and so do not need to
         ! abort here.
         call Logfile_stampMessage("Invalid Grid tree metadata after reading checkpoint;&
              & you may need to recompile with the Bittree feature if you want to restart&
              & from this checkpoint with a PARAMESH Grid.")
         call Driver_abort("The Grid tree metadata in the checkpoint is insufficient&
              & for restarting with this Grid implementation.")
      end if
      call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                    &  ylboundary, yrboundary, &
                                    &  zlboundary, zrboundary)
      call amr_morton_process()
   else
      call gr_initParameshDomainBboxes(xlboundary, xrboundary, &
                                    &  ylboundary, yrboundary, &
                                    &  zlboundary, zrboundary)
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

