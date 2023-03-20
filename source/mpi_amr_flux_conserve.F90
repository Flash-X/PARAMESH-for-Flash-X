!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_flux_conserve
!! NAME
!!
!!   amr_flux_conserve
!!
!! SYNOPSIS
!!
!!   Call amr_flux_conserve (mype, nsub)
!!   Call amr_flux_conserve (mype, nsub, flux_dir)
!!
!!   Call amr_flux_conserve (integer, integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The Calling processor.
!!
!!   integer, intent(in) :: nsub          
!!     The current time subcycle. If this is 1 then this info is used to 
!!     reset the temporary boundary flux arrays to 0. This argument only has
!!     an effect if variable time steps are being used.
!!
!!   optional, integer, intent(in) :: flux_dir
!!     Option integer which selects which coordinate direction to apply
!!     the flux conservation operation to:
!!     If flux_dir = 1 -> x direction
!!        flux_dir = 2 -> y direction
!!        flux_dir = 3 -> z direction
!!     If this argument is not specified, then the default behaviour is
!!     to operate on all directions.  Using this argument can be useful
!!     for Strang-split schemes to improve performance.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!
!! CALLS
!! 
!!   amr_flux_conserve_udt
!!   amr_flux_conserve_vdt
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit the fluxes stored in the arrays
!!   flux_x, flux_y, or flux_z are corrected at jumps in refinement.  This
!!   is either an averaging proceedure or a sum as selected by the user
!!   by adjusting the logical variables which control this behaviour
!!   that are read at runtime from the file 'amr_runtime_parameters'.
!!
!! DESCRIPTION
!!
!!   This is a wrapper routine which makes the appropriate Call to the
!!   routines which manage flux conservation at the boundaries between
!!   grid blocks of different refinement level.
!! 
!!   These routines get block boundary data from neighbors who are
!!   parents of leaf blocks. This is required in flux conserving schemes
!!   where the coarser block needs to use the same fluxes and mean pressures
!!   as will be used on the finer blocks across their shared boundary.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson for 
!!   directional guardcell filling.
!!   Klaus Weide (2021)    modified for pdg stuff
!!***

#include "paramesh_preprocessor.fh"

Subroutine amr_flux_conserve(mype,nsub,flux_dir,pdgNo)

!-----Use statements.
  Use paramesh_dimensions, only: ndim
  Use physicaldata, only: gr_thePdgs
  Use physicaldata, only: advance_all_levels, var_dt
  Use tree, only: lnblocks, nodetype
  Use paramesh_interfaces, only : amr_flux_conserve_udt,           & 
                                      amr_flux_conserve_vdt

  Implicit None

!-----Input/Output arguments.
  Integer, intent(in)  ::  mype,nsub
  Integer, optional, intent(in) :: flux_dir
  Integer, optional, intent(in) :: pdgNo

!-----Local variables
  Integer :: lb
  integer :: npdg, ig,sg,eg

!-----Begin executable code.
  npdg = 1
  if (present(pdgNo)) npdg = pdgNo
  if (npdg == -1) then
     sg = 1; eg = NUM_PDGS
  else
     sg = npdg; eg = npdg
  end if

  do ig = sg,eg

     If (lnblocks > 0) Then
        Do lb = 1,lnblocks

           If (nodetype(lb) == 1 .or. advance_all_levels) Then

!------Store fluxes in temporary storage
              if (allocated(gr_thePdgs(ig)%tflux_x)) gr_thePdgs(ig)%tflux_x(:,:,:,:,lb) = &
                                                      gr_thePdgs(ig)%flux_x(:,:,:,:,lb)
              if (ndim >= 2 .AND. allocated(gr_thePdgs(ig)%tflux_y)) then
                 gr_thePdgs(ig)%tflux_y(:,:,:,:,lb) = &
                  gr_thePdgs(ig)%flux_y(:,:,:,:,lb)
              end if
              if (ndim == 3 .AND. allocated(gr_thePdgs(ig)%tflux_z)) then
                 gr_thePdgs(ig)%tflux_z(:,:,:,:,lb) = &
                  gr_thePdgs(ig)%flux_z(:,:,:,:,lb)
              end if

           End If

        End Do
     End If

     If (var_dt) Then
        Call amr_flux_conserve_vdt(mype,nsub,gr_thePdgs(ig),ig) ! Called if variable dt
     Else
        Call amr_flux_conserve_udt(mype,gr_thePdgs(ig),ig,flux_dir) ! Called if uniform dt
     End If

  end do

  Return
End Subroutine amr_flux_conserve
