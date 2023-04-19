!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_restrict_unk_dg
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright 2023 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!   amr_restrict_unk_dg
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_dg (datain, dataout,        ioff, joff, koff,         pdg,        ig)
!!   Call amr_restrict_unk_dg (real array, real array, integer, integer, integer,TYPE(pdg_t),integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   ioff, joff, koff : offsets of the restricted data that is returned within the coarse
!!                      block; each of these numbers should be either 0 or half
!!                      the number of interior cells in the block for the
!!                      relevant direction.
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   DG quadrature-point data stored in 'unk'.
!!
!!   Data is passed in in the array 'datain' and returned in the array
!!   'dataout'.
!!   The restriction operates on a subset of variables in UNK selected
!!   by corresponding entries in the array interp_mask_unk_res.
!!
!!  This is a stub that needs to be overridden to be useful.
!!
!!  This stub version just aborts with a message.
!!
!! AUTHORS
!!
!!  Stub version created  -  Klaus Weide 2022-04-28
!!  Stub version updated  -  Klaus Weide 2023-03-16
!!  2023-03-20 Added pdg,ig arguments               - Klaus Weide
!!  2023-04-07 Eliminated ivar arg                 -  Klaus Weide
!!
!!***

Subroutine amr_restrict_unk_dg(datain,dataout,ioff,joff,koff,pdg,ig)

  !-----Use Statements
  Use Driver_interface, ONLY: Driver_abort
  use gr_pmPdgDecl, ONLY : pdg_t

  Implicit None

!-----Input/Output arguments.
  Real,    Intent(in)    :: datain(:,:,:,:)
  Real,    Intent(inout) :: dataout(:,:,:,:)
  Integer, Intent(in)    :: ioff,joff,koff
  type(pdg_t), intent(INOUT) :: pdg
  Integer, Intent(in)    :: ig

  call Driver_abort(&
       'An implementation of amr_restrict_unk_dg needs to be provided!')

End Subroutine amr_restrict_unk_dg
