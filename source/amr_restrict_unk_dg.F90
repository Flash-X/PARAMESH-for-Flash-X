!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_restrict_unk_dg
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
!!   amr_restrict_unk_dg
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_dg (datain, dataout, ivar)
!!   Call amr_restrict_unk_dg (real array, real array, integer, integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Real,    Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real,    Intent(inout) :: dataout(:,:,:,:)  data which is restricted and returned
!!   Integer, Intent(in)    :: ivar  variable number in unk to restrict
!!
!! DESCRIPTION
!!
!!   This routine performs interpolation for the restriction operation on
!!   DG quadrature-point data stored in 'unk'.
!!
!!   Data is passed in in the array 'datain' and returned in the array
!!   'dataout'.
!!   The last argument 'ivar' specifies which variable in 'unk' to apply
!!   the interpolation to.
!!
!!  This is a stub that needs to be overridden to be useful.
!!
!!  This stub version just aborts with a message.
!!
!! AUTHORS
!!
!!  Stub version created  -  Klaus Weide 2022-04-28
!!
!!***

Subroutine amr_restrict_unk_dg(datain,dataout,ivar)

  !-----Use Statements
  Use Driver_interface, ONLY: Driver_abort

  Implicit None

!-----Input/Output arguments.
  Real,    Intent(in)    :: datain(:,:,:,:)
  Real,    Intent(inout) :: dataout(:,:,:,:)
  Integer, Intent(in)    :: ivar

  call Driver_abort(&
       'An implementation of amr_restrict_unk_dg needs to be provided!')

End Subroutine amr_restrict_unk_dg
