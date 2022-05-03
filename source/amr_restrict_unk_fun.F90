!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_restrict_unk_fun
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
!!   amr_restrict_unk_fun
!!
!! SYNOPSIS
!!
!!   Call amr_restrict_unk_fun(datain, dataout)
!!   Call amr_restrict_unk_fun(real array, real array)
!!
!! ARGUMENTS
!!
!!   Real, Intent(in)    :: datain(:,:,:,:)  data to restrict
!!   Real, Intent(inout) :: dataout(:,:,:,:) restricted data to return
!!   
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_restrict_unk_genorder 
!!   amr_restrict_unk_user
!!   amr_restrict_unk_dg
!!
!! RETURNS
!!
!!   Restricted data returned in array 'dataout'.
!!
!! DESCRIPTION
!!   
!!   This routine performs restriction on the array datain and
!!   returns the result in dataout. Note that this does not update
!!   guard cell elements of dataout.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          January 1997
!!   Modified by Kevin Olson for high order restriction, 2004.
!!   Call amr_restrict_unk_dg for Thornado - Austin Harris, K. Weide 2022-04-28
!!***

Subroutine amr_restrict_unk_fun(datain,dataout)

!-----Use statements.
  Use paramesh_dimensions
  Use physicaldata
  Use paramesh_interfaces, only : amr_restrict_unk_genorder,       &
                                      amr_restrict_unk_user,           &
                                      amr_restrict_unk_dg

  Implicit None

!-----Input/Output arguments.
  Real, Intent(in)    :: datain(:,:,:,:)
  Real, Intent(inout) :: dataout(:,:,:,:)

!-----Local variables.
  Integer :: ivar, order

!-----Begin Executable code.

  Do ivar = 1, nvar

     If (int_gcell_on_cc(ivar)) Then

         If (interp_mask_unk_res(ivar) < 20) Then

!-----------Call the default interpolation routine for interpolation 
            order = interp_mask_unk_res(ivar)
            If (order <=0 .or. order > 5) order = 1
            Call amr_restrict_unk_genorder(datain,dataout,order,ivar)

         Elseif (interp_mask_unk_res(ivar) == 40) Then
!--------User defined interpolation to be used for
!prolongation/restriction from Thornado

            Call amr_restrict_unk_dg(datain,dataout,ivar)

         ElseIf (interp_mask_unk_res(ivar) >= 20) Then

!-----------Call a user defined routine for restriction
            Call amr_restrict_unk_user()

         End If

     End If

  End Do  ! End Do ivar = 1, nvar

End Subroutine amr_restrict_unk_fun




