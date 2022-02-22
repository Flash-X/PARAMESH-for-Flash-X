!----------------------------------------------------------------------
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2021 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
Module gr_pmDivpres_mod

  public :: prol_fc_divpres
  logical, PARAMETER :: prol_fc_divpres = .false.
  ! No DIVPRES variables

Contains

  subroutine prol_fc_divpres_init(n,tf,i_divf_fc_vars)
    use Driver_interface, ONLY: Driver_abort
    implicit none
      
    integer, intent(in) :: n, tf, i_divf_fc_vars(tf,n)

    call Driver_abort(&
           "Attempting to initialize divpres module, not comfigured in!")

  end subroutine prol_fc_divpres_init


  subroutine gr_pmDivpresApply(recvfx,recvfy,recvfz, &
       ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
       lb)
    use paramesh_dimensions, ONLY: il_bnd1,jl_bnd1,kl_bnd1

    implicit none
    Real,intent(inout) :: recvfx(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfy(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfz(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    integer, intent(in) :: ia,ib,ja,jb,ka,kb
    integer, intent(in) :: idest,ioff,joff,koff
    integer, intent(in) :: lb


  end subroutine gr_pmDivpresApply

end Module gr_pmDivpres_mod
