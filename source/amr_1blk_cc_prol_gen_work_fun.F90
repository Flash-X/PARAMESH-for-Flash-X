!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_1blk_cc_prol_gen_work_fun
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2016 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!   amr_1blk_cc_prol_gen_work_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_gen_work_fun (recvt, ia, ib, ja, jb, ka, kb,
!!                                      ioff, joff, koff, mype,
!!                                      lb, pe_p, lb_p, interp)
!!   Call amr_1blk_cc_prol_gen_work_fun (real, 
!!                                      integer, integer, integer, integer,
!!                                      integer, integer, integer, integer,
!!                                      integer, integer. integer)
!!  call amr_1blk_cc_prol_gen_work_fun(real(inout) :: recvt,
!!                                     integer(in) :: ia,
!!                                     integer(in) :: ib,
!!                                     integer(in) :: ja,
!!                                     integer(in) :: jb,
!!                                     integer(in) :: ka,
!!                                     integer(in) :: kb,
!!                                     integer(in) :: idest,
!!                                     integer(in) :: ioff,
!!                                     integer(in) :: joff,
!!                                     integer(in) :: koff,
!!                                     integer(in) :: mype,
!!                                     integer(in) :: lb,
!!                                     integer(in) :: pe_p,
!!                                     integer(in) :: lb_p,
!!                                     integer(in) :: interp)
!!
!! DESCRIPTION
!!
!!  This routine is a wrapper routine which calls the functions
!!  which prolong data for WORK. The argument idest selects the layer
!!  within WORK on which prolongation is required.
!!  The argument interp should select the actual interpolation function to be
!!  called, but currently only the values 30..32 are specially recognized.
!!
!! ARGUMENTS
!!
!!   recvt(:,:,:,:)   - Real,    Intent(inout)
!!     Data to prolong.
!!
!!   ia,ib,ja,jb,ka,kb,idest   - Integer, Intent(in)
!!     Indices in work1 array to place prolonged data
!!
!!   ioff,joff,koff   - Integer, Intent(in)
!!     Offsets
!!
!!   mype : my Processor Number
!!
!!   lb             - Integer, Intent(in)
!!     ID of local block to prolong
!!
!!   lb_p,pe_p      - Integer, Intent(in)
!!     ID and pe of parent block  (unused in this inplementation)
!!
!!   interp   - Integer, Intent(in)
!!     order of polynomial to use
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
!!   workspace
!!   prolong_arrays
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_1blk_cc_prol_work_mg
!!   amr_prolong_gen_work1_fun
!!   amr_1blk_cc_prol_work_inject
!!   amr_1blk_cc_prol_work_linear
!!   amr_1blk_cc_prol_work_user
!!   amr_1blk_cc_prol_work_genorder
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in work1.
!!
!! AUTHORS
!!
!! Written by Peter MacNeice January 2002.
!!
!! HISTORY
!!
!!  Written by Peter MacNeice January 2002.
!!  Modified to call amr_1blk_cc_prol_work_mg - Marcos Vanella 2010
!!  Modified to dispatch interp=30..32        - Klaus Weide Sep 2016
!!  Modified for GRID_WITH_MONOTONIC variant - Klaus Weide 2022-02-20
!!***

#include "paramesh_preprocessor.fh"


Subroutine amr_1blk_cc_prol_gen_work_fun(recvt,             &
             ia,ib,ja,jb,ka,kb,                             & 
             idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)


!-----Use Statements
#ifdef GRID_WITH_MONOTONIC
  Use paramesh_interfaces, Only :                  &
                        amr_1blk_cc_prol_work_mg
#else
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use workspace
  Use prolong_arrays

  Use paramesh_interfaces, Only :                  &
                        amr_1blk_cc_prol_work_inject,  & 
                        amr_1blk_cc_prol_work_linear,  & 
                        amr_1blk_cc_prol_work_user,    & 
                        amr_1blk_cc_prol_work_genorder,&
                        amr_1blk_cc_prol_work_mg
#endif


  implicit none


!-----Input/Output variables
  real,    intent(inout) :: recvt(:,:,:)
  integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
  integer, intent(in)    :: ioff,joff,koff,mype
  integer, intent(in)    :: lb,lb_p,pe_p
  integer, intent(in)    :: interp

!-----Begin Executable Code

  if (interp .GE. 30 .AND. interp .LE. 32) then
! Call special version developed for MC multigrid
     call amr_1blk_cc_prol_work_mg(recvt,       &
             ia,ib,ja,jb,ka,kb,                     &
             idest,ioff,joff,koff,lb,interp-30)

  else

#ifdef GRID_WITH_MONOTONIC
! Call the minimally changed subroutine from Paramesh2
     call amr_prolong_gen_work1_fun &
     &     (recvt,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &     mype,lb)
#else
! Apply the interpolation " native" to PARAMESH
     if (interp < 20) Then

!-----Simple Injection
        If (interp.eq.0)                            &
         Call amr_1blk_cc_prol_work_inject(recvt,   &
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype)

!-----Linear interpolation
        If (interp.eq.1)                            &
         Call amr_1blk_cc_prol_work_linear(recvt,   & 
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype)

!-----High order Lagrange polynomial interpolation
        If (interp > 1)                             &
         Call amr_1blk_cc_prol_work_genorder(recvt, & 
             ia,ib,ja,jb,ka,kb,                     & 
             idest,ioff,joff,koff,mype,interp)

     else if (interp == 20) Then

        Call amr_1blk_cc_prol_work_user()

     end If  ! End If (interp < 20)
#endif
  end if !  if (interp .GE. 30 .AND. interp .LE. 32)

  Return
End Subroutine amr_1blk_cc_prol_gen_work_fun
