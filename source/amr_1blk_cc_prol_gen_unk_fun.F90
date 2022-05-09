!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_1blk_cc_prol_gen_unk_fun
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
!!   amr_1blk_cc_prol_gen_unk_fun
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_gen_unk_fun (recv, ia, ib, ja, jb, ka, kb,
!!                                      ioff, joff, koff, mype,
!!                                      lb, pe_p, lb_p)
!!   Call amr_1blk_cc_prol_gen_unk_fun (real, 
!!                                      integer, integer, integer, integer,
!!                                      integer, integer, integer, integer,
!!                                      integer, integer)
!!
!! DESCRIPTION
!!
!!   This routine is a wrapper routine which calls the functions
!!   which prolong data for UNK. The logical array interp_mask_unk can
!!   be used to control which routine is actually operating on
!!   each variable stored within UNK.
!!
!! ARGUMENTS
!!
!!   recv(:,:,:,:)   - Real,    Intent(inout)
!!     Data to prolong.
!!
!!   ia,ib,ja,jb,ka,kb,idest   - Integer, Intent(in)
!!     Indices in unk1 array to place prolonged data
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
!! INCLUDES
!!   
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   prolong_arrays
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_prolong_gen_unk1_fun
!!   amr_1blk_cc_prol_inject
!!   amr_1blk_cc_prol_linear
!!   amr_1blk_cc_prol_genorder
!!   amr_1blk_cc_prol_dg
!!   amr_1blk_cc_prol_user
!!
!! RETURNS
!!
!!   Upon return prolonged data is placed in unk1.
!!
!! AUTHORS
!!
!!  Written by Peter MacNeice January 2002.
!!  Modified for GRID_WITH_MONOTONIC variant - Klaus Weide 2022-02-20
!!  Changes to call amr_1blk_cc_prol_dg for Thornado - Austin Harris 2021-12-06
!!***

#include "paramesh_preprocessor.fh"


subroutine amr_1blk_cc_prol_gen_unk_fun                &
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  & 
         mype,lb,pe_p,lb_p)

!-----Use Statements
  use timings, ONLY: timing_mpi, timer_amr_1blk_cc_prol_gen_unk
#ifndef GRID_WITH_MONOTONIC
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use prolong_arrays

  Use paramesh_interfaces, only :                  &
                       amr_1blk_cc_prol_inject,    & 
                       amr_1blk_cc_prol_linear,    & 
                       amr_1blk_cc_prol_genorder,  & 
                       amr_1blk_cc_prol_dg,    &
                       amr_1blk_cc_prol_user
#endif

  implicit none

!-----Include Statements
  Include 'mpif.h'

!-----Input/Output variables
  real,    intent(inout) :: recv(:,:,:,:)
  integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
  integer, intent(in)    :: ioff,joff,koff,mype
  integer, intent(in)    :: lb,lb_p,pe_p

!-----Local variables
  double precision :: time1
  integer :: ivar

!-----Begin Executable code

  if (timing_mpi) then
     time1 = mpi_wtime()
  end if  ! End If (timing_mpi)

#ifdef GRID_WITH_MONOTONIC
! Call the minimally changed subroutine from Paramesh2
  call amr_prolong_gen_unk1_fun &
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &     mype,lb)
#else
  Do ivar = 1, nvar
     If (int_gcell_on_cc(ivar)) Then

        If (interp_mask_unk(ivar) < 20) Then

           If (interp_mask_unk(ivar) == 0) then
!-----Simple Injection 
              Call amr_1blk_cc_prol_inject               &
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
            mype,ivar)
      
           Elseif (interp_mask_unk(ivar) == 1) Then
!-----Default multi-linear interpolation  
              Call amr_1blk_cc_prol_linear               &
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
           mype,ivar)

           Elseif (interp_mask_unk(ivar) > 1) Then
!-----High order Largrange ploynomial interpolation
              Call amr_1blk_cc_prol_genorder             &
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
            mype,ivar,interp_mask_unk(ivar))

           End If  ! End If (interp_mask_unk(ivar) == 0)

        Elseif (interp_mask_unk(ivar) == 20) Then

!--------User defined interpolation to be used for prolocation

           Call amr_1blk_cc_prol_user()

        Elseif (interp_mask_unk(ivar) == 40) Then

!--------User defined interpolation to be used for
!prolongation/restriction from Thornado

           Call amr_1blk_cc_prol_dg                      &
           (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
           mype,ivar)

        End If  ! End If (interp_mask_unk(ivar) < 20

     End If  ! Enf If (int_gcell_on_cc(ivar))
  End Do  ! End Do ivar = 1, nvar
#endif

  if (timing_mpi) then
     timer_amr_1blk_cc_prol_gen_unk =                 &
                   timer_amr_1blk_cc_prol_gen_unk     & 
                                + mpi_wtime() - time1
  end if  ! End If (timing_mpi)

  return
end subroutine amr_1blk_cc_prol_gen_unk_fun
