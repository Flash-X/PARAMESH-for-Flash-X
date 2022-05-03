!!****if* source/Grid/GridMain/AMR/Paramesh4/PM4_package/source/amr_1blk_cc_prol_dg
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
!!   amr_1blk_cc_prol_dg
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_cc_prol_dg (recv,ia,ib,ja,jb,ka,kb,idest,
!!                                 ioff,joff,koff,mype,ivar)
!!   Call amr_1blk_cc_prol_dg (real,
!!                                 integer, integer, integer, integer,
!!                                 integer, integer, integer, integer,
!!                                 integer, integer, integer, integer)
!!
!! DESCRIPTION
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

Subroutine amr_1blk_cc_prol_dg               &
        (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,  &
         mype,ivar)

  !-----Use Statements
  Use Driver_interface, ONLY: Driver_abort

  Implicit None

!-----Input/Output Variables
  Real,    Intent(inout) :: recv(:,:,:,:)
  Integer, Intent(in)    :: ia,ib,ja,jb,ka,kb
  Integer, Intent(in)    :: idest,ioff,joff,koff,mype
  Integer, Intent(in)    :: ivar

  call Driver_abort(&
       'An implementation of amr_1blk_cc_prol_dg needs to be provided!')

End Subroutine amr_1blk_cc_prol_dg
