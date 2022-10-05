!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_prolong_fun_init
!! NAME
!!
!!   amr_prolong_fun_init
!!
!! SYNOPSIS
!!
!!   Call amr_prolong_fun_init()
!!
!! ARGUMENTS
!!
!!   No arguments.
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_prolong_cc_fun_init
!!   amr_prolong_face_fun_init
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine calls the routines which compute the values of dx,dy and 
!!   dz and some index vectors used during the interpolation process. 
!!   These are used inside the prolongation routines
!!   saving needless repetitive computation at the cost of minimal storage
!!   space.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          July 1997
!!   Modified for pdg:  Klaus Weide        Sept 2022
!!
!!***

#include "Simulation.h"
      Subroutine amr_prolong_fun_init

      use gr_pmPdgDecl, ONLY : pdg_t
      use paramesh_interfaces, only : amr_prolong_cc_fun_init,         & 
                                      amr_prolong_face_fun_init
      use physicaldata, only: gr_thePdgs
      Implicit None

      integer :: ig

!-----Begin executable code.

      do ig = 1, NUM_PDGS
         Call amr_prolong_cc_fun_init(gr_thePdgs(ig),ig)
      end do
      Call amr_prolong_face_fun_init

      Return
    End Subroutine amr_prolong_fun_init
