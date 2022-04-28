!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_guardcell
!! NAME
!!
!!   amr_guardcell
!!
!! SYNOPSIS
!!
!!   call gr_amrGuardcellTail(mype, iopt)
!!
!!   call gr_amrGuardcellTail(integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: iopt  
!!     Selects whether to fill the guardcells for the arrays unk, 
!!     facevarx, facevary, facevarz, unk_e_x, unk_e_y, unk_e_z, and unk_n 
!!     (if iopt = 1) or work (if iopt 2).
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!! 
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!! 
!!   amr_1blk_guardcell_reset
!!   amr_restrict
!!   amr_1blk_guardcell
!!   mpi_amr_comm_setup
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells are filled with data for 
!!   all blocks.
!!
!! DESCRIPTION
!!
!!
!! SEE ALSO
!!   gr_amrGuardcellHead
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
#include "constants.h"

Module gr_amrGuardcellTail_mod

contains
  Subroutine gr_amrGuardcellTail(mype,iopt)

!-----Use Statements
   use Timers_interface, ONLY : Timers_start, Timers_stop
   Use paramesh_dimensions
   Use physicaldata
   Use workspace
   Use tree
   use paramesh_comm_data

   Use paramesh_interfaces, Only : amr_1blk_guardcell_reset

   Implicit none

!-----Input/Output Arguments
   Integer, intent(in) :: mype,iopt

!-----Local variables
   Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
   Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag,l_force_consist
   Integer :: lb,icoord
   Integer :: id,jd,kd
   Integer :: ilays,jlays,klays
   Integer :: i,j,k,ivar                                  
   Integer :: ip1,ip2,jp1,jp2,kp1,kp2
   Integer :: ilp,iup,jlp,jup,klp,kup
   Integer :: nprocs, ierr, tag_offset, iempty, iu, ju, ku, iopt0

!------------------------------------
!-----Begin Executable code section
!------------------------------------


   If (no_permanent_guardcells) Then

       If (mype == 0) Then
         Write(*,*) 'gr_amrGuardcellTail call ignored!'
         Write(*,*) 'NO_PERMANENT_GUARDCELLS is defined'
       End if
       Return

   Else  ! no_permanent_guardcells

!-----reinitialize addresses of cached parent blocks
      Call amr_1blk_guardcell_reset

!-----reset selections of guardcell variables to default
      int_gcell_on_cc(:) = .True.
      int_gcell_on_fc(:,:) = .True.
      int_gcell_on_ec(:,:) = .True.
      int_gcell_on_nc(:) = .True.

   Endif ! If (no_permanent_guardcells)

   Return
 End Subroutine gr_amrGuardcellTail
end Module gr_amrGuardcellTail_mod
