!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!------------------
! Guardcell filling

Subroutine mpi_amr_Write_guard_comm(nprocs)
  implicit none
  Integer, intent(in) :: nprocs
  return
end Subroutine mpi_amr_Write_guard_comm


Subroutine mpi_amr_read_guard_comm(nprocs)
  use gr_pmCommDataTypes, ONLY: GRID_PAT_GC, GRID_PAT_FCORR, &
                                    GRID_PAT_PROLONG, GRID_PAT_RESTRICT
  use gr_pmCommPatternData, ONLY: gr_pmActivateCommPattern, &
                                  gr_theActiveCommPattern
  !-----Use Statements
  Use physicaldata, ONLY : mpi_pattern_id
  Use tree, ONLY: grid_analysed_mpi, strt_buffer
  Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
  implicit none
  Integer, intent(in) :: nprocs
  If (grid_analysed_mpi.ne.1) Then
     Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
     Call amr_morton_process()
  End If
  call gr_pmActivateCommPattern(GRID_PAT_GC)
  strt_buffer = gr_theActiveCommPattern % strt_buffer
  mpi_pattern_id = 10

  Call amr_1blk_guardcell_reset
end Subroutine mpi_amr_read_guard_comm

Subroutine mpi_amr_write_prol_comm(nprocs)
!-----Use Statements
  Implicit None
!-----Input/Output Arguments
  Integer, intent(in) :: nprocs
End Subroutine mpi_amr_write_prol_comm





Subroutine mpi_amr_read_prol_comm(nprocs)

!-----Use Statements
      Use physicaldata, ONLY : mpi_pattern_id
      Use tree, ONLY: grid_analysed_mpi, strt_buffer
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
  use gr_pmCommDataTypes, ONLY: GRID_PAT_GC, GRID_PAT_FCORR, &
                                    GRID_PAT_PROLONG, GRID_PAT_RESTRICT
  use gr_pmCommPatternData, ONLY: gr_pmActivateCommPattern, &
                                  gr_theActiveCommPattern

      Implicit None

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs


!-----Begin Executable code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If
      call gr_pmActivateCommPattern(GRID_PAT_PROLONG)
      strt_buffer = gr_theActiveCommPattern % strt_buffer

      mpi_pattern_id = 20

      Return
      End Subroutine mpi_amr_read_prol_comm
 

 
!------------------
! Flux Conservation
  
Subroutine mpi_amr_write_flux_comm(nprocs)
  implicit none
  Integer, intent(in) :: nprocs
end Subroutine mpi_amr_write_flux_comm

Subroutine mpi_amr_read_flux_comm(nprocs)

!-----Use Statements
      Use physicaldata, ONLY : mpi_pattern_id
      Use tree, ONLY: grid_analysed_mpi, strt_buffer
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
  use gr_pmCommDataTypes, ONLY: GRID_PAT_GC, GRID_PAT_FCORR, &
                                    GRID_PAT_PROLONG, GRID_PAT_RESTRICT
  use gr_pmCommPatternData, ONLY: gr_pmActivateCommPattern, &
                                  gr_theActiveCommPattern

      Implicit None

!-----Input/Output Arguments
      Integer, intent(in) :: nprocs


!-----Begin Executable code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      call gr_pmActivateCommPattern(GRID_PAT_FCORR)
      strt_buffer = gr_theActiveCommPattern % strt_buffer
      mpi_pattern_id = 30

      Return
      End Subroutine mpi_amr_read_flux_comm
 


!------------------
! Restriction
Subroutine mpi_amr_write_restrict_comm(nprocs)
  implicit none
  Integer, intent(in) :: nprocs
end Subroutine mpi_amr_write_restrict_comm


Subroutine mpi_amr_read_restrict_comm(nprocs)

!-----Use Statements
  Use physicaldata, ONLY : mpi_pattern_id
  Use tree, ONLY: grid_analysed_mpi, strt_buffer
  Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
  use gr_pmCommDataTypes, ONLY: GRID_PAT_GC, GRID_PAT_FCORR, &
                                    GRID_PAT_PROLONG, GRID_PAT_RESTRICT
  use gr_pmCommPatternData, ONLY: gr_pmActivateCommPattern, &
                                  gr_theActiveCommPattern

  Implicit None

!-----Input/Output Statements
      Integer, intent(in) :: nprocs

!-----Begin Executable Code
      If (grid_analysed_mpi.ne.1) Then
        Write(*,*) 'PARAMESH ERROR: communication control info is ',   & 
             'being read, but it was never set up. You are ',          & 
             'probably missing a Call to amr_checkpoint_re or ',       & 
             'amr_refine_derefine.',                                   & 
             'These Call amr_morton_process, which is the Call',       & 
             ' that is actually missing. We will Call this for you',   & 
             ' now in the hope that this corrects your problem.',      & 
             ' Please review this before relying on results.'
        Call amr_morton_process()
      End If

      call gr_pmActivateCommPattern(GRID_PAT_RESTRICT)
      strt_buffer = gr_theActiveCommPattern % strt_buffer

      If (gr_theActiveCommPattern % num_recipient_pes > 0) Then

         Call amr_1blk_guardcell_reset

      End If

      mpi_pattern_id = 40

      Return
      End Subroutine mpi_amr_read_restrict_comm
