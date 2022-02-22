!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_plotfile_chombo
!! NAME
!!
!!   amr_plotfile_chombo
!!
!! SYNOPSIS
!!
!!   call amr_plotfile_chombo (file_num)
!!   call amr_checkpoint_wr_hdf5(integer)
!!
!! ARGUMENTS
!!
!!   integer, intent(in) :: file_num
!!     An integer number which will be appended to the end of the file name.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!! CALLS
!!
!! RETURNS
!!
!!   Does not return anything. 
!!
!! DESCRIPTION
!! 
!!  This is a dummy, placeholer routine for chombovis plotfile output.  
!!  This routine will return an error message if the chombovis plotfile 
!!  capability has not been installed.  The routine which actually writes 
!!  a chombovis file can be found in utilities/io/plotting/chombo.
!!
!! AUTHORS
!!
!!   Kevin Olson (2004)
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_plotfile_chombo (file_num)
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output arguments.
      integer, intent(in) :: file_num

!-----Local arrays and variables.
      integer :: mype, ierr

      Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

      If (mype == 0) then
         Print *,' WARNING: you are calling amr_plotfile_chombo '
         Print *,'          but your version of paramesh is not '
         Print *,'          yet configured to do this.          '
         Print *,'          Go to utilities/io/plotting/chombovis '
         Print *,'          in the main paramesh directory, run  '
         Print *,'          the INSTALL script, and recompile !!!'
      End If

      Return
      End Subroutine amr_plotfile_chombo

