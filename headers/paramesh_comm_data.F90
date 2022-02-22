!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! This module is meant to hold communication data items needed by
! various PARAMESH routines.

! Written: K.Olson January 2007

      Module paramesh_comm_data

      Integer :: amr_mpi_real  ! This a variable which holds either the 
                               ! MPI_REAL or MPI_DOUBLE PRECISION
                               ! type identifiers.  It was added to support
                               ! the g95 compiler.
      integer :: amr_mpi_meshComm

      !The following derived datatype is used in amr_set_runtime_parameters
      !to communicate the runtime parameters from the master to all
      !other processes.
      Integer, parameter :: NumTypes = 3
      Integer, parameter :: NumInt = 19
      Integer, parameter :: NumLog = 23
      Integer, parameter :: NumChar = 80
      !NumChar must be consistent with PM4_package/headers/io.F90
      !which has Character (Len=80) :: output_dir

      type ParameshParms
         integer, dimension(NumInt) :: IntParms
         integer, dimension(NumLog) :: LogParms ! Store 0 for FALSE, 1 for TRUE
         character, dimension(NumChar) :: CharParms
      end type ParameshParms

      End Module paramesh_comm_data
