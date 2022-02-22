!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
! Copyright (C) 2017 The University of Chicago
! Copyright (C) 2022 UChicago Argonne, LLC and contributors
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_set_runtime_parameters
!! NAME
!!
!!   amr_set_runtime_parameters
!!
!! SYNOPSIS
!!
!!   Call amr_set_runtime_parameters()
!!
!! ARGUMENTS
!!
!!   None
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
!!   io
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing returned
!!
!! DESCRIPTION
!!
!!   This routine reads in the runtime parameters which are used by
!!   PARAMESH to set array sizes.  The runtime parameters are stored
!!   in the file 'amr_runtime_parameters' which you must create.
!!   A copy of the file 'amr_runtime_parameters' must be available for 
!!   each processor to open and then to read from.
!!
!! AUTHORS
!!
!!   Kevin Olson
!!
!!   Modified by Chris Daley so that there is 1 file read instead of
!!   NPROC file reads.
!!
!!***

#include "paramesh_preprocessor.fh"
! Undefine these symbols if defined by Simulation.h, we want the ones from constants.h! - KW
#ifdef CARTESIAN
#undef CARTESIAN
#endif
#ifdef SPHERICAL
#undef SPHERICAL
#endif
#ifdef CYLINDRICAL
#undef CYLINDRICAL
#endif
#ifdef POLAR
#undef POLAR
#endif
#include "constants.h"

        Subroutine amr_set_runtime_parameters()

!-------Use statements.
        Use paramesh_dimensions
        Use physicaldata
        Use tree,    ONLY: mflags, nboundaries
        Use timings, ONLY: timing_mpi, timing_mpix
        Use io,      ONLY: output_dir, amr_log_file
        Use Paramesh_comm_data
        use Grid_data, ONLY: gr_geometry, gr_globalMe
        use Logfile_interface, ONLY: Logfile_stamp
#ifndef LIBRARY
        use Driver_interface, ONLY : Driver_abort
#endif
#ifndef USE_AMR_RUNTIME_PARAMETERS_FILE
        use RuntimeParameters_interface, ONLY: RuntimeParameters_get
#endif
        Implicit None

!-------Include statments
        Include 'mpif.h'

!-------Local variables
        logical :: fileExists, logicalValue
        type(ParameshParms) :: pmParms
        Integer :: mype, i, ierr

        !Unfortunately MPI_VERSION is an integer in mpif.h so we can't use it
        !to conditionally compile the appropriate MPI functions.  We provide a
        !custom macro named FLASH_MPI2 which the user can define to make use
        !of the better MPI-2 version.
#ifdef FLASH_MPI2
        Integer(MPI_ADDRESS_KIND), dimension(NumTypes) :: disp
        Integer(MPI_ADDRESS_KIND) :: base
        Integer, dimension(NumTypes) :: blocklen, type
        Integer :: newtype
#else
        Integer, dimension(NumInt) :: IntParms
        Integer, dimension(NumLog) :: LogParms
        Character, dimension(NumChar) :: CharParms
#endif

!-------Begin executable code.

        pmParms % LogParms(:) = 0
        do i = 1, NumChar
           !Insert a null in each element.  We do not want to send
           !uninitialized data.
           pmParms % CharParms (i) = char(0)
        end do

#ifdef USE_AMR_RUNTIME_PARAMETERS_FILE
        Call MPI_COMM_RANK(amr_mpi_meshComm,mype,ierr)

        !The master processor reads the amr_runtime_parameters file and
        !stores the data into a single derived datatype.
        If (mype == 0) Then
           Open (unit=35,                                                 & 
                file='amr_runtime_parameters',                           & 
                status='old',                                            & 
                action='READ',                                           & 
                form='formatted')

           !-------Integers
           do i = 1, NumInt
              Read (35,*) pmParms % IntParms (i)
           end do

           !-------Logicals
           !------ Read as logical, store 1 or 0 as integer
           do i = 1, NumLog
              Read (35,*) logicalValue
              if (logicalValue) then
                 pmParms % LogParms(i) =  1
              else
                 pmParms % LogParms(i) =  0
              end if
           end do

           !-------characters
           !A Fortran string type is different to an array of characters.
           !MPI understands an array of characters only, and so we copy
           !individual string characters into increasing elements of
           !our character array.

           Read (35,*) output_dir
           do i = 1, len_trim(output_dir)
              pmParms % CharParms (i) = output_dir(i:i)
           end do

           Close(35)
        End if

        !The master processor broadcasts the derived datatype to everyone and
        !then everyone (including the master) initializes Paramesh variables
        !from fields of the derived datatype.
#ifdef FLASH_MPI2
        !MPI_Get_address and MPI_Type_create_struct are only available in MPI-2
        !implementations.  The functions MPI_Address and MPI_Type_struct are
        !available in MPI-1 implementations, but have known issues.
        type(1) = MPI_INTEGER
        type(2) = MPI_INTEGER
        type(3) = MPI_CHARACTER

        blocklen(1) = NumInt
        blocklen(2) = NumLog
        blocklen(3) = NumChar

        call MPI_Get_address(pmParms % IntParms, disp(1), ierr)
        call MPI_Get_address(pmParms % LogParms, disp(2), ierr)
        call MPI_Get_address(pmParms % CharParms, disp(3), ierr)

        base = disp(1)
        disp(1) = disp(1) - base
        disp(2) = disp(2) - base
        disp(3) = disp(3) - base

        call MPI_Type_create_struct(NumTypes, blocklen, disp, type, newtype, ierr)
        call MPI_Type_commit(newtype, ierr)
        call MPI_Bcast(pmParms, 1, newtype, 0, amr_mpi_meshComm, ierr)
        call MPI_Type_free(newtype, ierr)
#else
        !We cannot use an elegant solution like above and maintain MPI-1
        !compatibility because we must replace MPI_Get_address with
        !MPI_Address and this has has known issues.  When we use MPI_Address
        !we get the following error on 64-bit platforms:
        !"Fatal error in MPI_Address: Invalid argument, error stack:
        !MPI_Address(210): An address does not fit into a Fortran INTEGER.
        !Use MPI_Get_address instead"
        !Our (ugly) solution is to broadcast 3 primitive type messages.
        If (mype == 0) Then
           IntParms(:) = pmParms % IntParms(:)
           LogParms(:) = pmParms % LogParms(:)
           CharParms(:) = pmParms % CharParms(:)
        End If

        call MPI_Bcast(IntParms, NumInt, MPI_INTEGER, 0, amr_mpi_meshComm, ierr)
        call MPI_Bcast(LogParms, NumLog, MPI_INTEGER, 0, amr_mpi_meshComm, ierr)
        call MPI_Bcast(CharParms, NumChar, MPI_CHARACTER, 0, amr_mpi_meshComm, ierr)

        pmParms % IntParms(:) = IntParms(:)
        pmParms % LogParms(:) = LogParms(:)
        pmParms % CharParms(:) = CharParms(:)
#endif

#else
        !The master processor reads the amr_runtime_parameters file and
        !stores the data into a single derived datatype.
        if (gr_globalMe == MASTER_PE) Then
           inquire (file = 'amr_runtime_parameters', exist = fileExists)
           if (fileExists) then
              print*,'NOTE: File amr_runtime_parameters exists, but will be ignored.'
              print*,'FLASH gr_pmrp* runtime parameters will be used to initialize PARAMESH runtime parameters.'
              call Logfile_stamp('File amr_runtime_parameters exists, but will be ignored.',&
                   '[amr_set_runtime_parameters] NOTE')
           end if
        end if


        call RuntimeParameters_get("gr_pmrpMaxblocks", pmParms%IntParms(1) )
        call RuntimeParameters_get("gr_pmrpNdim",      pmParms%IntParms(2) )
        call RuntimeParameters_get("gr_pmrpL2p5d",     pmParms%IntParms(3) )
        call RuntimeParameters_get("gr_pmrpNxb",       pmParms%IntParms(4) )
        call RuntimeParameters_get("gr_pmrpNyb",       pmParms%IntParms(5) )
        call RuntimeParameters_get("gr_pmrpNzb",       pmParms%IntParms(6) )
        call RuntimeParameters_get("gr_pmrpNvar",      pmParms%IntParms(7) )
        call RuntimeParameters_get("gr_pmrpNfacevar",  pmParms%IntParms(8) )
        call RuntimeParameters_get("gr_pmrpNvaredge",  pmParms%IntParms(9) )
        call RuntimeParameters_get("gr_pmrpNvarcorn",  pmParms%IntParms(10) )
        call RuntimeParameters_get("gr_pmrpNvarWork",  pmParms%IntParms(11) )
        call RuntimeParameters_get("gr_pmrpNguard",    pmParms%IntParms(12) )
        call RuntimeParameters_get("gr_pmrpNguardWork",pmParms%IntParms(13) )
        call RuntimeParameters_get("gr_pmrpNfluxvar",  pmParms%IntParms(14) )
        call RuntimeParameters_get("gr_pmrpNedgevar1", pmParms%IntParms(15) )
        call RuntimeParameters_get("gr_pmrpIfaceOff",  pmParms%IntParms(16) )
        call RuntimeParameters_get("gr_pmrpMflags",    pmParms%IntParms(17) )
        call RuntimeParameters_get("gr_pmrpNfieldDivf",pmParms%IntParms(18) )
        call RuntimeParameters_get("gr_pmrpNboundaries",pmParms%IntParms(19) )

        call rpget_LogAsInt       ("gr_pmrpDiagonals",            pmParms%LogParms(1) )
        call rpget_LogAsInt       ("gr_pmrpAmrErrorChecking",     pmParms%LogParms(2) )
        call rpget_LogAsInt       ("gr_pmrpNoPermanentGuardcells",pmParms%LogParms(3) )
        call rpget_LogAsInt       ("gr_pmrpAdvanceAllLevels", pmParms%LogParms(4) )
        call rpget_LogAsInt       ("gr_pmrpForceConsistency", pmParms%LogParms(5) )
        call rpget_LogAsInt       ("gr_pmrpConsvFluxes", pmParms%LogParms(6) )
        call rpget_LogAsInt       ("gr_pmrpConsvFluxDensities", pmParms%LogParms(7) )
        call rpget_LogAsInt       ("gr_pmrpEdgeValue", pmParms%LogParms(8) )
        call rpget_LogAsInt       ("gr_pmrpEdgeValueInteg", pmParms%LogParms(9) )
        call rpget_LogAsInt       ("gr_pmrpVarDt", pmParms%LogParms(10) )
        call rpget_LogAsInt       ("gr_pmrpPredCorr", pmParms%LogParms(11) )
        call rpget_LogAsInt       ("gr_pmrpEmptyCells", pmParms%LogParms(12) )
        call rpget_LogAsInt       ("gr_pmrpConserve", pmParms%LogParms(13) )
        call RuntimeParameters_get("gr_pmrpDivergenceFree", pmParms%LogParms(14) )
        call rpget_LogAsInt       ("gr_pmrpCurvilinear", pmParms%LogParms(15) )
        call rpget_LogAsInt       ("gr_pmrpCurvilinearConserve", pmParms%LogParms(16) )
        call rpget_LogAsInt       ("gr_pmrpCartesianPm", pmParms%LogParms(17) )
        call rpget_LogAsInt       ("gr_pmrpCylindricalPm", pmParms%LogParms(18) )
        call rpget_LogAsInt       ("gr_pmrpSphericalPm", pmParms%LogParms(19) )
        call rpget_LogAsInt       ("gr_pmrpPolarPm", pmParms%LogParms(20) )
        call rpget_LogAsInt       ("gr_pmrpLsingularLine", pmParms%LogParms(21) )
        call rpget_LogAsInt       ("gr_pmrpTimingMpi", pmParms%LogParms(22) )
        call rpget_LogAsInt       ("gr_pmrpTimingMpix", pmParms%LogParms(23) )

        if (pmParms%IntParms(1)==-1)  pmParms%IntParms(1) = MAXBLOCKS
        if (pmParms%IntParms(2)==-1)  pmParms%IntParms(2) = NDIM
        if (pmParms%IntParms(3)==-1) then
#ifdef FLASH_2P5DIM
           pmParms%IntParms(3) = 1
#else
           pmParms%IntParms(3) = 0
#endif
        end if
        if (pmParms%IntParms(4)==-1)  pmParms%IntParms(4) = NXB
        if (pmParms%IntParms(5)==-1)  pmParms%IntParms(5) = NYB
        if (pmParms%IntParms(6)==-1)  pmParms%IntParms(6) = NZB
        if (pmParms%IntParms(7)==-1)  pmParms%IntParms(7) = NUNK_VARS
        if (pmParms%IntParms(8)==-1)  pmParms%IntParms(8) = NFACE_VARS
        if (pmParms%IntParms(12)==-1) pmParms%IntParms(12) = NGUARD
        if (pmParms%IntParms(13)==-1) pmParms%IntParms(13) = NGUARD
        if (pmParms%IntParms(14)==-1) pmParms%IntParms(14) = NFLUXES
        if (pmParms%IntParms(15)==-1) then
#ifdef FLASH_NEDGE_VAR
           pmParms%IntParms(15) = FLASH_NEDGE_VAR
#else
           pmParms%IntParms(15) = 0
#endif
        end if
        if (pmParms%IntParms(18)==-1) then
#ifdef FLASH_NFIELD_DIVF
           pmParms%IntParms(18) = FLASH_NFIELD_DIVF
#else
           pmParms%IntParms(18) = 0
#endif
        end if
        if (pmParms%LogParms(14)==-1) then
           pmParms%LogParms(14) = 0
#ifdef DIVERGENCE_FREE
           pmParms%LogParms(14) = 1
#endif
        end if



#endif

#ifdef LIBRARY
#define SET_OR_CHECK(VAR,FIELD)  VAR = pmParms % FIELD
#else
#define SET_OR_CHECK(VAR,FIELD) if(VAR.NE.pmParms%FIELD) call die('VAR',pmParms%FIELD,VAR)
#endif
        
SET_OR_CHECK(maxblocks,IntParms(1))
SET_OR_CHECK(ndim,IntParms(2))
        l2p5d                   = pmParms % IntParms (3)
SET_OR_CHECK(nxb,IntParms(4))
SET_OR_CHECK(nyb,IntParms(5))
SET_OR_CHECK(nzb,IntParms(6))
SET_OR_CHECK(nvar,IntParms(7))
SET_OR_CHECK(nfacevar,IntParms(8))
        nvaredge                = pmParms % IntParms (9)
        nvarcorn                = pmParms % IntParms (10)
        nvar_work               = pmParms % IntParms (11)
SET_OR_CHECK(nguard,IntParms(12))
SET_OR_CHECK(nguard_work,IntParms(13))
        nfluxvar                = pmParms % IntParms (14)
        nedgevar1               = pmParms % IntParms (15)
        iface_off               = pmParms % IntParms (16)
        mflags                  = pmParms % IntParms (17)
        nfield_divf             = pmParms % IntParms (18)
        nboundaries             = pmParms % IntParms (19)

        diagonals               = (pmParms % LogParms (1) == 1)
        amr_error_checking      = (pmParms % LogParms (2) == 1)
        no_permanent_guardcells = (pmParms % LogParms (3) == 1)
        advance_all_levels      = (pmParms % LogParms (4) == 1)
        force_consistency       = (pmParms % LogParms (5) == 1)
        consv_fluxes            = (pmParms % LogParms (6) == 1)
        consv_flux_densities    = (pmParms % LogParms (7) == 1)
        edge_value              = (pmParms % LogParms (8) == 1)
        edge_value_integ        = (pmParms % LogParms (9) == 1)
        var_dt                  = (pmParms % LogParms (10) == 1)
        pred_corr               = (pmParms % LogParms (11) == 1)
        empty_cells             = (pmParms % LogParms (12) == 1)
        conserve                = (pmParms % LogParms (13) == 1)
        divergence_free         = (pmParms % LogParms (14) == 1)
        curvilinear             = (pmParms % LogParms (15) == 1)
        curvilinear_conserve    = (pmParms % LogParms (16) == 1)
        cartesian_pm            = (pmParms % LogParms (17) == 1)
        cylindrical_pm          = (pmParms % LogParms (18) == 1)
        spherical_pm            = (pmParms % LogParms (19) == 1)
        polar_pm                = (pmParms % LogParms (20) == 1)
        lsingular_line          = (pmParms % LogParms (21) == 1)
        timing_mpi              = (pmParms % LogParms (22) == 1)
        timing_mpix             = (pmParms % LogParms (23) == 1)

#ifdef USE_AMR_RUNTIME_PARAMETERS_FILE
        output_dir = ""
        do i = 1, NumChar
           if (pmParms % CharParms (i) == char(0)) then
              exit
           end if
           output_dir(i:i) = pmParms % CharParms (i)
        end do
#else
        call RuntimeParameters_get("gr_pmrpOutputDir", output_dir)
#endif

        amr_log_file = trim(output_dir) // 'amr.log'

        if (gr_geometry .NE. CARTESIAN .OR. cylindrical_pm .OR. spherical_pm .OR. polar_pm) then
           if (.not. curvilinear) then
              if (gr_globalMe==0) &
                   print*,'NOTE: Enabling curvilinear, cartesian_pm/cylindrical_pm/spherical_pm/polar_pm  so far was',&
                   cartesian_pm,cylindrical_pm,spherical_pm,polar_pm
              call Logfile_stamp('Enabling curvilinear suppport because of the selected geometry',&
                   '[amr_set_runtime_parameters] NOTE')
              curvilinear = .TRUE.
              curvilinear_conserve = .TRUE.
           end if
        end if

        Return

        contains
        subroutine die(var,filevalue,codevalue)
          character(len=*),intent(in) :: var
          integer,intent(in)          :: filevalue, codevalue
9         format('[amr_set_runitme_parameters]: The value of ',a,' read from the file amr_runtime_parameters &
               &is',I11,';')
10        format('the value compiled into the executable is',I11,', but they are required to be the same!')
          
          print  9,var,filevalue
          print 10,codevalue
          call Driver_abort("Bad value of "//var)
        end subroutine die

#ifndef USE_AMR_RUNTIME_PARAMETERS_FILE
        subroutine rpget_LogAsInt(name, value)
          character(len=*), intent(in)          :: name
          integer, intent(out)                  :: value
          logical :: logicalValue
          call RuntimeParameters_get(name, logicalValue)
          if(logicalValue) then
             value = 1
          else
             value = 0
          end if
        end subroutine rpget_LogAsInt
#endif

        End Subroutine amr_set_runtime_parameters
