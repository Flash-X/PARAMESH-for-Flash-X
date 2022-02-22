!!****if* source/Grid/GridMain/paramesh/gr_amr_dump_runtime_parameters
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2017 The University of Chicago
!!  Copyright (C) 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!  gr_amr_dump_runtime_parameters
!!
!! SYNOPSIS
!!
!!  call gr_amr_dump_runtime_parameters()
!!
!! DESCRIPTION
!!
!!  Dump AMR runtime parameters into a file amr_runtime_parameters.dump,
!!  in a form that can be read by PARAMESH subroutine
!!  amr_set_runtime_parameters in subsequent runs when PARAMESH has been
!!  compiled in LIBRARY mode.
!!
!!
!! ARGUMENTS
!!
!!   NONE 
!!
!! NOTES
!!
!!   AMR runtime parameters are not the same as FLASH runtime parameters.
!!   AMR runtime parameters are meaningful to PARAMESH3 (or later) in
!!   LIBRARY mode; they specify dimensions etc. which are otherwise
!!   hardcoded as F90 parameters.
!!
!!   PARAMESH need not be compiled in LIBRARY mode in ordere to call
!!   this subroutine and dump AMR runtime parameters; this is only
!!   required in order to use the resulting file (after renaming)
!!   for initialization.
!!
!!   Subroutine added to FLASH by KW, modeled on amr_set_runtime_parameters
!!   from PARAMESH 3.4. This file should be kept in synch with
!!   the PARAMESH subroutine amr_set_runtime_parameters.
!!
!!***

#include "paramesh_preprocessor.fh"
        subroutine gr_amr_dump_runtime_parameters()

        use paramesh_dimensions
        use physicaldata
        use tree
        use timings
        use io
        use Grid_data, ONLY : gr_meshMe, gr_meshNumProcs

        implicit none


        ! If there are "many" processors and this is not the master processor,
        ! do not write to file at all but return immediately. - KW
        if (gr_meshMe .NE. 0 .AND. gr_meshNumProcs .GT. 8)  return

        open (unit=35,                        &
     &        file='amr_runtime_parameters.dump',  &
     &        status='UNKNOWN',               &
     &        action='WRITE',                 &
     &        form='formatted')
! integers
        write (35,*) maxblocks        ,', maxblocks'
        write (35,*) ndim        ,', ndim'
        write (35,*) l2p5d        ,', l2p5d'
        write (35,*) nxb        ,', nxb'
        write (35,*) nyb        ,', nyb'
        write (35,*) nzb        ,', nzb'
        write (35,*) nvar        ,', nvar'
        write (35,*) nfacevar        ,', nfacevar'
        write (35,*) nvaredge        ,', nvaredge'
        write (35,*) nvarcorn        ,', nvarcorn'
        write (35,*) nvar_work        ,', nvar_work'
        write (35,*) nguard        ,', nguard'
        write (35,*) nguard_work        ,', nguard_work'
        write (35,*) nfluxvar        ,', nfluxvar'
        write (35,*) nedgevar1        ,', nedgevar1'
        write (35,*) iface_off        ,', iface_off'
        write (35,*) mflags        ,', mflags'
        write (35,*) nfield_divf        ,', nfield_divf'
        write (35,*) nboundaries        ,', nboundaries'
! logicals
        write (35,*) diagonals        ,', diagonals'
        write (35,*) amr_error_checking        ,', amr_error_checking'
        write (35,*) no_permanent_guardcells        ,', no_permanent_guardcells'
        write (35,*) advance_all_levels        ,', advance_all_levels'
        write (35,*) force_consistency        ,', force_consistency'
        write (35,*) consv_fluxes        ,', consv_fluxes'
        write (35,*) consv_flux_densities        ,', consv_flux_densities'
        write (35,*) edge_value        ,', edge_value'
        write (35,*) edge_value_integ        ,', edge_value_integ'
        write (35,*) var_dt        ,', var_dt'
        write (35,*) pred_corr        ,', pred_corr'
        write (35,*) empty_cells        ,', empty_cells'
        write (35,*) conserve        ,', conserve'
        write (35,*) divergence_free        ,', divergence_free'
        write (35,*) curvilinear        ,', curvilinear'
        write (35,*) curvilinear_conserve        ,', curvilinear_conserve'
        write (35,*) cartesian_pm        ,', cartesian_pm'
        write (35,*) cylindrical_pm        ,', cylindrical_pm'
        write (35,*) spherical_pm        ,', spherical_pm'
        write (35,*) polar_pm        ,', polar_pm'
        write (35,*) lsingular_line        ,', lsingular_line'
        write (35,*) timing_mpi        ,', timing_mpi'
        write (35,*) timing_mpix        ,', timing_mpix'
! characters
!!        write (35,*) "'",output_dir,"'"        ,', output_dir'
        write (35,'(A1,2A)') "'",trim(output_dir),"' , output_dir"
        write (35,*) '====='
#ifdef LIBRARY
        write (35,*) 'LIBRARY was defined'
#else
        write (35,*) 'LIBRARY was undefined'
#endif
#ifdef CURVILINEAR
        write (35,*) 'CURVILINEAR was defined'
#else
        write (35,*) 'CURVILINEAR was undefined'
#endif
#ifdef CURVILINEAR_CONSERVE
        write (35,*) 'CURVILINEAR_CONSERVE was defined'
#else
        write (35,*) 'CURVILINEAR_CONSERVE was undefined'
#endif
#ifdef SPHERICAL
        write (35,*) 'SPHERICAL was defined'
#else
        write (35,*) 'SPHERICAL was undefined'
#endif
#ifdef CYLINDRICAL
        write (35,*) 'CYLINDRICAL was defined'
#else
        write (35,*) 'CYLINDRICAL was undefined'
#endif

        close(35)

        return
        end subroutine gr_amr_dump_runtime_parameters
