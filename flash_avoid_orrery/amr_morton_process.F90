!!****if* source/amr_morton_process
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
!!   amr_morton_process
!!
!! SYNOPSIS
!!
!!   Call amr_morton_process()
!!
!! ARGUMENTS
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
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_morton_bnd
!!   mpi_morton_bnd_prolong
!!   mpi_morton_bnd_fluxcon
!!   mpi_morton_bnd_restrict
!!   mpi_setup
!!
!! RETURNS
!!
!!   Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine manages the construction of list of blocks which are
!!   communicated for guardcell filling, prolongation, restriction, and
!!   flux and edge fixups.  It does this by simply calling the routines
!!   mpi_morton_bnd, mpi_morton_bnd_prolong, mpi_morton_bnd_restrict and
!!   mpi_morton_bnd_fluxcon.  
!!
!! AUTHORS
!!
!!    Kevin Olson
!!
!!***

#include "constants.h"
#include "paramesh_preprocessor.fh"

      Subroutine amr_morton_process()

!-----Use Statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use timings
      Use io
      Use paramesh_mpi_interfaces, Only : mpi_morton_bnd,              & 
                                          mpi_morton_bnd_prolong,      & 
                                          mpi_morton_bnd_fluxcon,      & 
                                          mpi_morton_bnd_restrict,     & 
                                          mpi_setup
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Include Statements
      Include 'mpif.h'

!-----Local Variables
      Integer :: nprocs,mype,tag_offset,ierr

!-----Begin Executable code.

      CustomFlashVersion: If (use_flash_surr_blks_fill) then

         call amr_morton_process_flash ()

      Else

      Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
      Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

!-----call setup routines in preparation for calling 
!-----all the mpi_morton_bnd_XXX routines.

!-----Find the coordinate ranges
      Call mpi_amr_global_domain_limits
      Call mpi_setup(mype,nprocs)

!-----Set up surrounding blocks of all local blocks (must not precede
!-----setting of grid_xmin,... etc)

      Call find_surrblks()
      gsurrblks_set = +1

!-----Create guardcell filling communications information
      tag_offset = 100
      Call mpi_morton_bnd(mype,nprocs,tag_offset)
!-----Create prolongation communications information
      tag_offset = 100
      Call mpi_morton_bnd_prolong(mype,nprocs,tag_offset)
!-----Create flux fix and edge fix communications information
      tag_offset = 100
      Call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
!-----Create restriction communications information
      tag_offset = 100
      Call mpi_morton_bnd_restrict(mype,nprocs,tag_offset)


      End If CustomFlashVersion

      Return
      End Subroutine amr_morton_process


!FLASH's custom version of the same subroutine.  There are sufficient
!difference between this and the original that it is cleaner just
!to create a new subroutine.

!We can check flash surr_blks against orrery surr_blks by defining
!DEBUG_SURR_BLKS.
Subroutine amr_morton_process_flash()

  !-----Use Statements
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use timings
  Use io

  Use paramesh_mpi_interfaces, Only : mpi_morton_bnd,              & 
       mpi_morton_bnd_prolong,      & 
       mpi_morton_bnd_fluxcon,      & 
       mpi_morton_bnd_restrict,     & 
       mpi_setup

      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

  Implicit None

  !-----Include Statements
  Include 'mpif.h'

  !-----Local Variables
  Integer :: nprocs,mype,tag_offset,ierr
#ifdef DEBUG_SURR_BLKS
  Integer :: surr_blkst(3,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
  integer :: li, lj, lk, i, errorcode
  integer, dimension(mdim) :: gCell
  logical :: error
#endif
  !-----Begin Executable code.

  Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
  Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

  !-----call setup routines in preparation for calling 
  !-----all the mpi_morton_bnd_XXX routines.

  !-----Find the coordinate ranges
  Call mpi_amr_global_domain_limits
  Call mpi_setup(mype,nprocs)

  !-----Set up surrounding blocks of all local blocks (must not precede
  !-----setting of grid_xmin,... etc)
  if (.not.surr_blks_valid) then
     if (myPE == 0) then
        print *, "[amr_morton_process]: Initializing surr_blks using "//&
             "standard orrery implementation"
     end if
     Call find_surrblks()
  else

     !Do nothing because surr_blks is in a valid state.

#ifdef DEBUG_SURR_BLKS
     !In debug mode we compare our surr_blks calculation to orrery.

     if (myPE == 0) then
        print *, "[amr_morton_process]: Checking Flash surr_blks against "//&
             "standard orrery implementation"
     end if

     surr_blkst(1:3,:,:,:,1:lnblocks) = surr_blks(1:3,:,:,:,1:lnblocks)
     Call find_surrblks()
     Do i = 1,lnblocks
        do lk = 1, 1+2*k3d
           do lj = 1, 1+2*k2d
              do li = 1, 1+2*k1d
                 gCell = (/li,lj,lk/)
                 if ( surr_blks(1,li,lj,lk,i) /= surr_blkst(1,li,lj,lk,i) .or. &
                      surr_blks(2,li,lj,lk,i) /= surr_blkst(2,li,lj,lk,i) .or. &
                      surr_blks(3,li,lj,lk,i) /= surr_blkst(3,li,lj,lk,i) ) then

#ifdef DEBUG_SURR_BLKS_IGNORE_BOUNDARIES
                    !It takes two iterations of amr_morton_process before the
                    !physical boundaries are stored.  We need to define this
                    !macro if checking the surr_blks read from checkpoint file.
                    error = &
                          surr_blks(1,li,lj,lk,i) > PARAMESH_PHYSICAL_BOUNDARY .and. &
                          surr_blkst(1,li,lj,lk,i) > PARAMESH_PHYSICAL_BOUNDARY
#else
                    error = .true.
#endif

                    if (error) then
                       print *, "MISMATCH NEIGH! blk/proc", i, myPE, &
                            ", element", gCell(1:ndim), &
                            ". Mine:", surr_blkst(1:3,li,lj,lk,i), &
                            ", Actual:", surr_blks(1:3,li,lj,lk,i)
                       call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
                    end if
                 end if
              end do
           end do
        end do
     end Do

     call MPI_Barrier(amr_mpi_meshComm, ierr)
     if (myPE == 0) then
        print *, "[amr_morton_process]: Flash surr_blks is correct"
     end if

#endif
  end if
  gsurrblks_set = +1

  !-----Create guardcell filling communications information
  tag_offset = 100
  Call mpi_morton_bnd(mype,nprocs,tag_offset)
  !-----Create prolongation communications information
  tag_offset = 100
  Call mpi_morton_bnd_prolong(mype,nprocs,tag_offset)
  !-----Create flux fix and edge fix communications information
  tag_offset = 100
  Call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
  !-----Create restriction communications information
  tag_offset = 100
  Call mpi_morton_bnd_restrict(mype,nprocs,tag_offset)

  Return
End Subroutine amr_morton_process_flash
