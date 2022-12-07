!!****if* source/mpi_morton_bnd_restrict
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
!!   mpi_morton_bnd_restrict
!!
!! SYNOPSIS
!!
!!   Call mpi_morton_bnd_restrict(mype, nprocs, tag_offset, subPatNo)
!!   Call mpi_morton_bnd_restrict(integer, integer, integer, integer)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in)    :: mype       Local processor id.
!!   Integer, Intent(in)    :: nprocs     Number of processors.
!!   Integer, Intent(inout) :: tag_offset A unique id used in marking messages.
!!   subPatNo - Integer(in), OPTIONAL ::  request computation of non-default
!!                                        variant of the communication pattern
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   Flashx_mpi_implicitNone.fh
!!
!! USES
!!
!!   gr_pmCommDataTypes
!!   gr_pmCommPatternData
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   timings
!!   mpi_morton
!!   constants
!!   gr_pmCommDataTypes
!!   gr_pmCommPatternData
!!
!! CALLS
!!
!!    mpi_amr_write_restrict_comm
!!    process_fetch_list
!!    gr_pmCommPatternPtr
!!
!! RETURNS
!!
!!    Does not return anything.
!!
!! DESCRIPTION
!!
!!   This routine does a communications analysis for data restriction and
!!   constructs and stores lists of off-processor blocks which are to be 
!!   communicated during restriction.  Also stored are which sections of 
!!   the blocks to be fetched.
!!
!! AUTHORS
!!
!!    Written by Peter MacNeice  and Michael Gehmeyr, February 2000.
!!    Major simplification and rewrite by Kevin Olson, August 2007.
!!    Optional arg subPatNo                 Klaus Weide May 2022
!!
!! MODIFICATIONS
!!  2022-05-13 K. Weide  Use local pattern pointer to access comm pattern
!!  2022-05-20 K. Weide  Variant pattern for subPatNo=GRID_SUBPAT_RESTRICT_ANC
!!  2022-05-25 K. Weide  Variant pattern for GRID_SUBPAT_RESTRICT_FOR_FCORR
!!  2022-11-08 Klaus Weide  Added ONLY to USE physicaldata
!!  2022-11-30 Klaus Weide  USE paramesh_dimensions with ONLY
!!***

#include "paramesh_preprocessor.fh"
#include "Simulation.h"

      Subroutine mpi_morton_bnd_restrict (mype,                        &
                                          nprocs,                      &
                                          tag_offset,                  &
                                          subPatNo)

!-----Use Statements
      use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t, &
           GRID_PAT_RESTRICT, &
           GRID_SUBPAT_RESTRICT_DEFAULT, GRID_SUBPAT_RESTRICT_ANC, &
           GRID_SUBPAT_RESTRICT_FOR_FCORR
      use gr_pmCommPatternData, ONLY: gr_pmCommPatternPtr, &
           gr_pmPrintCommPattern
      Use paramesh_dimensions, ONLY: k2d, k3d
      Use physicaldata, ONLY: advance_all_levels
      Use tree
      Use timings
      Use mpi_morton, ONLY: npts_neigh
      Use constants

      Use paramesh_mpi_interfaces, only : mpi_amr_write_restrict_comm, & 
                                          process_fetch_list
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

!-----Include Statements
#include "Flashx_mpi_implicitNone.fh"

!-----Input/Output Variables
      Integer, Intent(in)    ::  mype,nprocs
      Integer, Intent(inout) ::  tag_offset
      Integer,OPTIONAL, intent(in) ::  subPatNo

!-----Local variables
      Integer :: lb,i,j,k,j00
      Integer :: ierror
      Integer :: istack
      Integer :: iproc
      Integer :: npts_neigh1,npts_neigh2
      Integer,Dimension (:),  Allocatable :: n_to_left
      Integer,Dimension (:,:),Allocatable :: fetch_list
      Integer,Dimension (:,:),Allocatable :: tfetch_list
      TYPE(gr_pmCommPattern_t),pointer :: pattern
      integer :: subPatLoc
      integer :: jf,dtype,dir0,b0,bc0

!-----Begin executable code.
      npts_neigh1 = npts_neigh
      npts_neigh2 = npts_neigh+100
      Allocate(fetch_list(3,npts_neigh2))
      Allocate(n_to_left(0:nprocs-1))

!----COMPUTE the number of blocks to the 'left' (ie. stored on processors with
!----smaller process ids) of every other processor
      Call MPI_ALLGATHER(lnblocks,                                     &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         n_to_left,                                    &
                         1,                                            &
                         MPI_INTEGER,                                  &
                         amr_mpi_meshComm,                               &
                         ierror)
                        
      Do iproc = nprocs-1, 1, -1
         n_to_left(iproc) = n_to_left(iproc-1)
      End Do
      n_to_left(iproc) = 0
      Do iproc = 2, nprocs-1
         n_to_left(iproc) = n_to_left(iproc) + n_to_left(iproc-1)
      End Do

      subPatLoc = GRID_SUBPAT_RESTRICT_DEFAULT
      if(present(subPatNo)) subPatLoc = subPatNo
      pattern => gr_pmCommPatternPtr(GRID_PAT_RESTRICT,subPatNo)

!-----Initializations
      pattern % commatrix_send(:) = 0
      pattern % commatrix_recv(:) = 0

!-----Construct a list of potential neighbors of all blocks on this
!-----processor, and potential neighbors of their parents.
!-----Exclude any which are on this processor.

      istack = 0
      Do lb = 1, lnblocks

      If (nodetype(lb) == 2 .or.                                       &
           (advance_all_levels .and. nodetype(lb) == 3)) Then
        if (advance_all_levels .OR. &
            (subPatLoc == GRID_SUBPAT_RESTRICT_DEFAULT) .OR.        &
            (subPatLoc == GRID_SUBPAT_RESTRICT_ANC .AND.        &
             any(surr_blks(1,1:3,1:1+2*k2d,1:1+2*k3d,lb) > 0 .and.      &
                 surr_blks(3,1:3,1:1+2*k2d,1:1+2*k3d,lb) == 1)) ) then

!------ADD OFF PROCESSOR CHILDREN OF BLOCK 'lb' TO FETCH LIST
        Do i = 1,nchild
        If (child(1,i,lb) > 0 .and.                            & 
            child(2,i,lb) .ne. mype) Then

            istack = istack + 1
            If (istack > npts_neigh1) Call expand_fetch_list
            fetch_list(1,istack) = child(1,i,lb)
            fetch_list(2,istack) = child(2,i,lb)
!-----------Fetch entire block
            fetch_list(3,istack) = 14

         End If  ! End If child
         End Do  ! End Do i = 1,nchild

        else if ((subPatLoc == GRID_SUBPAT_RESTRICT_FOR_FCORR .AND.        &
             hasLeafFaceNeighs(surr_blks(:,1:3,1:1+2*k2d,1:1+2*k3d,lb))) ) then
           childLoop:Do i = 1,nchild
              If (child(1,i,lb) > 0 .and. &
                  child(2,i,lb) .ne. mype) Then
                 dtype = 0
                 do jf = 1,nfaces
                    if (isLeafFaceNeigh(surr_blks(:,1:3,1:1+2*k2d,1:1+2*k3d,lb),jf)) then
                       dir0 = (jf-1) / 2
                       b0 = mod(jf-1,2)
                       bc0 = mod((i-1)/(2**dir0),2)
!!$                       print*,'i,jf,dir0,b0,bc0,dtype:',i,jf,dir0,b0,bc0,dtype
                       if (bc0 == b0) then
                          if (dtype == 0) then
                             dtype = 14 + (2*b0-1)*3**dir0
                          else
                             dtype = 14
                          end if

                          if (istack == 0 .OR. dtype .NE. 14) istack = istack + 1
                          If (istack > npts_neigh1) Call expand_fetch_list
                          fetch_list(1,istack) = child(1,i,lb)
                          fetch_list(2,istack) = child(2,i,lb)
                          fetch_list(3,istack) = dtype
                          if(dtype == 14) cycle childLoop
                       end if
                    end if
                 end do
              End If  ! End If child
           End Do childLoop  ! End Do i = 1,nchild

      end if  ! End if (advance_all_levels .OR. ...
      End If  ! End If (nodetype(lb) <= 2 .or. (advance_all_levels ...))

      End Do  ! End Do lb = 1, lnblocks

      Call process_fetch_list(pattern,                                 &
                              fetch_list,                              &
                              istack,                                  &
                              mype,                                    &
                              nprocs,                                  &
                              n_to_left,                               &
                              tag_offset)

!------Store communication info for future use
!!$       Call mpi_amr_write_restrict_comm(nprocs) !This is a no-op now.

!------Deallocate any memory which was dynamically allocated for local 
!------use in this routine.
       If (Allocated(fetch_list)) deallocate(fetch_list)
       If (Allocated(n_to_left)) deallocate(n_to_left)
#ifdef DEBUG
       call gr_pmPrintCommPattern(pattern,'mmbr:pattern',mype)
#endif
      Return

Contains
  logical function hasLeafFaceNeighs(surr)
    implicit none
    integer,intent(IN) :: surr(3,3,1+2*K2D,1+2*K3D)

    hasLeafFaceNeighs = .FALSE.
    if     (surr(1,1,1+K2D,1+K3D)>0 .AND. surr(3,1,1+K2D,1+K3D)==1) then
       hasLeafFaceNeighs = .TRUE.
    elseif (surr(1,3,1+K2D,1+K3D)>0 .AND. surr(3,3,1+K2D,1+K3D)==1) then
       hasLeafFaceNeighs = .TRUE.
#if NDIM > 1
    elseif (surr(1,2,1,1+K3D)>0 .AND. surr(3,2,1,1+K3D)==1) then
       hasLeafFaceNeighs = .TRUE.
    elseif (surr(1,2,3,1+K3D)>0 .AND. surr(3,2,3,1+K3D)==1) then
       hasLeafFaceNeighs = .TRUE.
#if NDIM == 3
    elseif (surr(1,2,2,1)>0 .AND. surr(3,2,2,1)==1) then
       hasLeafFaceNeighs = .TRUE.
    elseif (surr(1,2,2,3)>0 .AND. surr(3,2,2,3)==1) then
       hasLeafFaceNeighs = .TRUE.
#endif
#endif
    end if
  end function hasLeafFaceNeighs

  logical function isLeafFaceNeigh(surr,jf)
    implicit none
    integer,intent(IN) :: surr(3,3,1+2*K2D,1+2*K3D)
    integer,intent(IN) :: jf

    isLeafFaceNeigh = .FALSE.
    select case (jf)
    case(1)
       isLeafFaceNeigh = (surr(1,1,1+K2D,1+K3D)>0 .AND. surr(3,1,1+K2D,1+K3D)==1)
    case(2)
       isLeafFaceNeigh = (surr(1,3,1+K2D,1+K3D)>0 .AND. surr(3,3,1+K2D,1+K3D)==1)
#if NDIM > 1
    case(3)
       isLeafFaceNeigh = (surr(1,2,1,1+K3D)>0 .AND. surr(3,2,1,1+K3D)==1)
    case(4)
       isLeafFaceNeigh = (surr(1,2,3,1+K3D)>0 .AND. surr(3,2,3,1+K3D)==1)
#if NDIM == 3
    case(5)
       isLeafFaceNeigh = (surr(1,2,2,1)>0 .AND. surr(3,2,2,1)==1)
    case(6)
       isLeafFaceNeigh = (surr(1,2,2,3)>0 .AND. surr(3,2,2,3)==1)
#endif
#endif
    end select
  end function isLeafFaceNeigh

        Subroutine expand_fetch_list

               If (Allocated(tfetch_list)) deallocate(tfetch_list)
               Allocate(tfetch_list(3,npts_neigh2))
               tfetch_list(:,:istack-1) = fetch_list(:,:istack-1)
               npts_neigh1 = npts_neigh1 + 3000
               npts_neigh2 = npts_neigh2 + 3000
               deallocate(fetch_list)
               Allocate(fetch_list(3,npts_neigh2))
               fetch_list(:,:istack-1) = tfetch_list(:,:istack-1)
               Deallocate(tfetch_list)

        End Subroutine expand_fetch_list
      End Subroutine mpi_morton_bnd_restrict
