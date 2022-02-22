!!****fi source/mpi_amr_exchange_nodetype
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2010 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!   mpi_amr_exchange_nodetype
!!
!! SYNOPSIS
!!
!!   call mpi_amr_exchange_nodetype (myPE, localNumBlocks)
!!   call mpi_amr_exchange_nodetype (integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: myPE
!!     My processor ID.
!!   integer, intent(in) :: localNumBlocks
!!     Number of blocks on this processor.  This subroutine will 
!!     exchange nodetype information for block IDs 1 to localNumBlocks.
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
!! 
!! RETURNS
!!
!!   Nothing returned
!!
!! DESCRIPTION
!!
!!   This subroutine exchanges the surr_blks nodetype field between 
!!   blocks.  On entry it is assumed that nodetype array and surr_blks 
!!   block and processor field are up to date.  This subroutine has
!!   been introduced by the Flash center.
!!
!! AUTHORS
!!
!!   Chris Daley (2010).
!!
!!***

#include "paramesh_preprocessor.fh"
  
Subroutine mpi_amr_exchange_nodetype (myPE, localNumBlocks)
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

  Implicit None
  !-----Include statements.
  Include 'mpif.h'

  !-----Input/Output statements.
  Integer, Intent(in)    :: myPE, localNumBlocks

  !-----Local array and variables.
  Integer :: i, ineigh, ineigh_proc, nsend, nrecv, ierr, li, lj, lk
  Integer :: reqr(maxblocks), reqs(maxblocks)
  Integer :: statr(MPI_STATUS_SIZE,maxblocks), stats(MPI_STATUS_SIZE,maxblocks)
  integer, dimension(mdim) :: gCell, gs, gr
  Integer :: errorcode

  !Assert that the communication arrays are large enough (in case we use
  !this subroutine differently in future).
  if (localNumBlocks > maxblocks) then
     print *, "[mpi_amr_exchange_nodetype] Stack arrays are too small! "//&
          "Local blocks", localNumBlocks, "Max blocks:", maxblocks
     call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
  end if


  !NOTE: If the communication is changed so that receives are 
  !posted for *all* directions initially then req[rs] and stat[rs] size
  !must increase to:    
  !Integer :: stat[rs](MPI_STATUS_SIZE,maxblocks * (3**ndim - 1))
  !Integer :: req[rs](maxblocks * (3**ndim - 1))


  !-----Exchange node types between neighbors one direction at a time.
  kAxis: do lk = -k3d, k3d
     jAxis: do lj = -k2d, k2d
        iAxis: do li = -k1d, k1d
           gCell = (/li,lj,lk/)
           validRegion: if ( (0 /= sum(abs(gCell))) ) then

              gs = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)
              gr = (/(-li+k1d+1),(-lj+k2d+1),(-lk+k3d+1)/)


              nrecv = 0
              Do i = 1, localNumBlocks
                 If (surr_blks(1,gr(1),gr(2),gr(3),i) > 0) Then
                    If (surr_blks(2,gr(1),gr(2),gr(3),i).ne.myPE) Then
                       nrecv = nrecv + 1
                       Call MPI_IRECV(surr_blks(3,gr(1),gr(2),gr(3),i), & 
                            1, MPI_INTEGER, &
                            surr_blks(2,gr(1),gr(2),gr(3),i), &
                            surr_blks(1,gr(1),gr(2),gr(3),i), &
                            amr_mpi_meshComm, reqr(nrecv), ierr)
                    End If
                 End If
              End Do


              nsend = 0
              Do i = 1, localNumBlocks
                 ineigh = surr_blks(1,gs(1),gs(2),gs(3),i)
                 ineigh_proc = surr_blks(2,gs(1),gs(2),gs(3),i)
                 If (ineigh > 0) Then
                    If (surr_blks(2,gs(1),gs(2),gs(3),i).ne.myPE) Then
                       nsend = nsend + 1
                       Call MPI_ISEND(nodetype(i), 1, MPI_INTEGER, &
                            ineigh_proc, i, amr_mpi_meshComm, reqs(nsend), ierr)
                    Else
                       surr_blks(3,gr(1),gr(2),gr(3),ineigh) = nodetype(i)
                    End If
                 End If
              End Do


              If (nrecv > 0) Then
                 Call MPI_WAITALL(nrecv,reqr,statr,ierr)
              End If
              If (nsend > 0) Then
                 Call MPI_WAITALL(nsend,reqs,stats,ierr)
              End If


           else
              Do i = 1, localNumBlocks
                 surr_blks(3,1+k1d,1+k2d,1+k3d,i) = nodetype(i)
              End Do
           end if validRegion
        end do iAxis
     end do jAxis
  end do kAxis


  !Ensure that the nodetype slot in surr_blks is set to -1 if the
  !block slot in surr_blks is negative.  This gives us the same result
  !as the orrery calculation.
  do i = 1, localNumBlocks
     do lk = 1, 1+2*k3d
        do lj = 1, 1+2*k2d
           do li = 1, 1+2*k1d
              if (surr_blks(1,li,lj,lk,i) < 0) then
                 surr_blks(3,li,lj,lk,i) = -1
              end if
           end do
        end do
     end do
  end do

End Subroutine mpi_amr_exchange_nodetype
