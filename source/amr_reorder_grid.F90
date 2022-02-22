!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_reorder_grid
!! NAME
!!
!!   amr_reorder_grid
!!
!! SYNOPSIS
!!
!!   Call amr_reorder_grid()
!!
!! ARGUMENTS
!!
!!   No arguments.
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
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_morton_order
!!   amr_morton_process
!!
!! RETURNS
!!
!!   Nothing.
!!
!! DESCRIPTION
!!
!!   Calling this routine reorders the blocks for load balancing
!!   WITHOUT moving the block data (such as is in the 'unk' array).
!!
!! AUTHORS
!!
!!   Peter MacNeice and Kevin Olson.
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!#define DEBUG_AMR

      Subroutine amr_reorder_grid

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use paramesh_interfaces, only : amr_morton_order
      Use paramesh_mpi_interfaces, only : amr_morton_process
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Local variables and arrays.
      Integer :: nprocs,mype
      Integer :: lnblocks_old
      Integer :: ierr, i
      Integer :: lnblocks2, tot_blocks, tot_blocksa
      Integer :: max_blocks, min_blocks
      Logical :: l_move_solution, l_reorder_grid

!-----Begin executable code.
      lnblocks_old = lnblocks

      Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
      Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

      newchild(:) = .False.

      Call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER,           & 
                          MPI_SUM,amr_mpi_meshComm,ierr)

#ifdef DEBUG_AMR
      Write(*,*) 'refderef: tot_blocks ',tot_blocks,mype
      If (mype == 0) Then
         Print *,' tot_blocks before ',tot_blocks
      End If

!-----I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      Call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,          & 
                          MPI_MAX,amr_mpi_meshComm,ierr)
      Call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,          & 
                          MPI_MIN,amr_mpi_meshComm,ierr)

      If (mype == 0) Then
         Print *, ' max_blocks 1',max_blocks
         Print *, ' min_blocks 1',min_blocks
      End If
#endif

!-----set work per block values
      if(gr_btSortByWork.AND..NOT.gr_btCustomWork) then
        Do i = 1,lnblocks
         if (nodetype(i).eq.1) work_block(i) = gr_btWorkDefaultLeaf
         if (nodetype(i).gt.1) work_block(i) = gr_btWorkDefaultPar
        end do
      end if
      if(.not.gr_btSortByWork) work_block(1:lnblocks) = 1.

      l_move_solution = .False.
      l_reorder_grid = .True.
      Call amr_morton_order (lnblocks_old,nprocs,mype,                 & 
                             l_move_solution,                          & 
                             l_reorder_grid)

      ! Reset work_block if not exchanging work
      if(.NOT.gr_btExchangeWork) then
        do i=1,lnblocks
          if(nodetype(i).eq.1) work_block(i) = gr_btWorkDefaultLeaf
          if(nodetype(i).gt.1) work_block(i) = gr_btWorkDefaultPar
        end do
      end if
      
!-----Copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      Call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER,         & 
                          MPI_SUM,amr_mpi_meshComm,ierr)
      Call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,          & 
                          MPI_MAX,amr_mpi_meshComm,ierr)
      Call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,          & 
                          MPI_MIN,amr_mpi_meshComm,ierr)

      Call amr_morton_process()

!-----set grid modification flag
      grid_changed = 1
      grid_analysed_mpi = 1

      Return
      End Subroutine amr_reorder_grid
