!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_morton_order
!! NAME
!!
!!   amr_morton_order
!!
!! SYNOPSIS
!!
!!   call amr_morton_order(lnblocks_old, nprocs, mype, l_move_solution, reorder_grid)
!!
!!   call amr_morton_order(integer, integer, integer, logical, optional logical)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: lnblocks_old
!!     The number of block on the calling processor before any reordering is done.
!!
!!   integer, intent(in) :: nprocs
!!     The number for processors used.
!!
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   logical, intent(in) :: l_move_solution
!!     Logical switch.  If true then the data on each block (e.g. unk, facevarx,...etc.)
!!     is moved to the new location in the morton ordered list of blocks.  If it is
!!     false the data is not moved. 
!!
!!   logical, intent(in), optional :: reorder_grid
!!     An optional logical switch.  If this is set to true the blocks are only being
!!     reordered and are NOT assumed to be in any particular order.  If only reordering
!!     is being done, a specialized sorting routine is called.  This argument
!!     should only be set to true by the routine 'amr_reorder_grid' and user.
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
!!   io
!!   paramesh_interfaces
!!
!! CALLS
!! 
!!   amr_compute_morton
!!   amr_sort_morton
!!   amr_sort_morton_reorder_grid
!!   amr_sort_by_work
!!   amr_migrate_tree_data
!!   amr_redist_blk
!!    
!! RETURNS
!!
!!  Nothing returned.  Upon exit, the blocks are Morton ordered according to the
!!  morton space filling curve.
!!
!! DESCRIPTION
!!
!!  The is the wrapper routine which controls the process whereby the blocks are
!!  morton order and distributed across processors.  In general, the user should
!!  not need to call this subroutine.  Rather, it is called by either 
!!  'amr_refine_derefine', 'amr_reorder_grid', or 'amr_checkpoint_re' and the user
!!  should only call those routines to organize blocks for their particular 
!!  problem.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2001).
!!
!!***

#include "paramesh_preprocessor.fh"

      Subroutine amr_morton_order (lnblocks_old,nprocs,mype,           & 
                                   l_move_solution,                    & 
                                   reorder_grid)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_interfaces, only : amr_compute_morton,              & 
                                      amr_migrate_tree_data,           & 
                                      amr_redist_blk,                  & 
                                      amr_sort_by_work,                & 
                                      amr_sort_morton,                 & 
                                      amr_sort_morton_reorder_grid

      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Include statements.
      Include 'mpif.h'

!-----Input/Output statements.
      Integer, Intent(in)    :: lnblocks_old,nprocs,mype
      Logical, Intent(in)    :: l_move_solution
      Logical, Intent(in), Optional :: reorder_grid

!-----local variables and arrays.
      Integer :: new_loc(2,maxblocks_tr)
      Integer :: tot_blocks
      Integer :: ierr,lb
      Integer :: mort_no(6,2*maxblocks_tr)
      Integer :: ireduce_datain(1),ireduce_dataout(1)
      Logical :: ltemp
      Logical,Save :: first = .TRUE.

!-----Begin executable code.

!-----compute morton numbers for each cell
      Call amr_compute_morton (mort_no)
 
!-----Sort these morton numbers into order. The subroutine amr_sort_morton
!-----returns the array new_loc which gives the new locations that 
!-----the cells are to move to (local address is the first arguement
!-----and processor number is the second).

      new_loc(:,:) = -1
      ltemp = .False.
      If (present(reorder_grid)) ltemp = reorder_grid

      If (.Not.ltemp) Then
         Call amr_sort_morton (mort_no,new_loc,nprocs)
      Else
         Call amr_sort_morton_reorder_grid (mort_no,new_loc,nprocs)
      End If

      first = .FALSE.

!-----The following Call to sort_by_work attempts to realign the 
!-----sorted list returned by sort_morton such that the work load is 
!-----balanced across processors.
      ireduce_datain(1) = lnblocks
      Call MPI_ALLREDUCE (ireduce_datain,ireduce_dataout,         &
                          1,MPI_INTEGER,MPI_SUM,amr_mpi_meshComm,ierr)
      tot_blocks = ireduce_dataout(1)

      If (tot_blocks.gt.2*nprocs) Then
         Call amr_sort_by_work (new_loc,nprocs,mype)
      End If

      Call amr_migrate_tree_data (new_loc,nprocs,mype)

!-----move blocks of data to new locations
      If (l_move_solution)                                             & 
        Call amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)
      lnblocks = new_lnblocks

      Return
      End Subroutine amr_morton_order




