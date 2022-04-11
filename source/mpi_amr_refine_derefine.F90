#define FLASH_DEBUG_AMR
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_refine_derefine
!! NAME
!!
!!   amr_refine_derefine
!!
!! SYNOPSIS
!!
!!   call amr_refine_derefine(OPTIONAL, logical (IN) : force_rebalance)
!!
!!
!! ARGUMENTS
!!
!!   force_rebalance : optional, indicates whether the grid should be updated even if no
!!                     blocks are to be refined/derefined. If True, the distribution of
!!                     blocks over procs will be rebalanced by the work array. Default=False.
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
!!   mpi_morton
!!   timings
!!   io
!!   paramesh_interfaces
!!
!! CALLS
!!
!!   amr_check_refine
!!   amr_check_derefine
!!   amr_refine_blocks
!!   amr_derefine_blocks
!!   amr_morton_order
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit blocks marked for refinement have 
!!   been created, blocks marked for derefinement have been eliminated, and 
!!   the blocks have been reordered to achieve load balance.
!!
!! DESCRIPTION
!! 
!!  Subroutine to refine or derefine blocks.  This routine is called by the 
!!  user who sets the refine or derefine flags to be true or false.  These 
!!  flags are logical variables called 'refine' and 'derefine' and are stored 
!!  for each block.  If a block is marked for refinement, amr_refine_derefine 
!!  will create its new child blocks.  Also, tests will be executed to see if 
!!  any other blocks need to also refine (by calling amr_check_refine) to 
!!  ensure that the a jump in refinement of more than one level is not created.
!!
!!  If a block is marked for derefinement, amr_refine_derefine first checks 
!!  to make that the block can derefine and not create a jump in refinement of
!!  more than one level.  A check is also run to check that all the siblings of 
!!  the derefining block's siblings are also marked for derefinement.  If these 
!!  tests succeed, the block is removed from the list of blocks.
!!
!!  Once these operations are completed, the routine 'amr_morton_order' is 
!!  called and the tree data structure is reorganized to acheive load balance 
!!  using a morton space filling curve.  After this routine is called, the 
!!  routine 'amr_redist_blk' is called, which actually moves the block data 
!!  into the correct positions in the morton order list of blocks.
!!
!!  Finally, the routine 'amr_morton_process' is called.  This routine 
!!  computes the communications patterns needed for guardcell filling, 
!!  restriction, and prologation and stores them for later use. 
!!
!! AUTHORS
!!
!!  Kevin Olson (1997)
!!
!! HISTORY
!!
!!  force_rebalance and other tweaks   Tom Klosterman, Klaus Weide   2019-2022
!!***

#include "paramesh_preprocessor.fh"
#include "Simulation.h"
#include "constants.h"

      Subroutine amr_refine_derefine(force_rebalance)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_interfaces, Only : amr_check_refine,                & 
                                      amr_check_derefine,              & 
                                      amr_refine_blocks,               & 
                                      amr_derefine_blocks,             & 
                                      amr_morton_order
      Use paramesh_mpi_interfaces, Only : mpi_amr_singular_line
      use Logfile_interface, ONLY:  Logfile_stamp, Logfile_stampMessage 
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

#include "Flashx_mpi_implicitNone.fh"

      logical,intent(in),optional :: force_rebalance

!-----Local variables and arrays.
      Integer :: lnblocks2,tot_blocks,tot_blocksa,icontinue
      Integer :: icontinue_ref,icontinue_deref
      Integer :: icontinue2,max_blocks
      Integer :: min_blocks
      Integer :: nprocs,mype
      Integer :: i,l
      Integer :: lnblocks_old
#ifdef FLASH_DEBUG_AMR
      integer :: lnblocks_leaf
      integer min_blocks_leaf, max_blocks_leaf, tot_blocksa_leaf
#endif
      Integer :: istrategy
      Integer :: ierrorcode, ierr
      Logical :: l_move_solution
      Logical :: refinet(maxblocks_tr)
#ifdef FLASH_DEBUG_AMR
      character(len=32), dimension(3,2) :: block_buff
      character(len=32)                 :: int_to_str
#endif
      Logical,save :: first_Call = .True.
      logical :: rebalance

!-----Begin executable code.

      Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
      Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

!-----error trap for lrefine_max and lrefine_min
      If (lrefine_max < 1.Or.lrefine_max > 100) Then
        Write(*,*) 'PARAMESH error : lrefine_max has a bad value'      & 
                   ,lrefine_max
        Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      End If
      If (lrefine_min < 1.Or.lrefine_min > 100.Or.                     & 
         lrefine_min > lrefine_max) Then
        Write(*,*) 'PARAMESH error : lrefine_min or lrefine_max ',     & 
         'has a bad value : lrefine_min= ',lrefine_min, & 
         ' lrefine_max= ',lrefine_max
        Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      End If

!-----enforce refinement level limits
!-----first upper limit
      Where (lrefine==lrefine_max) refine = .False.
!-----Then lower limit
      Where (lrefine==lrefine_min) derefine = .False.
!-----finally force grid to refine toward base level if too coarse
      Where ( (lrefine<lrefine_min) .And.                              & 
              (nodetype==1) ) refine = .True.

      refinet(1:lnblocks) = refine(1:lnblocks)
      newchild(:) = .FALSE.

!-------CHECK derefinements and refinements

!-------test to see if any refinements have been requested.
      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l = 1,lnblocks
            If (nodetype(l) == 1.And.refine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          & 
                          MPI_MAX,amr_mpi_meshComm,ierr)

      icontinue_ref = icontinue2

      If (spherical_pm) Then
         istrategy = 0
         If (ndim == 3.And.lsingular_line) istrategy = 1
         Call mpi_amr_singular_line(istrategy,nprocs)
      End If

      Call amr_check_refine (nprocs,mype,icontinue_ref)

!-----test to see if any derefinements have been requested
      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l=1,lnblocks
            If (nodetype(l) == 1 .and. derefine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          & 
                          MPI_MAX,amr_mpi_meshComm,ierr)
      icontinue_deref = icontinue2

      If (icontinue_deref > 0) Call amr_check_derefine (mype)

!-----test to see if any derefinements have been requested which passed 
!-----the tests in amr_check_derefine
      If (icontinue_deref > 0) Then

      icontinue=0
      icontinue_deref = 0
      If (lnblocks > 0) Then
         Do l=1,lnblocks
            If (nodetype(l) == 1.And.derefine(l)) Then
               icontinue=1
               Exit
            End If
         End Do
      End If
      Call MPI_ALLREDUCE (icontinue,icontinue2,1,MPI_INTEGER,          &
                          MPI_MAX,amr_mpi_meshComm,ierr)
      icontinue_deref = icontinue2

      End If

      rebalance = .false.
      if(present(force_rebalance)) rebalance = force_rebalance
      If (icontinue_ref   == 0 .And.                                   &
          icontinue_deref == 0 .And.                                   &
          .Not.first_Call      .And.                                   &
          .Not.rebalance ) return

      
!-----NOW Actually refine and derefine the mesh
#ifdef FLASH_DEBUG_AMR
      call Logfile_stamp( 'initiating refinement', &
     &                   '[GRID amr_refine_derefine]')
#endif


#ifdef DEBUG_AMR
      if (mype.eq.0) then         
         open (unit=30,file=amr_log_file,status='unknown', & 
     &        position='append')
         write (30,*)' starting REFINE_BLOCKS'
         close(30)
         print *,' starting REFINE_BLOCKS'
      end if
#endif

      lnblocks_old = lnblocks
      If (icontinue_ref > 0) Call amr_refine_blocks (nprocs,mype)

      If (icontinue_deref > 0.Or.icontinue_ref > 0)                    & 
                     Call amr_derefine_blocks(lnblocks_old,mype)

      Call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER,           & 
                          MPI_SUM,amr_mpi_meshComm,ierr)

      if (myPE.eq.0) then
         !This logfile stamp is useful because it tells us the global number
         !of blocks requested by the refinement criteria.  New blocks are not
         !yet created!  We log this information now because it will be
         !lost if Paramesh aborts in amr_sort_by_work (as happens if the
         !block request cannot be satisfied).
         call Logfile_stamp(tot_blocks, '[GRID amr_refine_derefine]: '//&
              'redist. phase.  tot blks requested')         
         if (tot_blocks > MAXBLOCKS * nprocs) then
            !Maybe we could dump a checkpoint file here?
            call Logfile_stamp('[GRID amr_refine_derefine]: '//&
                 'Refinement will create too many blocks!', 'ERROR')
            print *, "Too many blocks!  Increase MAXBLOCKS or use more processors."
            !I won't abort yet because there may be some code in Paramesh that
            !I'm unaware of that allows us to survive this situation.
         end if
      end if

!-----set work values
#ifdef FLASH_DEBUG_AMR
      lnblocks_leaf = 0
#endif
      Do i = 1,lnblocks
         if (nodetype(i).eq.1) then
            if(gr_btSortByWork.AND. &
                 (.NOT.gr_btCustomWork .OR. (first_Call .AND. work_block(i) < gr_btWorkBoundsLeaf(LOW)))) &
              work_block(i) = gr_btWorkDefaultLeaf
#ifdef FLASH_DEBUG_AMR
            lnblocks_leaf = lnblocks_leaf + 1
#endif
         endif
         if (nodetype(i) >= 2) then
            if(gr_btSortByWork.AND. &
                 (.NOT.gr_btCustomWork .OR. (first_Call .AND. work_block(i) < gr_btWorkBoundsPar(LOW)))) &
              work_block(i) = gr_btWorkDefaultPar
         endif
      end do
      if(.not.gr_btSortByWork) work_block(1:lnblocks) = 1.

      l_move_solution = .True.
      Call amr_morton_order (lnblocks_old,nprocs,mype,                 & 
                             l_move_solution)

#ifdef DEBUG_AMR
      if (mype == 0) print *,' exited amr_morton_order'
#endif

      ! Reset work_block if not exchanging work
      if(.NOT.gr_btExchangeWork) then
        do i=1,lnblocks
          if(nodetype(i).eq.1) work_block(i) = gr_btWorkDefaultLeaf
          if(nodetype(i).gt.1) work_block(i) = gr_btWorkDefaultPar
        end do
      end if

#ifdef FLASH_DEBUG_AMR
      ! The counting above gives nonsensical results for min and max leaf block
      ! counts, apparently this has to be done AFTER calling amr_morton_order.
      ! - KW
      lnblocks_leaf = 0
      do i = 1,lnblocks
         if (nodetype(i).eq.1) then
            lnblocks_leaf = lnblocks_leaf + 1
         endif
      end do
#endif

!-----I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      Call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER,         & 
                          MPI_SUM,amr_mpi_meshComm,ierr)
      Call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,          & 
                          MPI_MAX,amr_mpi_meshComm,ierr)
      Call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,          & 
                          MPI_MIN,amr_mpi_meshComm,ierr)
#ifdef FLASH_DEBUG_AMR
      call MPI_ALLREDUCE (lnblocks_leaf,tot_blocksa_leaf,1,MPI_INTEGER,&
     &                    MPI_SUM,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (lnblocks_leaf,max_blocks_leaf,1,MPI_INTEGER,&
     &                    MPI_MAX,amr_mpi_meshComm,ierr)
      call MPI_ALLREDUCE (lnblocks_leaf,min_blocks_leaf,1,MPI_INTEGER,&
     &                    MPI_MIN,amr_mpi_meshComm,ierr)

#endif

#ifdef DEBUG_AMR
      if (mype == 0) print *,' exited MPI_ALLREDUCE'
#endif

      if (mype.eq.0) then
#ifdef DEBUG_AMR
         open (unit=30,file=amr_log_file,status='unknown', & 
     &        position='append')
         write (30,*) ' tot_blocks after ',tot_blocksa
         write (30,*) ' max_blocks 2',max_blocks
         write (30,*) ' min_blocks 2',min_blocks
         print *, ' tot_blocks after ',tot_blocksa
         print *, ' max_blocks 2',max_blocks
         print *, ' min_blocks 2',min_blocks
         close (30)
#endif
#ifdef FLASH_DEBUG_AMR
         
         
         ! write to string array and pass it to logfile;
         ! array dimensions are hard-coded in declaration
         
         write (block_buff(1, 1), '(a)') ' min blks '
         write (int_to_str, '(i7)') min_blocks
         write (block_buff(1, 2), '(a)') trim(adjustl(int_to_str))
         
         write (block_buff(2, 1), '(a)') 'max blks '
         write (int_to_str, '(i7)') max_blocks
         write (block_buff(2, 2), '(a)') trim(adjustl(int_to_str))
         
         write (block_buff(3, 1), '(a)') 'tot blks '
         write (int_to_str, '(i7)') tot_blocksa
         write (block_buff(3, 2), '(a)') trim(adjustl(int_to_str))


         call Logfile_stampMessage( '[GRID amr_refine_derefine]' //&
     &    trim(block_buff(1,1)) //' '// &
     &    trim(block_buff(1,2)) //'    '// trim(block_buff(2,1)) // &
     &    ' '// trim(block_buff(2,2)) //'    '// trim(block_buff(3,1)) &
     &    //' '// trim(block_buff(3,2)))


         write (block_buff(1, 1), '(a)') ' min leaf blks '
         write (int_to_str, '(i7)') min_blocks_leaf

         write (block_buff(1, 2), '(a)') trim(adjustl(int_to_str))
         
         write (block_buff(2, 1), '(a)') 'max leaf blks '
         write (int_to_str, '(i7)') max_blocks_leaf

         write (block_buff(2, 2), '(a)') trim(adjustl(int_to_str))
         
         write (block_buff(3, 1), '(a)') 'tot leaf blks '
         write (int_to_str, '(i7)') tot_blocksa_leaf

         write (block_buff(3, 2), '(a)') trim(adjustl(int_to_str))
         print *, 'refined: total leaf blocks = ', tot_blocksa_leaf
         print *, 'refined: total blocks = ', tot_blocksa


         call Logfile_stampMessage( '[GRID amr_refine_derefine]' //&
     &    trim(block_buff(1,1)) //' '// &
     &    trim(block_buff(1,2)) //'    '// trim(block_buff(2,1)) // &
     &    ' '// trim(block_buff(2,2)) //'    '// trim(block_buff(3,1)) &
     &    //' '// trim(block_buff(3,2)))
         
#endif
         
      end if      









#ifdef DEBUG_AMR
      if (tot_blocksa.ne.tot_blocks) then
         print *,' ERROR: tot_blocksa.ne.tot_blocks ', & 
     &        tot_blocksa,tot_blocks
         call MPI_ABORT(amr_mpi_meshComm,ierrorcode,ierr)
      end if
#endif

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown', & 
     &        position='append')
         write (30,*) ' done REFINE_DEREFINE '
         write (30,*) ' '
         print *, ' done REFINE_DEREFINE '
         print *,' '
         close(30)
      end if
#endif

#ifdef FLASH_DEBUG_AMR
      if (mype.eq.0) then
      call Logfile_stamp( 'refinement complete', &
     &                   '[GRID amr_refine_derefine]')
      end if
#endif

      Call amr_morton_process()

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file=amr_log_file,status='unknown', & 
     &        position='append')
         write (30,*) ' done MORTON_PROCESS '
         write (30,*) ' '
         print *, ' done MORTON_PROCESS '
         print *,' '
         close(30)
      end if
#endif

!-----Set up an array of cell sizes for each grid refinement level.
!-----These can be Used to minimize variation due to roundoff, but
!-----should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      If (ndim > 1)                                                    & 
        level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      If (ndim == 3)                                                   & 
        level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      Do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      End Do

!-----set grid modification flag
      grid_changed = 1
      grid_analysed_mpi = 1
#ifdef CALL_BOUNDARY_BLOCK_INFO
      Call mpi_amr_boundary_block_info(mype,nprocs)
#endif

      first_Call = .False.
      Return
      End Subroutine amr_refine_derefine

