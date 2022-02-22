!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* bittree/amr_refine_derefine
!! NAME
!!
!!   amr_refine_derefine
!!
!! SYNOPSIS
!!
!!   call amr_refine_derefine(OPTIONAL, logical(IN) : force_rebalance)
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
!!   amr_morton_order_bittree
!!   amr_morton_process
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit blocks marked for refinement have 
!!   been created, blocks marked for derefinement have been eliminated, and 
!!   the blocks have been reordered to acheive load balance.
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
!!  In this Bittree version, Bittree is also updated throughout. After the
!!  refinement pattern is settled, the routine 'amr_morton_order_bittree' is
!!  called, which generates tree data for the new local list of blocks. Then 
!!  it calls 'amr_redist_blk' is called, which actually moves the block data 
!!  into the correct positions in the morton order list of blocks.
!!
!!  Finally, the routine 'amr_morton_process' is called.  This routine 
!!  computes the communications patterns needed for guardcell filling, 
!!  restriction, and prologation and stores them for later use. 
!!
!! AUTHORS
!!
!!  Kevin Olson (1997)
!!  Tom Klosterman (2019)
!!
!!***

#include "paramesh_preprocessor.fh"
#include "Simulation.h"

      Subroutine amr_refine_derefine(force_rebalance)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_interfaces, Only : amr_check_refine,                &
                                      amr_check_derefine
      Use paramesh_mpi_interfaces, Only : mpi_amr_singular_line
      use Logfile_interface, ONLY:  Logfile_stamp, Logfile_stampMessage 
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
      use bittree
      use Driver_interface, only: Driver_abort
      use iso_c_binding, only: c_int,c_bool
      
      Implicit None

!-----Include statements.
      Include 'mpif.h'

      logical,intent(in),optional :: force_rebalance

!-----Local variables and arrays.
      integer, allocatable :: blks_per_proc(:)
      Integer :: lnblocks2,tot_blocks,tot_blocksa,icontinue
      Integer :: icontinue_ref,icontinue_deref
      Integer :: icontinue2,max_blocks
      Integer :: min_blocks
      Integer :: nprocs,mype
      Integer :: i,j
#ifdef FLASH_DEBUG_AMR
      integer,allocatable :: leaves_per_proc(:),bitids(:)
      integer :: lnblocks_leaf
      integer :: min_blocks_leaf, max_blocks_leaf, tot_blocksa_leaf
      logical(c_bool) :: is_par
#endif
      Integer :: istrategy
      Integer :: ierrorcode, ierr
      Logical :: l_move_solution
      Logical :: refinet(maxblocks_tr)

      Integer :: lcoord(3)
      integer(c_int) :: bcount,lcount
      integer(c_int) :: ref_count, ref_count2
      logical(c_bool) :: ctrue = .TRUE.
      logical(c_bool) :: cfalse = .FALSE.

      character(len=32), dimension(3,2) :: block_buff
      character(len=32)                 :: int_to_str
      Logical,save :: first_Call = .True.
      logical :: rebalance

!-----Begin executable code.
      Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
      Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' entered amr_refine_derefine'
#endif

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

      newchild(:) = .FALSE.

!---------------------------------------------------------------------
!-----Initialize bittree refinement and mark leaves to be refined
      call bittree_refine_init()

      do i=1,lnblocks
        if(refine(i).AND.nodetype(i).eq.1) then
          call gr_getIntCoords(i,lcoord)
          call gr_btRefineMark(lrefine(i), lcoord)
        end if
      end do !l=1,blocks
      call bittree_refine_reduce(amr_mpi_meshComm)

!-----If any blocks are marked for refinement, call check_refine
      call bittree_delta_count(ref_count)
      if (ref_count>0) &
        call amr_check_refine(nprocs,mype,ref_count)
      call bittree_delta_count(ref_count)

!-----Now mark leaves to be derefined
      do i=1,lnblocks
        if(derefine(i).AND.nodetype(i).eq.1) then
          call gr_getIntCoords(i,lcoord)
          call gr_btDerefineMark(lrefine(i), lcoord )
        end if
      end do !l=1,blocks 
      call bittree_refine_reduce(amr_mpi_meshComm) 

!-----If delta count has increased, (thus some blocks marked for
!-----derefinement), call check_derefine
      call bittree_delta_count(ref_count2)
      If (ref_count2 > ref_count) &
        call amr_check_derefine(mype)
      call bittree_delta_count(ref_count2)

!-----Generate updated Bittree
      call bittree_refine_update()

!-----If no blocks marked for nodetype change, apply refinement and return
      rebalance = .false.
      if(present(force_rebalance)) rebalance = force_rebalance
      if (ref_count2.eq.0.AND..NOT.first_call.AND..NOT.rebalance) then
        call bittree_refine_apply()
        return
      end if
      first_call = .FALSE.

!-----Make logfile stamp with updated block count
      if (myPE.eq.0) then
         call bittree_block_count(ctrue,bcount)
         tot_blocks = bcount
         !This logfile stamp is useful because it tells us the global number
         !of blocks requested by the refinement criteria.  New blocks are not
         !yet created!  We log this information now because it will be
         !lost if Paramesh aborts in amr_sort_by_work (as happens if the
         !block request cannot be satisfied).
         call Logfile_stamp(tot_blocks,'[GRID amr_refine_derefine]: '//&
              'redist. phase.  tot blks requested')         
         if (tot_blocks > MAXBLOCKS * nprocs) then
            !Maybe we could dump a checkpoint file here?
            call Logfile_stamp('[GRID amr_refine_derefine]: '//&
                 'Refinement will create too many blocks!', 'ERROR')
            print *, "Too many blocks!  Increase MAXBLOCKS or &
                     &use more processors."
            !I won't abort yet because there may be some code in Paramesh that
            !I'm unaware of that allows us to survive this situation.
         end if
      end if

#ifdef FLASH_DEBUG_AMR
      call Logfile_stamp( 'initiating amr_morton_order_bittree', &
                          '[GRID amr_refine_derefine]')
#endif

!---------------------------------------------------------------------
!-----Call amr_morton_order bittree, which updates tree data (which replaces
!-----the functionality of amr_refine_blocks and amr_derefine_blocks), then 
!-----moves physical data
      call amr_morton_order_bittree(nprocs,mype,ref_count2)


!-----Finalize Bittree update
      call bittree_refine_apply()

#ifdef FLASH_DEBUG_AMR
      call amr_verify_bittree()
#endif


!-----On first proc, run some debugging and logstamps
      if(mype.eq.0) then
!-------Use localMortUB to compute total, max, and min blocks
        allocate(blks_per_proc(nprocs))
        blks_per_proc(1) = localMortUB(1)
        do i=2,nprocs
          blks_per_proc(i) = localMortUB(i) - localMortUB(i-1)
        end do
        tot_blocksa = sum(blks_per_proc)
        max_blocks = maxval(blks_per_proc)
        min_blocks = minval(blks_per_proc)
        
        print *, 'Done with refinement: total blocks = ', tot_blocksa

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

        deallocate(blks_per_proc)

#ifdef DEBUG_AMR
        open (unit=30,file=amr_log_file,status='unknown', & 
              position='append')
        write (30,*) ' tot_blocks after ',tot_blocksa
        write (30,*) ' max_blocks 2',max_blocks
        write (30,*) ' min_blocks 2',min_blocks
        print *, ' tot_blocks after ',tot_blocksa
        print *, ' max_blocks 2',max_blocks
        print *, ' min_blocks 2',min_blocks
        close (30)
#endif
         


#ifdef DEBUG_AMR
      if (tot_blocksa.ne.tot_blocks) then
         print *,' ERROR: tot_blocks different from &
         &localMortUB and Bittree: ',tot_blocksa,tot_blocks
         call MPI_ABORT(amr_mpi_meshComm,ierrorcode,ierr)
      end if
#endif

#ifdef DEBUG_AMR
        open (unit=30,file=amr_log_file,status='unknown', & 
     &        position='append')
        write (30,*) ' done REFINE_DEREFINE '
        write (30,*) ' '
        print *, ' done REFINE_DEREFINE '
        print *,' '
        close(30)
#endif

#ifdef FLASH_DEBUG_AMR
        call Logfile_stamp('refinement complete', &
                           '[GRID amr_refine_derefine]')
#endif
      
      end if !mype=0

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

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'pe ',mype,' exiting amr_refine_derefine'
#endif
      
      Return
      End Subroutine amr_refine_derefine
