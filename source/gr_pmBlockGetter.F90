!!  An object for "getting" the next block that is "ready for processing".
!
!!  "Getting": the primary method is call "get()". It returns the blockID
!!  of the next block that is ready, or -1 if there is no such block,
!!  that is, when the set of blocks that need to be returned has been exhausted.
!!  Of course, each blockID will be returned at most once.
!
!!  "Ready for processing":
!!  Processing can consist of applying operations that (1) use previous
!!  cell values (including those in guard cells) and/or (2) update cell
!!  values in block interiors; this is the primary use. Processing can
!!  also include updating guard cell values (example: Eos_guardCells),
!!  copying data, etc.
!!  Ready for processing means that data to be used as input is
!!  available (in the case of block cell data, stored in UNK),
!!  including guard cells where needed, AND that data locations to be
!!  updated by processing can be written to without data conflicts (in the
!!  case of block cell data, this means that later invocations of the
!!  get() method on the same getter object in the same cycle will not
!!  need to use data that are in danger of being overwritten).

!!  2021-01-09 K. Weide  Created for asynchronous domain data communications

#include "constants.h"
#include "Simulation.h"

module gr_pmBlockGetter
    use Grid_data, ONLY: nprocs => gr_meshNumProcs
    use Grid_data, ONLY: mype => gr_meshMe
    use gr_pmCommDataTypes, ONLY: gr_pmCommCtl_t
    use gr_parameshInterface, ONLY: gr_blockMatch
    use gr_flashHook_interfaces, ONLY: i27b ! a kind type parameter
    use tree, ONLY : lnblocks, surr_blks, laddress, strt_buffer
    use tree, ONLY : nodetype, child, parent, which_child, newchild, lrefine
    use physicaldata,        ONLY: advance_all_levels
    use mpi_morton,          ONLY: commatrix_recv, ladd_strt, ladd_end
    use paramesh_dimensions, ONLY : ndim, nguard,nguard_work

  implicit none

  private
  public :: gr_pmBlockGetterBuild, gr_pmBlockGetterDestroy

  integer,parameter:: i8=selected_int_kind(17)

  integer,parameter,public:: GRID_PAT_GC=10, GRID_PAT_PROLONG=20, GRID_PAT_FCORR=30, GRID_PAT_RESTRICT=40

  integer,parameter::S_initial=1,S_error=2,S_checkQueue=3,S_checkLocal=4,S_checkReady=5,S_checkRemote=6,S_wait=7,S_gotData=8,&
                     S_shutdown=9

  type, public :: gr_pmBlockGetter_t
     integer :: nodetype = LEAF
     integer :: lev      = INVALID_LEVEL

     logical :: receivedAll = .FALSE.
     integer :: receivedCount = 0

     integer :: iopt = 1
     logical :: lcc  = .TRUE.
     logical :: lfc = .FALSE.,lec=.FALSE.,lnc=.FALSE.
     integer :: buf_dim = 0
     integer :: nlayersx,nlayersy,nlayersz

     integer :: patternFamily = GRID_PAT_GC

     type(gr_pmCommCtl_t) :: commCtl

     integer :: state = S_initial
     integer :: next_state = S_initial
     logical :: doDrainQueue = .FALSE.
     logical :: inputsExhausted = .FALSE.
     logical :: defRemDepScan = .FALSE.
     logical :: ldependentsDontMatter = .FALSE.
     logical :: holdLdepends = .FALSE.
     integer :: tot_metaMissing = 0
     integer :: tot_ldependents_count = 0
     integer :: ldependents_count(MAXBLOCKS)
     integer :: numReleased           = 0
     integer :: oldNumReleased        = 0
     integer :: numUnkReleased        = 0
     integer :: oldNumUnkReleased     = 0
     integer :: remainingToDeliver    = 0
     integer :: Qcur(-1:MAXBLOCKS)
     integer :: Qcur_inP = 1, Qcur_outP = 1
     integer :: tot_no_blocks_to_be_received
     integer :: readyInUnk(MAXBLOCKS)
     integer(i8) :: localAvailFlagMaskGc(MAXBLOCKS)
     integer(i8) :: remoteAvailFlagMaskGc(MAXBLOCKS)
     integer(i8) :: needFlagMaskGc(MAXBLOCKS)
     integer(i8) :: remoteNeedFlagMaskGc(MAXBLOCKS)
     logical,ALLOCATABLE :: section_processed(:)
     logical,ALLOCATABLE :: psections_processed(:)
!!$     logical :: section_processed(tot_no_blocks_to_be_received)
!!$     logical :: psections_processed(nprocs)

   contains
     procedure :: get     => pmBlockGet
  end type gr_pmBlockGetter_t

contains
!#define DEBUG_STATEMACHINE
  subroutine gr_pmBlockGetterBuild(getter, nodetype, level, &
                                   patFam,          &
                                   iopt,            &
                                   lcc,lfc,lec,lnc, & 
                                   buf_dim,         &
                                   nlayersx,nlayersy,nlayersz)
#include "Flashx_mpi_implicitNone.fh"

        type(gr_pmBlockGetter_t), intent(OUT) :: getter
        integer,                  intent(IN)  :: nodetype
        integer,                  intent(IN)  :: level
        integer,                  intent(IN)  :: patFam
        logical,                  intent(in)  :: lcc,lfc,lec,lnc
        integer,                  intent(IN)  :: iopt
        integer,                  intent(IN)  :: buf_dim
        integer,                  intent(IN)  :: nlayersx,nlayersy,nlayersz

        getter%nodetype = nodetype
        getter%lev = level
        getter%patternFamily = patFam
        if (patFam == GRID_PAT_GC) then
           getter%ldependentsDontMatter = .FALSE.
        else
           getter%ldependentsDontMatter = .TRUE. !trying something...
        end if

        getter%iopt = iopt
        getter%lcc  = lcc
        getter%lfc  = lfc
        getter%lec  = lec
        getter%lnc  = lnc
        getter%buf_dim = buf_dim
        getter%nlayersx = nlayersx
        getter%nlayersy = nlayersy
        getter%nlayersz = nlayersz

        associate(commCtl => getter % commCtl)
          commCtl          % numReq = 0
          commCtl          % allReceivedCount = 0
          allocate(commCtl % recvrequest(nprocs))
          allocate(commCtl % recvstatus (MPI_STATUS_SIZE,nprocs)) ! We could use MPI_STATUSES_IGNORE instead.
          allocate(commCtl % receivedIndices(nprocs))
          commCtl          % numSReq = 0
          commCtl          % allSentCount = 0
          allocate(commCtl % sendRequest(nprocs))
          allocate(commCtl % sendStatus (MPI_STATUS_SIZE,nprocs)) ! We could use MPI_STATUSES_IGNORE instead.
        end associate

        getter % tot_no_blocks_to_be_received = sum(commatrix_recv(:))
        allocate(getter % section_processed(getter % tot_no_blocks_to_be_received))
        allocate(getter % psections_processed(0:nprocs-1))
#ifdef DEBUG_STATEMACHINE
         print*,'@',mype,' gr_pmBlockGetter is built; section_processed,psections_processed allocate with sizes',&
              size(getter % section_processed),size(getter % psections_processed), &
              ', ntype,level:',nodetype,level
#endif

  end subroutine gr_pmBlockGetterBuild

  subroutine gr_pmBlockGetterDestroy(getter)
#include "Flashx_mpi_implicitNone.fh"

        type(gr_pmBlockGetter_t), intent(INOUT) :: getter

        integer :: ierrorcode
        logical :: testFlag
        integer :: sentCount
        integer :: receivedCount
        integer :: i

        getter%buf_dim = 0

        associate(commCtl => getter % commCtl,      &
                  patFam  => getter % patternFamily &
                  )
          associate(ji               => commCtl % numSReq,          &
                    sendrequest      => commCtl % sendrequest,     &
                    sendstatus       => commCtl % sendstatus,      &
                    receivedInd      => commCtl % receivedIndices, &
                    allSentCount => commCtl % allSentCount)
            if (allocated(commCtl % sendRequest)) then
               do while (allSentCount < ji)
#define DEBUG_ASYNC_COMM
180               format('GetterDestroy(',I2,') @',I5,1x,A,I3,' of',I3,1x,A,' requests',A)
#ifdef DEBUG_ASYNC_COMM_MORE
                  print 180,patFam,mype,'Waiting for',(ji-allSentCount),ji,'send','...'
!                  print*,'GetterDestroy @',mype,' Waiting for', (ji-allSentCount),' of',ji,' send requests...'
#endif
                  if (allSentCount < ji) then
                     call MPI_Waitsome(ji,sendrequest,sentCount, &
                          receivedInd, sendstatus, ierrorcode)
                     allSentCount = allSentCount + sentCount
#ifdef DEBUG_GRIDCOMM
                     print 96,mype,sentCount,ji,(receivedInd(i),i=1,sentCount)
96                   format('On',I2,': WaitSome sent',I2,'/',I2,',',20(1x,I4,:,','))
#endif
                  end if
               end do
               deallocate(commCtl % sendRequest)
               deallocate(commCtl % sendStatus )
               If (Allocated(commCtl % sendBuf)) Deallocate(commCtl % sendBuf)
            end if
            ji = 0
            allSentCount = 0
          end associate

          if (allocated(commCtl % recvrequest)) then
             associate(receivedAll      => getter % receivedAll,        &
                        ij               => commCtl % numReq,          &
                        recvrequest      => commCtl % recvrequest,     &
                        recvstatus       => commCtl % recvstatus,      &
                        receivedInd      => commCtl % receivedIndices, &
                        allReceivedCount => commCtl % allReceivedCount &
                        )
               if (allReceivedCount < ij) then
#ifdef DEBUG_ASYNC_COMM
                  print 180,patFam,mype,'Testing for',(ij-allReceivedCount),ij,'receive','!'
!                  print*,'GetterDestroy @',mype,' Testing for', (ij-allReceivedCount),' of',ij,' receive requests!'
#endif
                  call MPI_TestAll(ij,recvrequest,testFlag,&
                       recvstatus, ierrorcode)
                  if (.NOT. testFlag) then
!!$                  do i=1,ij
!!$                     print*,'@',mype,' state',state,'Cancelling request',i,' ...'
!!$                     call MPI_Cancel(recvrequest(i),ierrorcode)
!!$                  end do
!!$                  print*,'@',mype,' state',state,'Cleanup: MPI_Cancel of',ij,' requests appears successful.'
#ifdef DEBUG_ASYNC_COMM
                     print 180,patFam,mype,'Waiting for',(ij-allReceivedCount),ij,'receive','...'
!                     print*,'GetterDestroy @',mype,' Waiting for', (ij-allReceivedCount),' of',ij,' receive requests...'
#endif
                     do while (allReceivedCount < ij)
                        call MPI_WaitSome(ij,recvrequest,receivedCount, &
                          receivedInd, recvstatus, ierrorcode)
                        allReceivedCount = allReceivedCount + receivedCount
#ifdef DEBUG_ASYNC_COMM
399                     format('On',I2,': WaitSome rec.procs',I2,'/',I2,',',20(1x,I4,:,','))
                        print 399,mype,receivedCount,ij,(recvstatus(MPI_SOURCE,receivedInd(i)),i=1,receivedCount)
#endif
                     end do
#ifdef DEBUG_GRIDCOMM
                     print*,'@',mype,' state',state,'Cleanup: MPI_WaitAll of',ij,' requests was successful.'
#endif
                  else
#ifdef DEBUG_ASYNC_COMM
499                  format(18x,'@ ',I4,' TestAll rec.procs',I2,'/',I2,',',20(1x,I4,':',I4,:,','))
                     receivedInd(1:ij-allReceivedCount) = &
                          pack([(i,i=1,ij)], &
                               (recvstatus(MPI_TAG,1:ij).NE.MPI_ANY_TAG))
                     print 499,mype,ij-allReceivedCount,ij,(                      receivedInd(i),&
                                                            recvstatus(MPI_SOURCE,receivedInd(i)),&
                                                               i=1,ij-allReceivedCount)
#endif
#ifdef DEBUG_GRIDCOMM
                     print*,'@',mype,' state',state,'Cleanup: MPI_TestAll of',ij,' requests was successful.'
#endif
                     allReceivedCount = ij
                  end if
               end if
             end associate
             deallocate(commCtl % recvrequest)
             deallocate(commCtl % recvstatus )
             deallocate(commCtl % receivedIndices)
#ifdef DEBUG_STATEMACHINE
          else
             print*,'@',mype,' gr_pmBlockGetter destroy: commCtl was not properly initialized?'
#endif
          end if
          commCtl          % numReq = 0
          commCtl          % allReceivedCount = 0
        end associate
        if (allocated(getter % section_processed)) then
           deallocate(getter % section_processed)
           deallocate(getter % psections_processed)
#ifdef DEBUG_STATEMACHINE
        else
           print*,'@',mype,' gr_pmBlockGetter destroy: getter was not properly initialized!'
#endif
        end if
  end subroutine gr_pmBlockGetterDestroy

!!  Some remarks / ideas / plans:
!!  Re readyInUnk:
!!    The set of values should be extented, for tracking the progress
!!    of each local block:
!!      0 - initial (if block is wanted?)
!!      1 - block (with all its grid data) is ready in perm storage.
!!      2 - block has been queued.
!!          (Also, unwanted blocks could start with this?)
!!      3 - block has been delivered. (Currently only useful
!!          for debugging.)

  function pmBlockGet_(this, commCtl) result(blockID)
    use mpi_morton,           ONLY: temprecv_buf

#include "Flashx_mpi_implicitNone.fh"

    integer ::                     blockID
    type(gr_pmBlockGetter_t),intent(inout) :: this
    type(gr_pmCommCtl_t), intent(INOUT)     :: commCtl

    ASYNCHRONOUS :: temprecv_buf
    integer :: ierrorcode
!!$    integer :: receivedCount
    integer :: i,j,invDir,k
    integer :: isrc,clb
    integer :: prevNumReleased
    integer :: queue_length
    integer :: blk, destBlk, destProc, ppBlk, ppProc
    integer :: neighDir, neighBlk, neighProc, NDirs, NDirsReverse
    integer :: misNei, parBlk, pneighDir, pneighBlk, pneighProc
    integer :: childBlk, childProc
    integer :: ia,ib,ja,jb,ka,kb
    integer :: iD,jD,kD, iDa,iDb, jDa,jDb, kDa,kDb
    integer :: ioff,joff,koff,nguarda
    logical :: have_parent_surr, testFlag
    logical :: progress, skip, doGotoTop
    integer :: sentCount

    logical :: get_data_quickly = .FALSE. !to become runtime paramter
    logical :: get_data_very_quickly = .FALSE. !to become runtime paramter
!!$    logical :: push_data_out_quickly = .FALSE. !to become runtime paramter
    logical :: push_data_out_very_quickly = .TRUE. !to become runtime paramter
    integer,parameter :: QUEUE_PACKET_SIZE = 1 !to become runtime paramter?

    associate(receivedAll      => this % receivedAll,        &
              ij               => commCtl % numReq,          &
              recvrequest      => commCtl % recvrequest,     &
              recvstatus       => commCtl % recvstatus,      &
              receivedInd      => commCtl % receivedIndices, &
              allReceivedCount => commCtl % allReceivedCount,&
              pfam             => this % patternFamily,      &
              iopt             => this % iopt,               &
              lcc              => this % lcc,                &
              lfc              => this % lfc,                &
              lec              => this % lec,                &
              lnc              => this % lnc,                &
              buf_dim          => this % buf_dim,            &
              nlayersx         => this % nlayersx,           &
              nlayersy         => this % nlayersy,           &
              nlayersz         => this % nlayersz,           &
              ntype            => this % nodetype,           &
              lev              => this % lev,                &
              state            => this % state,              &
              next_state       => this % next_state,         &
              doDrainQueue     => this % doDrainQueue,        &
              inputsExhausted  => this % inputsExhausted,     &
              deferRemoteDependencyScan => this % defRemDepScan, &
              holdLdepends     => this % holdLdepends,        &
              tot_metaMissing  => this % tot_metaMissing,     &
              tot_ldependents_count => this % tot_ldependents_count,    & !may nor be needed
              ldependents_count     => this % ldependents_count,    &
              numReleased           => this % numReleased,              &
              oldNumReleased           => this % oldNumReleased,              &
              numUnkReleased           => this % numUnkReleased,              &
              oldNumUnkReleased           => this % oldNumUnkReleased,              &
              remainingToDeliver    => this % remainingToDeliver,       &
              readyInUnk            => this % readyInUnk,    &
              Qcur                  => this % Qcur,          &
              Qcur_inP              => this % Qcur_inP,      &
              Qcur_outP             => this % Qcur_outP,     &
              localAvailFlagMaskGc  => this % localAvailFlagMaskGc,  &
              remoteAvailFlagMaskGc => this % remoteAvailFlagMaskGc, &
              needFlagMaskGc        => this % needFlagMaskGc,        &
              remoteNeedFlagMaskGc  => this % remoteNeedFlagMaskGc,  & ! maybe UNNEEDED, use needFlagMaskGc?
              section_processed     => this % section_processed,     &
              psections_processed   => this % psections_processed,   &
              receivedCount         => this % receivedCount          &
              )

      select case (pFam)
      case(GRID_PAT_GC)
         NDirs = 3**NDIM
         NDirsReverse = 3**NDIM
      case(GRID_PAT_PROLONG)
         NDirs = 3**NDIM
         NDirsReverse = 3**NDIM
      case(GRID_PAT_FCORR)
         NDirs = 2 * NDIM
         NDirsReverse = 2 * NDIM
      case default !(GRID_PAT_RESTRICT)
         NDirs = 2**NDIM
         NDirsReverse = 1
      end select

1     continue
#ifdef DEBUG_STATEMACHINE
    print*,'@',mype,' state',state,' next_state',next_state,' at the Top!'
#endif

      do 

#ifdef DEBUG_STATEMACHINE
         print*,'@',mype,' state',state,' ->',next_state
#endif

         state = next_state

         select case (state)

         case(S_initial)
            nguarda = max(nguard,nguard_work)
            readyInUnk(:) = 0   !will use 1 for 'ready', 2 for 'ready and queued'
            Qcur_inP = 1
            Qcur_outP = 1
            ldependents_count(:) = 0
            section_processed(:) = .FALSE.
            psections_processed(:) = .FALSE.
            ! etc. etc. ?
            ! (first_call):
            !# Set per-gcregion bit flags for (a) "gcregion needed"; initialize another field
            do blk = 1,lnblocks
               localAvailFlagMaskGc(blk)  = int(z'00000000') ! local availability of Gc region data
               remoteAvailFlagMaskGc(blk) = int(z'00000000') ! initialize
               remoteNeedFlagMaskGc(blk)  = int(z'00000000') ! initialize
               needFlagMaskGc(blk)  = int(z'00000000') ! initialize
            end do

            !# Set per-gcregion bit flags for (b) "gcregion locally available"
            do blk = 1,lnblocks
               if (isWanted(blk,ntype,lev)) then      ! if this block is among the ones requested (by nodetype, etc.)
                  remainingToDeliver = remainingToDeliver + 1
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' Block',blk,' is wanted, type/lev:',&
                      nodetype(blk),lrefine(blk),', initially',ldependents_count(blk),' dependents.'
#endif
                  select case (pFam)
                  case(GRID_PAT_GC)
#if N_DIM == 1
                     needFlagMaskGc(blk) = int(o'0000000005') ! 3D # from comm pattern?
#endif
#if N_DIM == 2
                     needFlagMaskGc(blk) = int(o'0000000757') ! 2D # from comm pattern?
#endif
#if N_DIM == 3
                     needFlagMaskGc(blk) = int(o'0777757777') ! 3D # from comm pattern?
#endif
                  case(GRID_PAT_RESTRICT)
#if N_DIM == 1
                     needFlagMaskGc(blk) = int(o'00000003') ! 3D # from comm pattern?
#endif
#if N_DIM == 2
                     needFlagMaskGc(blk) = int(z'0000000f') ! 2D # from comm pattern?
#endif
#if N_DIM == 3
                     needFlagMaskGc(blk) = int(z'000000ff') ! 3D # from comm pattern?
#endif
                  end select
                  misNei = 0
                  ka = 1; kb = 1+2*K3D
                  ja = 1; jb = 1+2*K2D
                  ia = 1; ib = 3
                  have_parent_surr = (parent(2,blk) == mype)
!!#ifdef PM_OPTIMIZE_MORTONBND_FETCHLIST
                  if (have_parent_surr) then
                     ASSOCIATE(psurr_blk => surr_blks(:,:,:,:,parent(1,blk)))
                     ! Do this optimization only if parent does not touch a domain boundary
                     ! anywhere (otherwise boundary condition can be called with invalid
                     ! input data). Cf. mpi_morton_bnd.
                       if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
#if NDIM >= 3
                          if (nguarda .LE. 3NZB/2) then
                             koff = mod((which_child(blk)-1)/4,2)
                             ka = 1+koff; kb = 1+K3D+koff
                          end if
#endif
#if NDIM >= 2
                          if (nguarda .LE. NYB/2) then
                             joff = mod((which_child(blk)-1)/2,2)
                             ja = 1+joff; jb = 1+K2D+joff
                          end if
#endif
                          if (nguarda .LE. NXB/2) then
                             ioff = mod(which_child(blk)-1,2)
                             ia = 1+ioff; ib = 2+ioff
                          end if
                       end if
                     end ASSOCIATE
                  end if
!!#endif
                  iDa = ia; jDa = ja; kDa = ka
                  iDb = ib; jDb = jb; kDb = kb
                  do neighDir = 1,dependencyDirs(blk,pFam)!NDirs
                     if (neighDir > NDirs) EXIT
                     call neighInfo(blk,neighDir,neighBlk,neighProc,pFam)
#ifdef DEBUG_STATEMACHINE
                     print*,'@',mype,' state',state,' Block',blk, &
                           ', neighDir/neighBlk/neighProc:',&
                              neighDir,neighBlk,neighProc
#endif
                     if (neighProc == mype) then ! neighBlk is local:
                        if (neighBlk .NE. blk) then  ! condition for self-dependence at periodic bdry
                           ldependents_count(neighblk) = ldependents_count(neighblk) + 1
                           tot_ldependents_count = tot_ldependents_count + 1
                        end if
                        call bset(localAvailFlagMaskGc(blk), neighDir) ! set neighDir bit in localAvail bitmask
                     else if (neighProc .GE. 0) then
                        call bset(remoteNeedFlagMaskGc(blk), neighDir) ! set neighDir bit in remoteNeed bitmask
                     else if (neighBlk .LE. PARAMESH_PHYSICAL_BOUNDARY) then
                        call bclr(needFlagMaskGc(blk), neighDir) ! unset neighDir bit in need bitmask
                     else if ((neighBlk == -1 .OR. neighProc == -1) .AND. &
                              nodetype(blk) == 1) then !should only happen for surr_blks?
                        misNei = misNei + 1                        ! a missing surr_blks neighbor!
#ifdef DEBUG_STATEMACHINE
                        print*,'@',mype,' state',state,' Block',blk,' has',misNei,' missing surr_blks neighbors!'
                        print 30, mype,state,blk,needFlagMaskGc(blk),localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
                        if (pFam == GRID_PAT_GC) then
                           call bclr(needFlagMaskGc(blk), neighDir) ! unset neighDir bit in need bitmask

                           if (have_parent_surr) then
                              iD  = mod( neighDir-1,   3)+1
                              jD  = mod((neighDir-1)/3,3)+1
                              kD  =     (neighDir-1)/9   +1
                              ASSOCIATE(psurr_blk => surr_blks(:,:,:,:,parent(1,blk)))
                                if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
                                   kDa = min(kDa,kD-1) ; kDb = max(kDb,kD+1)
                                   jDa = min(jDa,jD-1) ; jDb = max(jDb,jD+1)
                                   iDa = min(iDa,iD-1) ; iDb = max(iDb,iD+1)
                                end if
                              end ASSOCIATE
                           end if  !(parent(2,blk) == mype)
                        end if     !if (pFam
                     end if     ! if (neighProc == mype
                  end do        !do neighDir

                  if (misNei > 0) then

                     if (have_parent_surr .AND. pFam == GRID_PAT_GC) then

                        ASSOCIATE(psurr_blk => surr_blks(:,:,:,:,parent(1,blk)))
#ifdef PM_OPTIMIZE_MORTONBND_FETCHLIST
                                ! Do this optimization only if parent does not touch a domain boundary
                                ! anywhere (otherwise boundary condition can be called with invalid
                                ! input data). Cf. mpi_morton_bnd.
                          if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
                             nguarda = max(nguard,nguard_work)
#if NDIM >= 3
                             if (nguarda .LE. NZB/2) then
                                ka = max(ka,kDa); kb = min(kb,kDb)
                             end if
#endif
#if NDIM >= 2
                             if (nguarda .LE. NYB/2) then
                                ja = max(ja,jDa); jb = min(jb,jDb)
                             end if
#endif
                             if (nguarda .LE. NXB/2) then
                                ia = max(ia,iDa); ib = min(ib,iDb)
                             end if
                          end if
#endif
!!$                                print*,'@',mype,' state',state,' Block',blk,'in ASSOCIATE...'
                          do k = ka,kb
                             do j = ja,jb
                                do i = ia,ib
                                   pneighBlk  = psurr_blk(1,i,j,k)
                                   pneighProc = psurr_blk(2,i,j,k)
#ifdef DEBUG_STATEMACHINE
                                   print*,'@',mype,' state',state,' Block',blk,' i,j,k:',i,j,k,', pneighBlk/pneighProc:',&
                                                pneighBlk,pneighProc
#endif
                                   if (pNeighBlk > 0 .AND. pneighProc > -1) then !.AND. (psurr_blk(3,ia,ja,ka) < 3) ?
                                      pneighDir = 3**NDIM + i + K2D*(j-1)*3 + K3D*(k-1)*9
#ifdef DEBUG_STATEMACHINE
                                      print*,'@',mype,' state',state,' Block',blk,' i,j,k:',i,j,k, &
                                                 ', pneighDir/pneighBlk/pneighProc:',&
                                                    pneighDir,pneighBlk,pneighProc
#endif
                                      if (pneighProc == mype) then ! pneighBlk is local:
                                         ldependents_count(pneighBlk) = ldependents_count(pneighBlk) + 1
                                         tot_ldependents_count = tot_ldependents_count + 1
                                         call bset(localAvailFlagMaskGc(blk), pneighDir) ! set pneighDir bit in localAvail bitmask
                                      end if
                                      call bset(needFlagMaskGc(blk), pneighDir) ! set neighDir bit in remoteNeed bitmask
                                      if (pneighProc .NE. mype) then
                                         call bset(remoteNeedFlagMaskGc(blk), pneighDir) ! set neighDir bit in remoteNeed bitmask
                                      end if
                                   end if !(pNeighBlk > 0
                                end do
                             end do
                          end do
! DEV: No special handling of (lsingular_line), (cf. mpi_morton_bnd)
                        end ASSOCIATE
                     else
#define PARENTMETA_BIT (2*3**NDIM+1)
                        call bset(needFlagMaskGc(blk), PARENTMETA_BIT)
                        tot_metaMissing = tot_metaMissing + 1
#ifdef DEBUG_STATEMACHINE
                        print*,'@',mype,' state',state,' Block',blk,' - TRANSITION holdLdepends to T from',holdLdepends,&
                                ', missing meta',tot_metaMissing
#endif
                        holdLdepends = .TRUE.
                     end if  !(parent(2,blk) == mype)
                  end if     !if (misNei > 0)
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' Block',blk,' is wanted, type/lev:',&
                      nodetype(blk),lrefine(blk),', now',ldependents_count(blk),' dependents.'
                  print 30, mype,state,blk,needFlagMaskGc(blk),localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
               end if           !if (isWanted
            end do              !do blk

            next_state = S_checkLocal ! or S_wait or S_checkRemote if(get_data_very_quickly)
            if (remainingToDeliver == 0) next_state = S_shutdown
#ifdef DEBUG_STATEMACHINE
            print*,'@',mype,' state',state,' Queue has',Qcur_inP-Qcur_outP,' entries:',&
                 Qcur(Qcur_outP:Qcur_inP-1),' (',Qcur_outP,' :',Qcur_inP-1,')'
#endif

         case(S_checkQueue)

            queue_length = Qcur_inP - Qcur_outP
#ifdef DEBUG_STATEMACHINE
            if (queue_length <= 10) then
               print*,'@',mype,' state',state,' Queue has',queue_length,' entries:',&
                 Qcur(Qcur_outP:Qcur_inP-1),' (',Qcur_outP,' :',Qcur_inP-1,'), remain',remainingToDeliver
            else
               print*,'@',mype,' state',state,' Queue has',queue_length,' entries:',&
                 Qcur(Qcur_outP:Qcur_outP+3),' ...',Qcur(Qcur_inP-4:Qcur_inP-1),' (',Qcur_outP,' :',Qcur_inP-1,')'
            end if
#endif
            if(Qcur_outP < Qcur_inP .AND. (inputsExhausted .OR. doDrainQueue .OR. queue_length >= QUEUE_PACKET_SIZE)) then
!            if(Q_for_current_step not empty .AND. (inputsExhausted .OR. doDrainQueue .OR. queue_length >= QUEUE_PACKET_SIZE)) then
               next_state = S_checkQueue ! again, unneccesary assignment
               if (queue_length >= QUEUE_PACKET_SIZE) then
                  doDrainQueue = .TRUE.
               else if (queue_length <= 1) then
                  doDrainQueue = .FALSE.
               end if
               remainingToDeliver = remainingToDeliver - 1
               if (remainingToDeliver == 0) next_state = S_shutdown

               blockID = Qcur(Qcur_outP)
               Qcur_outP = Qcur_outP + 1
#ifdef DEBUG_STATEMACHINE
               print*,'@',mype,' state',state,' RETURN block',blockID
#endif
               readyInUnk(blockID) = 3
               RETURN
            else
               if (numUnkReleased .NE. oldNumUnkReleased .AND. .NOT. get_data_very_quickly) then
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' oldNumUnkReleased,numUnkReleased,oldNumReleased,numReleased=',&
                             oldNumUnkReleased,numUnkReleased,oldNumReleased,numReleased
#endif
                  next_state = S_checkReady
               else if (numReleased .NE. oldNumReleased .AND. .NOT. get_data_quickly) then
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' oldNumUnkReleased,numUnkReleased,oldNumReleased,numReleased=',&
                             oldNumUnkReleased,numUnkReleased,oldNumReleased,numReleased
#endif
                  next_state = S_checkLocal
               else if (.NOT. inputsExhausted) then
                  next_state = S_wait
               else
                  next_state = S_checkRemote
               end if
            end if


         case(S_checkReady)

            next_state = S_wait
            do blk = 1,lnblocks
               if (isWanted(blk,ntype,lev)) then
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' Block',blk,' is wanted, type/lev:',&
                      nodetype(blk),lrefine(blk),', has',ldependents_count(blk),' dependents.'
#endif
!                  if (readyInUnk(blk) == 1 .AND. ldependents_count(blk) == 0) then
                  if (readyInUnk(blk) == 1 .AND. ldependents_ok(blk)) then
                     Qcur(Qcur_inP) = blk
                     Qcur_inP = Qcur_inp + 1
                     readyInUnk(blk) = 2
#ifdef DEBUG_STATEMACHINE
                     print*,'@',mype,' state',state,' Block',blk,' ready in Unk, queue has',Qcur_inP-Qcur_outP,' entries:',&
                          Qcur(Qcur_outP:Qcur_inP-1),' (',Qcur_outP,' :',Qcur_inP-1,')'
#endif
                     next_state = S_checkQueue
                  end if
               end if
            end do
            oldNumUnkReleased = numUnkReleased


         case(S_checkLocal)

            next_state = S_wait
            prevNumReleased = numReleased
            do blk = 1,lnblocks
               if (isWanted(blk,ntype,lev) .AND. readyInUnk(blk)==0) then
#ifdef DEBUG_STATEMACHINE
#   if N_DIM == 3
#      define OCTFLAG_FORMAT o19.19
#   else
#      define OCTFLAG_FORMAT o10.10
#   endif
                  print*,'@',mype,' state',state,' Block',blk,' is wanted, type/lev:',&
                      nodetype(blk),lrefine(blk),', not readyInUnk,',ldependents_count(blk),' dependents.'
30                format(1x,'@',I12,'  state',I12,'  Block',I12,&
                       ' ... needed: o',OCTFLAG_FORMAT,', local/remote avail: o',OCTFLAG_FORMAT,',o',OCTFLAG_FORMAT)
35                format(1x,'+',I12,'  state',I12,'  Block',I12,&
                       ' ... needed: o',OCTFLAG_FORMAT,', local/remote avail: o',OCTFLAG_FORMAT,',o',OCTFLAG_FORMAT)
                  print 30, mype,state,blk,needFlagMaskGc(blk),localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
                  if (neededGcsAreAvailable(this,blk)) then
#ifdef DEBUG_STATEMACHINE
                     print*,'@',mype,' state',state,' Block',blk,' ... needed guardcells are available!'
#endif
                     do neighDir = 1,dependencyDirs(blk,pFam)!NDirs
                        call neighInfo(blk,neighDir,neighBlk,neighProc,pFam,skip)
                        if (neighProc == mype .AND. neighBlk .NE. blk .AND. .NOT. skip) then
                           ldependents_count(neighblk) = ldependents_count(neighblk) - 1
#ifdef DEBUG_STATEMACHINE
                           print*,'@',mype,' state',state,' Block',neighBlk,' (dir',neighDir,')', &
                                                  ' had ldependents decreased to',ldependents_count(neighBlk),&
                                                  ' because we release local blk',blk
#endif
                           if (ldependents_ok(neighblk)) then
                              if (readyInUnk(neighblk) > 0) then
                                 numUnkReleased = numUnkReleased + 1
                              else
                                 numReleased    = numReleased + 1
                              end if
                           end if
                           tot_ldependents_count = tot_ldependents_count - 1
                           if (.NOT. get_data_quickly) next_state = S_checkLocal ! again
                        end if
                     end do
                     call update_perm_griddata(this,blk) ! also applies domain boundary conditions
                     
                     readyInUnk(blk) = 1
#ifdef DEBUG_STATEMACHINE
                     print*,'@',mype,' state',state,' Block',blk,' became ready in Unk.'
#endif
!                     if (ldependents_count(blk) == 0) then
                     if (ldependents_ok(blk)) then
                        Qcur(Qcur_inP) = blk
                        Qcur_inP = Qcur_inp + 1
                        readyInUnk(blk) = 2
                        next_state = S_checkQueue
                        if (push_data_out_very_quickly) then
#ifdef DEBUG_STATEMACHINE
                           print*,'@',mype,' state',state,' CYCLE to Top'
#endif
                           goto 1 
                        end if
                     else
                        next_state = S_checkReady
                     end if
                  end if        !if (neededGcsAreAvailable
               end if
            end do

            if (numReleased == prevNumReleased) oldNumReleased = numReleased
            if (numReleased == oldNumReleased) then
               next_state = S_wait
            end if

         case(S_wait)

            ! We must only go into this state if we are sure that there is not anything else to do.
            next_state = S_gotData

            if (receivedAll) then
               next_state = S_error
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' - receivedAll is True!'
#endif
            else if ((allReceivedCount .GT. ij) .OR. &
                     (allReceivedCount == ij .AND. receivedCount == 0)) then
               next_state = S_error
#ifdef DEBUG_STATEMACHINE
               print*,'@',mype,' state',state,' - allReceivedCount is',allReceivedCount,' of',ij,', ',&
                    'receivedCount',receivedCount,', inputsExh',inputsExhausted
#endif
               if (.NOT. inputsExhausted) then
                  inputsExhausted = .TRUE.
                  next_state = S_checkQueue
               end if
            else if (receivedCount == 0) then

!#define DEBUG_GRIDCOMM
#define HAVE_WAIT
               ! We can wait without having to cycle
#ifdef HAVE_WAIT
               ! Wait for new data from any remote
               call MPI_Waitsome(ij,recvrequest,receivedCount, &
                       receivedInd, recvstatus, ierrorcode)
               allReceivedCount = allReceivedCount + receivedCount
#ifdef DEBUG_GRIDCOMM
198             format('Gtr S_wait',I2,': WaitSome received',I2,'/',I2,',',20(1x,I4,:,','))
               print 198,mype,receivedCount,ij,(receivedInd(i),i=1,receivedCount)
199             format('Gtr S_wait',I2,': WaitSome rec.procs',I2,'/',I2,',',20(1x,I5,:,','))
               print 199,mype,receivedCount,ij,(recvstatus(MPI_SOURCE,i),i=1,receivedCount)
#endif
               do i=1,receivedCount
                  isrc = recvstatus(MPI_SOURCE,i)
                  call pmMpiUnpackBlksFromProc(isrc, iopt, &
                          lcc,lfc,lec,lnc, & 
                          buf_dim,temprecv_buf, & 
                          nlayersx,nlayersy,nlayersz)
               end do

               ! receivedAll = .TRUE.
               ! this%curBlk = 0
#else
!            poll for data
!            if (no new data from any remote) then
!               next_state = S_wait ! again
!               sleep some time or yield to another thread
!            end if
               call MPI_Testsome(ij,recvrequest,receivedCount, &
                                receivedInd, recvstatus, ierrorcode)
#ifdef DEBUG_GRIDCOMM
298             format('Gtr S_wait',I2,': TestSome received',I2,'/',I2,',',20(1x,I4,:,','))
               print 298,mype,receivedCount,ij,(receivedInd(i),i=1,receivedCount)
299             format('Gtr S_wait',I2,': TestSome rec.procs',I2,'/',I2,',',20(1x,I5,:,','))
               print 299,mype,receivedCount,ij,(recvstatus(MPI_SOURCE,i),i=1,receivedCount)
#endif
               if (receivedCount == 0) then
                  next_state = S_wait
                     ! ... we could sleep here for some time ...
               else if (receivedCount == MPI_UNDEFINED) then
                  next_state = S_error
               else
                  allReceivedCount = allReceivedCount + receivedCount
                  do i=1,receivedCount
                     isrc = recvstatus(MPI_SOURCE,i)
                     call pmMpiUnpackBlksFromProc(isrc, iopt, &
                          lcc,lfc,lec,lnc, & 
                          buf_dim,temprecv_buf, & 
                          nlayersx,nlayersy,nlayersz)
                  end do

               end if
               ! receivedAll = .TRUE.
               ! this%curBlk = 0
#endif
            end if
#ifdef DEBUG_STATEMACHINE
            print*,'E',mype,' state',state,' - allReceivedCount is',allReceivedCount,' of',ij,', ',&
                    'receivedCount',receivedCount,', inputsExh',inputsExhausted
#endif
            ! Maybe invent some better way to pass the location of new data (if available) to the next state.
            ! For now, we just pass on directly to the state S_checkRemote whatever has been set by the
            ! MPI_Waitsome/MPI_Testsome calls above:
            ! (receivedCount, recvstatus).)

         case(S_gotData)

            next_state = S_checkRemote


         case(S_checkRemote)

            next_state = S_wait

            prevNumReleased = numReleased
            progress = .FALSE.
#ifdef DEBUG_STATEMACHINE
            print*,'@',mype,' state',state,' - allReceivedCount is',allReceivedCount,' of',ij,', ',&
                    'receivedCount',receivedCount,', inputsExh',inputsExhausted
#endif
            ! Start a scan of the receive buffer
            ! for every_section(recv_buffer_for_current_step): ! each section represents data from a pp block
            do i=1,receivedCount
               isrc = recvstatus(MPI_SOURCE,i)
               if (.NOT. psections_processed(isrc)) then
                  ! Iterate over the pp blocks received from isrc.
                  do clb=ladd_strt(isrc),ladd_end(isrc) !'cached local block', an index in the high region
                     if (.NOT. section_processed(clb-strt_buffer+1)) then

                        ! get ppBlk from section_data - (ppBlk,isrc) identify the remote block.
                        ppBlk = laddress(1,clb)
#ifdef DEBUG_STATEMACHINE
                        print*,'@',mype,' state',state,'clb',clb,' now holds (ppBlk,isrc)',ppBlk,isrc,&
                                 ' nodetype',nodetype(clb),' meta,def:',tot_metaMissing,deferRemoteDependencyScan
#endif
                        if (pFam == GRID_PAT_GC) then
                           if (tot_metaMissing > 0 .AND. nodetype(clb) == 2) then
                              Cloop0: do k=1,2**NDIM
                                 childProc = child(2,k,clb)
                                 if (childProc == mype) then
                                    childBlk = child(1,k,clb)
!!!                                 if (btst(needFlagMaskGc(childBlk),PARENTMETA_BIT)) then
                                    if (btstClr(needFlagMaskGc(childBlk),PARENTMETA_BIT)) then
                                       tot_metaMissing = tot_metaMissing - 1
#ifdef DEBUG_STATEMACHINE
                                       print*,'@',mype,' state',state,'CBlock',childBlk,' found parent,',&
                                         ' tot_metaMissing reduced to',tot_metaMissing
#endif
                                       call deferredSetDepends(this, childBlk, clb)
#ifdef DEBUG_STATEMACHINE
                                       if (.NOT. deferRemoteDependencyScan) then
                                          print*,'@',mype,' state',state,'clb',clb,' found potential Uncle,',&
                                            ' setting deferRemoteDependencyScan = .TRUE. !'
                                       else
                                          print*,'@',mype,' state',state,'clb',clb,' found potential Uncle,',&
                                            ' confirming deferRemoteDependencyScan == .TRUE.'
                                       end if
#endif
                                       deferRemoteDependencyScan = .TRUE.

                                       if (tot_metaMissing == 0) then
                                          holdLdepends = .FALSE.
                                          exit Cloop0
                                       end if
                                    end if
                                 end if
                              end do Cloop0
                           end if
                           if (ANY(surr_blks(3,1:3,1:1+2*K2D,1:1+2*K3D,clb) == 2)) then
#ifdef DEBUG_STATEMACHINE
                              if (.NOT. deferRemoteDependencyScan) then
                                 print*,'@',mype,' state',state,'clb',clb,' found potential Uncle,',&
                                         ' setting deferRemoteDependencyScan = .TRUE. !'
                              else
                                 print*,'@',mype,' state',state,'clb',clb,' found potential Uncle,',&
                                         ' confirming deferRemoteDependencyScan == .TRUE.'
                              end if
#endif
                              deferRemoteDependencyScan = .TRUE.
                           end if
                        end if  !if(pFam
                        ! neighDir = direction of ppBlk as seen from destBlk
                        ! invDir   = direction of destBlk as seen from ppBlock
                        do invDir = NDirsReverse,1,-1
                           call neighInfoReverse(clb,invDir,destBlk,destProc,pFam)
                           if (destProc == mype) then
                              if (isWanted(destBlk,ntype,lev) .AND. readyInUnk(destBlk)==0) then
#ifdef DEBUG_STATEMACHINE
                                 print*,'@',mype,' state',state,' Block',destBlk,' is wanted, type/lev:',&
                      nodetype(destBlk),lrefine(destBlk),', not readyInUnk,',&
                                      ldependents_count(destBlk),' dependents, invDir',invDir,&
                                      ' from (',ppBlk,',',isrc,')'
#endif
                                 neighDir = revDir2Dir(invDir,clb,pFam)
                                 if (needs(destBlk,neighDir,ppBlk,isrc,pFam)) then

                                    call bset(remoteAvailFlagMaskGc(destBlk), neighDir) ! set neighDir bit in remoteAvail bitmask
#ifdef DEBUG_STATEMACHINE
                                    print 30, mype,state,destBlk,needFlagMaskGc(destBlk),&
                                            localAvailFlagMaskGc(destBlk), remoteAvailFlagMaskGc(destBlk)
#endif
                                    call testDepsAndProcessBlock(this,destBlk,checkingRemote=.TRUE., &
                                         prevNumReleased=prevNumReleased, &
                                         doExitLoop=doGotoTop,outNextState=next_state)
                                    if (doGotoTop) goto 1
                                 end if    !needs(...)
                              end if       !isWanted(...)
                           end if          !(destProc == mype)
                        end do             !do invDir = NDirsReverse,1,-1
                        section_processed(clb-strt_buffer+1) = .TRUE.
                        progress = .TRUE.
                     end if
                  end do
                  psections_processed(isrc) = .TRUE.
               end if           !if (.NOT. psections_processed(isrc))
            end do              !do i=1,receivedCount

            if (receivedCount > 0) then !unnecessary or even dysfunctional test?
               if (deferRemoteDependencyScan) then !Scan for remote 'uncle' dependencies
#ifdef DEBUG_STATEMACHINE
                  print*,'>',mype,' state',state,&
                                      ' Starting deferred RemoteDependencyScan'
#endif
                  do blk = 1,lnblocks
#ifdef DEBUG_STATEMACHINE
                                 print*,'>',mype,' state',state,' Block',blk, &
                                      ', deferred RemoteDependencyScan:',&
                                      'isWanted,readyInUnk:',isWanted(blk,ntype,lev),readyInUnk(blk)
#endif
                     if (isWanted(blk,ntype,lev)) then      ! if this block is among the ones requested (by nodetype, etc.)
                        if (readyInUnk(blk) == 0) then
#ifdef DEBUG_STATEMACHINE
                           print 35, mype,state,blk,needFlagMaskGc(blk),localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
                           if (.NOT.btst(needFlagMaskGc(blk),PARENTMETA_BIT)) then
                              DLoop0: do neighDir = 3**NDIM+1, NDirs*2
                                 call neighInfo(blk,neighDir,ppBlk,ppProc,pFam)
#ifdef DEBUG_STATEMACHINE
                                 print*,'>',mype,' state',state,' Block',blk, &
                                      ', neighDir/neighBlk/neighProc:',&
                                      neighDir,ppBlk,ppProc,' (deferRemoteDependencyScan)'
#endif
                                 if (ppProc .NE. mype .AND. ppProc .GE. 0) then
                                    if (btst(needFlagMaskGc(blk),neighDir)) then
                                    if (needs(blk,neighDir,ppBlk,ppProc,pFam)) then
                                       if (haveReceivedBlkProc(ppBlk,ppProc)) then
                                          call bset(remoteAvailFlagMaskGc(blk), neighDir)
#ifdef DEBUG_STATEMACHINE
                                          print 35, mype,state,blk,needFlagMaskGc(blk),&
                                               localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
                                          call testDepsAndProcessBlock(this,blk,checkingRemote=.TRUE., &
                                               prevNumReleased=prevNumReleased, &
                                               doExitLoop=doGotoTop,outNextState=next_state)
                                          if (doGotoTop) goto 1
                                          if (readyInUnk(blk) > 0) EXIT DLoop0
                                       end if !if (haveReceivedBlkProc
                                    end if !if (needs
                                    end if !if (btst(needFlag,neighDir
                                 end if    !if (ppProc
                              end do DLoop0  ! xneighDir
                           end if            ! if (.NOT.btst
                        end if          !if (readyInUnk
                     end if             !if (isWanted
                  end do                !do blk
                  deferRemoteDependencyScan = .FALSE. !Just did it and completed it.
               end if                                 !if (deferRemoteDependencyScan
            end if                                    !if (receivedCount > 0)

            if (next_state == S_wait .and. inputsExhausted) then
               if (progress) then
                  next_state = S_checkQueue
               else
                  next_state = S_checkReady
               end if
            end if
#ifdef DEBUG_STATEMACHINE
            print*,'E',mype,' state',state,' - allReceivedCount is',allReceivedCount,' of',ij,', ',&
                    'receivedCount',receivedCount,' becomes 0, inputsExh',inputsExhausted
#endif
            receivedCount = 0   !We can MPI_Wait/MPI_Test again in S_wait!

         case(S_shutdown)

#ifdef DEBUG_STATEMACHINE
            if (.NOT. receivedAll) then
               print*,'@',mype,' state',state,'BlockGetter reached S_shutdown but receivedAll is FALSE!'
            end if
            if (allReceivedCount < ij) then
               print*,'@',mype,' state',state,'BlockGetter reached S_shutdown but (allReceivedCount < ij):',&
               '(',allReceivedCount,' <',ij,')!'
            end if
            if (.NOT. inputsExhausted) then
               print*,'@',mype,' state',state,'BlockGetter reached S_shutdown but inputsExhausted is FALSE!'
            end if
#endif
180         format('Getter_shutdn(',I2,') @',I5,1x,A,I3,' of',I3,1x,A,' requests',A)
#define DEFER_WAITS_UNTIL_DESTRUCTION
#ifndef DEFER_WAITS_UNTIL_DESTRUCTION
            associate(ji               => commCtl % numSReq,          &
                      sendrequest      => commCtl % sendrequest,     &
                      sendstatus       => commCtl % sendstatus,      &
                      allSentCount => commCtl % allSentCount)
              do while (allSentCount < ji)
#ifdef DEBUG_ASYNC_COMM
                 print 180,pfam,mype,'Waiting for',(ji-allSentCount),ji,'send','...'
!                 print*,'Getter_shutdn @',mype,' Waiting for', (ji-allSentCount),' of',ji,' send requests...'
#endif
                 if (allSentCount < ji) then
                    call MPI_Waitsome(ji,sendrequest,sentCount, &
                         receivedInd, sendstatus, ierrorcode)
                    allSentCount = allSentCount + sentCount
#ifdef DEBUG_GRIDCOMM
                    print 96,mype,sentCount,ji,(receivedInd(i),i=1,sentCount)
96                  format('On',I2,': WaitSome sent',I2,'/',I2,',',20(1x,I4,:,','))
#endif
                 end if
              end do
            end associate

            if (allReceivedCount < ij) then
#ifdef DEBUG_ASYNC_COMM
               print 180,pfam,mype,'Testing for',(ij-allReceivedCount),ij,'receive','!'
!               print*,'Getter_shutdn @',mype,' Testing for', (ij-allReceivedCount),' of',ij,' receive requests!'
#endif
               call MPI_TestAll(ij,recvrequest,testFlag,&
                    recvstatus, ierrorcode)
               if (.NOT. testFlag) then
!!$                  do i=1,ij
!!$                     print*,'@',mype,' state',state,'Cancelling request',i,' ...'
!!$                     call MPI_Cancel(recvrequest(i),ierrorcode)
!!$                  end do
!!$                  print*,'@',mype,' state',state,'Cleanup: MPI_Cancel of',ij,' requests appears successful.'
#ifdef DEBUG_ASYNC_COMM
                  print 180,pfam,mype,'Waiting for',(ij-allReceivedCount),ij,'receive','...'
!                  print*,'Getter_shutdn @',mype,' Waiting for', (ij-allReceivedCount),' of',ij,' receive requests...'
#endif
                  call MPI_WaitAll(ij,recvrequest,&
                       MPI_STATUSES_IGNORE, ierrorcode)
#ifdef DEBUG_GRIDCOMM
                  print*,'@',mype,' state',state,'Cleanup: MPI_WaitAll of',ij,' requests was successful.'
#endif
                  allReceivedCount = ij
               else
#ifdef DEBUG_GRIDCOMM
                  print*,'@',mype,' state',state,'Cleanup: MPI_TestAll of',ij,' requests was successful.'
#endif
                  allReceivedCount = ij
               end if
            end if

            if (.NOT. receivedAll) then
!               call Driver_abort('BlockGetter reached S_shutdown but receivedAll is FALSE!')
            end if
            if (allReceivedCount < ij) then
               call Driver_abort('BlockGetter reached S_shutdown but (allReceivedCount < ij)!')
            end if
            if (.NOT. inputsExhausted) then
!               call Driver_abort('BlockGetter reached S_shutdown but inputsExhausted is FALSE!')
            end if
#elif defined(DEBUG_ASYNC_COMM)
!!$            associate(ji               => commCtl % numSReq,          &
!!$                      allSentCount => commCtl % allSentCount)
!!$              if (allSentCount < ji) &
!!$                 print*,'Getter_shutdn@',mype,' Deferring Wait for', allSentCount,' of',ji,' send requests.'
!!$            end associate
            if (allReceivedCount < ij) then
               print 180,pfam,mype,'Deferring Test for',(ij-allReceivedCount),ij,'receive','.'
!               print*,'Getter_shutdn @',mype,' Deferring Test for', (ij-allReceivedCount),' of',ij,' receive requests.'
            end if
#endif

            next_state = S_error
            !!TODO!! cleanup if necessary
            blockID = -1
#ifdef DEBUG_STATEMACHINE
            print*,'@',mype,' state',state,' RETURN block',blockID
#endif
            RETURN ! (-1,-1)          ! End of stream indication if needed; or go to S_error directly.


         case(S_error)

            call Driver_abort('gr_pmBlockGetter: State machine reached state S_error')

         end select
      end do


    end associate

    ! Unreachable code:
    blockID = -1

  contains

    PURE logical function ldependents_ok(lb)
      integer, intent(IN) :: lb
      if (this % ldependentsDontMatter) then
         ldependents_ok = .TRUE.
      else if (this % ldependents_count(lb) .NE. 0) then
         ldependents_ok = .FALSE.
      else if (.NOT. this % holdLdepends) then
         ldependents_ok = .TRUE.
      else if (ANY(surr_blks(3,1:3,1:1+2*K2D,1:1+2*K3D,lb) == 2)) then
         ldependents_ok = .FALSE.
      else
         ldependents_ok = .TRUE.
      end if
    end function ldependents_ok

    PURE integer function revDir2Dir(invDir,lb,pFam) result(dir)
      integer, intent(in) :: invDir
      integer, intent(in) :: lb
      integer, intent(in) :: pFam

      select case (pFam)
      case(GRID_PAT_GC)
         dir = NDirs + 1 - invDir
      case(GRID_PAT_PROLONG)
         dir = NDirs + 1 - invDir !??
      case(GRID_PAT_FCORR)
         dir = mod(6-invDir,2) + 1 + ((invDir-1)/2)*2
      case(GRID_PAT_RESTRICT)
         dir = which_child(lb)
      end select
    end function revDir2Dir

  ! Does the (local) lb need (remote_blk,remote_pe) for the given direction?
    logical function needs(lb,nDir,remote_blk,remote_pe,pFam)
      integer,intent(in) :: lb,nDir,remote_blk,remote_pe
      integer,intent(in) :: pFam

      logical :: needAny
      integer :: iDir,neighBlk,neighProc

      needs = .FALSE.              !!DEV: TODO for other comm patterns


      ASSOCIATE(need        => this % needFlagMaskGc(lb))
        needAny = btst(need,nDir)
      end ASSOCIATE


      if (needAny) then
         select case (pfam)
         case(GRID_PAT_GC)
            call neighInfo(lb,nDir,neighBlk,neighProc)
            if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
               needs = .TRUE.
               return
            end if
!!$         if (nodetype(lb) == 1    .AND. &
!!$              parent(2,lb) == mype .AND. &
!!$              any(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,lb) > -20 .and.    & 
!!$              surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,lb) < 0)) then
!!$            call neighInfo(parent(1,lb),nDir,neighBlk,neighProc)
!!$            do iDir = 1,3**NDIM
!!$               if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
!!$                  needs = .TRUE.
!!$                  return
!!$               end if
!!$            end do
!!$         end if
         case(GRID_PAT_RESTRICT)
            neighBlk  = child(1,nDir,lb)
            neighProc = child(2,nDir,lb)
            if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
               needs = .TRUE.
               return
            end if
         end select
      end if

  end function needs

  logical function haveReceivedBlkProc(remBlk,remProc)
    integer,intent(in) :: remBlk,remProc
    ! Note: Currently make no used of remProc.

    associate(              psections_processed   => this % psections_processed)

      if (psections_processed(remProc)) then
         haveReceivedBlkProc = .TRUE.
      else
         haveReceivedBlkProc = .FALSE.
      end if
      ! DEV: make any used of ( section_processed(clb-strt_buffer+1) = .TRUE. ) ?
    end associate

  end function haveReceivedBlkProc

  end function pmBlockGet_


  function pmBlockGet(this) result(blockID)
    integer ::                     blockID
    class(gr_pmBlockGetter_t),intent(inout) :: this

    blockID = pmBlockGet_(this, this % commCtl)
  end function pmBlockGet

  ! Test whether the dependencies of a given block are now satisfied according to the bit masks and,
  ! if yes, process the block towards delivery.
  !   destBlk should identify a local block that is
  !     * wanted,
  !     * not yet ready in UNK.
  subroutine testDepsAndProcessBlock(this,destBlk,checkingRemote,prevNumReleased,doExitLoop,outNextState)
    type(gr_pmBlockGetter_t),intent(inout) :: this
    integer,intent(in)    :: destBlk
    logical,intent(in)    :: checkingRemote
    integer,intent(in)    :: prevNumReleased
    LOGICAL,INTENT(OUT)   :: doExitLoop
    INTEGER,INTENT(INOUT) :: outNextState

    integer :: inDir
    integer :: neighBlk, neighProc, NDirs
    logical :: skip

    logical :: get_data_quickly = .FALSE. !to become runtime paramter
    logical :: push_data_out_quickly = .FALSE. !to become runtime paramter
    logical :: push_data_out_very_quickly = .FALSE. !to become runtime paramter

    doExitLoop = .FALSE.
    associate(pfam             => this % patternFamily,      &
              state            => this % state,              &
              tot_ldependents_count => this % tot_ldependents_count,    & !may nor be needed
              ldependents_count     => this % ldependents_count,    &
              ldependentsDontMatter     => this % ldependentsDontMatter,    &
              numReleased           => this % numReleased,              &
              numUnkReleased           => this % numUnkReleased,              &
              readyInUnk            => this % readyInUnk,    &
              Qcur                  => this % Qcur,          &
              Qcur_inP              => this % Qcur_inP          &
              )

      select case (pFam)
      case(GRID_PAT_GC)
         NDirs = 3**NDIM
      case(GRID_PAT_PROLONG)
         NDirs = 3**NDIM
      case(GRID_PAT_FCORR)
         NDirs = 2 * NDIM
      case default !(GRID_PAT_RESTRICT)
         NDirs = 2**NDIM
      end select

!!$           ! was this the last missing piece?
      if (neededGcsAreAvailable(this,destBlk)) then

         do inDir = 1,dependencyDirs(destBlk,pFam)!NDirs
            call neighInfo(destBlk,inDir,neighBlk,neighProc,pFam,skip)
            if (neighProc == mype .AND. neighBlk .NE. destBlk .AND. .NOT. skip) then
               ldependents_count(neighblk) = ldependents_count(neighblk) - 1
#ifdef DEBUG_STATEMACHINE
               print*,'p',mype,' state',state,' Block',neighBlk,&
                  ' had ldependents decreased to',ldependents_count(neighBlk),&
                  ' because we release local destBlk',destBlk
#endif
               ! if (ldependents_ok(neighblk)) then ! if(.NOT.checkingRemote) ?
               if (ldependentsDontMatter .OR. ldependents_count(neighblk) == 0) then
                  if (readyInUnk(neighblk) > 0) then
                     numUnkReleased = numUnkReleased + 1
                  else
                     numReleased    = numReleased + 1
                  end if
               end if
               tot_ldependents_count = tot_ldependents_count - 1
               if (.NOT. get_data_quickly) then
                  if (numReleased > prevNumReleased) outNextState = S_checkLocal
               end if
            end if
         end do !inDir

         if (readyInUnk(destBlk) < 1) then
               ! This does update destBlk in UNK, makes it ready to be delivered.
            call update_perm_griddata(this,destBlk) ! also applies domain boundary conditions
            readyInUnk(destBlk) = 1
#ifdef DEBUG_STATEMACHINE
            print*,'p',mype,' state',state,' Block',destBlk,' became ready in Unk.'
#endif
         end if
         if (checkingRemote) outNextState = S_checkReady

         if (ldependents_ok(destBlk)) then ! if destBlk has no further local dependents
            if (readyInUnk(destBlk) < 2) then !ALWAYS TRUE when we get here (?)
               Qcur(Qcur_inP) = destBlk
               Qcur_inP = Qcur_inp + 1
               readyInUnk(destBlk) = 2
#ifdef DEBUG_STATEMACHINE
               print*,'p',mype,' state',state,' Block',destBlk,' was queued.'
#endif
               if ((checkingRemote .AND. push_data_out_quickly) .OR. &
                      (.NOT.checkingRemote)                      ) then
                  outNextState = S_checkQueue
               end if
               if ((     checkingRemote .AND. push_data_out_quickly     ) .OR. &
                   (.NOT.checkingRemote .AND. push_data_out_very_quickly)) then
#ifdef DEBUG_STATEMACHINE
                  print*,'p',mype,' state',state,' DoExitLoop to Top'
#endif
                  goto 11
               end if
            end if !readyInUnk(destBlk)
         else
            if (.NOT. checkingRemote) outNextState = S_checkReady
         end if !ldependents_ok(destBlk)

      end if !neededGcsAreAvailable(...)
    end associate
    RETURN

11  doExitLoop = .TRUE.
  contains

    PURE logical function ldependents_ok(lb)
      integer, intent(IN) :: lb
      if (this % ldependentsDontMatter) then
         ldependents_ok = .TRUE.
      else if (this % ldependents_count(lb) .NE. 0) then
         ldependents_ok = .FALSE.
      else if (.NOT. this % holdLdepends) then
         ldependents_ok = .TRUE.
      else if (ANY(surr_blks(3,1:3,1:1+2*K2D,1:1+2*K3D,lb) == 2)) then
         ldependents_ok = .FALSE.
      else
         ldependents_ok = .TRUE.
      end if
    end function ldependents_ok

  end subroutine testDepsAndProcessBlock

  SUBROUTINE deferredSetDepends(this, blk, pclb)
    type(gr_pmBlockGetter_t),intent(inout) :: this
    integer, intent(in) :: blk, pclb !'parent cached local block'

    integer :: i,j,k
    integer :: neighDir, neighBlk, neighProc
    integer :: misNei, pneighDir, pneighBlk, pneighProc
    integer :: ia,ib,ja,jb,ka,kb
    integer :: iD,jD,kD, iDa,iDb, jDa,jDb, kDa,kDb
    integer :: ioff,joff,koff,nguarda
    integer :: pFam
    integer :: state = S_checkRemote

    pFam = this % patternFamily
    ka = 1; kb = 1+2*K3D
    ja = 1; jb = 1+2*K2D
    ia = 1; ib = 3
    ASSOCIATE(psurr_blk => surr_blks(:,:,:,:,pclb), &
              tot_ldependents_count => this % tot_ldependents_count,    & ! needed?
              ldependents_count     => this % ldependents_count,    &
              localAvailFlagMaskGc  => this % localAvailFlagMaskGc,  &
              remoteAvailFlagMaskGc => this % remoteAvailFlagMaskGc, &
              needFlagMaskGc        => this % needFlagMaskGc,        &
              remoteNeedFlagMaskGc  => this % remoteNeedFlagMaskGc  & ! maybe UNNEEDED, use needFlagMaskGc?
            )
!!#ifdef PM_OPTIMIZE_MORTONBND_FETCHLIST
         ! Do this optimization only if parent does not touch a domain boundary
         ! anywhere (otherwise boundary condition can be called with invalid
         ! input data). Cf. mpi_morton_bnd.
      if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
         nguarda = max(nguard,nguard_work)
#if NDIM >= 3
         if (nguarda .LE. 3NZB/2) then
            koff = mod((which_child(blk)-1)/4,2)
            ka = 1+koff; kb = 1+K3D+koff
         end if
#endif
#if NDIM >= 2
         if (nguarda .LE. NYB/2) then
            joff = mod((which_child(blk)-1)/2,2)
            ja = 1+joff; jb = 1+K2D+joff
         end if
#endif
         if (nguarda .LE. NXB/2) then
            ioff = mod(which_child(blk)-1,2)
            ia = 1+ioff; ib = 2+ioff
         end if
      end if
!!#endif

      misNei = 0
      iDa = ia; jDa = ja; kDa = ka
      iDb = ib; jDb = jb; kDb = kb
      do neighDir = 1,3**NDIM
         call neighInfo(blk,neighDir,neighBlk,neighProc,pFam)
#ifdef DEBUG_STATEMACHINE
30       format(1x,'D',I12,'  state',I12,'  Block',I12,&
                       ' ... needed: o',OCTFLAG_FORMAT,', local/remote avail: o',OCTFLAG_FORMAT,',o',OCTFLAG_FORMAT)
         print*,'D',mype,' state',state,' Block',blk, &
              ', neighDir/neighBlk/neighProc:',&
              neighDir,neighBlk,neighProc,' (was deferred)'
#endif
         if (neighProc == mype) then ! neighBlk is local:
         else if (neighProc .GE. 0) then
         else if (neighBlk .LE. PARAMESH_PHYSICAL_BOUNDARY) then
         else if ((neighBlk == -1 .OR. neighProc == -1) .AND. &
            nodetype(blk) == 1) then !should only happen for surr_blks?
            misNei = misNei + 1                        ! a missing surr_blks neighbor!
#ifdef DEBUG_STATEMACHINE
            print*,'@',mype,' state',state,' Block',blk,' has',misNei,' missing surr_blks neighbors!'
            print 30, mype,state,blk,needFlagMaskGc(blk),localAvailFlagMaskGc(blk), remoteAvailFlagMaskGc(blk)
#endif
            iD  = mod( neighDir-1,   3)+1
            jD  = mod((neighDir-1)/3,3)+1
            kD  =     (neighDir-1)/9   +1
            if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
               kDa = min(kDa,kD-1) ; kDb = max(kDb,kD+1)
               jDa = min(jDa,jD-1) ; jDb = max(jDb,jD+1)
               iDa = min(iDa,iD-1) ; iDb = max(iDb,iD+1)
            end if
         end if     ! if (neighProc == mype
      end do        !do neighDir

      if (misNei == 0) then
         print*,'@',mype,' state',state,' Block',blk,' having',misNei,&
              ' missing surr_blks neighbors SHOULD NOT HAVE BEEN DEFERRED!'
         RETURN
      end if
     
#ifdef PM_OPTIMIZE_MORTONBND_FETCHLIST
             ! Do this optimization only if parent does not touch a domain boundary
             ! anywhere (otherwise boundary condition can be called with invalid
             ! input data). Cf. mpi_morton_bnd.
      if (All(psurr_blk(1,1:3,1:1+2*K2D,1:1+2*K3D) > -20)) then
         nguarda = max(nguard,nguard_work)
#if NDIM >= 3
         if (nguarda .LE. NZB/2) then
            ka = max(ka,kDa); kb = min(kb,kDb)
         end if
#endif
#if NDIM >= 2
         if (nguarda .LE. NYB/2) then
            ja = max(ja,jDa); jb = min(jb,jDb)
         end if
#endif
         if (nguarda .LE. NXB/2) then
            ia = max(ia,iDa); ib = min(ib,iDb)
         end if
      end if
#endif
!!$                                print*,'@',mype,' state',state,' Block',blk,'in ASSOCIATE...'
      do k = ka,kb
         do j = ja,jb
            do i = ia,ib
               pneighBlk  = psurr_blk(1,i,j,k)
               pneighProc = psurr_blk(2,i,j,k)
#ifdef DEBUG_STATEMACHINE
               print*,'@',mype,' state',state,' Block',blk,' i,j,k:',i,j,k,', pneighBlk/pneighProc:',&
                         pneighBlk,pneighProc,' (deferred)'
#endif
               if (pNeighBlk > 0 .AND. pneighProc > -1) then !.AND. (psurr_blk(3,ia,ja,ka) < 3) ?
                  pneighDir = 3**NDIM + i + K2D*(j-1)*3 + K3D*(k-1)*9
#ifdef DEBUG_STATEMACHINE
                  print*,'@',mype,' state',state,' Block',blk,' i,j,k:',i,j,k, &
                            ', pneighDir/pneighBlk/pneighProc:',&
                            pneighDir,pneighBlk,pneighProc,' (deferred)'
#endif
                  if (pneighProc == mype) then ! pneighBlk is local:
                     ldependents_count(pneighBlk) = ldependents_count(pneighBlk) + 1
                     tot_ldependents_count = tot_ldependents_count + 1
                     call bset(localAvailFlagMaskGc(blk), pneighDir) ! set pneighDir bit in localAvail bitmask
                  end if
                  call bset(needFlagMaskGc(blk), pneighDir) ! set neighDir bit in remoteNeed bitmask
                  if (pneighProc .NE. mype) then
                     call bset(remoteNeedFlagMaskGc(blk), pneighDir) ! set neighDir bit in remoteNeed bitmask
                  end if
               end if !(pNeighBlk > 0
            end do
         end do
      end do
           ! DEV: No special handling of (lsingular_line), (cf. mpi_morton_bnd)
    end ASSOCIATE


  END SUBROUTINE deferredSetDepends


  subroutine bset(i,pos)
    integer(i8),intent(INOUT) :: i
    integer,intent(in) ::    pos
    i = ibset(i,pos-1)
  end subroutine bset
  subroutine bclr(i,pos)
    integer(i8),intent(INOUT) :: i
    integer,intent(in) ::    pos
    i = ibclr(i,pos-1)
  end subroutine bclr
  logical function btst(i,pos)
    integer(i8),intent(IN) :: i
    integer,intent(in) ::    pos
    btst = (ibits(i,pos-1,1) == 1_i8)
  end function btst
  logical function btstClr(i,pos)
    integer(i8),intent(INOUT) :: i
    integer,intent(in) ::    pos
    btstClr = (ibits(i,pos-1,1) == 1_i8)
    if (btstClr) call bclr(i,pos)
  end function btstClr

  ! return the number of "directions" in which lb can (potentially) have dependencies
  PURE integer function dependencyDirs(lb,pFam) result(n)
    integer,intent(in)  :: lb,pFam

    n = 0
    select case (pFam)
    case(GRID_PAT_GC)
       n = 3**NDIM
       if (lb > 0) then
          if (nodetype(lb) == 1    .AND. &
              parent(1,lb) > 0     .AND. &
              any(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,lb) > -20 .and.    &
                  surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,lb) < 0)) then
             n = 2*n
          end if
       end if
    case(GRID_PAT_PROLONG)
       if (nodetype(lb) == 1    .AND. newchild(lb)) then
          n = 3**NDIM
       end if
    case(GRID_PAT_FCORR)
       if (nodetype(lb) == 1    .AND. &
           parent(1,lb) > 0     .AND. &
            any(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,lb) > 0 .and.    &
                surr_blks(3,1:3,1:1+2*K2D,1:1+2*K3D,lb) == 2)) then
          n = 2 * NDIM
       end if
    case default !(GRID_PAT_RESTRICT)
       if (nodetype(lb) == 2 .or.                                       &
           (advance_all_levels .and. nodetype(lb) == 3)) then
          n = 2**NDIM
       end If
    end select
  end function dependencyDirs

  subroutine neighInfo_surr(surr,in,neighBlk,neighProc)
    integer,intent(in)  :: surr(3,3**NDIM)
    integer,intent(in)  :: in
    integer,intent(OUT) :: neighBlk,neighProc
    neighBlk  = surr(1,in)
    neighProc = surr(2,in)
  end subroutine neighInfo_surr
  subroutine neighInfo_2surr(lb,in,neighBlk,neighProc,skip)
    integer,intent(in)  :: lb
    integer,intent(in)  :: in
    integer,intent(OUT) :: neighBlk,neighProc
    logical,intent(OUT),OPTIONAL :: skip
    integer,parameter :: n1=3**NDIM, n2=2*n1!, n3=n2+2**NDIM
    integer :: parentBlk, parentProc, neighDir
    integer :: iblk, remote_block, remote_pe
    integer :: iD,jD,kD, iC,jC,kC
    logical :: lfound
    neighBlk  = -1
    neighProc = -1
    if (present(skip)) skip = .FALSE.
    if (in .LE. n1) then
       call neighinfo_surr(surr_blks(1,1,1,1,lb),in,neighBlk,neighProc)
    else if (in .LE. n2) then
       if (parent(1,lb) .GE. 1 .AND. parent(2,lb) .GE. 0) then
          parentBlk  = parent(1,lb)
          parentProc = parent(2,lb)
          neighDir   = in-n1
          if (parentProc == mype) then
             call neighinfo_surr(surr_blks(1,1,1,1,parentBlk),in-n1,neighBlk,neighProc)
             if (present(skip)) then
                if (neighProc == mype) then
                   if (All(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,parentBlk) > -20)) then
                      iC = mod( which_child(lb)-1   ,2)*2; iD  = mod( neighDir-1,   3)
                      jC = mod((which_child(lb)-1)/2,2)*2; jD  = mod((neighDir-1)/3,3)
                      kC = mod((which_child(lb)-1)/4,2)*2; kD  =     (neighDir-1)/9
                      if (abs(iD-iC)>1) skip = .TRUE.
                      if (abs(jD-jC)>1) skip = .TRUE.
                      if (abs(kD-kC)>1) skip = .TRUE.
#ifdef DEBUG_STATEMACHINE
                      print*,'  neighInfo_2surr: Skip =',skip,',iC,iD,jC,jD,kC,kD:',iC,iD,jC,jD,kC,kD
#endif
                   end if
                end if
             end if
          else if (n2-in == in-(n1+1)) then !the parent itself
             neighBlk  = parent(1,lb)
             neighProc = parent(2,lb)
          else if (parentBlk > 0) then ! .AND. psections_processed(parentProc)) then
             lfound = .False.
             iblk = ladd_strt(parent(2,lb))
             Do While(.Not.lfound.And.                                  &
                      iblk <= ladd_end(parent(2,lb)))
                If (parent(1,lb) == laddress(1,iblk).And.                &
                    parent(2,lb) == laddress(2,iblk) ) Then
                   remote_block = iblk
                   remote_pe    = mype
                   lfound = .True.
                Else
                   iblk = iblk+1
                End If
             End Do
             If (.Not.lfound) Then
                Write(*,*) 'Error in surr_blks_lkup for parent : ',       &
                     'remote block ',parent(:,lb),                     &
                     ' not located on pe ',mype,                       &
                     ' while processing blk ',lb,mype
                Call Driver_abort('gr_pmBlockGetter: surr_blks_lkup for parent failed')
             End If

             call neighinfo_surr(surr_blks(1,1,1,1,remote_block),in-n1,neighBlk,neighProc)
             if (present(skip)) then
                if (neighProc == mype) then
                   if (All(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,remote_block) > -20)) then
                      iC = mod( which_child(lb)-1   ,2)*2; iD  = mod( neighDir-1,   3)
                      jC = mod((which_child(lb)-1)/2,2)*2; jD  = mod((neighDir-1)/3,3)
                      kC = mod((which_child(lb)-1)/4,2)*2; kD  =     (neighDir-1)/9
                      if (abs(iD-iC)>1) skip = .TRUE.
                      if (abs(jD-jC)>1) skip = .TRUE.
                      if (abs(kD-kC)>1) skip = .TRUE.
#ifdef DEBUG_STATEMACHINE
                      print*,'  neighInfo_2surr: Skip =',skip,',iC,iD,jC,jD,kC,kD:',iC,iD,jC,jD,kC,kD
#endif
                   end if
                end if
             end if

          else
             ! .....
          end if
       end if
!!$    else if (in .LE. n3) then
!!$       neighBlk  = child(1,in-n2,lb)
!!$       neighProc = child(2,in-n2,lb)
    end if
  end subroutine neighInfo_2surr
  subroutine neighInfo(lb,in,neighBlk,neighProc,pFam,skip)
    integer,intent(in)  :: lb
    integer,intent(in)  :: in
    integer,intent(OUT) :: neighBlk,neighProc
    integer,intent(in),OPTIONAL :: pFam
    logical,intent(OUT),OPTIONAL :: skip
    integer,parameter :: n1=3**NDIM, n2=2*n1, n3=n2+2**NDIM
    integer :: myPFam
    neighBlk  = -1
    neighProc = -1
    myPFam = GRID_PAT_GC
    if (present(pFam)) myPFam = pFam

    if (myPFam == GRID_PAT_GC) then
       call neighInfo_2surr(lb,in,neighBlk,neighProc,skip)
    else
       if (present(skip)) skip = .FALSE.
       neighBlk  = child(1,in,lb)
       neighProc = child(2,in,lb)
    end if
  end subroutine neighInfo
  subroutine neighInfoReverse(lb,in,neighBlk,neighProc,pFam,skip)
    integer,intent(in)  :: lb
    integer,intent(in)  :: in
    integer,intent(OUT) :: neighBlk,neighProc
    integer,intent(in),OPTIONAL :: pFam
    logical,intent(OUT),OPTIONAL :: skip
    integer,parameter :: n1=3**NDIM, n2=2*n1, n3=n2+2**NDIM
    integer :: myPFam
    neighBlk  = -1
    neighProc = -1
    myPFam = GRID_PAT_GC
    if (present(pFam)) myPFam = pFam

    if (myPFam == GRID_PAT_GC) then
       call neighInfo_2surr(lb,in,neighBlk,neighProc,skip)
    else
       if (present(skip)) skip = .FALSE.
       neighBlk  = parent(1,lb)
       neighProc = parent(2,lb)
    end if
  end subroutine neighInfoReverse


  logical function isWanted(blk,nodetype,lev)
    integer,intent(in) :: blk
    integer,intent(in) :: nodetype
    integer,intent(in) :: lev
    if (lev == UNSPEC_LEVEL) then
       ! No level given at creation
       isWanted = gr_blockMatch(blk, nodetype)
    else
       isWanted = gr_blockMatch(blk, nodetype, lev)
    end if
  end function isWanted
  logical function neededGcsAreAvailable(getter,blk)
    type(gr_pmBlockGetter_t), intent(IN) :: getter
    integer, intent(in)                  :: blk

    ASSOCIATE(locAvail    => getter % localAvailFlagMaskGc(blk)    , &
              remoteAvail => getter % remoteAvailFlagMaskGc(blk)   , &
              need        => getter % needFlagMaskGc(blk))
      neededGcsAreAvailable = (iand(ior(locAvail, remoteAvail), need) == need)
    end ASSOCIATE
  end function neededGcsAreAvailable

  subroutine update_perm_griddata(this,blk)
    use paramesh_interfaces, ONLY: amr_1blk_guardcell, &
                                   gr_amr1blkGcToPerm, &
                                   gr_amr1blkRestrict
    type(gr_pmBlockGetter_t),intent(inout) :: this
    integer, intent(in) :: blk  !ID of the local block about to be shipped out

    integer,parameter :: nlayers = NGUARD !should be ignored by called subroutines
    logical,parameter :: l_srl_only = .FALSE. !for now
    integer,parameter :: icoord = 0 !for now
    logical,parameter :: ldiag = .TRUE. !for now

#if N_DIM == 1
    integer(kind=i27b) :: REGIONMASK = int(o'0000000007',i27b)
#endif
#if N_DIM == 2
    integer(kind=i27b) :: REGIONMASK = int(o'0000000777',i27b)
#endif
#if N_DIM == 3
    integer(kind=i27b) :: REGIONMASK = int(o'0777777777',i27b)
#endif

    integer(kind=i27b) :: parentPresentRegions

    associate(pfam             => this % patternFamily,      &
              iopt             => this % iopt,               &
              lcc              => this % lcc,                &
              lfc              => this % lfc,                &
              lec              => this % lec,                &
              lnc              => this % lnc,                &
              nlayersx         => this % nlayersx,           &
              nlayersy         => this % nlayersy,           &
              nlayersz         => this % nlayersz,           &
              locAvail         => this % localAvailFlagMaskGc(blk), &
              remoteAvail      => this % remoteAvailFlagMaskGc(blk) &
              )
      select case (pfam)
      case(GRID_PAT_GC)
         parentPresentRegions=iand(int(shiftr(ior(locAvail,remoteAvail),3**NDIM), &
                                       kind=i27b),&
                                   REGIONMASK)
#ifdef DEBUG_STATEMACHINE
40       format(1x,A1,I12,'  Block',I12,' passing to amr_1blk_guardcell:',&
              ' parentPresentRegions = o',o9.9)
         print 40,'@',mype,blk,parentPresentRegions
#endif
         ! also applies domain boundary conditions
         call amr_1blk_guardcell(mype,iopt,nlayers,blk,mype, &
              lcc,lfc,lec,lnc,l_srl_only,                    &
              icoord,ldiag,                                  &
              nlayersx,nlayersy,nlayersz,                    &
              parentPresentRegions)
          call gr_amr1blkGcToPerm(mype,iopt,nlayers,blk,        &
                                 lcc,lfc,lec,lnc,              &
                                 nlayersx,nlayersy,nlayersz)

      case(GRID_PAT_RESTRICT)
         call gr_amr1blkRestrict(mype,blk,iopt,lcc,lfc,lec,lnc)
      case default
         call Driver_abort("gr_pmBlockGetter: unsupported pattern family!")
      end select
    end associate
  end subroutine update_perm_griddata

  ! Does the (local) destBlk need (remote_blk,remote_pe) for ANY direction?
  logical function needsForAnyDir(this,destBlk,remote_blk,remote_pe) result(needs)
    type(gr_pmBlockGetter_t),intent(in) :: this
    integer,intent(in) :: destBlk,remote_blk,remote_pe

    integer :: iDir,neighBlk,neighProc

    needs = .FALSE.              !!DEV: TODO for other comm patterns

    associate(pfam             => this % patternFamily      &
              )
      select case (pfam)
      case(GRID_PAT_GC)
         do iDir = 1,3**NDIM
            call neighInfo(destBlk,iDir,neighBlk,neighProc)
            if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
               needs = .TRUE.
               return
            end if
         end do
         if (nodetype(destBlk) == 1    .AND. &
             parent(2,destBlk) == mype .AND. &
             any(surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,destBlk) > -20 .and.    &
                 surr_blks(1,1:3,1:1+2*K2D,1:1+2*K3D,destBlk) < 0)) then
            call neighInfo(parent(1,destBlk),iDir,neighBlk,neighProc)
            do iDir = 1,3**NDIM
               if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
                  needs = .TRUE.
                  return
               end if
            end do
         end if
      case(GRID_PAT_RESTRICT)
         do iDir = 1,2**NDIM
            neighBlk  = child(1,iDir,destBlk)
            neighProc = child(2,iDir,destBlk)
            if (neighBlk == remote_blk .AND. neighProc == remote_pe) then
               needs = .TRUE.
               return
            end if
         end do
      end select
    end associate

  end function needsForAnyDir


  ! An auxiliary routine that updates metainformation for
  ! pp blocks that have just been received.
  ! Variant of, and derived from, mpi_unpack_blocks.
  subroutine pmMpiUnpackBlksFromProc(sproc,iopt, & 
     &                             lcc,lfc,lec,lnc, & 
     &                             buf_dim,R_buffer, & 
     &                             nlayersx,nlayersy,nlayersz)
    !------------------------------------------------------------------------
    !
    ! This subroutine unpacks all blocks which have been received on mype
    ! from one specifc sending proc.
    ! "Unpacking" here just means extracting Grid metainformation
    ! from a region in R_buffer (which should normally be a region in
    ! temprecv_buf into which data have been received) and copying
    ! that metainformation into the proper PARAMESH arrays; actually
    ! block payload data is not copied here but is left in
    ! temprecv_buf for subroutines like mpi_1blk_guardcell to access.
    ! Arrays like laddress, ladd_strt, ladd_end are NOT modified here,
    ! they should already have been updated according to the new
    ! communication pattern when this subroutine is called.
    !
    !
    ! Written: mpi_unpack_blocks   Maharaj Bhat & Michael Gehmeyr  March 2000
    ! Adapted: pmMpiUnpackBlksFromProc  Klaus Weide                  Jan 2021
    !------------------------------------------------------------------------
    !
    ! Arguments:
    !      sproc          sending processor id
    !      iopt           option setting for work array
    !      lcc            if true include unk data in buffer
    !      lfc            if true include facevar data in buffer
    !      lec            if true include unk_e_? data in buffer
    !      lnc            if true include unk_n data in buffer
    !      buf_dim        dimension of buffer
    !      R_buffer       receive buffer
    !
    !------------------------------------------------------------------------
    use paramesh_dimensions, ONLY: maxblocks
    use mpi_morton, ONLY: ir_buf
    use paramesh_comm_data, ONLY: amr_mpi_meshComm

    use paramesh_mpi_interfaces, only : mpi_put_buffer

#include "Flashx_mpi_implicitNone.fh"

    integer, intent(in) :: sproc,buf_dim,iopt
    logical, intent(in) :: lcc,lfc,lec,lnc
    real,    intent(IN),ASYNCHRONOUS ::  R_buffer(buf_dim)
    integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! local variables

    integer :: lblk, lnumb, lb
    integer :: index
    integer :: ierrorcode,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    lnumb = commatrix_recv(sproc+1)
    if(lnumb.gt.maxblocks) then
       call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
    endif
    index = ir_buf(1,sproc+1)


    do lblk=1,lnumb
       lb = lblk + ladd_strt(sproc) - 1

       ! unpack all arrays from buffer into lb
#ifdef DEBUG
       write(*,*) 'pe ',mype,' lblk ',lblk,' unpacking starting ', &
            &        ' at index ',index,' buf_dim ',buf_dim
#endif /* DEBUG */
       call mpi_put_buffer( &
            &         lb,iopt,index,lcc,lfc,lec,lnc,buf_dim,R_buffer, &
            &         nlayersx,nlayersy,nlayersz)
#ifdef DEBUG
       write(*,*) 'pe ',mype,' lblk ',lblk,' unpacked into ',lb
#endif /* DEBUG */

    enddo

#ifdef DEBUG
    if (index .gt. buf_dim) then
       print *,' '
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' Rbuffer is too small by ',index-size(R_buffer,dim=1)
       print *,' index = ',index,' buf_dim = ',buf_dim,' @',mype
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' '
    end if
    if (index .eq. buf_dim) then
       print *,' '
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' End of Rbuffer reached exactly by pmMpiUnpackBlksFromProc!'
       print *,' index = ',index,' buf_dim = ',buf_dim,' @',mype
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' '
    end if
    if (index .ne. ir_buf(2,sproc+1)+1) then
       print *,' '
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' Rbuffer mismatch by ',&
            index - (ir_buf(2,sproc+1)+1),' @',mype,'  !'
       print *,' index = ',index,' buf_dim = ',buf_dim,' last', &
            ir_buf(2,sproc+1)
       print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       print *,' '
    end if
#endif
  end subroutine pmMpiUnpackBlksFromProc

end module gr_pmBlockGetter
