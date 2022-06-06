#if N_DIM > 1
#define DEBUG_DAT
#endif
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/gr_mpiAmrComm
!! NAME
!!
!!   gr_mpiAmrComm
!! 
!! SYNOPSIS
!!
!!   Call gr_mpiAmrComm(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset,
!!                           stages,
!!                           nlayersx,nlayersy,nlayersz,
!!                           flux_dir)
!!   Call gr_mpiAmrComm(mype,nprocs,
!!                           lguard,lprolong,
!!                           lflux,ledge,lrestrict,lfulltree,
!!                           iopt,lcc,lfc,lec,lnc,tag_offset)
!!
!!   Call gr_mpiAmrComm(integer, integer,
!!                           logical, logical,
!!                           logical, logical, logical, logical,
!!                           integer, logical, logical, logical, logical, integer
!!                           optional integer,
!!                           optional integer, optional integer, optional, integer,
!!                           optional integer)
!!
!! ARGUMENTS      
!!
!!    integer, intent(in) :: mype,nprocs,iopt
!!      mype   -> The local processor number.
!!      nprocs -> The total number of processors being used.
!!       
!!   logical, intent(in)  :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
!!      lguard    -> If true, info for guardcell filling will be
!!                    communicated.
!!      lprolong  -> If true, info for prolongation will be
!!                    communicated.
!!      lflux     -> If true, info for flux conservation will be
!!                    communicated.
!!      ledge     -> If true, info for edge averaging will be
!!                    communicated.
!!      lrestrict -> If true, info for restriction will be
!!                    communicated.
!!      lfulltree -> If true, info for restricting data up the entire
!!                    tree will be communicated.
!!
!!   integer, intent(in) :: iopt
!!      iopt   -> Controls which arrays are operated on:
!!                If iopt = 1 then 'unk', 'facevarx(y,z)', 'unk_e_x(y,z)'
!!                and 'unk_n'
!!                If iopt >= 2 then 'work'.
!!
!!   logical, intent(in) :: lcc,lfc,lec,lnc
!!      lcc -> A logical switch controlling whether unk or work data
!!              is filled.
!!      lfc -> A logical switch controlling whether facevar data
!!              is filled.
!!      lec -> A logical switch controlling whether unk_e_? data
!!              is filled.
!!      lnc -> A logical switch controlling whether unk_n data
!!              is filled.
!!
!!   integer, intent(inout) :: tag_offset
!!      tag_offset -> Defines the last tag number used for an mpi message.
!!                    This can be almost anything, but 100 is usually a 
!!                    good choice.
!!
!!   integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
!!     Optional integer arguments specifying how many guardcells to 
!!      exchange in each coordinate direction.  If these are not 
!!      specified then all guardcells are filled.
!!
!!   integer, intent(in), optional :: flux_dir
!!     Optional argument specifying which direction to operate on 
!!       for flux conservation. 
!!       If flux_dir = 1, operate on 'x' direction.
!!       If flux_dir = 2, operate on 'y' direction.
!!       If flux_dir = 3, operate on 'z' direction.
!!     If not specified, all directions are operated on.
!!
!!   stages - a mask of bitflags.
!!            Tells this routine which parts of the work of
!!            communcating it shoudl actually do.
!!            Bitwise OR of the following integer values:
!!               1 (sendBit)
!!               2 (postRecvBit)
!!               4 waitBit
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   mpi_morton
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   mpi_pack_blocks
!!   mpi_Sbuffer_size
!!   mpi_unpack_blocks
!!   mpi_Rbuffer_size
!!   mpi_pack_edges
!!   mpi_unpack_edges
!!   mpi_pack_fluxes
!!   mpi_unpack_fluxes
!!   mpi_amr_read_guard_comm
!!   mpi_amr_read_prol_comm
!!   mpi_amr_read_flux_comm
!!   mpi_amr_read_restrict_comm
!!   mpi_set_message_sizes
!!   mpi_xchange_blocks
!!
!! RETURNS
!!
!!   Upon exit, all necessary block data has been communicated to the local
!!   Calling processor for the operation selected through the specified
!!   arguments.
!!
!! DESCRIPTION
!!
!!   This routine is the pre-analysis operation which precedes a code block 
!!   in which interprocessor communications may be required.  What is does
!!   is set up and communicate the necessary block information for performing
!!   various communications tasks.  Those tasks are:
!!     1) guardcell filling
!!     2) prolongation
!!     3) restriction
!!     4) flux conservation
!!     5) edge averaging
!!
!! AUTHORS
!!
!!   Peter MacNeice (June 2000) with modifications by Kevin Olson for
!!   directional guardcell filling and flux conservation.
!!
!! MODIFICATIONS
!!  2020-12-11 K. Weide  Created from mpi_amr_comm_setup for async dd comms
!!***

#include "paramesh_preprocessor.fh"
#include "constants.h"
Module gr_mpiAmrComm_mod

contains
      Subroutine gr_mpiAmrComm(mype,nprocs,                       & 
                                    lguard,lprolong,                   & 
                                    lflux,ledge,lrestrict,lfulltree,   & 
                                    iopt,lcc,lfc,lec,lnc,tag_offset,   &
                                    getter,                            &
!!$                                    stages,                            &
                                    ntype,level,                       &
                                    nlayersx,nlayersy,nlayersz,        & 
                                    flux_dir)
!#define DEBUG
!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree
      Use mpi_morton, ONLY: temprecv_buf,                              &
                            ladd_strt, ladd_end
#ifdef DEBUG
      Use mpi_morton, ONLY: commatrix_send, commatrix_recv, max_no_to_send
#endif
      Use paramesh_mpi_interfaces, Only :                              & 
                                      mpi_pack_blocks,                 & 
                                      mpi_Sbuffer_size,                & 
                                      mpi_unpack_blocks,               & 
                                      mpi_Rbuffer_size,                & 
                                      mpi_pack_edges,                  & 
                                      mpi_unpack_edges,                & 
                                      mpi_pack_fluxes,                 & 
                                      mpi_unpack_fluxes,               & 
                                      mpi_amr_read_guard_comm,         & 
                                      mpi_amr_read_prol_comm,          & 
                                      mpi_amr_read_flux_comm,          & 
                                      mpi_amr_read_restrict_comm,      & 
                                      mpi_set_message_sizes,           & 
                                      mpi_xchange_blocks
      use gr_parameshInterface, ONLY : gr_mpiXchangeBlocks
#ifdef DEBUG
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
#endif
      use gr_pmCommDataTypes, ONLY: GRID_PAT_GC,GRID_PAT_PROLONG,GRID_PAT_FCORR,GRID_PAT_RESTRICT
      use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t, gr_pmCommCtl_t
      use gr_pmBlockGetter,   ONLY: gr_pmBlockGetter_t, gr_pmBlockGetterBuild

#ifdef DEBUG_DAT
      use Logfile_interface, ONLY: Logfile_stamp
#endif

      Implicit None

!-----Include statements.
! #ifdef DEBUG  ! currently also used for MPI_STATUS_SIZE
      Include 'mpif.h'
! #endif

!-----Input/Output arguments.
      Integer, Intent(in)    :: mype,nprocs,iopt
      Integer, Intent(inout) :: tag_offset
      Logical, Intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      Logical, Intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      type(gr_pmBlockGetter_t), intent(OUT), OPTIONAL, TARGET :: getter
!!$      integer, intent(in), optional :: stages
      integer, intent(IN), optional :: ntype
      integer, intent(IN), optional :: level
      Integer, Intent(in), Optional :: nlayersx,nlayersy,nlayersz
      Integer, Intent(in), Optional :: flux_dir

      integer,parameter :: sendBit=1, postRecvBit=2, waitBit=4

!-----Local variables and arrays.
      Integer :: buffer_dim_send, buffer_dim_recv
      Integer :: len_surr_blks
      Integer :: offset_tree
#ifdef DEBUG
      Integer :: max_blks_sent
      Integer :: itemp
      Integer :: ierror
#endif
!!$      integer :: myStages
      Integer :: nlayerstx, nlayersty, nlayerstz
      integer :: ntypeLoc, lev, patFam
      Integer :: flux_dirt
      Integer :: ii,jj
      integer :: lb

      logical :: doSend, doPostRecv, doWait

      type(gr_pmCommPattern_t),SAVE :: commPat !NOT YET USED
      type(gr_pmCommCtl_t),SAVE,TARGET :: commCtl     ! This one is used.
      type(gr_pmCommCtl_t),POINTER :: pCommCtl
!!$      type(gr_pmBlockGetter_t) :: myGetter
!!$      type(gr_pmBlockGetter_t),POINTER :: pGetter

      ASYNCHRONOUS :: temprecv_buf

#ifdef DEBUG_DAT
      character(len=32), dimension(2,2) :: block_buff
      character(len=32)                 :: int_to_str
      integer, save :: lastWritten_buffer_dim_send = -1
      integer, save :: lastWritten_buffer_dim_recv = -1
#endif

!-----Begin executable code.

! The actions done here can be summarized into the following stages, in order:
! (omitting debugging stuff, minor things like handling optional arguments)
! * Guard cell mask handling (using gcell_on_?c if (lguard))
! * Set message sizes
! * Swap in the right communcation pattern
! * Compute sizes of send_buf and temprecv_buf (rank 0: sometimes log them)
! * Allocate temprecv_buf
! * CASE A (using getter)
!   * Initialize/allocate local commCt object
! * CASE B (using getter)
!   * Build the GETTER (which includes a commCtl subobject)
! * Allocate commCtl % sendBuf
! * Pack data into commCtl % sendBuf
! * Call gr_mpiXchangeBlocks(stages=3): Communicate, part1: POST sends & receives
! * Set ladd_strt, ladd_end from comm pattern's strt_buff, laddress
! * CASE A (traditional)
!   * Call gr_mpiXchangeBlocks(stages=4): Finish Communicating!
!   * Unpack metainfo for blocks (call mpi_unpack_{blocks,fluxes,edges})
!   * Cleanup/deallocate commCtl object
! * CASE B (using getter)
!   * Return the GETTER to caller
! In CASE B, caller will do something like
!     * DO UNTIL EXHAUSTED:
!       * GETTER % GET()    -  note: unpacking metainfo happens inside
!     * Destroy GETTER (This should include cleanup/deallocate of commCtl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$      if (present(stages)) then
!!$         myStages = stages
!!$      else
!!$         myStages = 7           !all in
!!$      end if
!!$      doSend     = IAND(myStages, sendBit    ) /= 0
!!$      doPostRecv = IAND(myStages, postRecvBit) /= 0
!!$      doWait     = IAND(myStages, waitBit    ) /= 0


#ifdef DEBUG
      write(*,*) 'pe ',mype,' entered gr_mpiAmrComm: ' & 
     &           ,' max_no_to_send ', &
     &           max_no_to_send,' tag_offset ',tag_offset, &
     &           ' nprocs ',nprocs, & 
     &           '  gcell_on_cc ', gcell_on_cc,' iopt ',iopt
      print*,'@',mype,'lguard,lrestrict,lfulltree,lprolong,lflux,ledge:',&
                       lguard,lrestrict,lfulltree,lprolong,lflux,ledge

#endif /* DEBUG */
      If (iopt == 1) Then



!-----install user selection for guardcell variables on reset defaults.
      If (lguard) Then
        int_gcell_on_cc = gcell_on_cc
        int_gcell_on_fc = gcell_on_fc
        int_gcell_on_ec = gcell_on_ec
        int_gcell_on_nc = gcell_on_nc
        lguard_in_progress = .True.

        If (nvar > 0) Then
          ngcell_on_cc = 0
          Do ii = 1,nvar
            If (int_gcell_on_cc(ii)) Then
              ngcell_on_cc =  ngcell_on_cc + 1
              gcell_on_cc_pointer(ngcell_on_cc) = ii
            End If
          End Do
        End If

        If (nfacevar > 0) Then
          ngcell_on_fc = 0
          Do ii = 1,nfacevar
          Do jj = 1,3
            If (int_gcell_on_fc(jj,ii)) Then
              ngcell_on_fc(jj) =  ngcell_on_fc(jj) + 1
              gcell_on_fc_pointer(jj,ngcell_on_fc(jj)) = ii
            End If
          End Do
          End Do
        End If

        If (nvaredge > 0) Then
          ngcell_on_ec = 0
          Do ii = 1,nvaredge
          Do jj = 1,3
            If (int_gcell_on_ec(jj,ii)) Then
              ngcell_on_ec(jj) =  ngcell_on_ec(jj) + 1
              gcell_on_ec_pointer(jj,ngcell_on_ec(jj)) = ii
            End If
          End Do
          End Do
        End If

        If (nvarcorn > 0) Then
          ngcell_on_nc = 0
          Do ii = 1,nvarcorn
            If (int_gcell_on_nc(ii)) Then
              ngcell_on_nc =  ngcell_on_nc + 1
              gcell_on_nc_pointer(ngcell_on_nc) = ii
            End If
          End Do
        End If

      Else                                

        int_gcell_on_cc(:) = .True.
        int_gcell_on_fc(:,:) = .True.
        int_gcell_on_ec(:,:) = .True.
        int_gcell_on_nc(:) = .True.
        lguard_in_progress = .False.
        ngcell_on_cc = nvar
        Do ii=1,nvar
          gcell_on_cc_pointer(ii) = ii
        End Do
        ngcell_on_fc = nfacevar
        Do ii=1,nfacevar
          gcell_on_fc_pointer(:,ii) = ii
        End Do
        ngcell_on_ec = nvaredge
        Do ii=1,nvaredge
          gcell_on_ec_pointer(:,ii) = ii
        End Do
        ngcell_on_nc = nvarcorn
        Do ii=1,nvarcorn
          gcell_on_nc_pointer(ii) = ii
        End Do

      End If  ! End If (lguard)

      If (Present(nlayersx)) Then
         nlayerstx = nlayersx
      Else
         nlayerstx = nguard
      End If

      If (Present(nlayersy)) Then
         nlayersty = nlayersy
      Else
         nlayersty = nguard
      End If

      If (Present(nlayersz)) Then
         nlayerstz = nlayersz
      Else
         nlayerstz = nguard
      End If

      Else

      If (Present(nlayersx)) Then
         nlayerstx = nlayersx
      Else
         nlayerstx = nguard_work
      End If

      If (Present(nlayersy)) Then
         nlayersty = nlayersy
      Else
         nlayersty = nguard_work
      End If

      If (Present(nlayersz)) Then
         nlayerstz = nlayersz
      Else
         nlayerstz = nguard_work
      End If

      End If  ! End If (iopt == 1)

      If (lguard) Then
         If (nxb/nguard < 2) nlayerstx = min(nlayerstx+1,   nguard)
         If (nyb/nguard < 2) nlayersty = min(nlayersty+k2d, nguard)
         If (nzb/nguard < 2) nlayerstz = min(nlayerstz+k3d, nguard)
      End If
      
      If (Present(flux_dir)) Then
         flux_dirt = flux_dir
      Else
         flux_dirt = 0
      End If

      Call mpi_set_message_sizes(iopt,nlayerstx,nlayersty,nlayerstz)


      If (lguard.and.(.not.lrestrict) .or. lfulltree ) Then

#ifdef DEBUG
         print*,'@',mype,' will now Call mpi_amr_read_guard_comm(nprocs) from gr_mpiAmrComm.'
#endif
         patFam = GRID_PAT_GC
         Call mpi_amr_read_guard_comm(nprocs)

      ElseIf (lprolong) Then

         patFam = GRID_PAT_PROLONG
         Call mpi_amr_read_prol_comm(nprocs)

      ElseIf ((lflux.or.ledge).and.(.not.lrestrict)) Then

         patFam = GRID_PAT_FCORR
        Call mpi_amr_read_flux_comm(nprocs)

      ElseIf (lrestrict) Then

#ifdef DEBUG
         print*,'@',mype,' will now Call mpi_amr_read_restrict_comm(nprocs) from gr_mpiAmrComm.'
#endif
         patFam = GRID_PAT_RESTRICT
        Call mpi_amr_read_restrict_comm(nprocs)

      End If

#ifdef DEBUG
      print*,'@',mype,' has commatrix_send=',commatrix_send
      print*,'@',mype,' has commatrix_recv=',commatrix_recv
      itemp = max(sum(commatrix_send), sum(commatrix_recv))
      Call MPI_ALLREDUCE (itemp,                                       & 
                          max_blks_sent,                               & 
                          1,                                           & 
                          MPI_INTEGER,                                 & 
                          MPI_MAX,                                     & 
                          amr_mpi_meshComm,                              & 
                          ierror)
      write(*,*) 'pe ',mype,' comm pattern in gr_mpiAmrComm: ' &
                 ,' family ',patFam                            &
     &           ,' max_blks_sent ',max_blks_sent
#endif

!-----Buffer space is always needed for data and tree information.

      offset_tree = -1          ! will not be used in following calls
      If (lguard.or.lprolong) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .True., .False., .False., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .True.,.False., .False., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .False., .True., .False., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .False.,.True., .False., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      ElseIf (ledge) Then

        Call mpi_Sbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_send,offset_tree,             & 
                              .False., .False., .True., flux_dir,      &
                              nlayerstx,nlayersty,nlayerstz)

        Call mpi_Rbuffer_size(mype,nprocs,iopt,lcc,lfc,lec,lnc,        & 
                              buffer_dim_recv,                         & 
                              .False.,.False., .True., flux_dir,       &
                              nlayerstx,nlayersty,nlayerstz)

      End If  ! End If (lguard.or.lprolong)

#ifdef DEBUG_DAT
      if (lguard .and. .not. lrestrict) then
         if (mype .eq. 0) then  ! Only do this for MASTER_PE
            if ((lastWritten_buffer_dim_send .NE. buffer_dim_send) .OR. & 
     &           (lastWritten_buffer_dim_recv   .NE. buffer_dim_recv)) then
               write (block_buff(1,1), '(a)') 'buffer_dim_send'
               write (int_to_str, '(i11,a1)') buffer_dim_send, ','
               write (block_buff(1,2),'(a)')  trim(adjustl(int_to_str))

               write (block_buff(2,1), '(a)') 'buffer_dim_recv'
               write (int_to_str, '(i11)') buffer_dim_recv
               write (block_buff(2,2), '(a)') trim(adjustl(int_to_str))
               
               call Logfile_stamp( block_buff, 2, 2, & 
     &              '[gr_mpiAmrComm]')
               lastWritten_buffer_dim_send = buffer_dim_send
               lastWritten_buffer_dim_recv = buffer_dim_recv
            end if
         end if
      end if
#endif

!-----Set up buffer size info, including space for tree data
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
      offset_tree = 32+16+len_surr_blks

      If (Allocated(temprecv_buf)) Deallocate(temprecv_buf)
      Allocate( temprecv_buf(buffer_dim_recv))

!!$      Call mpi_xchange_blocks(mype,nprocs,tag_offset,                  &
!!$                              buffer_dim_send,send_buf,                &
!!$                              buffer_dim_recv,temprecv_buf)
      !! DEMOSTRATION how completion of communication can be deferred:
      IF(PRESENT(getter)) THEN
         if (present(ntype)) then
            ntypeLoc = ntype
         else if (lguard.and.(.not.lrestrict) .or. lfulltree ) Then
            ntypeLoc = ACTIVE_BLKS
         ElseIf (lprolong) Then
            ntypeLoc = ACTIVE_BLKS
         ElseIf ((lflux.or.ledge).and.(.not.lrestrict)) Then
            ntypeLoc = LEAF
         ElseIf (lrestrict) Then
            ntypeLoc = PARENT_BLK
         else
            ntypeLoc = ACTIVE_BLKS
         end if
         if (present(level)) then
            lev = level
         else if (lguard.and.(.not.lrestrict) .or. lfulltree ) Then
            lev = UNSPEC_LEVEL
         ElseIf (lprolong) Then
            lev = UNSPEC_LEVEL
         ElseIf ((lflux.or.ledge).and.(.not.lrestrict)) Then
            lev = UNSPEC_LEVEL
         ElseIf (lrestrict) Then
            lev = UNSPEC_LEVEL
         else
            lev = UNSPEC_LEVEL
         end if
         ! Build a getter and leave completion of communications to it
         call gr_pmBlockGetterBuild(getter,ntypeLoc,lev,&
                                    patFam,          &
                                    iopt,            &
     &                              lcc,lfc,lec,lnc, &
     &                              buffer_dim_recv, &
     &                              nlayerstx,nlayersty,nlayerstz)
         pCommCtl => getter % commCtl
!!$         print*,'@',mype,' has pCommCtl % numReq =',pCommCtl % numReq,' after call gr_pmBlockGetterBuild',ntypeLoc,lev &
!!$                 ,' family ',patFam
      ELSE
         commCtl          % numReq = 0
         commCtl          % allReceivedCount = 0
         allocate(commCtl % recvrequest(nprocs))
         allocate(commCtl % recvstatus (MPI_STATUS_SIZE,nprocs)) ! We could use MPI_STATUSES_IGNORE instead.
         allocate(commCtl % receivedIndices(nprocs))
         commCtl          % numSReq = 0
         commCtl          % allSentCount = 0
         allocate(commCtl % sendRequest(nprocs))
         allocate(commCtl % sendStatus (MPI_STATUS_SIZE,nprocs)) ! We could use MPI_STATUSES_IGNORE instead.
         pCommCtl => commCtl
      END IF

      If (Allocated(pCommCtl % sendBuf)) Deallocate(pCommCtl % sendBuf)
      Allocate( pCommCtl % sendBuf(buffer_dim_send))

      If (lguard.or.lprolong) Then

        Call mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc,         &
                             buffer_dim_send,pCommCtl % sendBuf,offset_tree,     &
                             nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then

        Call mpi_pack_fluxes(mype,nprocs,buffer_dim_send,pCommCtl % sendBuf,     &
                             offset_tree,flux_dir)

      ElseIf (ledge) Then

        Call mpi_pack_edges(mype,nprocs,buffer_dim_send,pCommCtl % sendBuf,      &
                            offset_tree)

      End If  ! End If (lguard.or.lprolong)

      ! Post receives and sends:
      Call gr_mpiXchangeBlocks(mype,nprocs,tag_offset,                 &
                              buffer_dim_send,pCommCtl % sendBuf,                &
                              buffer_dim_recv,temprecv_buf,            &
                              3,                                       &
                              pCommCtl)

!-----Store laddress pe limits to help optimize searching
      If (last_buffer > strt_buffer.and.nprocs > 1) Then

         ladd_strt(1:nprocs-1) = last_buffer
         ladd_strt(0         ) = strt_buffer
         Do jj= last_buffer,strt_buffer,-1
           ladd_strt(laddress(2,jj)) = jj
         End Do
         Do jj= nprocs-2,0,-1
           ladd_strt(jj) = minval(ladd_strt(jj:jj+1))
         End Do

         ladd_end(0:nprocs-2) = strt_buffer
         ladd_end(nprocs-1  ) = last_buffer
         Do jj= strt_buffer,last_buffer
           ladd_end(laddress(2,jj)) = jj
         End Do
         Do jj= 1,nprocs-1
           ladd_end(jj) = maxval(ladd_end(jj-1:jj))
         End Do

      Else

        ladd_strt = strt_buffer
        ladd_end  = strt_buffer

      End If

      IF(PRESENT(getter)) THEN
         ! A getter was built above.
         ! We shall leave completion of communications to it.
!!$         print*,'@',mype,' has pCommCtl % numReq =',pCommCtl % numReq,' after call gr_mpiXchangeBlocks'

         if (.NOT.(lguard.or.lprolong)) then
            PRINT*,'GETTER USE NOT YET SUPPORTED!'
         end if
         GOTO 999
!         RETURN                 ! RETURN NOW!
      else
!!$         print*,'@',mype,' has pCommCtl % numReq =',pCommCtl % numReq,' after call gr_mpiXchangeBlocks w/o getter'
      END IF
         
      ! Wait for completions:
      Call gr_mpiXchangeBlocks(mype,nprocs,tag_offset,                 &
                                  buffer_dim_send,commCtl % sendBuf,       &
                                  buffer_dim_recv,temprecv_buf,            &
                                  4,                                       &
                                  commCtl)


      If (lguard.or.lprolong) Then

         Call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc,             &
                                buffer_dim_recv,temprecv_buf,          &
                                nlayerstx,nlayersty,nlayerstz)

      ElseIf (lflux) Then
         
         Call mpi_unpack_fluxes(mype,buffer_dim_recv,temprecv_buf,     & 
                                flux_dir)

      ElseIf (ledge) Then

         Call mpi_unpack_edges(mype,buffer_dim_recv,temprecv_buf)

      End If  ! End If (lguard.or.lprolong)

      deallocate(commCtl % recvrequest)
      deallocate(commCtl % recvstatus )
      deallocate(commCtl % receivedIndices)

      deallocate(commCtl % sendRequest)
      deallocate(commCtl % sendStatus )

      If (Allocated(commCtl % sendBuf)) Deallocate(commCtl % sendBuf)

      ! NOTE: label 999 was moved HERE since this routine was changed to NOT wait for all
      ! outstanding MPI_Isend's to complete before returning. Otherwise send_buf may get
      ! deallocated while still being in use!
999   continue


      Return
      End Subroutine gr_mpiAmrComm

End Module gr_mpiAmrComm_mod

