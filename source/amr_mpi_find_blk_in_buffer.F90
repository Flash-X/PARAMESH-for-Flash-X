#define DEBUG

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_mpi_find_blk_in_buffer
!! NAME
!!
!!   amr_mpi_find_blk_in_buffer
!!
!! SYNOPSIS
!!
!!    Call amr_mpi_find_blk_in_buffer(mype,remote_block,remote_pe,idest,
!!                                    dtype,index0,lfound)
!!    Call amr_mpi_find_blk_in_buffer(integer,integer,integer,integer,
!!                                    integer,integer,logical)
!!
!! ARGUMENTS
!!
!!   Integer, Intent(in)  :: mype           local processor
!!   Integer, Intent(in)  :: remote_block   remote block 
!!   Integer, Intent(in)  :: remote_pe      remote processor
!!   Integer, Intent(in)  :: idest
!!   Integer, Intent(out) :: dtype          message type - an integer between 1 and 27
!!                                          indicating the section of a blck contained in
!!                                          the message segment from (remote_block,remote_pe).
!!   Integer, Intent(out) :: index0         the address index0+1 is where you should start
!!                                           reading the data from this message.
!!   Logical, Intent(out) :: lfound         logical indicating if the block is found or not
!!    
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   Flashx_mpi_implicitNone.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   mpi_morton
!!
!! CALLS
!!
!!   amr_abort
!!
!! RETURNS
!!
!!   Returns message type in dtype.
!!   Returns starting address in buffer where to find the requested block 
!!    in the variable index0+1.
!!   Returns if the block is found or not in lfound.
!!
!! DESCRIPTION
!!
!!   This routine finds where data for a remote block with address
!!   (remote_block,remote_pe) is in the recv buffer. It returns
!!   the message type, dtype, and the address in recv of the data word
!!   preceeding the first word of real data.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          May 2001
!!
!! MODIFICATIONS
!!
!!  2021 Dec   K. Weide  minor tweaks
!!  2022-05-23 K. Weide  error and corner case handling
!!  2022-11-08 K. Weide  Added more output for debugging, by default disabled
!!
!!***

#include "paramesh_preprocessor.fh"

      Subroutine amr_mpi_find_blk_in_buffer(                           & 
              mype,remote_block,remote_pe,idest,dtype,index0,lfound)

!-----Use statements.
#ifdef DEBUG_LITE
      use Grid_data, ONLY: nprocs => gr_meshNumProcs
      use gr_pmCommPatternData, ONLY: gr_theActiveCommPattern
      Use physicaldata, ONLY: mpi_pattern_id
#endif
      Use tree, only : lnblocks, strt_buffer, laddress
      Use mpi_morton
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm
#ifdef DEBUG
      Use physicaldata, ONLY : mpi_pattern_id
#endif

!-----Implicit and Include statements.
#include "Flashx_mpi_implicitNone.fh"


!-----Input/Output arguments.
      Integer, Intent(in)  :: mype,remote_pe,remote_block,idest
      Integer, Intent(out) :: dtype,index0
      Logical, Intent(out) :: lfound

!-----Local arrays and variables.
      integer :: jseg,seg_no,iaddress,no_of_comms,jpe,jpe0
      integer :: ierrorcode,ierr,no_of_segments
      integer :: rem_pe,rem_blk,seg_offset
      integer :: iseg_no,jj
      logical :: llfound

      real :: dtypeReal

#if defined(DEBUG) || defined(DEBUG_LITE)
!!$      integer :: pe_source(nprocs)
      integer,ALLOCATABLE :: pe_source(:)
      intrinsic pack
#endif

!-----Begin executable code.

      If (remote_pe.ne.mype) Then
#ifdef DEBUG
        print*,'@',mype,' amr_mpi_find_blk_in_buffer: Case 1 for',remote_block, remote_pe
#endif
        rem_blk = remote_block
        rem_pe  = remote_pe

        llfound = .False.
        jj = ladd_strt(rem_pe)
        iseg_no = 0
        Do While(.Not.llfound.and.jj <= ladd_end(rem_pe))
          If (rem_blk == laddress(1,jj).and.                           & 
             rem_pe == laddress(2,jj)) Then
            llfound = .True.
            iseg_no = jj - strt_buffer + 1
          Else
            jj = jj+1
          End If
        End Do

      ElseIf (remote_pe == mype.and.remote_block > lnblocks) Then
#ifdef DEBUG
        print*,'@',mype,' amr_mpi_find_blk_in_buffer: Case 2 for',remote_block, remote_pe
#endif
        rem_blk = laddress(1,remote_block)
        rem_pe  = laddress(2,remote_block)
        iseg_no = remote_block - strt_buffer + 1
        llfound = .True.

      else    ! avoid compiler warnings
#ifdef DEBUG
        print*,'@',mype,' amr_mpi_find_blk_in_buffer: Case 3 for',remote_block, remote_pe
#endif
        rem_pe = mype
!!$        dtype  = 14
!!$        index0 = -1
!!$        lfound = .FALSE.

      End If  ! End If (remote_pe.ne.mype)

      If (rem_pe.ne.mype) Then

#ifdef DEBUG_LITE
      ASSOCIATE(p => gr_theActiveCommPattern)
!-------locate rem_pe in the list of sending processors
!!$        pe_source(:) = -1
        no_of_comms = count(p%commatrix_recv /= 0)
        pe_source = pack(p%commatrix_recv, p%commatrix_recv /= 0)
        jpe0 = 0
        do jj = 1,nprocs
           if (p%commatrix_recv(jj) > 0) then
              jpe0 = jpe0 + 1
              pe_source(jpe0) = jj
           end if
        end do
        lfound = .False.
        jpe  = 0
        jpe0 = 0
        Do While((.Not.lfound).and.(jpe0 < no_of_comms))
          jpe0 = jpe0+1
          If (rem_pe == pe_source(jpe0)-1) Then
            lfound = .True.
            jpe = jpe0
          End If
        End Do
!-------If rem_pe is not located stop with error message
        If (jpe == 0) Then
          If (idest == 2) return
          Write(*,*) 'Paramesh error : pe ',mype,                      &  
           ' pe address of required data is not in the list of ',      & 
           'communicating pes. ',                                      & 
           ' remote_block ',remote_block,                              & 
           ' remote_pe ',remote_pe,                                    & 
           ' rem_blk ',rem_blk,                                          & 
           ' rem_pe ',rem_pe,                                          & 
           ' laddress ',laddress
          Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
        End If
        no_of_segments = size(p%to_be_received,2)
        lfound = .False.
        jseg = 0
        seg_no = jseg
        Do While((.Not.lfound).and.(jseg < no_of_segments))
          jseg = jseg+1
          If (p%to_be_received(1,jseg,jpe) == rem_blk .and.              &
              (p%to_be_received(2,jseg,jpe) == -1     .OR.               &
               p%to_be_received(2,jseg,jpe) == rem_pe+1   ) ) Then
            lfound = .True.
            seg_no = jseg
          End If
        End Do
!-------Now compute where the list of segments from proc rem_pe actually begins
!-------in the complete list of message segments received on this processor
        seg_offset = 0
        If (rem_pe > 0) seg_offset = sum(p%commatrix_recv(1:rem_pe))
        seg_no = seg_no + seg_offset
#ifdef DEBUGX
        If (seg_no.ne.iseg_no) Then
          Write(*,*) 'seg_no and iseg_no are different ',              & 
                      seg_no,iseg_no
          Call amr_abort()
        End If
#endif /* DEBUGX */
      end ASSOCIATE
#endif /* DEBUG_LITE */

        seg_no = iseg_no
        lfound = llfound

!-------If the requested segment is not located stop with error message
        If (seg_no == 0) Then
#ifndef DEBUG
          If (idest == 2) return
#endif
#ifdef DEBUG_LITE
          ASSOCIATE(p => gr_theActiveCommPattern)
          Write(*,*) 'Paramesh error : ',                              & 
           'message segment required is not in the list of ',          & 
           'segments received.: proc ',mype,' to_be_received ',        & 
           p%to_be_received(:,:,jpe),' mpi_pattern_id ',mpi_pattern_id
          end ASSOCIATE
#else
          Write(*,*) 'Paramesh error : ',                              & 
           'message segment required is not in the list of ',          & 
           'segments received.'
#endif
          Call MPI_ABORT(amr_mpi_meshComm,ierrorcode,ierr)
        End If  ! End If (seg_no == 0)

!-------set start address for this segment in R_buffer
        iaddress = mess_segment_loc(seg_no)

!-------Read out message into appropriate part of unk1 or work1
        dtypeReal = temprecv_buf(iaddress+2)
        if (.NOT.(dtypeReal .GE. 14.0 - 0.5*(Real(3**(N_DIM)-1.0)) .AND. &
                  dtypeReal .LE. 14.0 + 0.5*(Real(3**(N_DIM)-1.0)))) then
#ifdef DEBUG_LITE
          ASSOCIATE(p => gr_theActiveCommPattern)
           Write(*,*) 'amr_mpi_find_blk_in_buffer: pe ',mype,          &
           ' message segment with unexpected dtype ',dtypeReal,        &
           ' remote_block ',remote_block,                              &
           ' remote_pe ',remote_pe,                                    &
           ' rem_blk ',rem_blk,                                        &
           ' rem_pe ',rem_pe,                                          &
           ' seg_no ',seg_no,' idest',idest,                           &
           ' pe_source ',pe_source,                                    &
           ' to_be_received ',                                         &
           p%to_be_received(:,:,:),' mpi_pattern_id ',mpi_pattern_id
           end ASSOCIATE
#endif
           dtype = 0
           index0 = iaddress+2
           lfound = .FALSE.
           If (idest == 2) return
           ierrorcode = 12345
           dtype = anint(dtypeReal)
           call MPI_ABORT(amr_mpi_meshComm,ierrorcode,ierr)
        end if
        dtype = anint(temprecv_buf(iaddress+2))

!-------We must Write the message into recv first in case it is larger
!-------than the range we actually need.

        index0 = iaddress+2

!-------NOTE: start address of data is now index0+1

      End If  ! If (rem_pe.ne.mype)

#ifdef DEBUG
      400 format(1x,'@ ',I0,50x,' remote: ',I3,'@',I0,' idest ',I0,' -- found ',&
                 L1,' dtype ',I0,' seg_no ',I0,' iaddr ',I0)!,3(1x,1P,G11.4))
      print 400,mype,rem_blk,rem_pe,idest,lfound,dtype,iaddress!,temprecv_buf(iaddress:iaddress+2)
      print*,temprecv_buf(iaddress:iaddress+2)
#endif
      Return
      End Subroutine amr_mpi_find_blk_in_buffer
! Local Variables:
! f90-do-indent: 2
! f90-if-indent: 2
! indent-tabs-mode: nil
! End:
