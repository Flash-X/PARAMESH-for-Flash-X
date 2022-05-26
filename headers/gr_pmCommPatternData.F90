!!****ih* headers/gr_pmCommPatternData
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
!!   gr_pmCommPatternData
!!
!! SYNOPSIS
!!
!!   use gr_pmCommPatternData
!!
!! USES
!!
!!   gr_pmCommDataTypes ,ONLY: gr_pmCommPattern_t
!!   tree
!!
!! DESCRIPTION
!!
!!   Fortran module which holds data structures for managing
!!   Paramesh communication patterns.
!!
!! HISTORY
!!
!!  2022-05-13 K. Weide  created
!!  2022-05-20 K. Weide  added space for alternative comm patterns
!!  2022-05-25 K. Weide  added space for an additional restriction subpattern
!!***

module gr_pmCommPatternData
  use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t

  implicit none
  integer,parameter :: NUM_COMM_PATTERN_FAMILIES = 4

  ! At least the following variants could be distinguished based on
  ! node type (not specific refinement level):
  !  DEFAULT; LEAF; ACTIVE_BLKS; PARENT; ALL_BLKS
  integer,parameter,dimension(NUM_COMM_PATTERN_FAMILIES) :: &
       NUM_COMM_PATTERN_VARIANTS = (/2,1,1,3/)
  ! Unused now:
  integer,parameter :: MAX_COMM_PATTERN_VARIANTS = 3 ! maxval(NUM_COMM_PATTERN_VARIANTS)

  ! There is going to be one of the following for each pattern family:
  type gr_pmCommPattFamSet_t
     TYPE(gr_pmCommPattern_t),allocatable,dimension(:) :: pat
  end type gr_pmCommPattFamSet_t


  ! The global holder for all the comm patterns:
  TYPE(gr_pmCommPattFamSet_t),dimension(NUM_COMM_PATTERN_FAMILIES), &
       Target :: gr_pmTheCommPatterns

  TYPE(gr_pmCommPattern_t),POINTER :: gr_theActiveCommPattern => NULL()

contains
  subroutine gr_pmInvalidateCommPatterns() !Just resets 'valid' flags
    integer :: f,p

    do f = 1,NUM_COMM_PATTERN_FAMILIES
       if(allocated(gr_pmTheCommPatterns(f) % pat)) then
          do p = 1,size(gr_pmTheCommPatterns(f) % pat, 1)

             gr_pmTheCommPatterns(f) % pat(p) % valid = .FALSE.
          end do
       end if
    end do
  end subroutine gr_pmInvalidateCommPatterns

  subroutine gr_pmResetCommPatterns()
    use tree, ONLY: laddress
    integer :: f,p

    do f = 1,NUM_COMM_PATTERN_FAMILIES
       if(allocated(gr_pmTheCommPatterns(f) % pat)) then
          do p = 1,size(gr_pmTheCommPatterns(f) % pat, 1)

             if (associated(gr_theActiveCommPattern, &
                            gr_pmTheCommPatterns(f) % pat(p))) then
                nullify(gr_theActiveCommPattern)
             end if
             if(allocated(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_recv)) &
                  deallocate(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_recv)
             if(allocated(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_send)) &
                  deallocate(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_send)
             if(allocated(gr_pmTheCommPatterns(f) % pat(p) %  to_be_received)) &
                  deallocate(gr_pmTheCommPatterns(f) % pat(p) %  to_be_received)
             if(allocated(gr_pmTheCommPatterns(f) % pat(p) %  to_be_sent)) &
                  deallocate(gr_pmTheCommPatterns(f) % pat(p) %  to_be_sent)
             gr_pmTheCommPatterns(f) % pat(p) %  strt_buffer = -1
             if(associated(laddress, &
                  gr_pmTheCommPatterns(f) % pat(p) %  laddress)) then
                nullify(laddress)
             end if
             if(allocated(gr_pmTheCommPatterns(f) % pat(p) %  laddress)) &
                  deallocate(gr_pmTheCommPatterns(f) % pat(p) %  laddress)
             gr_pmTheCommPatterns(f) % pat(p) % valid = .FALSE.
          end do
       end if
    end do
  end subroutine gr_pmResetCommPatterns

  subroutine gr_pmDeallocateCommPatterns()
    integer :: f

    call gr_pmResetCommPatterns()
    do f = 1,NUM_COMM_PATTERN_FAMILIES
       if(allocated(gr_pmTheCommPatterns(f) % pat)) then
          deallocate(gr_pmTheCommPatterns(f) % pat)
       end if
    end do
  end subroutine gr_pmDeallocateCommPatterns

  subroutine gr_pmInitCommPatterns(nprocs,maxBlocksAlloc)
    integer,intent(in) :: nprocs, maxblocksAlloc
    integer :: f,p

    call gr_pmResetCommPatterns()
    do f = 1,NUM_COMM_PATTERN_FAMILIES
       allocate(gr_pmTheCommPatterns(f) % pat(NUM_COMM_PATTERN_VARIANTS(f)))
       do p = 1,size(gr_pmTheCommPatterns(f) % pat, 1)
          allocate(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_recv(nprocs))
          allocate(gr_pmTheCommPatterns(f) % pat(p) %  commatrix_send(nprocs))
          gr_pmTheCommPatterns(f) % pat(p) %  commatrix_recv(:) = 0 !- 999
          gr_pmTheCommPatterns(f) % pat(p) %  commatrix_send(:) = 0 !- 888
          ! to_be_received,to_be_received not initially allocated!
!!$          allocate(gr_pmTheCommPatterns(f) % pat(p) %  to_be_received)
!!$          allocate(gr_pmTheCommPatterns(f) % pat(p) %  to_be_sent)
          allocate(gr_pmTheCommPatterns(f) % pat(p) %  laddress(1:2,1:maxblocksAlloc))
          gr_pmTheCommPatterns(f) % pat(p) %  strt_buffer = -1
          gr_pmTheCommPatterns(f) % pat(p) % id = f * 10 + (p-1)
          gr_pmTheCommPatterns(f) % pat(p) % valid = .FALSE.
       end do
    end do
  end subroutine gr_pmInitCommPatterns

  function gr_pmCommPatternPtr(family,variant) result(ptr)
    TYPE(gr_pmCommPattern_t),POINTER :: ptr
    integer, intent(in) :: family
    integer, intent(in),OPTIONAL :: variant

    integer :: f,p
    nullify(ptr)

    if (family > 0 .AND. family <= NUM_COMM_PATTERN_FAMILIES) then
       f = family               !use small integer directly as is
    else
       f = family / 10          !need error checking...
    end if
    if (.NOT. present(variant)) then
       p = 1                    !DEFAULT variant; initially, that's
       !all there is.
    else
       p = variant              !need error checking...
    end if
    ptr => gr_pmTheCommPatterns(f) % pat(p)
  end function gr_pmCommPatternPtr

  subroutine gr_pmActivateCommPattern(family,variant)
    use tree, ONLY: laddress
    integer, intent(in) :: family
    integer, intent(in),OPTIONAL :: variant

    gr_theActiveCommPattern => gr_pmCommPatternPtr(family,variant)
    laddress => gr_theActiveCommPattern % laddress
  end subroutine gr_pmActivateCommPattern

#include "FortranLangFeatures.fh"
  subroutine gr_pmPrintCommPattern(pattern,str,i)
    TYPE(gr_pmCommPattern_t),POINTER_INTENT_IN :: pattern
    character(len=*),intent(in) :: str
    integer,intent(in) :: i
    print*,'> pe', i,' ',str
    if (.NOT. associated(pattern)) then
       print*,'(NULL pattern)'
       return
    end if

#define MAXPRINT 22
    if(allocated(pattern % commatrix_send)) then
       if (size(pattern % commatrix_send,1) <= MAXPRINT) then
          print*,'| commatrix_send:',pattern % commatrix_send
       else
          print*,'| commatrix_send:', &
               pattern % commatrix_send(1:MAXPRINT),&
               ' ... #',MAXPRINT,' of',&
               size(pattern % commatrix_send,1)
       end if
    else
       print*,'| (commatrix_send unallocated)'
    end if
    if(allocated(pattern % commatrix_recv)) then
       if (size(pattern % commatrix_recv,1) <= MAXPRINT) then
          print*,'| commatrix_recv:',pattern % commatrix_recv
       else
          print*,'| commatrix_recv:', &
               pattern % commatrix_recv(1:MAXPRINT),&
               ' ... #',MAXPRINT,' of',&
               size(pattern % commatrix_recv,1)
       end if
    else
       print*,'| (commatrix_recv unallocated)'
    end if

    if(allocated(pattern % to_be_sent)) then
       if (size(pattern % to_be_sent,3) <= MAXPRINT) then
          if (size(pattern % to_be_sent,2) <= MAXPRINT) then
             print*,'| 2be_s',pattern % to_be_sent
          else
             print*,'| 2be_s', &
                  pattern % to_be_sent(:,1:MAXPRINT,:),&
                  ' ... # (1..3,1..',&
                  MAXPRINT,',1..',size(pattern % to_be_sent,3),') of (',&
                  size(pattern % to_be_sent,1),&
                  size(pattern % to_be_sent,2),&
                  size(pattern % to_be_sent,3),')'
          end if
       else
          if (size(pattern % to_be_sent,2) <= MAXPRINT) then
             print*,'| 2be_s', &
                  pattern % to_be_sent(:,:,1:MAXPRINT),&
                  ' ... # (1..3,1..',size(pattern % to_be_sent,2),&
                  ',1..',MAXPRINT,') of (',&
                  size(pattern % to_be_sent,1),&
                  size(pattern % to_be_sent,2),&
                  size(pattern % to_be_sent,3),')'
          end if

       end if
    else
       print*,'| (to_be_sent unallocated)'
    end if
    if(allocated(pattern % to_be_received)) then
       if (size(pattern % to_be_received,3) <= MAXPRINT) then
          if (size(pattern % to_be_received,2) <= MAXPRINT) then
             print*,'| 2be_r',pattern % to_be_received
          else
             print*,'| 2be_r', &
                  pattern % to_be_received(:,1:MAXPRINT,:),&
                  ' ... # (1..3,1..',&
                  MAXPRINT,',1..',size(pattern % to_be_received,3),') of (',&
                  size(pattern % to_be_received,1),&
                  size(pattern % to_be_received,2),&
                  size(pattern % to_be_received,3),')'
          end if
       else
          if (size(pattern % to_be_received,2) <= MAXPRINT) then
             print*,'| 2be_r', &
                  pattern % to_be_received(:,:,1:MAXPRINT),&
                  ' ... # (1..3,1..',size(pattern % to_be_received,2),&
                  ',1..',MAXPRINT,') of (',&
                  size(pattern % to_be_received,1),&
                  size(pattern % to_be_received,2),&
                  size(pattern % to_be_received,3),')'
          end if

       end if
    else
       print*,'| (to_be_received unallocated)'
    end if

    if(allocated(pattern % laddress)) then
       if (size(pattern % laddress,2) <= MAXPRINT) then
          print*,'| ladrs',pattern % laddress
       else
          print*,'| ladrs', &
               pattern % laddress(:,&
               size(pattern % laddress,2)-(MAXPRINT)+1 &
              :size(pattern % laddress,2)),&
               ' ... last #',MAXPRINT,' of',&
               size(pattern % laddress,2)
       end if
    else
       print*,'| (laddress unallocated)'
    end if

    print*,'| strt_buffer', pattern % strt_buffer
    print*,'| num_recipient_pes', pattern % num_recipient_pes
    print*,'| id', pattern % id
  end subroutine gr_pmPrintCommPattern
end module gr_pmCommPatternData
