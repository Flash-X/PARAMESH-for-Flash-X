module gr_pmCommPatternData
  use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t

  implicit none
  integer,parameter :: NUM_COMM_PATTERN_FAMILIES = 4

  ! At least the following variants could be distinguished based on
  ! node type (not specific refinement level):
  !  DEFAULT; LEAF; ACTIVE_BLKS; PARENT; ALL_BLKS
  integer,parameter,dimension(NUM_COMM_PATTERN_FAMILIES) :: &
       NUM_COMM_PATTERN_VARIANTS = (/5,5,5,5/)
  integer,parameter :: MAX_COMM_PATTERN_VARIANTS = maxval(NUM_COMM_PATTERN_VARIANTS)

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
    integer :: f,p

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
end module gr_pmCommPatternData
