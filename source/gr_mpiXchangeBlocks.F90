subroutine gr_mpiXchangeBlocks(mype,nprocs, tag_offset, &
                              buf_dim_send, S_buffer, &
                              buf_dim_recv, R_buffer, &
                              stages,                 &
                              commCtl,                &
                              commPat,                &
                              commShaped)

!------------------------------------------------------------------------
!
! This routine uses the general hand-shaking scheme for message passing 
! between pe's. 
! sends and receives buffers containing all requested block data.
!
!
! Written: (mpi_xchange_blocks) Maharaj Bhat & Michael Gehmeyr   March 2000
! Modified:  Klaus Weide                                         December 2020
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      buf_dim_send        dimension of send buffers
!      buf_dim_recv        dimension of recv buffers
!      tag_offset     offset for MPI tag
!      S_buffer       send buffer
!      R_buffer       recv buffer
!      stages         controls which stage(s) to do
!
!------------------------------------------------------------------------

  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  use mpi_morton, ONLY: commatrix_send, commatrix_recv, &
                        is_buf,         ir_buf
  Use paramesh_comm_data, ONLY: amr_mpi_real, amr_mpi_meshComm
  use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t, gr_pmCommCtl_t, &
                                gr_pmCommShaped_t

#include "Flashx_mpi_implicitNone.fh"

  integer, intent(in)    :: mype,nprocs,buf_dim_send,buf_dim_recv
  integer, intent(inout) :: tag_offset
  real,    intent(in)    :: S_buffer(buf_dim_send)
  real,    intent(inout),ASYNCHRONOUS :: R_buffer(buf_dim_recv)
  integer, intent(in)    :: stages
  type(gr_pmCommPattern_t), intent(in),OPTIONAL :: commPat
  type(gr_pmCommShaped_t),  intent(in),OPTIONAL :: commShaped
  type(gr_pmCommCtl_t), intent(INOUT)     :: commCtl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! local variables

  integer :: i,j
  integer :: istrt,ilast,isize
  integer :: isrc,idest, itag
  integer :: recvrequest(nprocs), sendRequest(nprocs)
  integer :: recvstatus(MPI_STATUS_SIZE,nprocs)
      integer :: sendStatus(MPI_STATUS_SIZE,nprocs)
  integer :: ierrorcode,ierr
  integer :: ij, ji
  integer :: receivedCount, sentCount
  logical :: doSend, doPostRecv, doWait

  integer,parameter :: sendBit=1, postRecvBit=2, waitBit=4

  intrinsic IAND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!      R_buffer = 0

#ifdef USEBARS
  Call MPI_BARRIER(amr_mpi_meshComm, ierr)
#endif

  doSend     = IAND(stages, sendBit    ) /= 0
  doPostRecv = IAND(stages, postRecvBit) /= 0
  doWait     = IAND(stages, waitBit    ) /= 0

  ij = commCtl % numReq
  ji = commCtl % numSReq

!!$  associate(commatrix_send => commPat % commatrix_send, &
!!$            commatrix_recv => commPat % commatrix_recv, &
!!$            is_buf         => commShaped % is_buf,         &
!!$            ir_buf         => commShaped % ir_buf          )
  associate(recvrequest => commCtl % recvrequest, &
            recvstatus  => commCtl % recvstatus,  &
            receivedInd => commCtl % receivedIndices, &
            allReceivedCount => commCtl % allReceivedCount, &
            sendRequest => commCtl % sendRequest, &
            sendStatus  => commCtl % sendStatus,  &
            allSentCount => commCtl % allSentCount &
            )

    if (commatrix_send(mype+1) > 0  &
         &    .or.  &
         &    commatrix_recv(mype+1) > 0) then
       write(*,*) 'Paramesh error :  error in xchange : pe ',mype, &
            &      ' diagonal element of commatrix is non-zero ', &
            &      commatrix_recv(mype+1), &
            &      commatrix_send(mype+1)
       call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
    endif

    if (doPostRecv) then
       ij = 0

       do i = 1,nprocs
          isrc = i-1
          idest= mype
#ifdef PM_UNIQUE_MPI_TAGS
          itag = isrc*nprocs + idest+1 + tag_offset
#else
          itag = tag_offset !OK because 1 msg per process.
#endif
          ! receive to pe=j
          if(commatrix_recv(i).gt.0) then
             ij = ij+1
             istrt = ir_buf(1,i)
             ilast = ir_buf(2,i)
             isize = ilast-istrt+1

             if (amr_error_checking) then
                if(isize.gt.buf_dim_recv) then
                   write(*,*) 'PARAMESH ERROR : gr_mpiXchangeBlocks 1:', &
                        &                    ' message is bigger than buffer space  ' &
                        &           ,' isrc ',isrc,' idest ',idest,' i ',i, &
                        &           ' istrt ',istrt,' ilast ',ilast,' isize ',isize, &
                        &           ' buf_dim_recv ',buf_dim_recv
                   call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
                endif
             endif
             R_buffer(istrt:istrt+2) = 0.0 ! Zero out the first few entries.
             call Mpi_Irecv(R_buffer(istrt),isize, &
                  &                     amr_mpi_real, &
                  &                     isrc ,itag,amr_mpi_meshComm, &
                  &                     recvrequest(ij),ierr)
          endif
       enddo
       commCtl % numReq = ij
    end if                        !doPostRecv

    if (doSend) then
       ji = 0
       do j = mype+1,nprocs
          isrc = mype
          idest= j-1
#ifdef PM_UNIQUE_MPI_TAGS
          itag = isrc*nprocs + idest+1 + tag_offset
#else
          itag = tag_offset !OK because 1 msg per process.
#endif
          ! send from mype=i
          if(commatrix_send(j).gt.0) then
            ji = ji+1
             istrt = is_buf(1,j)
             ilast = is_buf(2,j)
             isize = ilast-istrt+1

             if (amr_error_checking) then
                if(isize.gt.buf_dim_send) then
                   write(*,*) 'PARAMESH ERROR : ', &
                        &                   ' message is bigger than buffer space  '
                   call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
                endif
             endif
             call MPI_Isend(S_buffer(istrt),isize, &
                  &                      amr_mpi_real, &
                  &                      idest,itag,amr_mpi_meshComm, &
                            sendRequest(ji),ierr)
          endif
       enddo

       do j = 1,mype
          isrc = mype
          idest= j-1
#ifdef PM_UNIQUE_MPI_TAGS
          itag = isrc*nprocs + idest+1 + tag_offset
#else
          itag = tag_offset !OK because 1 msg per process.
#endif
          ! send from mype=i
          if(commatrix_send(j).gt.0) then
            ji = ji+1
             istrt = is_buf(1,j)
             ilast = is_buf(2,j)
             isize = ilast-istrt+1

             if (amr_error_checking) then
                if(isize.gt.buf_dim_send) then
                   write(*,*) 'PARAMESH ERROR : gr_mpiXchangeBlocks 2:', &
                        &                   ' message is bigger than buffer space  '
                   call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
                endif
             endif
             call MPI_Isend(S_buffer(istrt),isize, &
                  &                      amr_mpi_real, &
                  &                      idest,itag,amr_mpi_meshComm, &
                            sendRequest(ji),ierr)
          endif
       enddo
       commCtl % numSReq = ji
#ifdef EARLY_ISEND_WAIT
       do while (allSentCount < ji)
          if (allSentCount < ji) then
             call MPI_Waitsome(ji,sendrequest,sentCount, &
                            receivedInd, sendstatus, ierrorcode)
             allSentCount = allSentCount + sentCount
#ifdef DEBUG_GRIDCOMM
             print 96,mype,sentCount,ji,(receivedInd(i),i=1,sentCount)
#endif
          end if
       end do
#endif
    end if                        !doSend

    if (doWait) then
!!$       if(ij.gt.0) &
!!$            &   call MPI_Waitall(ij,recvrequest,recvstatus, &
!!$            &                    ierrorcode)
       do while (allReceivedCount < ij .OR. allSentCount < ji)
          if (allSentCount < ji) then
             call MPI_Waitsome(ji,sendrequest,sentCount, &
                            receivedInd, sendstatus, ierrorcode)
             allSentCount = allSentCount + sentCount
#ifdef DEBUG_GRIDCOMM
96           format('On',I2,': WaitSome sent',I2,'/',I2,',',20(1x,I4,:,','))
             print 96,mype,sentCount,ji,(receivedInd(i),i=1,sentCount)
#endif
          end if
          if (allReceivedCount < ij) then
             call MPI_Waitsome(ij,recvrequest,receivedCount, &
                            receivedInd, recvstatus, ierrorcode)
             allReceivedCount = allReceivedCount + receivedCount
#ifdef DEBUG_GRIDCOMM
98           format('On',I2,': WaitSome received',I2,'/',I2,',',20(1x,I4,:,','))
             print 98,mype,receivedCount,ij,(receivedInd(i),i=1,receivedCount)
99           format('On',I2,': WaitSome rec.procs',I2,'/',I2,',',20(1x,I4,:,','))
             print 99,mype,receivedCount,ij,(recvstatus(MPI_SOURCE,receivedInd(i)),i=1,receivedCount)
#endif
          end if
       end do

  ! reset offset with largest tag
#ifdef PM_UNIQUE_MPI_TAGS
       tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset
#endif

    end if                        !doWait
  end associate


end subroutine gr_mpiXchangeBlocks
