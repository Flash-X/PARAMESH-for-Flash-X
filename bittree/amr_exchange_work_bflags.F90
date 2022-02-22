!!****if* source/Grid/GridMain/AMR/Paramesh4/bittree/amr_exchange_work_bflags
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
!!   amr_exchange_work_bflags
!!
!! HISTORY
!!
!!   Tom Klosterman
!!
!!***

#include "paramesh_preprocessor.fh"
  
      subroutine amr_exchange_work_bflags(nprocs,mype,lnblocks_old, &
                                      new_loc,old_loc,new_child)

!-----Use statements
      use tree
      use paramesh_dimensions, only : maxblocks
      use paramesh_comm_data, ONLY : amr_mpi_meshComm, amr_mpi_real

      Implicit None

!-----Include statements
      Include 'mpif.h'

!-----Input/Output arguments.
      integer, intent(in) :: nprocs,mype
      integer, intent(in) :: lnblocks_old
      integer, intent(in) :: new_loc(2,maxblocks_tr)
      integer, intent(in) :: old_loc(2,maxblocks_tr)
      integer, intent(in) :: new_child(2,mchild,maxblocks_tr)

!-----Local variables. and arrays.
      Integer :: i,j
      real :: work_blockt(maxblocks_tr), work_scaled(maxblocks_tr)
      integer :: bflagst(mflags,maxblocks)
      integer :: nrecv, nsend, ierr
      integer :: reqr(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE, maxblocks_tr)

!-----Initialize temp arrays
      work_blockt(:) = -1.
      bflagst(:,:) = -1

!-----Scale work for parents to send to new children
      work_scaled = work_block * gr_btWorkChildScaling

!-----Post receives into temp arrays
      nrecv = 0
      do i=1,new_lnblocks

!-------Receive work_block from old_loc.
!-------New children inherit value from parent.
        if(gr_btExchangeWork.AND.(old_loc(2,i).ne.mype)) then
          nrecv = nrecv + 1
          call MPI_IRECV(work_blockt(i),1,amr_mpi_real, &
                         old_loc(2,i),i,amr_mpi_meshComm, &
                         reqr(nrecv),ierr)
        end if

!-------Receive bflags from old_loc.
!-------New children inherit value from parent.
        if(gr_btExchangeBflags.AND.(old_loc(2,i).ne.mype)) then
          nrecv = nrecv + 1
          call MPI_IRECV(bflagst(1,i),mflags,MPI_INTEGER, &
                         old_loc(2,i),i+maxblocks_tr,amr_mpi_meshComm, &
                         reqr(nrecv),ierr)
        end if
      end do !1,new_lnblocks

!-----Post sends from work_block and bflags.
      nsend=0
      do i=1,lnblocks_old
        if(new_loc(1,i).gt.0) then

          if(new_loc(2,i).ne.mype) then
            if(gr_btExchangeWork) then
              nsend=nsend+1
              call MPI_SSEND(work_block(i),1,amr_mpi_real, &
                             new_loc(2,i),new_loc(1,i), &
                             amr_mpi_meshComm, ierr)
            end if
            if(gr_btExchangeBflags) then
              nsend=nsend+1
              call MPI_SSEND(bflags(1,i),mflags,MPI_INTEGER, &
                             new_loc(2,i),new_loc(1,i)+maxblocks_tr, &
                             amr_mpi_meshComm, ierr)
            end if
          else
            if(gr_btExchangeWork) work_blockt(new_loc(1,i)) = &
                                                       work_block(i)
            if(gr_btExchangeBflags) bflagst(:,new_loc(1,i)) = &
                                                       bflags(:,i)
          end if !new_loc(2,i).ne.mype

!---------If new_child has nonzero entries, the block was refined.
!---------Post sends to all of its new children.
          if(any(new_child(:,:,i).gt.0)) then
            do j=1,nchild
              if(new_child(2,j,i).ne.mype) then
                if(gr_btExchangeWork) then
                  nsend = nsend+1
                  call MPI_SSEND(work_scaled(i),1,amr_mpi_real, &
                                 new_child(2,j,i),new_child(1,j,i), &
                                 amr_mpi_meshComm, ierr)
                end if
                if(gr_btExchangeBflags) then
                  nsend = nsend+1
                  call MPI_SSEND(bflags(1,i),mflags,MPI_INTEGER, &
                                 new_child(2,j,i), &
                                 new_child(1,j,i)+maxblocks_tr, &
                                 amr_mpi_meshComm, ierr)
                end if
              else
                if(gr_btExchangeWork) work_blockt(new_child(1,j,i)) = &
                              work_scaled(i)
                if(gr_btExchangeBflags) bflagst(:,new_child(1,j,i)) = &
                              bflags(:,i)
              end if !new_child(2,j,i).ne.mype
            end do !j=1,nchild
          end if !any(new_child(:,:,i).gt.0)

        end if !new_loc(1,i).gt.0
      end do !1,lnblocks_old

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

!-----Copy from temp arrays into proper arrays.
      do i=1,new_lnblocks
        if(gr_btExchangeWork) work_block(i) = work_blockt(i)
        if(gr_btExchangeBflags) bflags(:,i) = bflagst(:,i)
      end do !1,new_lnblocks

!-----Fill with default values if not exchanging.
      if(.NOT.gr_btExchangeWork) then
        where(nodetype.eq.1) work_block = gr_btWorkDefaultLeaf
        where(nodetype.gt.1) work_block = gr_btWorkDefaultPar
      end if
      if(.NOT.gr_btExchangeBflags) bflags(:,:) = -1

      end subroutine amr_exchange_work_bflags
