!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------


#include "paramesh_preprocessor.fh"
!#define DEBUG
      subroutine mpi_unpack_blocks(commatrixRecv,mype,iopt, & 
     &                             lcc,lfc,lec,lnc, & 
     &                             buf_dim,R_buffer,ig, &
     &                             nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine unpacks all blocks which are to be received on mype.
! It further stores the local (receiving) block id, the neighboring remote 
! (sending) block id, and the local guard block id into the array laddress 
! which is to be used in the subroutine mpi_1blk_guardcell.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!! MODIFICATIONS
!!  Modified for commatrixRecv arg : Klaus Weide          May 2022
!!  2022-11-08 Klaus Weide  Use paramesh_dimensions with ONLY
!------------------------------------------------------------------------
!
! Arguments:
!      commatrixRecv  a component from the current communication pattern
!      mype           current processor id
!      iopt           option setting for work array
!      lcc            if true include unk data in buffer
!      lfc            if true include facevar data in buffer
!      lec            if true include unk_e_? data in buffer
!      lnc            if true include unk_n data in buffer
!      buf_dim        dimension of buffer
!      R_buffer       receive buffer 
!
!------------------------------------------------------------------------
      use paramesh_dimensions, ONLY: maxblocks_alloc
      use physicaldata
      use tree
      use paramesh_comm_data

      use paramesh_mpi_interfaces, only : mpi_put_buffer

#include "Flashx_mpi_implicitNone.fh"
#include "FortranLangFeatures.fh"

      integer, CONTIGUOUS_INTENT(in) :: commatrixRecv(:)
      integer, intent(in) :: mype,buf_dim,iopt
      logical, intent(in) :: lcc,lfc,lec,lnc
      real,    intent(IN) ::  R_buffer(buf_dim)
      integer, intent(in) :: ig
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lblk, lnumb, lb
      integer :: index
      integer :: ierrorcode,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      lnumb = sum(commatrixRecv(:))
      if(lnumb.gt.maxblocks_alloc) then
            call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      endif
      index = 1

 
      do lblk=1,lnumb
        lb = lblk + strt_buffer - 1

                                  ! unpack all arrays from buffer into lb
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacking starting ', & 
     &        ' at index ',index,' buf_dim ',buf_dim
#endif /* DEBUG */
        call mpi_put_buffer( & 
     &         lb,iopt,index,lcc,lfc,lec,lnc,buf_dim,R_buffer,ig, &
     &         nlayersx,nlayersy,nlayersz)
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacked into ',lb
#endif /* DEBUG */

      enddo

#ifdef DEBUG
      if (index .ne. buf_dim) then
      print *,' '
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' Rbuffer is too big by ',size(R_buffer,dim=1)-index,mype
      print *,' index = ',index,' buf_dim = ',buf_dim,mype
      print *,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
      print *,' '
      end if
#endif

      return
      end subroutine mpi_unpack_blocks

!------------------------------------------------------------------------

      subroutine mpi_Rbuffer_size(pattern,mype,nprocs,iopt, & 
     &                           lcc,lfc,lec,lnc, & 
     &                           buf_dim, & 
     &                           block_sections, fluxes, edges, pdg,ig, flux_dir, &
     &                           nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine computes the exact size of the receive buffer
!
!
! Written :     Kevin Olson January 2007
!! MODIFICATIONS
!!  Modified :    Klaus Weide          May 2022
!!  2021-12-15 Klaus Weide  Made PDG-aware
!!  2022-11-08 Klaus Weide  Removed unnecessary Use paramesh_dimensions
!------------------------------------------------------------------------
!
! Arguments:
!      pattern        current communication pattern
!      mype           current processor id
!      nprocs         No. of processors
!      iopt           option setting for work array
!      lcc            if true include unk data in buffer
!      lfc            if true include facevar data in buffer
!      lec            if true include unk_e_? data in buffer
!      lnc            if true include unk_n data in buffer
!      buf_dim        dimension of buffer
!
!------------------------------------------------------------------------
      use gr_pmPdgDecl, ONLY : pdg_t
      use gr_pmCommDataTypes, ONLY: gr_pmCommPattern_t
      use physicaldata
      use tree
      use paramesh_comm_data

      use paramesh_mpi_interfaces, only : mpi_get_Rbuffer_size, &
                                          mpi_get_Rbuffer_size_fluxes, &
                                          mpi_get_Rbuffer_size_edges

      implicit none

      TYPE(gr_pmCommPattern_t),intent(in) :: pattern
      integer, intent(in)  :: mype, nprocs, iopt
      integer, intent(out) :: buf_dim
      logical, intent(in)  :: lcc,lfc,lec,lnc
      logical, intent(in)  :: block_sections, fluxes, edges
      type(pdg_t), intent(IN) :: pdg
      integer, intent(in)  :: ig
      integer, intent(in), optional :: flux_dir
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb
      integer :: index
      integer :: ierrorcode,ierr
      integer :: next_pe, irpe, dtype, iblk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      index = 1

      next_pe = 0
      do irpe = 1, nprocs
         if (pattern%commatrix_recv(irpe) > 0) then
            next_pe = next_pe + 1
            do iblk = 1, pattern%commatrix_recv(irpe)
               if (pattern%to_be_received(1,iblk,next_pe) > 0) then
                  lb = pattern%to_be_received(1,iblk,next_pe)
                  dtype = pattern%to_be_received(3,iblk,next_pe)

                                   ! unpack all arrays from buffer into lb
#ifdef DEBUG
        write(*,*) 'pe ',mype,' iblk ',iblk,' unpacking starting ', & 
     &        ' at index ',index,' buf_dim ',buf_dim
#endif /* DEBUG */

                  if (block_sections) then

                  call mpi_get_Rbuffer_size(  & 
     &                   lb,dtype,iopt,index, &
     &                   lcc,lfc,lec,lnc,ig,  &
     &                   nlayersx,nlayersy,nlayersz)

                  elseif (fluxes) then

                  call mpi_get_Rbuffer_size_fluxes(  & 
     &                   lb,dtype,index, pdg,ig, flux_dir)

                  elseif (edges) then

                  call mpi_get_Rbuffer_size_edges(  & 
     &                   lb,dtype,index)

                  end if
#ifdef DEBUG
        write(*,*) 'pe ',mype,' iblk ',iblk,' unpacked into ',lb
#endif /* DEBUG */
                end if
             end do
          end if
      enddo

      if (index > 0) then
         buf_dim = index
      else
         buf_dim = 1
      end if

      return
      end subroutine mpi_Rbuffer_size
