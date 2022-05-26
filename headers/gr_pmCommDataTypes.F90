!!****ih* headers/gr_pmCommDataTypes
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
!!   gr_pmCommDataTypes
!!
!! SYNOPSIS
!!
!!   use gr_pmCommDataTypes
!!
!! DESCRIPTION
!!
!!   Fortran module which holds data types and some constants
!!   related to Paramesh data communications.
!!
!! HISTORY
!!
!!  2020 Dec.  K. Weide  created from information collected elsewhere
!!  2021       K. Weide  reorg & additions for supporting PM async comm
!!
!! MODIFICATIONS
!!  2022-05-13 K. Weide  tweaked components of gr_pmCommPattern_t for use
!!  2022-05-20 K. Weide  added id omponent, constants like GRID_SUBPAT_GC_OPT(2)
!!  2022-05-25 K. Weide  added constant GRID_SUBPAT_RESTRICT_FOR_FCORR(3)
!!***

module gr_pmCommDataTypes
  implicit none

  type gr_pmCommPattern_t
     ! COMMUNICATION PATTERN PROPER BEGIN !
     ! For each other pe, how many pp blocks will I send to it?  (pp = "potentially partial")
     integer,dimension(:),allocatable :: commatrix_recv !set: process_fetch_list, SWAPINOUT (commatrix_guard(:,1)(etc.))
     ! For each other pe, how many pp blocks will I receive from it?
     integer,dimension(:),allocatable :: commatrix_send !set: process_fetch_list, SWAPINOUT (commatrix_guard(:,2)(etc.))

     ! For each remote pe from which I receive data:
     !     For each pp block to receive:
     !         global identifier (plus extra) (blockID,pe,dtype,[...]) of the fragment to receive
     Integer,dimension(:,:,:),allocatable :: to_be_received !set: process_fetch_list, SWAPINOUT, use: mpi_pack_blocks(etc.)
     ! For each remote pe to which I send data:
     !     For each pp block to send:
     !         global identifier (plus extra) (blockID,MYPE,dtype{,[#prior fragments in recipient's buffer for all 27 dtypes]})
     !         of the fragment to send
     Integer,dimension(:,:,:),allocatable :: to_be_sent !set: process_fetch_list, SWAPINOUT, use: mpi_pack_blocks(etc.)

     ! For each pp block (to be) received, its global identifier (blockID,pe):
     Integer, Allocatable :: laddress(:,:) ! tree.F90, set: process_fetch_list, use: mpi_amr_1blk_guardcell etc., SWAPINOUT
     ! COMMUNICATION PATTERN PROPER END !

     ! TO HOW MANY procs will I send? (minimum 0?) ! Can be easily recomputed bu counting nonzeros in  commatrix_send
     Integer                          :: num_recipient_pes !set: process_fetch_list, SWAPINOUT
     ! Guard block starting index, MUST SATISFY lnblocks < strt_buffer <= maxblocks_alloc :
     Integer                          :: strt_buffer !set: process_fetch_list, SWAPINOUT (strt_guard(etc.))
     integer :: id = -999
     logical :: valid = .FALSE.
  end type gr_pmCommPattern_t


  type gr_pmCommShaped_t
     ! For each pp block to receive from a remote pe:
     !     (starting index in recipient's global temprecv_buf,
     !      ending   index in recipient's global temprecv_buf), both in terms of FLASH_REALs.
     Integer,dimension(:,:),allocatable :: ir_buf       ! set: mpi_pack_blocks <- mpi_amr_comm_setup
     ! For each pp block to send to a remote pe:
     !     (starting index in sender's global send buffer,
     !      ending   index in sender's global send buffer,
     !      starting index in recipient's global temprecv_buf), all in terms of FLASH_REALs.
     Integer,dimension(:,:),allocatable :: is_buf       ! set: mpi_pack_blocks <- mpi_amr_comm_setup

     Integer,dimension(2*27) :: message_size_cc         ! set: mpi_set_message_sizes <- mpi_amr_comm_setup
     Integer,dimension(2*27) :: message_size_fcx,message_size_fcy,message_size_fcz ! etc...; as previous
     Integer,dimension(:),allocatable :: mess_segment_loc ! set: mpi_pack_blocks(etc.), use: amr_mpi_find_blk_in_buffer
     Integer,dimension(:),allocatable :: ladd_strt,ladd_end !set: mpi_amr_comm_setup from strt_buffer
  end type gr_pmCommShaped_t

  type gr_pmCommExtra_t
     Integer,dimension(:),allocatable :: pe_remote ! unused
     Integer,dimension(:),allocatable :: pe_source ! unused, was set: mpi_amr_store_comm_info, used: mpi_unpack_tree_info
     Integer,dimension(:),allocatable :: pe_destination ! unused


     Integer  :: largest_no_of_blocks ! unused

  end type gr_pmCommExtra_t

  type gr_pmCommCtl_t
     integer,allocatable :: recvrequest(:), sendRequest(:)
     integer,allocatable :: recvstatus(:,:),sendStatus(:,:)
     integer,allocatable :: receivedIndices(:)
     integer :: numReq, numSReq
     integer :: allReceivedCount, allSentCount
     Real, Dimension (:), allocatable :: sendBuf
  end type gr_pmCommCtl_t

  type gr_pmCommBuf_t
  end type gr_pmCommBuf_t

  integer,parameter,public:: GRID_PAT_GC=10, GRID_PAT_PROLONG=20, GRID_PAT_FCORR=30, GRID_PAT_RESTRICT=40

  integer,parameter,public:: GRID_SUBPAT_GC_DEFAULT=1, GRID_SUBPAT_GC_OPT=2

  integer,parameter,public:: GRID_SUBPAT_RESTRICT_DEFAULT=1, GRID_SUBPAT_RESTRICT_ANC=2
  integer,parameter,public:: GRID_SUBPAT_RESTRICT_FOR_FCORR=3

end module gr_pmCommDataTypes
