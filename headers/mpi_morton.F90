!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------------
! mpi_morton Module



      Module mpi_morton

      Use paramesh_dimensions

      Private

      Public :: npts_neigh
      Integer, Parameter :: npts_neigh = 3000

! 
! variables for storing the morton environment
      Public :: ir_buf,is_buf


      Integer, Save,dimension(:,:),allocatable :: ir_buf
      Integer, Save,dimension(:,:),allocatable :: is_buf


! list of block edges which need diagonal info during edge averaging
      Public :: edge_mark,no_of_diagonal_edges
      Integer, Save :: edge_mark(6:6,4,npts_neigh)
      Integer, Save :: no_of_diagonal_edges

! Used to make searching of laddress more efficient
      Public :: ladd_strt,ladd_end
      Integer,Save,dimension(:),allocatable :: ladd_strt,ladd_end

!new code

      Public :: message_size_cc
      Public :: message_size_fcx, message_size_fcy, message_size_fcz
      Public :: message_size_ec
      Public :: message_size_nc
      Public :: message_size_wk
      Public :: mess_segment_loc
      Integer,Save,dimension(2*27) :: message_size_cc
      Integer,Save,dimension(2*27) :: message_size_fcx
      Integer,Save,dimension(2*27) :: message_size_fcy
      Integer,Save,dimension(2*27) :: message_size_fcz
      Integer,Save,dimension(2*27) :: message_size_ec
      Integer,Save,dimension(2*27) :: message_size_nc
      Integer,Save,dimension(2*27) :: message_size_wk
      Integer,Save,dimension(:),allocatable :: mess_segment_loc

      Public :: temprecv_buf
      Real, Save,dimension(:),allocatable :: temprecv_buf

      Public :: l_datapacked
      Logical,Save,dimension(5) :: l_datapacked

!new code end

      Public :: lperiodicx,lperiodicy,lperiodicz
      Logical, Save :: lperiodicx
      Logical, Save :: lperiodicy
      Logical, Save :: lperiodicz

      Public :: treeinfo

      Type treeinfo
        Real coord(3)
        Real bsize(3)
        Real bnd_box(2,3)
        Integer parent(2)
        Integer which_child
        Logical newchild
        Integer neigh(2,6)
        Integer lrefine
        Integer nodetype
        Integer empty
      End Type treeinfo


      End Module mpi_morton
!-----------------------------------------------------------------------
