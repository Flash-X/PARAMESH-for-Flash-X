!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! paramesh_dimensions Module
!------------------------------------------------------------------------------


      Module paramesh_dimensions

#include "paramesh_preprocessor.fh"

!------------------
! Model dimensions
!------------------

! set physical dimension of model and number of edges on each grid block
#ifndef LIBRARY
      integer, parameter ::  ndim = N_DIM
      integer, parameter ::  nbedges = ndim*2**(ndim-1)
#else
      Integer :: ndim
      Integer :: nbedges
#endif

      Integer ::  l2p5d

! an increment variable for the z dimension to enable the same code to
! work for 2D or 3D models.
#ifndef LIBRARY
      integer, parameter :: k3d=(ndim-1)/2
      integer, parameter :: k2d=ndim/2
      integer, parameter :: k1d=1
#else
      Integer :: k3d, k2d, k1d
#endif

!----------------------
! Grid Block Properties
!----------------------

! set size of grid blocks
! The spatial relationship between grid points in parents
! and their children is subtley different depending on whether
! the block dimensions (ie nxb,nyb,nzb) are even or odd. This has 
! significant consequences when defining interpolation within
! restriction or prolongation operations.
#ifndef LIBRARY
      integer, parameter :: nxb = NXB
      integer, parameter :: nyb = max(1,(NYB) * k2d)
      integer, parameter :: nzb = max(1,(NZB) * k3d)
      integer, parameter :: maxdim=max(nxb,nyb,nzb)
#else
      Integer :: nxb
      Integer :: nyb
      Integer :: nzb
      Integer :: maxdim
#endif

! these guard cell offsets are required to accomodate differences
! in cases when block dimensions are odd or even
      Integer :: gc_off_x
      Integer :: gc_off_y
      Integer :: gc_off_z

! set the maximum number of blocks per processor
#ifndef LIBRARY
      integer, parameter :: maxblocks = MAXBLOCKS
      integer, parameter :: maxblocks_alloc = maxblocks*10
#else
      Integer :: maxblocks
      Integer :: maxblocks_alloc
#endif

! set the number of guard cell layers at each boundary
#ifndef LIBRARY
      integer, parameter :: nguard = NGUARD
#else
      Integer, Save :: nguard
#endif
#ifndef LIBRARY
#ifdef FL_NON_PERMANENT_GUARDCELLS
      integer, parameter :: npgs = 0
#else
      integer, parameter :: npgs = 1
#endif
#else
      Integer, Save :: npgs 
#endif

!----------------------
! Solution Variables
!----------------------

! Set number of unknowns associated with each grid cell.
! If you are using a multi-step timestep integration algorithm 
! (eg predictor-corrector) then the recommended value for nvar is 
!             nvarp*(nphase + 1)
! where nvarp denotes the number of physical variables (ie. 1D hydro would
! have nvarp=3, for mass, momentum and energy), nphase is the number of 
! stages in an integration timestep (ie predictor-corrector would have 
! nphase=2). Similar considerations apply for nfacevar.

#ifndef LIBRARY
      integer, parameter :: nvar = NUNK_VARS

! The number of data words needed on a cell face is set by nfacevar.
      integer, parameter :: nfacevar = NFACE_VARS

#else
      Integer :: nvar

! The number of data words needed on a cell face is set by nfacevar.
      Integer :: nfacevar
#endif

! The number of data words needed on cell edges is set by nvaredge.
      Integer :: nvaredge

! The number of data words needed at cell corners is set by nvarcorn.
      Integer :: nvarcorn

! The convention for relating variables associated with cell faces to the
! variables defined at cell centers is as follows:

! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell.
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell.

      Integer, Save      :: iface_off

!------------------------------------------------
! ! Declare dimensions for the solution variables
!------------------------------------------------

! Cell centered data bounds
#ifndef LIBRARY
      integer, parameter :: il_bnd=1, iu_bnd=nxb+2*nguard*npgs
      integer, parameter :: jl_bnd=1, ju_bnd=nyb+2*nguard*npgs*k2d
      integer, parameter :: kl_bnd=1, ku_bnd=nzb+2*nguard*npgs*k3d
#else
      Integer :: il_bnd, iu_bnd
      Integer :: jl_bnd, ju_bnd
      Integer :: kl_bnd, ku_bnd
#endif

! Cell centered data bounds defining the block interior (i.e. the part of
! the block NOT including guardcells
#ifndef LIBRARY
      integer, parameter :: il_bndi=nguard*npgs+1
      integer, parameter :: iu_bndi=nguard*npgs+nxb
      integer, parameter :: jl_bndi=nguard*npgs*k2d+1 
      integer, parameter :: ju_bndi=nguard*npgs*k2d+nyb
      integer, parameter :: kl_bndi=nguard*npgs*k3d+1
      integer, parameter :: ku_bndi=nguard*npgs*k3d+nzb
#else
      Integer :: il_bndi, iu_bndi
      Integer :: jl_bndi, ju_bndi
      Integer :: kl_bndi, ku_bndi
#endif

#ifndef LIBRARY
      integer, parameter :: nbndvar=max(1,nfacevar)
#else
      Integer :: nbndvar
#endif
      Integer :: nbndvare
      Integer :: nbndvarc
      Integer :: maxblocksf
      Integer :: maxblocksue
      Integer :: maxblocksn

! set data length of grid blocks
! cell center
      Integer :: len_block
! cell face centers
      Integer :: len_blockfx 
      Integer :: len_blockfy
      Integer :: len_blockfz
! cell edge centers
      Integer :: len_blockex
      Integer :: len_blockey
      Integer :: len_blockez
! cell corner
      Integer :: len_blockn
! cell face centers for exchange of block faces only
      Integer :: len_blockfxf
      Integer :: len_blockfyf
      Integer :: len_blockfzf

! Set the number of padded blocks required for the 1blk guardcell routines.
! This should be 2 in almost all circumstances, one block for the current
! working leaf node, and one for its parent.
      Integer, parameter :: npblks=2

! Set index bounds with guardcells included.
#ifndef LIBRARY
      integer, parameter :: il_bnd1=1,iu_bnd1=nxb+2*nguard
      integer, parameter :: jl_bnd1=1,ju_bnd1=nyb+2*nguard*k2d
      integer, parameter :: kl_bnd1=1,ku_bnd1=nzb+2*nguard*k3d
#else
      Integer :: il_bnd1,iu_bnd1
      Integer :: jl_bnd1,ju_bnd1
      Integer :: kl_bnd1,ku_bnd1
#endif

! Set length of messages required when working blocks are to be passed.
      Integer :: len_block1
      Integer :: len_blockfx1 
      Integer :: len_blockfy1
      Integer :: len_blockfz1 
      Integer :: len_blockex1
      Integer :: len_blockey1
      Integer :: len_blockez1 
      Integer :: len_blockn1
      Integer :: maxblocks_gt
      Integer :: maxblocksf_gt
      Integer :: maxblocksue_gt
      Integer :: maxblocksn_gt

!---------------------------------------
! Declare dimensions for the WORK arrays
!---------------------------------------

! Set number of guard cells associated with the workspace array.
#ifndef LIBRARY
      integer, parameter :: nguard_work = NGUARD
      integer, parameter :: ngw2=2*nguard_work
#else
      Integer :: nguard_work 
      Integer :: ngw2
#endif

! Set number of variables which the workspace array can handle.
      Integer :: nvar_work
      Integer :: maxblocksw

#ifndef LIBRARY
      integer, parameter :: ilw=1,iuw=nxb+ngw2*npgs
      integer, parameter :: jlw=1,juw=nyb+ngw2*npgs*k2d
      integer, parameter :: klw=1,kuw=nzb+ngw2*npgs*k3d
#else
      Integer :: ilw,iuw
      Integer :: jlw,juw
      Integer :: klw,kuw
#endif
      Integer :: len_wblock
! Set index bounds with guardcells included.
#ifndef LIBRARY
      integer, parameter :: ilw1=1,iuw1=nxb+ngw2
      integer, parameter :: jlw1=1,juw1=nyb+ngw2*k2d
      integer, parameter :: klw1=1,kuw1=nzb+ngw2*k3d
#else
      Integer :: ilw1,iuw1
      Integer :: jlw1,juw1
      Integer :: klw1,kuw1
#endif
      Integer :: len_wblock1

!-----------------
! Timestep control
!-----------------

! common block for timestep control
      Integer, parameter :: maxlevels=50

!---------------------------------
! Conservation at refinement jumps
!---------------------------------
      Real, Save :: red_f

!---------------------------------
! mpi message type control
!---------------------------------
      Integer, Save :: nmax_lays

!---------------------------------
! No. of fields for which divergence free prolongation
! will be required
      Integer, Save      :: nfield_divf 

!---------------------------------


      End Module paramesh_dimensions
!-----------------------------------------------------------------
