!! MODIFICATIONS
!!  2022-10-10 Klaus Weide  Added new logical doRedist flag
!!  2022-10-13 Klaus Weide  Added 7 slots densVar...velzVar to pdgConst_t
!!  2022-10-14 Klaus Weide  Corrected type of densVar...velzVar to integer
module gr_pmPdgDecl
  implicit none

  type pdgConst_t                    ! pdg == physicaldata group
     ! from paramesh_dimensions.F90
     Integer :: nxb
     Integer :: nyb
     Integer :: nzb

     integer :: nguard

     integer :: nvar

     Integer :: il_bnd, iu_bnd
     Integer :: jl_bnd, ju_bnd
     Integer :: kl_bnd, ku_bnd
     Integer :: il_bndi, iu_bndi
     Integer :: jl_bndi, ju_bndi
     Integer :: kl_bndi, ku_bndi
     Integer :: il_bnd1,iu_bnd1
     Integer :: jl_bnd1,ju_bnd1
     Integer :: kl_bnd1,ku_bnd1
     
     ! Various variables that are treated specially in various places
     ! in Grid implementation, based on the VARIABLE names given them:
     integer :: densVar = -1
     integer :: presVar = -1
     integer :: enerVar = -1
     integer :: eintVar = -1
     integer :: velxVar = -1
     integer :: velyVar = -1
     integer :: velzVar = -1
  end type pdgConst_t

  type pdg_t                    ! pdg == physicaldata group
     ! from physicaldata.F90
     Real,Allocatable ::  unk(:,:,:,:,:)

     Real, Allocatable :: unk1(:,:,:,:,:)

     Real, Allocatable ::  flux_x(:,:,:,:,:)
     Real, Allocatable ::  flux_y(:,:,:,:,:)
     Real, Allocatable ::  flux_z(:,:,:,:,:)

     Real, Allocatable :: tflux_x(:,:,:,:,:)
     Real, Allocatable :: tflux_y(:,:,:,:,:)
     Real, Allocatable :: tflux_z(:,:,:,:,:)

     Real, Allocatable :: recvarxf(:,:,:,:)
     Real, Allocatable :: recvaryf(:,:,:,:)
     Real, Allocatable :: recvarzf(:,:,:,:)
     Real, Allocatable :: bndtempx1(:,:,:,:)
     Real, Allocatable :: bndtempy1(:,:,:,:)
     Real, Allocatable :: bndtempz1(:,:,:,:)

     ! from prolong_arrays.F90
     Real, Allocatable :: prol_dx(:)
     Real, Allocatable :: prol_dy(:) 
     Real, Allocatable :: prol_dz(:)
     Integer, Allocatable :: prol_indexx(:,:,:)
     Integer, Allocatable :: prol_indexy(:,:,:)
     Integer, Allocatable :: prol_indexz(:,:,:)

     ! integer variables from physicaldata.F90
     integer :: nfluxvar
     integer :: nfluxes
     integer :: maxblocksfl

     logical :: doRedist = .TRUE.
  end type pdg_t
end module gr_pmPdgDecl
