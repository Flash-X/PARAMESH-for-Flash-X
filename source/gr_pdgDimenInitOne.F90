
#include "Simulation.h"

subroutine gr_pdgDimenInitOne(pdgDimen, nvar,nguard,nx,ny,nz, npgsArg,k2dArg,k3dArg)
  use gr_pmPdgDecl, ONLY : pdgConst_t

  implicit none

  type(pdgConst_t), intent(OUT) :: pdgDimen
  integer, intent(IN) :: nvar
  integer, intent(IN) :: nguard
  integer, intent(IN) :: nx,ny,nz
  integer, intent(IN),OPTIONAL :: npgsArg,k2dArg,k3dArg

  integer :: npgs,k2d,k3d

  if (present(k2dArg)) then
     k2d = k2dArg
  else
     k2d = K2D
  end if
  if (present(k3dArg)) then
     k3d = k3dArg
  else
     k3d = K3D
  end if
  if (present(npgsArg)) then
     npgs = npgsArg
  else
     npgs = 1
  end if

  pdgDimen % nvar   = nvar
  pdgDimen % nguard = nguard
  pdgDimen % nxb    = nx
  pdgDimen % nyb    = ny
  pdgDimen % nzb    = nz

  pdgDimen % il_bnd       = 1
  pdgDimen % jl_bnd       = 1
  pdgDimen % kl_bnd       = 1
  pdgDimen % iu_bnd       = nx+2*nguard*npgs
  pdgDimen % ju_bnd       = ny+2*nguard*npgs*k2d
  pdgDimen % ku_bnd       = nz+2*nguard*npgs*k3d
  pdgDimen % il_bndi      = nguard*npgs+1
  pdgDimen % iu_bndi      = nguard*npgs+nx
  pdgDimen % jl_bndi      = nguard*npgs*k2d+1
  pdgDimen % ju_bndi      = nguard*npgs*k2d+ny
  pdgDimen % kl_bndi      = nguard*npgs*k3d+1
  pdgDimen % ku_bndi      = nguard*npgs*k3d+nz
  pdgDimen % il_bnd1      = 1
  pdgDimen % jl_bnd1      = 1
  pdgDimen % kl_bnd1      = 1
  pdgDimen % iu_bnd1      = nx+2*nguard
  pdgDimen % ju_bnd1      = ny+2*nguard*k2d
  pdgDimen % ku_bnd1      = nz+2*nguard*k3d

end subroutine gr_pdgDimenInitOne
