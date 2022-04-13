
!!REORDER(5): unk, unk1, flux_[xyz], tflux_[xyz]
!!REORDER(4): recvar[xyz]f

subroutine gr_pdgInitOne(pdg,pdgDimen, nfluxvar,nfluxes,maxblocksfl)
  use gr_pmPdgDecl, ONLY : pdg_t, pdgConst_t
  use paramesh_dimensions, ONLY: maxblocks, npblks

  implicit none

  type(pdg_t), intent(OUT) :: pdg
  type(pdgConst_t), intent(IN) :: pdgDimen
  integer, intent(IN) :: nfluxvar,nfluxes,maxblocksfl

  ASSOCIATE(il_bnd      => pdgDimen % il_bnd,  &
            iu_bnd      => pdgDimen % iu_bnd,  &
            jl_bnd      => pdgDimen % jl_bnd,  &
            ju_bnd      => pdgDimen % ju_bnd,  &
            kl_bnd      => pdgDimen % kl_bnd,  &
            ku_bnd      => pdgDimen % ku_bnd,  &
            il_bnd1     => pdgDimen % il_bnd1,  &
            iu_bnd1     => pdgDimen % iu_bnd1,  &
            jl_bnd1     => pdgDimen % jl_bnd1,  &
            ju_bnd1     => pdgDimen % ju_bnd1,  &
            kl_bnd1     => pdgDimen % kl_bnd1,  &
            ku_bnd1     => pdgDimen % ku_bnd1,  &
            il_bndi     => pdgDimen % il_bndi,  &
            iu_bndi     => pdgDimen % iu_bndi,  &
            jl_bndi     => pdgDimen % jl_bndi,  &
            ju_bndi     => pdgDimen % ju_bndi,  &
            kl_bndi     => pdgDimen % kl_bndi,  &
            ku_bndi     => pdgDimen % ku_bndi,  &
            nvar        => pdgDimen % nvar      &
            )
    Allocate(           &
        pdg % unk(nvar,          & 
            il_bnd:iu_bnd, & 
            jl_bnd:ju_bnd, & 
            kl_bnd:ku_bnd, & 
            maxblocks))

    Allocate(pdg % unk1(nvar,            &
                     il_bnd1:iu_bnd1, &
                     jl_bnd1:ju_bnd1, & 
                     kl_bnd1:ku_bnd1, &
                     npblks))

!-----Allocate arrays for flux fix-up at refinement jumps

    Allocate(                & 
       pdg % flux_x(nfluxes,         &
              1:2,             & 
              jl_bndi:ju_bndi, &
              kl_bndi:ku_bndi, &
              maxblocksfl))
    Allocate(                & 
       pdg % flux_y(nfluxes,         &
              il_bndi:iu_bndi, & 
              1:2,             &
              kl_bndi:ku_bndi, &
              maxblocksfl))
    Allocate(                & 
       pdg % flux_z(nfluxes,         &
              il_bndi:iu_bndi, & 
              jl_bndi:ju_bndi, &
              1:2,             &
              maxblocksfl))

    Allocate(                & 
       pdg % tflux_x(nfluxes,         &
              1:2,             & 
              jl_bndi:ju_bndi, &
              kl_bndi:ku_bndi, &
              maxblocksfl))
    Allocate(                & 
       pdg % tflux_y(nfluxes,         &
              il_bndi:iu_bndi, & 
              1:2,             &
              kl_bndi:ku_bndi, &
              maxblocksfl))
    Allocate(                & 
       pdg % tflux_z(nfluxes,         &
              il_bndi:iu_bndi, & 
              jl_bndi:ju_bndi, &
              1:2,             &
              maxblocksfl))

    Allocate(pdg % recvarxf(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi))
    Allocate(pdg % recvaryf(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi))
    Allocate(pdg % recvarzf(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2))
    Allocate(pdg % bndtempx1(nfluxes,1:2,jl_bndi:ju_bndi,kl_bndi:ku_bndi))
    Allocate(pdg % bndtempy1(nfluxes,il_bndi:iu_bndi,1:2,kl_bndi:ku_bndi))
    Allocate(pdg % bndtempz1(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi,1:2))

    Allocate(pdg % prol_dx(il_bnd1:iu_bnd1))
    Allocate(pdg % prol_dy(jl_bnd1:ju_bnd1))
    Allocate(pdg % prol_dz(kl_bnd1:ku_bnd1))
    Allocate(pdg % prol_indexx(2,il_bnd1:iu_bnd1,2))
    Allocate(pdg % prol_indexy(2,jl_bnd1:ju_bnd1,2))
    Allocate(pdg % prol_indexz(2,kl_bnd1:ku_bnd1,2))

    pdg % nfluxvar    = nfluxvar
    pdg % nfluxes     = nfluxes
    pdg % maxblocksfl = maxblocksfl

  end ASSOCIATE
end subroutine gr_pdgInitOne
