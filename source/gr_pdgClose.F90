subroutine gr_pdgCloseOne(pdg)
  use gr_pmPdgDecl, ONLY : pdg_t
  use paramesh_dimensions, ONLY: maxblocks, npblks

  implicit none

  type(pdg_t), intent(INOUT) :: pdg

    deallocate(           &
        pdg % unk)

    deallocate(pdg % unk1)

!-----Deallocate arrays for flux fix-up at refinement jumps

    deallocate(                & 
       pdg % flux_x)
    deallocate(                & 
       pdg % flux_y)
    deallocate(                & 
       pdg % flux_z)

    deallocate(pdg % recvarxf)
    deallocate(pdg % recvaryf)
    deallocate(pdg % recvarzf)
    deallocate(pdg % bndtempx1)
    deallocate(pdg % bndtempy1)
    deallocate(pdg % bndtempz1)

    deallocate(pdg % prol_dx)
    deallocate(pdg % prol_dy)
    deallocate(pdg % prol_dz)
    deallocate(pdg % prol_indexx)
    deallocate(pdg % prol_indexy)
    deallocate(pdg % prol_indexz)

    pdg % nfluxvar    = 0
    pdg % nfluxes     = 0
    pdg % maxblocksfl = 0

end subroutine gr_pdgCloseOne
