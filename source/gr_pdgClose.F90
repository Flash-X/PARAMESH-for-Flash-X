!! AUTHORS
!!
!!   Klaus Weide
!!
!! MODIFICATIONS
!!
!!  2022-10-31 Klaus Weide  moved gcell_on_cc flag array into pdg_t
!!  2022-11-08 K. Weide  moved cell_ geometry arrays from physicaldata to pdg_t
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

    Deallocate(pdg % cell_vol)
    Deallocate(pdg % cell_area1)
    Deallocate(pdg % cell_area2)
    Deallocate(pdg % cell_area3)
    Deallocate(pdg % cell_leng1)
    Deallocate(pdg % cell_leng2)
    Deallocate(pdg % cell_leng3)
    Deallocate(pdg % cell_face_coord1)
    Deallocate(pdg % cell_face_coord2)
    Deallocate(pdg % cell_face_coord3)

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

    deallocate(pdg % gcell_on_cc)

    pdg % nfluxvar    = 0
    pdg % nfluxes     = 0
    pdg % maxblocksfl = 0

end subroutine gr_pdgCloseOne
