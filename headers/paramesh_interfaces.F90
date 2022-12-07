!!****ih* headers/paramesh_interfaces
!! NOTICE
!!  This file is from PARAMESH - an adaptive mesh library.
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
!!   paramesh_interfaces
!!
!! SYNOPSIS
!!
!!   use paramesh_interfaces
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!
!!  2022-05-20 K. Weide  added optional presentRegions to amr_1blk_guardcell_srl
!!  2022-05-23 K. Weide  added parentPresentRegions to amr_1blk_guardcell
!!  2022-05-27 K. Weide  Added pdgNo,pdg,ig args; variant routines; gr_pmPdgDecl
!!  2022-10-07 K. Weide  Added pdg,ig args to more routines
!!  2022-10-10 K. Weide  Added interface for send_block_data;
!!                       adjusted intent for some new_loc arguments to IN.
!!  2022-10-26 K. Weide  added amr_prolong_gen_unk1_fun interface
!!                Changed intent for some recv arguments to IN
!!  2022-11-02 K. Weide  added ig to amr_restrict_unk_fun interface
!!  2022-11-08 K. Weide  moved cell_ geometry arrays from physicaldata to pdg_t
!!  2022-11-08 K. Weide  added pdg,ig arguments to amr_1blk_cc_prol_gen_unk_fun
!!  2022-11-08 K. Weide  added pdg,ig arguments to amr_1blk_to_perm
!!  2022-11-08 K. Weide  added ig argument to amr_restrict_unk_genorder
!!  2022-12-03 K. Weide  amr_1blk_cc_prol_inject interface with pdg,ig args
!!***

!#ifdef HAVE_CONFIG_H
!#include <config.h>
!#endif

#include "FortranLangFeatures.fh"
#include "Simulation.h"
#include "paramesh_preprocessor.fh"

      module paramesh_interfaces

      interface
      subroutine amr_1blk_bcset(mype,ibc,lb,pe,                          & 
     &    idest,iopt,ibnd,jbnd,kbnd,surrblks,ig)
      integer, intent(in) :: mype,ibc,lb,pe
      integer, intent(in) :: idest,iopt,ibnd,jbnd,kbnd
      integer, intent(in) :: surrblks(:,:,:,:)
      integer, intent(in) :: ig
      end subroutine amr_1blk_bcset
      end interface

      interface
      subroutine amr_1blk_cc_cp_remote(mype,remote_pe,remote_block,      & 
     &    idest,iopt,id,jd,kd,is,js,ks,ilays,jlays,klays,nblk_ind, & 
     &    ipolar,pdg,ig)
        use gr_pmPdgDecl, ONLY : pdg_t
        implicit none
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,iopt,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays,nblk_ind
      integer, intent(in) :: ipolar(:)
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_1blk_cc_cp_remote
      end interface

      interface
      subroutine amr_1blk_nc_cp_remote(mype,remote_pe,remote_block,      & 
     &    idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,         & 
     &    ip3,jp3,kp3,nblk_ind,ig)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip3,jp3,kp3
      integer, intent(in) :: nblk_ind
      integer, intent(in) :: ig
      end subroutine amr_1blk_nc_cp_remote
      end interface

      interface
         subroutine amr_1blk_cc_prol_gen_unk_fun(recv,ia,ib,ja,jb,ka,kb, &
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p, pdg,ig)
           use gr_pmPdgDecl, ONLY : pdg_t
           implicit none
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)
      type(pdg_t),intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_1blk_cc_prol_gen_unk_fun
      end interface

      interface
      subroutine amr_1blk_cc_prol_inject(recv,ia,ib,ja,jb,ka,kb,         & 
     &       idest,ioff,joff,koff,mype,ivar,pdg,ig)
           use gr_pmPdgDecl, ONLY : pdg_t
           implicit none
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(in) :: recv(:,:,:,:)
      type(pdg_t),intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_1blk_cc_prol_inject
      end interface

      interface
      subroutine amr_1blk_cc_prol_linear(recv,ia,ib,ja,jb,ka,kb,          & 
     &       idest,ioff,joff,koff,mype,ivar,pdg)
        use gr_pmPdgDecl, ONLY : pdg_t
        implicit none
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(IN) :: recv(:,:,:,:)
      type(pdg_t),intent(INOUT) :: pdg
      end subroutine amr_1blk_cc_prol_linear
      end interface

      interface
      subroutine amr_1blk_cc_prol_genorder(recv,ia,ib,ja,jb,ka,kb,        & 
     &       idest,ioff,joff,koff,mype,ivar,order,pdg,ig)
        use gr_pmPdgDecl, ONLY : pdg_t
        implicit none
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,order
      real,    intent(IN) :: recv(:,:,:,:)
      type(pdg_t),intent(INOUT) :: pdg
      Integer, Intent(in) :: ig
      end subroutine amr_1blk_cc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_cc_prol_user()
      end subroutine amr_1blk_cc_prol_user
      end interface

      interface
      subroutine amr_1blk_cc_prol_dg(recv,ia,ib,ja,jb,ka,kb,          &
     &       idest,ioff,joff,koff,mype,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_cc_prol_dg
      end interface

      interface amr_prolong_gen_unk1_fun
         subroutine amr_prolong_gen_unk1_fun &
     &       (recv, ia, ib, ja, jb, ka, kb, idest, &
     &        ioff, joff, koff, mype, isg, pdg,ig)
           use gr_pmPdgDecl, ONLY : pdg_t
           implicit none
           real, CONTIGUOUS_INTENT(IN), dimension(:,:,:,:) :: recv
           integer, intent(IN) :: idest, isg, mype, ia, ib, ja, jb, ka, kb,ioff,joff,koff
           type(pdg_t),intent(INOUT) :: pdg
           Integer, Intent(in) :: ig
         end subroutine amr_prolong_gen_unk1_fun
      end interface

      interface
      subroutine amr_1blk_cc_prol_gen_work_fun(recv,                     & 
     &       ia,ib,ja,jb,ka,kb,                                          & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p,interp
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_gen_work_fun
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_inject(recv,                     & 
     &       ia,ib,ja,jb,ka,kb,                                        & 
     &       idest,ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_inject
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_linear(recv,                   & 
     &       ia,ib,ja,jb,ka,kb,                                      & 
     &       idest,ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_linear
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_genorder(recv,             & 
     &       ia,ib,ja,jb,ka,kb,                                  & 
     &       idest,ioff,joff,koff,mype,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: order
      real,    intent(inout) :: recv(:,:,:)
      end subroutine amr_1blk_cc_prol_work_genorder
      end interface

      interface
         subroutine amr_1blk_cc_prol_work_mg(recv, &
                ia,ib,ja,jb,ka,kb,                 &
                idest,ioff,joff,koff,lb,order)
           integer, intent(in) :: ia,ib,ja,jb,ka,kb
           integer, intent(in) :: idest,ioff,joff,koff,lb
           integer, intent(in) :: order
           real,    intent(inout) :: recv(:,:,:)
         end subroutine amr_1blk_cc_prol_work_mg
      end interface

      interface
      subroutine amr_1blk_cc_prol_work_user()
      end subroutine amr_1blk_cc_prol_work_user
      end interface

      interface
      subroutine amr_1blk_copy_soln(level)
      integer, intent(in) :: level
      end subroutine amr_1blk_copy_soln
      end interface

      interface
      subroutine amr_1blk_ec_cp_remote(mype,remote_pe,remote_block, & 
     &   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,    & 
     &    ip2,jp2,kp2,ip3,jp3,kp3,iface,nblk_ind,ig)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3,iface
      integer, intent(in) :: nblk_ind
      integer, intent(in) :: ig
      end subroutine amr_1blk_ec_cp_remote 
      end interface

      interface
      subroutine amr_1blk_fc_cp_remote(mype,remote_pe,remote_block, & 
     &   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,    & 
     &    ip2,jp2,kp2,iface,nblk_ind,ipolar,                       &
          ig,                                                      &
          curBlock,ibnd,jbnd,kbnd,surrblks)
      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,iface
      integer, intent(in) :: nblk_ind
      integer, intent(in) :: ipolar(:)
      integer, intent(in) :: ig
      integer,OPTIONAL, intent(in) :: curBlock
      integer,OPTIONAL, intent(in) :: ibnd,jbnd,kbnd
      integer,OPTIONAL, intent(in) :: surrblks(:,:,:,:)
      end subroutine amr_1blk_fc_cp_remote 
      end interface

      interface
      subroutine amr_1blk_ec_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype,iface)
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype,iface
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_gen_fun
      end interface

      interface
      subroutine amr_1blk_ec_prol_linear & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &      mype,ivar,iedge_dir)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar, iedge_dir
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_linear
      end interface

      interface
      subroutine amr_1blk_ec_prol_genorder & 
     &     (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, & 
     &      mype,ivar,iedge_dir,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,iedge_dir,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_ec_prol_genorder
      end interface

      interface
      subroutine amr_1blk_ec_prol_user()
      end subroutine amr_1blk_ec_prol_user
      end interface

      interface
      subroutine amr_1blk_nc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest, & 
     &       ioff,joff,koff,mype)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_gen_fun
      end interface

      interface
      subroutine amr_1blk_nc_prol_linear(recv,ia,ib,ja,jb,ka,kb,idest,   & 
     &       ioff,joff,koff,mype,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_linear
      end interface

      interface
      subroutine amr_1blk_nc_prol_genorder(                              & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,ivar,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_nc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_nc_prol_user()
      end subroutine amr_1blk_nc_prol_user
      end interface

      interface
      subroutine amr_1blk_fc_prol_gen_fun(recv,ia,ib,ja,jb,ka,kb,idest,  & 
     &       ioff,joff,koff,mype,lb,pe_p,lb_p,iface)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: lb,lb_p,pe_p
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_gen_fun
      end interface

      interface
      subroutine prol_fc_dbz_init(n,i_fc_vars)
      integer, intent(in) :: n, i_fc_vars(:,:)
      end subroutine prol_fc_dbz_init
      end interface

      interface
      function prol_fc_dbz_varp(ivar, iface) result(ldbz)
      integer, intent(in) :: ivar, iface
      logical :: ldbz
      end function prol_fc_dbz_varp
      end interface
      
      interface
      subroutine amr_1blk_fc_prol_dbz(                                   & 
     &        recvfx, recvfy, recvfz,                                    & 
     &        nfacevar, iv1, iv2, iv3                                    & 
     &        ,ia,ib,ja,jb,ka,kb,                                        & 
     &        idest,ioff,joff,koff,                                      & 
     &        mype,lb,parent_pe,parent_blk)
      real, intent(inout), dimension(:,:,:,:) :: recvfx,recvfy,recvfz
      integer, intent(in)    :: nfacevar, iv1, iv2, iv3
      integer, intent(in)    :: ia,ib,ja,jb,ka,kb
      integer, intent(in)    :: idest,ioff,joff,koff
      integer, intent(in)    :: mype,lb,parent_pe,parent_blk
      end subroutine amr_1blk_fc_prol_dbz
      end interface

      interface
      subroutine amr_1blk_fc_prol_inject(                                & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_inject
      end interface


      interface
      subroutine amr_1blk_fc_prol_linear(                                & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype,iface
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_linear
      end interface

      interface
      subroutine amr_1blk_fc_prol_genorder(                              & 
     &       recv,ia,ib,ja,jb,ka,kb,idest,                               & 
     &       ioff,joff,koff,mype,iface,ivar,order)
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,iface,order
      real,    intent(inout) :: recv(:,:,:,:)
      end subroutine amr_1blk_fc_prol_genorder
      end interface

      interface
      subroutine amr_1blk_fc_prol_user()
      end subroutine amr_1blk_fc_prol_user
      end interface

      interface
      subroutine amr_1blk_guardcell(mype,iopt,nlayers,lb,pe,             & 
     &                              lcc,lfc,lec,lnc,                     & 
     &                              l_srl_only,icoord,ldiag,             & 
                                    pdg,ig,                              &
     &                              nlayersx, nlayersy, nlayersz,        &
     &                              parentPresentRegions)
      use gr_pmPdgDecl, ONLY : pdg_t
      implicit none
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      integer(kind=selected_int_kind(9)), intent(in), optional :: parentPresentRegions
      end subroutine amr_1blk_guardcell
      end interface


      interface
      subroutine amr_1blk_guardcell_c_to_f(mype,lb,pe,iopt,nlayers,      & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayersx,nlayersy,nlayersz)
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx, nlayersy, nlayersz
      end subroutine amr_1blk_guardcell_c_to_f
      end interface


      interface
      subroutine amr_1blk_guardcell_reset
      end subroutine amr_1blk_guardcell_reset
      end interface


      interface
      subroutine amr_1blk_guardcell_srl(mype,pe,lb,iblock,iopt,nlayers,  & 
     &                         surrblks,lcc,lfc,lec,lnc,icoord,ldiag,    & 
     &                         nlayers0x,nlayers0y,nlayers0z, & 
     &                         ipolar,pdg,ig,curBlock,presentRegions)
      use gr_pmPdgDecl, ONLY : pdg_t
      use gr_flashHook_interfaces, ONLY : i27b
      implicit none
      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayers0x, nlayers0y, nlayers0z
      integer, intent(in) :: ipolar(2)
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      integer,OPTIONAL, intent(in) :: curBlock
      integer(kind=i27b),OPTIONAL, intent(in) :: presentRegions
      end subroutine amr_1blk_guardcell_srl
      end interface


      interface
      subroutine amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc)
      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc
      end subroutine amr_1blk_restrict
      end interface


      interface
      subroutine amr_1blk_save_soln
      end subroutine amr_1blk_save_soln
      end interface

      interface
         subroutine amr_1blk_to_perm(lcc,lfc,lec,lnc,lb,iopt,idest,pdg,ig)
           use gr_pmPdgDecl, ONLY : pdg_t
           implicit none
      integer, intent(in) :: lb,iopt,idest
      logical, intent(in) :: lcc,lfc,lec,lnc
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_1blk_to_perm
      end interface


      interface
      subroutine amr_abort()
      end subroutine amr_abort
      end interface

      interface
      subroutine amr_bi_sort(list,gid,npp)
      integer,intent(inout) :: list(:)
      integer,intent(inout) :: gid(:)
      integer,intent(in) :: npp
      end subroutine amr_bi_sort
      end interface

      interface
      subroutine amr_bc_block(jface,ibc,iopt,l,mype)
      integer, intent(in) :: jface,ibc,iopt,l,mype
      end subroutine amr_bc_block
      end interface


      interface
         subroutine amr_block_geometry(lb,pe,pdg,ig)
           use gr_pmPdgDecl, ONLY: pdg_t
           implicit none
           integer, intent(in) :: lb,pe
           type(pdg_t), intent(INOUT) :: pdg
           integer, intent(in) :: ig
         end subroutine amr_block_geometry
      end interface

      interface
      subroutine user_coord_transfm(x,y)
      integer, intent(inout) :: x,y
      end subroutine user_coord_transfm
      end interface

      interface
      subroutine amr_checkpoint_re(iunit1, & 
     &                             l_with_guardcells, & 
     &                             check_format, & 
     &                             user_attr_1, & 
     &                             user_attr_2, & 
     &                             user_attr_3, & 
     &                             user_attr_4, & 
     &                             user_attr_5)

      integer, intent(in)                      :: iunit1
      logical, optional, intent(in)            :: l_with_guardcells
      character (len=80), optional, intent(in) :: check_format
      real, optional, intent(out)              :: user_attr_1, & 
     &                                            user_attr_2, & 
     &                                            user_attr_3, & 
     &                                            user_attr_4, & 
     &                                            user_attr_5
      end subroutine amr_checkpoint_re
      end interface

      interface
      subroutine amr_checkpoint_re_default(iunit1,l_with_guardcells, & 
     &  user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)
      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(out) :: & 
     &  user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5
      end subroutine amr_checkpoint_re_default
      end interface

      interface
      subroutine amr_checkpoint_re_hdf5(iunit1, & 
     &                                  l_with_guardcells, & 
     &                                  user_attr_1, & 
     &                                  user_attr_2, & 
     &                                  user_attr_3, & 
     &                                  user_attr_4, & 
     &                                  user_attr_5)

      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(out)   :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5
      end subroutine amr_checkpoint_re_hdf5
      end interface


      interface
      subroutine amr_checkpoint_re_mpiio(iunit1, & 
     &                                  l_with_guardcells, & 
     &                                  user_attr_1, & 
     &                                  user_attr_2, & 
     &                                  user_attr_3, & 
     &                                  user_attr_4, & 
     &                                  user_attr_5)

      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(out)   :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5
      end subroutine amr_checkpoint_re_mpiio
      end interface


      interface
      subroutine amr_checkpoint_wr(iunit1, & 
     &                             l_with_guardcells, & 
     &                             check_format,      & 
     &                             user_attr_1, & 
     &                             user_attr_2, & 
     &                             user_attr_3, & 
     &                             user_attr_4, & 
     &                             user_attr_5)

      integer, intent(in)                      :: iunit1
      logical, optional, intent(in)            :: l_with_guardcells
      character (len=80), optional, intent(in) :: check_format
      real, optional, intent(in)               :: user_attr_1, & 
     &                                            user_attr_2, & 
     &                                            user_attr_3, & 
     &                                            user_attr_4, & 
     &                                            user_attr_5

      end subroutine amr_checkpoint_wr
      end interface

      interface
      subroutine amr_checkpoint_wr_default(iunit1, & 
     &                                     l_with_guardcells, & 
     &  user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5)

      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(in) :: & 
     &  user_attr_1,user_attr_2,user_attr_3,user_attr_4,user_attr_5
      end subroutine amr_checkpoint_wr_default
      end interface

      interface
      subroutine amr_checkpoint_wr_hdf5(iunit1, & 
     &                                  l_with_guardcells, & 
     &                                  user_attr_1, & 
     &                                  user_attr_2, & 
     &                                  user_attr_3, & 
     &                                  user_attr_4, & 
     &                                  user_attr_5)

      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(in)    :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5
      end subroutine amr_checkpoint_wr_hdf5
      end interface

      interface
      subroutine amr_checkpoint_wr_mpiio(iunit1, & 
     &                                  l_with_guardcells, & 
     &                                  user_attr_1, & 
     &                                  user_attr_2, & 
     &                                  user_attr_3, & 
     &                                  user_attr_4, & 
     &                                  user_attr_5)

      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(in)    :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5
      end subroutine amr_checkpoint_wr_mpiio
      end interface

      interface
      subroutine amr_plotfile_chombo(iunit1)
      integer, intent(in) :: iunit1
      end subroutine amr_plotfile_chombo
      end interface

      interface
      subroutine amr_plotfile_matlab(iunit1, & 
     &                               l_with_guardcells, & 
     &                               user_attr_1, & 
     &                               user_attr_2, & 
     &                               user_attr_3, & 
     &                               user_attr_4, & 
     &                               user_attr_5)
      integer, intent(in) :: iunit1
      logical, optional, intent(in) :: l_with_guardcells
      real, optional, intent(in)    :: user_attr_1, & 
     &                                 user_attr_2, & 
     &                                 user_attr_3, & 
     &                                 user_attr_4, & 
     &                                 user_attr_5
      end subroutine amr_plotfile_matlab
      end interface

      interface
      subroutine amr_check_derefine(mype)
      integer, intent(in) :: mype
      end subroutine amr_check_derefine
      end interface

      interface
      subroutine amr_check_refine(nprocs,mype,icontinue)
      integer, intent(in) :: nprocs,mype
      integer, intent(in) :: icontinue
      end subroutine amr_check_refine
      end interface


      interface
      subroutine amr_close
      end subroutine amr_close
      end interface

      interface
      subroutine amr_derefine_blocks(lnblocks_old,mype)
      integer, intent(in)    :: mype
      integer, intent(inout) :: lnblocks_old
      end subroutine amr_derefine_blocks
      end interface

      interface
      subroutine amr_compute_morton (mort_no)
      integer, intent(out) ::  mort_no(:,:)
      end subroutine amr_compute_morton 
      end interface


      interface
      subroutine amr_edge_average(mype,lfullblock,nsub)
      integer, intent(in)  ::  mype,nsub
      logical, intent(in)  ::  lfullblock
      end subroutine amr_edge_average
      end interface


      interface
      subroutine amr_edge_average_udt(mype)
      integer, intent(in)  ::  mype
      end subroutine amr_edge_average_udt
      end interface


      interface
      subroutine amr_edge_average_vdt(mype,nsub)
      integer, intent(in)  ::  mype,nsub
      end subroutine amr_edge_average_vdt
      end interface

      interface
      subroutine amr_edge_diagonal_check(mype)
      integer, intent(in)  ::  mype
      end subroutine amr_edge_diagonal_check
      end interface

      interface
      subroutine amr_flush(iunit)
      integer, intent(in)  ::  iunit
      end subroutine amr_flush
      end interface

      interface
      subroutine amr_flux_conserve(mype,nsub,flux_dir,pdgNo)
        implicit none
        integer, intent(in)  ::  mype,nsub
        integer, optional, intent(in) :: flux_dir
        integer, optional, intent(in) :: pdgNo
      end subroutine amr_flux_conserve
      end interface

      interface
      subroutine amr_flux_conserve_udt(mype,pdg,ig,flux_dir)
      use gr_pmPdgDecl, ONLY : pdg_t
      implicit none
      integer, intent(in)  ::  mype
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      integer, optional, intent(in) :: flux_dir
      end subroutine amr_flux_conserve_udt
      end interface

      interface
      subroutine amr_flux_conserve_vdt(mype,nsub,pdg,ig)
      use gr_pmPdgDecl, ONLY : pdg_t
      implicit none
      integer, intent(in)  ::  mype,nsub
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_flux_conserve_vdt
      end interface

      interface
      subroutine amr_gsurrounding_blks(mype,ldiag)
      integer, intent(in)    ::  mype
      logical, intent(in)    ::  ldiag
      end subroutine amr_gsurrounding_blks
      end interface


      interface amr_guardcell
!!$      subroutine amr_guardcell(mype,iopt,nlayers, & 
!!$                               nlayersx,nlayersy,nlayersz,&
!!$                               maxNodetype_gcWanted)
!!$      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
!!$      integer, intent(in)  ::  mype,iopt,nlayers
!!$      integer,OPTIONAL, intent(in) :: maxNodetype_gcWanted
!!$      end subroutine amr_guardcell
      subroutine amr_guardcell_pdgNo(mype,iopt,nlayers, & 
                               nlayersx,nlayersy,nlayersz,&
                               maxNodetype_gcWanted,pdgNo)
        implicit none
        integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
        integer, intent(in)  ::  mype,iopt,nlayers
        integer,OPTIONAL, intent(in) :: maxNodetype_gcWanted
        integer, intent(in), optional :: pdgNo
      end subroutine amr_guardcell_pdgNo
      subroutine amr_guardcell_onePdg(mype,iopt,nlayers, pdg,ig, &
                               nlayersx,nlayersy,nlayersz,&
                               maxNodetype_gcWanted)
        use gr_pmPdgDecl, ONLY : pdg_t
        implicit none
        integer, intent(in)  ::  mype,iopt,nlayers
        type(pdg_t), intent(INOUT) :: pdg
        integer, intent(in) :: ig
        integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
        integer,OPTIONAL, intent(in) :: maxNodetype_gcWanted
      end subroutine amr_guardcell_onePdg
      end interface

      interface
      subroutine amr_initialize
      end subroutine amr_initialize
      end interface

      interface
         subroutine gr_pdgDimenInitOne(pdgDimen, nvar,nguard,nx,ny,nz, npgsArg,k2dArg,k3dArg)
           use gr_pmPdgDecl, ONLY : pdgConst_t
           implicit none
           type(pdgConst_t), intent(OUT) :: pdgDimen
           integer, intent(IN) :: nvar
           integer, intent(IN) :: nguard
           integer, intent(IN) :: nx,ny,nz
           integer, intent(IN),OPTIONAL :: npgsArg,k2dArg,k3dArg
         end subroutine gr_pdgDimenInitOne
      end interface


      interface
      subroutine amr_migrate_tree_data(new_loc,nprocs,mype)
      use tree, ONLY: maxblocks_tr
      integer, intent(in)    ::  mype
      integer, intent(in)    ::  nprocs
      integer, intent(inout) ::  new_loc(2,maxblocks_tr)
      end subroutine amr_migrate_tree_data
      end interface


      interface
      subroutine amr_morton_order (lnblocks_old,nprocs,mype, & 
     &                             l_move_solution, & 
     &                             reorder_grid)
      integer, intent(in) ::  mype
      integer, intent(in) ::  nprocs,lnblocks_old
      logical, intent(in) ::  l_move_solution
      logical, intent(in), optional ::  reorder_grid
      end subroutine amr_morton_order
      end interface


      interface
      subroutine amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt,idest,pdg,ig)
      use gr_pmPdgDecl, ONLY : pdg_t
      implicit none
      integer, intent(in) ::  lb,pe,iopt,idest
      logical, intent(in) ::  lcc,lfc,lec,lnc
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      end subroutine amr_perm_to_1blk
      end interface

      interface
      subroutine amr_mpi_find_blk_in_buffer( & 
     &       mype,remote_block,remote_pe,idest,dtype,index,lfound)
      integer, intent(in)  :: mype,remote_pe,remote_block,idest
      integer, intent(out) :: dtype,index
      logical, intent(out) :: lfound
      end subroutine amr_mpi_find_blk_in_buffer
      end interface


      interface amr_prolong
      subroutine amr_prolong_pdgNo(mype,iopt,nlayers,pdgNo)
      integer, intent(in) ::  mype,iopt,nlayers
      integer, intent(in),optional ::  pdgNo
      end subroutine amr_prolong_pdgNo
      end interface

      interface
         subroutine amr_prolong_cc_fun_init(pdg,ig)
           use gr_pmPdgDecl, ONLY : pdg_t
           implicit none
           type(pdg_t), intent(INOUT) :: pdg
           integer, intent(in)    :: ig
         end subroutine amr_prolong_cc_fun_init
      end interface

      interface
      subroutine amr_prolong_face_fun_init
      end subroutine amr_prolong_face_fun_init
      end interface

      interface
!      subroutine amr_prolong_fc_divbconsist(mype)
      subroutine amr_prolong_fc_divbconsist(mype,level,nfield,ig)
      integer, intent(in) ::  mype
      integer, intent(in) ::  level
      integer, intent(in) ::  nfield
      integer, intent(in) ::  ig
      end subroutine amr_prolong_fc_divbconsist
      end interface


      interface
      subroutine amr_prolong_fun_init
      end subroutine amr_prolong_fun_init
      end interface

      interface
      subroutine amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)
      integer, intent(in)    ::  nprocs
      integer, CONTIGUOUS_INTENT(in) :: new_loc(:,:)
      integer, intent(in)    ::  lnblocks_old
      integer, intent(in)    ::  mype
      end subroutine amr_redist_blk
      end interface

      interface
         Subroutine send_block_data (lb, new_loc, old_loc, free,       &
                                  moved, sent,                         &
                                  lnblocks_old, mype, nmoved,          &
                                  test, point_to,                      &
                                  reqs, nsend, unk_int_types,          &
                                  facex_int_type, facey_int_type,      &
                                  facez_int_type, edgex_int_type,      &
                                  edgey_int_type, edgez_int_type,      &
                                  unkn_int_type)
           use paramesh_dimensions, ONLY: maxblocks
           use tree, ONLY: maxblocks_tr
           implicit none
           Integer,intent(in) :: new_loc(2,maxblocks_tr), old_loc(2,maxblocks_tr)
           Logical,intent(INOUT) :: free(maxblocks), moved(maxblocks), sent(maxblocks)
           Integer,intent(in) :: lb, lnblocks_old, mype
           Integer,intent(INOUT) :: reqs(maxblocks_tr), nsend
           Integer,intent(INOUT) :: nmoved
           Integer,intent(INOUT) :: point_to(maxblocks),test(maxblocks)
           Integer,intent(in) :: unk_int_types(1:NUM_PDGS)
           Integer,intent(in) :: facex_int_type, facey_int_type, facez_int_type
           Integer,intent(in) :: edgex_int_type, edgey_int_type, edgez_int_type
           Integer,intent(in) :: unkn_int_type
         end Subroutine send_block_data
      end interface


      interface
      subroutine amr_refine_blocks (nprocs,mype)
      integer, intent(in)    :: nprocs,mype
      end subroutine amr_refine_blocks
      end interface


      interface
      subroutine amr_refine_derefine(force_rebalance)
      logical,intent(in),optional :: force_rebalance
      end subroutine amr_refine_derefine
      end interface


      interface amr_restrict
!!$      subroutine amr_restrict(mype,iopt,iempty,filling_guardcells)
!!$      integer, intent(in)    :: mype,iopt,iempty
!!$      logical, optional, intent(in) :: filling_guardcells
!!$      end subroutine amr_restrict
         subroutine amr_restrict_pdgNo(mype,iopt,iempty,filling_guardcells,pdgNo)
           implicit none
           integer, intent(in)    :: mype,iopt,iempty
           logical, optional, intent(in) :: filling_guardcells
           integer, intent(in), optional :: pdgNo
         end subroutine amr_restrict_pdgNo
      end interface

      interface
      subroutine amr_restrict_bnd_data(mype,flux_dir,pdg,ig)
      use gr_pmPdgDecl, ONLY : pdg_t
      implicit none
      integer, intent(in)    :: flux_dir
      integer, intent(in)    :: mype
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in)    :: ig
      end subroutine amr_restrict_bnd_data
      end interface

      interface
      subroutine amr_restrict_bnd_data_vdt(mype,ig)
      integer, intent(in)    :: mype
      integer, intent(in)    :: ig
      end subroutine amr_restrict_bnd_data_vdt
      end interface

      interface
      subroutine amr_restrict_edge(icoord)
      integer, intent(in)    :: icoord
      end subroutine amr_restrict_edge
      end interface

      interface
      subroutine amr_restrict_edge_data(mype)
      integer, intent(in)    :: mype
      end subroutine amr_restrict_edge_data
      end interface

      interface
      subroutine amr_restrict_edge_data_vdt(mype)
      integer, intent(in)    :: mype
      end subroutine amr_restrict_edge_data_vdt
      end interface

      interface
      subroutine amr_restrict_ec_fun(recv,temp,icoord)
      integer, intent(in)    :: icoord
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_ec_fun
      end interface

      interface
      subroutine amr_restrict_ec_genorder(recv,temp,icoord,order,ivar)
      integer, intent(in)    :: icoord, order, ivar
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_ec_genorder
      end interface

      interface
      subroutine amr_restrict_ec_user()
      end subroutine amr_restrict_ec_user
      end interface

      interface
      subroutine amr_restrict_fc_fun(recv,temp,icoord)
      integer, intent(in)    :: icoord
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_fc_fun
      end interface

      interface
      subroutine amr_restrict_fc_genorder(recv,temp,icoord,order,ivar)
      integer, intent(in)    :: icoord, order, ivar
      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      end subroutine amr_restrict_fc_genorder
      end interface

      interface
         subroutine amr_restrict_fc_ins(recv,temp,icoord,order,ivar)
           integer, intent(in)    :: icoord, order, ivar
           real,    intent(in)    :: recv(:,:,:,:)
           real,    intent(inout) :: temp(:,:,:,:)
         end subroutine amr_restrict_fc_ins
      end interface

      interface
      subroutine amr_restrict_fc_user()
      end subroutine amr_restrict_fc_user
      end interface

      interface
      subroutine amr_restrict_red(icoord)
      integer, intent(in)    :: icoord
      end subroutine amr_restrict_red
      end interface

      interface
      subroutine amr_restrict_unk_fun(datain,dataout, ig)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ig
      end subroutine amr_restrict_unk_fun
      end interface

      interface
      subroutine amr_restrict_unk_genorder(datain,dataout,order,ivar,ig)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: order, ivar
      integer, intent(in) :: ig
      end subroutine amr_restrict_unk_genorder
      end interface

      interface
      subroutine amr_restrict_unk_user()
      end subroutine amr_restrict_unk_user
      end interface

      interface
      subroutine amr_restrict_unk_dg(datain,dataout,ivar)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ivar
      end subroutine amr_restrict_unk_dg
      end interface

      interface
      subroutine amr_restrict_nc_fun(datain,dataout)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      end subroutine amr_restrict_nc_fun
      end interface

      interface
      subroutine amr_restrict_nc_genorder(datain,dataout,ivar)
      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: ivar
      end subroutine amr_restrict_nc_genorder
      end interface

      interface
      subroutine amr_restrict_nc_user()
      end subroutine amr_restrict_nc_user
      end interface

      interface
      subroutine amr_restrict_work_fun(datain,dataout,iopt)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      integer, intent(in) :: iopt
      end subroutine amr_restrict_work_fun
      end interface

      interface
      subroutine amr_restrict_work_genorder(datain,dataout,iopt,order)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      integer, intent(in) :: iopt, order
      end subroutine amr_restrict_work_genorder
      end interface

      interface
      subroutine amr_restrict_work_user()
      end subroutine amr_restrict_work_user
      end interface

      interface
      subroutine amr_restrict_work_fun_recip(datain,dataout)
      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)
      end subroutine amr_restrict_work_fun_recip
      end interface


      interface
      subroutine amr_ser_distribute (nprocs,mype,lnblocks_old)
      integer, intent(in)  ::  nprocs,mype,lnblocks_old
      end subroutine amr_ser_distribute
      end interface


      interface
      subroutine amr_sort_by_work (new_loc,nprocs,mype)
      integer, intent(in)    :: mype
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  nprocs
      end subroutine amr_sort_by_work
      end interface

      interface
      subroutine morton_sort(mort_no,ix,iend)
      integer, intent(in) :: iend
      integer, intent(inout) :: mort_no(6,iend), ix(iend)
      end subroutine morton_sort
      end interface

      interface
      subroutine amr_sort_morton (mort_no,new_loc,nprocs)
      integer, intent(inout) ::  mort_no(:,:)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  nprocs
      end subroutine amr_sort_morton
      end interface

      interface
      subroutine amr_sort_morton_reorder_grid (mort_no,new_loc,nprocs)
      integer, intent(inout) ::  mort_no(:,:)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  nprocs
      end subroutine amr_sort_morton_reorder_grid
      end interface

      interface
      subroutine amr_surrounding_blks(mype,pe,lb,surrblks,ldiag)
      integer, intent(in)    ::  mype,pe,lb
      integer, intent(inout) ::  surrblks(2,3,3,3)
      logical, intent(in)    ::  ldiag
      end subroutine amr_surrounding_blks
      end interface

      interface
      subroutine amr_test_refinement(mype,lrefine_min,lrefine_max)
      integer, intent(in)    ::  mype,lrefine_min,lrefine_max
      end subroutine amr_test_refinement
      end interface

      interface
      subroutine comm_finish
      end subroutine comm_finish
      end interface

      interface
      subroutine comm_start(MaxProcs,nprocs,mype)
      integer, intent(out) :: nprocs,mype
      integer, intent(in)  :: MaxProcs
      end subroutine comm_start
      end interface

      interface
      subroutine comm_int_sum_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_sum_to_all
      end interface

      interface
      subroutine comm_int_min_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_min_to_all
      subroutine comm_int_min_to_all1(i)
        integer, intent(inout) :: i
      end subroutine comm_int_min_to_all1
      end interface

      interface
      subroutine comm_int_max_to_all(target,source)
      integer, intent(in) :: source
      integer, intent(out)  :: target
      end subroutine comm_int_max_to_all
      subroutine comm_int_max_to_all1(i)
        integer, intent(inout) :: i
      end subroutine comm_int_max_to_all1
      end interface


      interface
      subroutine comm_real_sum_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_sum_to_all
      end interface

      interface
      subroutine comm_real_min_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_min_to_all
      end interface

      interface
      subroutine comm_real_max_to_all(target,source)
      real, intent(in) :: source
      real, intent(out)  :: target
      end subroutine comm_real_max_to_all
      end interface

      interface
      subroutine test_neigh_data(mype,istep)
      integer, intent(in)    ::  mype,istep
      end subroutine test_neigh_data
      end interface

      interface
      subroutine fill_old_loc(new_loc,old_loc,nprocs,mype)
      integer, intent(in)    :: mype,nprocs
      integer, intent(in)    :: new_loc(:,:)
      integer, intent(out)   :: old_loc(:,:)
      end subroutine fill_old_loc
      end interface

      interface
      subroutine gtest_neigh_data(mype,istep,test_a)
      integer, intent(in)    ::  mype,istep
      real,    intent(in)    ::  test_a
      end subroutine gtest_neigh_data
      end interface

      interface
      subroutine mesh_test(mype)
      integer, intent(in)    ::  mype
      end subroutine mesh_test
      end interface

      interface
      subroutine guardcell_test(mype)
      integer, intent(in)    ::  mype
      end subroutine guardcell_test
      end interface

      interface 
      subroutine init_sparse_solver
      end subroutine init_sparse_solver
      end interface

      interface
      subroutine amr_1blk_fc_clean_divb(                                 & 
     &        nfacevar_in,                                               & 
     &        ia,ib,ja,jb,ka,kb,                                         & 
     &        ionea,ioneb,                                               & 
     &        jonea,joneb,                                               & 
     &        konea,koneb,                                               & 
     &        idest,ioff,joff,koff,                                      & 
     &        mype,lb,parent_pe,parent_blk )
      integer, intent(in) :: nfacevar_in
      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: ionea,ioneb
      integer, intent(in) :: jonea,joneb
      integer, intent(in) :: konea,koneb
      integer, intent(in) :: idest, ioff, joff, koff
      integer, intent(in) :: mype, lb, parent_pe, parent_blk
      end subroutine amr_1blk_fc_clean_divb 
      end interface

      interface
      subroutine prol_fc_clean_divb_test(flag)
      logical, intent(in) :: flag
      end subroutine prol_fc_clean_divb_test 
      end interface

      interface
      subroutine prol_fc_clean_divb_test_report(nerrors)
      integer, intent(inout) :: nerrors
      end subroutine prol_fc_clean_divb_test_report
      end interface

      interface 
      subroutine amr_q_sort (ix,n,ia,ib)
      integer, intent(in) :: n
      integer, dimension(n),  intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib
      end subroutine amr_q_sort
      end interface

      interface 
      subroutine amr_q_sort_real (ix,n,ia,ib)
      integer, intent(in) :: n
      real,    dimension(n),  intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib
      end subroutine amr_q_sort_real
      end interface

      end module paramesh_interfaces
