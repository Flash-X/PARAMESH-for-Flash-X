!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/gr_amr1blkGcToPerm
!! NAME
!!
!!   gr_amr1blkGcToPerm
!!
!! SYNOPSIS
!!
!!   call gr_amr1blkGcToPerm(mype, iopt, nlayers, lb, lcc,lfc,lec,lnc, pdg,ig)
!!   call gr_amr1blkGcToPerm(mype, iopt, nlayers, lb, lcc,lfc,lec,lnc, pdg,ig,nlayersx,nlayersy,nlayersz)
!!
!!   call gr_amr1blkGcToPerm(integer, integer, integer, integer,
!!                           logical, logical, logical, logical,
!!                           TYPE(pdg_t), integer,
!!                     optional integer, optional integer, optional integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
!!
!!   integer, intent(in) :: lb           
!!        The local ID of the block for which data should be copied..
!!
!!   integer, intent(in) :: iopt  
!!     Selects whether to fill the guardcells for the arrays unk, 
!!     facevarx, facevary, facevarz, unk_e_x, unk_e_y, unk_e_z, and unk_n 
!!     (if iopt = 1) or work (if iopt 2).
!!
!!   integer, intent(in) :: nlayers  
!!     Dummy variable which does nothing.  Included for consistency with older 
!!     PARAMESH versions.
!!
!!   logical, intent(in) :: lcc, lfc, lec, lnc
!!        Logical switches which indicate which data is to be copied.
!!         lcc -> cell centered
!!         lfc -> face centered
!!         lec -> edge centered     (FLASH: unused)
!!         lnc -> node centered     (FLASH: unused)
!!
!!   optional, integer, intent(in) :: nlayersx, nlayersy, nlayersz 
!!     Optional integers which select how many guardcell layers to fill in each
!!     direction.  If these variables are not passed in, then the default is to
!!     fill all allocated guardcell space.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   gr_pmPdgDecl
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   gr_pmFlashHookData
!!
!! CALLS
!! 
!!   (amr_abort)
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells data for block lb
!!   have been copied from 1blk storage to permanent storage.
!!
!! DESCRIPTION
!!
!!   Routine to copy guard cell data for one block from the oneBlk arrays
!!   (such as unk1, facevarx1, etc.) to permanent storage (such as unk,
!!   facevarx, etc.).
!!   This subroutine is meant to be called after amr_1blk_guardcell.
!!
!!   Some special handling for face variables is included, this has not
!!   been tested.
!!
!!   The implementation has been extracted from mpi_amr_guardcell.F90.
!!
!! AUTHORS
!!
!!         mpi_amr_guardcell Peter MacNeice (1997) with modifications by Kevin Olson
!!  Extracted: gr_amr1blkGcToPerm   Klaus Weide                   Jan 2021
!!
!!  2022-12-12 K. Weide  Consolidating PDG and PmAsync features
!!***

!!REORDER(5): unk1, facevar[xyz]1
!!REORDER(5): unk, facevar[xyz]
!!REORDER(5): unk_e_[xyz]1, unk_n1
!!REORDER(5): unk_e_[xyz], unk_n

#include "paramesh_preprocessor.fh"

Subroutine gr_amr1blkGcToPerm(mype,iopt,nlayers,lb,         &
                                    lcc,lfc,lec,lnc,                   & 
                               pdg,ig,                      &
                               nlayersx,nlayersy,nlayersz)

!-----Use Statements
  use gr_pmPdgDecl, ONLY : pdg_t
  use paramesh_dimensions, ONLY: gr_thePdgDimens
  use paramesh_dimensions, only: ndim,k2d,k3d,nguard_work,npgs, nfacevar,nvarcorn,nvaredge

      Use physicaldata, only: facevarx,   facevary,   facevarz, &
                              gt_facevarx,gt_facevary,gt_facevarz, &
                              facevarx1,  facevary1,  facevarz1
      Use physicaldata, only: int_gcell_on_cc,int_gcell_on_fc,int_gcell_on_ec,int_gcell_on_nc
      Use physicaldata, only:                     gcell_on_fc,    gcell_on_ec,    gcell_on_nc
      Use physicaldata, only: force_consistency
      Use workspace, ONLY: work1, work
      Use tree,      ONLY: lnblocks, neigh

      use gr_pmFlashHookData, ONLY : alwaysFcAtBC => gr_pmAlwaysFillFcGcAtDomainBC

      Implicit none

!-----Include statements
      Include 'mpif.h'

!-----Input/Output Arguments
      Integer, intent(in) :: mype,iopt,nlayers,lb
      Logical, VALUE      :: lcc,lfc,lec,lnc
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      Integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!-----Local variables
      Logical :: l_force_consist
      Integer :: id,jd,kd
      Integer :: ilays,jlays,klays
      Integer :: nlayers0x, nlayers0y, nlayers0z, nguard0
      Integer :: nguard_npgs
      Integer :: i,j,k,ivar                                  
      Integer :: ip1,ip2,jp1,jp2,kp1,kp2
      Integer :: ilp,iup,jlp,jup,klp,kup
      Integer :: iu, ju, ku, iopt0

!------------------------------------
!-----Begin Executable code section
!------------------------------------

  ASSOCIATE(nxb         => gr_thePdgDimens(ig) % nxb,      &
            nyb         => gr_thePdgDimens(ig) % nyb,      &
            nzb         => gr_thePdgDimens(ig) % nzb,      &
            nguard      => gr_thePdgDimens(ig) % nguard,   &
            nvar        => gr_thePdgDimens(ig) % nvar,     &
            unk         => pdg % unk,      &
            unk1        => pdg % unk1,      &
            gcell_on_cc => pdg % gcell_on_cc      &
            )

      If (iopt == 1) Then

!--------set users selections of guardcell variables
         int_gcell_on_cc = gcell_on_cc
         int_gcell_on_fc = gcell_on_fc
         int_gcell_on_ec = gcell_on_ec
         int_gcell_on_nc = gcell_on_nc

         If (.not.present(nlayersx)) Then
            nlayers0x = nguard
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard
         Else
            nlayers0z = nlayersz
         End if
      Else
         If (.not.present(nlayersx)) Then
            nlayers0x = nguard_work
         Else
            nlayers0x = nlayersx
         End if
         If (.not.present(nlayersy)) Then
            nlayers0y = nguard_work
         Else
            nlayers0y = nlayersy
         End if
         If (.not.present(nlayersz)) Then
            nlayers0z = nguard_work
         Else
            nlayers0z = nlayersz
         End if
      End if ! End If (iopt == 1)

      If (iopt == 1) Then
        nguard0 = nguard
      Else
        nguard0 = nguard_work
      End if

      nguard_npgs = nguard*npgs

!------make sure that nlayers and iopt are set consistently.
       If (iopt == 1.and.nlayers.ne.nguard) Then
         If (mype == 0) Then
           Write(*,*) 'PARAMESH ERROR !'
           Write(*,*) 'Error in guardcell - iopt and nlayers'
           Write(*,*) 'are not consistent. For iopt=1 you must'
           Write(*,*) 'set nlayers=nguard.'
         Endif
         Call amr_abort
        Else If (iopt >= 2.and.nlayers > nguard_work) Then
         If (mype == 0) Then
           Write(*,*) 'PARAMESH ERROR !'
           Write(*,*) 'Error in guardcell - iopt and nlayers'
           Write(*,*) 'are not consistent. For iopt>=2 you must'
           Write(*,*) 'set nlayers le nguard_work.'
         endif
         Call amr_abort
       Endif   ! end If (iopt == 1.and.nlayers.ne.nguard)

      If (iopt == 1) Then
        If (nvar > 0) lcc = .True.
        If (nfacevar > 0) lfc = .True.
        If (nvaredge > 0) lec = .True.
        If (nvarcorn > 0) lnc = .True.
      Else If (iopt >= 2) Then
        lcc = .True.
      Endif

      l_force_consist = .False.
      If (force_consistency) Then
       l_force_consist = .True.
       If (lfc) Then
           gt_facevarx(:,1,:,:,lb) = facevarx(:,1+nguard_npgs,:,:,lb)
           gt_facevarx(:,2,:,:,lb) =                                   & 
                facevarx(:,nxb+1+nguard_npgs,:,:,lb)
          If (ndim >= 2) Then
           gt_facevary(:,:,1,:,lb) =                                   & 
                   facevary(:,:,1+nguard_npgs*k2d,:,lb)
           gt_facevary(:,:,1+k2d,:,lb) =                               & 
                facevary(:,:,nyb+(1+nguard_npgs)*k2d,:,lb)
          End if
          If (ndim == 3) Then
           gt_facevarz(:,:,:,1,lb) =                                   & 
                   facevarz(:,:,:,1+nguard_npgs*k3d,lb)
           gt_facevarz(:,:,:,1+k3d,lb) =                               & 
                facevarz(:,:,:,nzb+(1+nguard_npgs)*k3d,lb)
          End if
       End if   ! end if (lfc)
      End if    ! end if force_consistency)


      If (lnblocks > 0) Then

        Do k = 1,1+2*k3d
         klp = 0
         kup = 0
         If (k == 1) Then
           klays = nlayers0z*k3d
           kd = nguard0*k3d+1 - klays
           kp1 = 0
           kp2 = 0
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,5,lb) .le. -20))) kup = k3d
         Else if (k == 2) Then
           klays = nzb*k3d
           kd = nguard0*k3d+1
           kp1 = 0
           kp2 = k3d
         Else !if (k == 3) Then
           klays = nlayers0z*k3d
           kd = (nguard0+nzb)*k3d + 1
           kp1 = k3d
           kp2 = k3d
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,6,lb) .le. -20))) klp = -k3d
         End if
         ku = kd + klays - k3d
        Do j = 1,1+2*k2d
         jlp = 0
         jup = 0
         If (j == 1) Then
           jlays = nlayers0y*k2d
           jd = nguard0*k2d+1 - jlays
           jp1 = 0
           jp2 = 0
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,3,lb) .le. -20))) jup = k2d
         Else if (j == 2) Then
           jlays = nyb*k2d
           jd = nguard0*k2d+1
           jp1 = 0
           jp2 = k2d
         Else !if (j == 3) Then
           jlays = nlayers0y*k2d
           jd = (nguard0+nyb)*k2d + 1
           jp1 = k2d
           jp2 = k2d
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,4,lb) .le. -20))) jlp = -k2d
         End if
         ju = jd + jlays - k2d
       Do i = 1,3
         ilp = 0
         iup = 0
         If (i == 1) Then
           ilays = nlayers0x
           id = nguard0+1 - ilays
           ip1 = 0
           ip2 = 0
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,1,lb) .le. -20))) iup = 1
         Else if (i == 2) Then
           ilays = nxb
           id = nguard0+1
           ip1 = 0
           ip2 = 1
         Else !if (i == 3) Then
           ilays = nlayers0x
           id = nguard0+nxb + 1
           ip1 = 1
           ip2 = 1
           If ((l_force_consist) .or. (alwaysFcAtBC .AND. (neigh(1,2,lb) .le. -20))) ilp = -1
         End if
         iu = id + ilays - 1

         If (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d) Then

         Else

          If (lcc) Then
           If (iopt == 1) Then
             Do ivar=1,nvar
               If (int_gcell_on_cc(ivar)) Then
                 unk(ivar,id:iu,jd:ju,kd:ku,lb) =                      & 
                  unk1(ivar,id:iu,jd:ju,kd:ku,1)
               Endif
             Enddo
           Else
            iopt0 = iopt-1
            work(id:iu,jd:ju,kd:ku,lb,iopt0) =                         & 
              work1(id:iu,jd:ju,kd:ku,1)
           End if   ! end if iopt == 1
          End if    ! end if (lcc)

          If (lfc) Then
           Do ivar = 1,nfacevar
            If (int_gcell_on_fc(1,ivar)) Then
             facevarx( ivar,id+ip1+ilp:iu+ip2+iup,                     & 
                            jd:ju,kd:ku,lb) =                          & 
             facevarx1( ivar,id+ip1+ilp:iu+ip2+iup,                    & 
                             jd:ju,kd:ku,1)
            End if
            If (ndim > 1) Then
             If (int_gcell_on_fc(2,ivar)) Then
              facevary( ivar,id:iu,jd+jp1+jlp:ju+jp2+jup,              & 
                             kd:ku,lb) =                               & 
              facevary1(ivar,id:iu,jd+jp1+jlp:ju+jp2+jup,              & 
                              kd:ku,1)
             End if
             If (ndim == 3) Then
              If (int_gcell_on_fc(3,ivar)) Then
               facevarz( ivar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,lb) =  & 
                facevarz1(ivar,id:iu,jd:ju,kd+kp1+klp:ku+kp2+kup,1)
              End if
             End if  ! end if (ndim == 3)
            End if   ! end if (ndim > 1)
           End do    ! end do ivar = 1, nfacevar
          End if     ! end if (lfc)

#ifdef FLASH_PMFEATURE_UNUSED
          If (lec) Then
           Do ivar = 1, nvaredge
            If (ndim > 1) Then
             If (int_gcell_on_ec(1,ivar)) Then
              unk_e_x( ivar,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,lb) =    & 
              unk_e_x1(ivar,id:iu,jd+jp1:ju+jp2,kd+kp1:ku+kp2,1)
             End if
             If (int_gcell_on_ec(2,ivar)) Then
              unk_e_y( ivar,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,lb) =    & 
              unk_e_y1(ivar,id+ip1:iu+ip2,jd:ju,kd+kp1:ku+kp2,1)
             End If
             If (ndim == 3) Then
              If (int_gcell_on_ec(3,ivar)) Then
               unk_e_z( ivar,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,lb) =   & 
               unk_e_z1(ivar,id+ip1:iu+ip2,jd+jp1:ju+jp2,kd:ku,1)
              End If
             End If ! End If (ndim == 3)
            End If  ! End If (ndim > 1)
           End Do   ! End Do ivar = 1, nvaredge
          End If    ! End If (lec)

          If (lnc) Then
           Do ivar = 1, nvarcorn
            If (int_gcell_on_nc(ivar)) Then
             unk_n( ivar,                                              & 
                 id+ip1:iu+ip2,                                        & 
                 jd+jp1:ju+jp2,                                        & 
                 kd+kp1:ku+kp2,lb) =                                   & 
                  unk_n1(ivar,                                         & 
                         id+ip1:iu+ip2,                                & 
                         jd+jp1:ju+jp2,                                & 
                         kd+kp1:ku+kp2,1)
            End If ! End If (int_gcell_on_nc(ivar))
           End Do  ! End Do ivar = 1, nvarcorn
          End If   ! Dnd If (lnc)
#endif

        End If     ! end if (i  ==  2 .and. j  ==  1+k2d .and. k  ==  1+k3d)

      End Do  ! End Do i = 1,3
      End Do  ! End Do j = 1,1+2*k2d
      End Do  ! End Do k = 1,1+2*k3d

      End if  ! End If (lnblocks >= 1)
    end ASSOCIATE
      Return
    End Subroutine gr_amr1blkGcToPerm


