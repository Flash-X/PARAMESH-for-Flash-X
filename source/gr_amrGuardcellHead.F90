!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_guardcell
!! NAME
!!
!!   amr_guardcell
!!
!! SYNOPSIS
!!
!!  use gr_itParam_mod, ONLY: gr_itParam_t
!!  use gr_pmBlockGetter, ONLY: gr_pmBlockGetter_t
!!
!!   call gr_amrGuardcellHead(mype, iopt, nlayers)
!!   call gr_amrGuardcellHead(mype, iopt, nlayers, nlayersx, nlayersy, nlayersz)
!!   call gr_amrGuardcellHead(mype, iopt, nlayers, nlayersx, nlayersy, nlayersz, maxNodetype_gcWanted)
!!
!!   call gr_amrGuardcellHead(integer, integer, integer, 
!!                     optional integer, optional integer, optional integer, optional integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(in) :: mype  
!!     The calling processor.
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
!!   paramesh_dimensions
!!   physicaldata
!!   workspace
!!   tree
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!! 
!!   amr_1blk_guardcell_reset
!!   amr_restrict
!!   amr_1blk_guardcell
!!   mpi_amr_comm_setup
!!    
!! RETURNS
!!
!!   Does not return anything.  Upon exit, guardcells are filled with data for 
!!   all blocks.
!!
!! DESCRIPTION
!!
!!   Routine to perform guardcell filling in the case that permanent guardcells
!!   are user.  Calling this routine will result in the guardcells being filled
!!   for all child blocks if 'advance_all_levels' is turned OFF or in ALL the blocks
!!   having their guardcells filled if 'advance_all_levels' is turned ON.
!!   The user can choose how many layers of guardcells are filled on each face
!!   of a block by specifying the optional arguments 'nlayersx', 'nlayersy', and
!!   'nlayersz'.
!!
!! AUTHORS
!!
!!   Peter MacNeice (1997) with modifications by Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
#include "constants.h"

Module gr_amrGuardcellHead_mod

contains
    Subroutine gr_amrGuardcellHead(itContext, getter,mype,iopt,nlayers,         & 
                               nlayersx,nlayersy,nlayersz, &
                               maxNodetype_gcWanted)

!-----Use Statements
      use Timers_interface, ONLY : Timers_start, Timers_stop
      Use paramesh_dimensions
      Use physicaldata
      Use workspace
      Use tree
      use paramesh_comm_data

      Use paramesh_interfaces, Only : amr_1blk_guardcell_reset, & 
                                      amr_restrict,             & 
                                      amr_1blk_guardcell

!!$      Use paramesh_mpi_interfaces, Only : mpi_amr_comm_setup
      Use gr_mpiAmrComm_mod,   only : gr_mpiAmrComm
      use gr_itParam_mod, ONLY: gr_itParam_t
      use gr_pmBlockGetter,   ONLY: gr_pmBlockGetter_t, gr_pmBlockGetterDestroy

      use gr_pmFlashHookData, ONLY : alwaysFcAtBC => gr_pmAlwaysFillFcGcAtDomainBC

      Implicit none

!-----Include statements
      Include 'mpif.h'

!-----Input/Output Arguments
      type(gr_itParam_t), intent(INOUT) :: itContext
      type(gr_pmBlockGetter_t),intent(INOUT), TARGET :: getter
      Integer, intent(in) :: mype,iopt,nlayers
      Integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      Integer, intent(in), optional :: maxNodetype_gcWanted

!-----Local variables
      Logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      Logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag,l_force_consist
      Integer :: lb,icoord
      Integer :: id,jd,kd
      Integer :: ilays,jlays,klays
      Integer :: nlayers0x, nlayers0y, nlayers0z, nguard0
      Integer :: nguard_npgs
      Integer :: i,j,k,ivar                                  
      Integer :: ip1,ip2,jp1,jp2,kp1,kp2
      Integer :: ilp,iup,jlp,jup,klp,kup
      Integer :: nprocs, ierr, tag_offset, iempty, iu, ju, ku, iopt0
      Integer :: maxNodetype_gcWanted_loc
      integer :: ntype,lev

!------------------------------------
!-----Begin Executable code section
!------------------------------------

      If ((.not.diagonals) .and. (iopt .ne. 2)) Then
         Write(*,*) 'gr_amrGuardcellHead:  diagonals off'
      End if

      If (present(maxNodetype_gcWanted)) Then
         maxNodetype_gcWanted_loc = maxNodetype_gcWanted
      Else
         maxNodetype_gcWanted_loc = -1
      End If

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

      call MPI_COMM_SIZE(amr_mpi_meshComm,nprocs,ierr)

      If (no_permanent_guardcells) Then

       If (mype == 0) Then
         Write(*,*) 'gr_amrGuardcellHead call ignored!'
         Write(*,*) 'NO_PERMANENT_GUARDCELLS is defined'
       End if
       Return

      Else  ! no_permanent_guardcells

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

!-----Reinitialize addresses of cached parent blocks
      Call amr_1blk_guardcell_reset

      lcc = .False.
      lfc = .False.
      lec = .False.
      lnc = .False.
      If (iopt == 1) Then
        If (nvar > 0) lcc = .True.
        If (nfacevar > 0) lfc = .True.
        If (nvaredge > 0) lec = .True.
        If (nvarcorn > 0) lnc = .True.
      Else If (iopt >= 2) Then
        lcc = .True.
      Endif

!-----Restrict solution to parent blocks
      If (.not.advance_all_levels) Then
         iempty = 0
         call amr_restrict(mype,iopt,iempty,.True.)
         call amr_1blk_guardcell_reset
         Call MPI_BARRIER(amr_mpi_meshComm, ierr)
      End if


      l_force_consist = .False.
      If (force_consistency) Then
       l_force_consist = .True.
       If (lfc) Then
        Do lb = 1,lnblocks
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
        End do  ! end do lb = 1, lnblocks
       End if   ! end if (lfc)
      End if    ! end if force_consistency)

      tag_offset = 100
      lguard    = .True.
      lprolong  = .False.
      lflux     = .False.
      ledge     = .False.
      lrestrict = .False.
      lfulltree = .False.
      lev = UNSPEC_LEVEL        ! or as requested??
      if (maxNodetype_gcWanted_loc < 0) then
         ntype = ACTIVE_BLKS
      else if (maxNodetype_gcWanted_loc == 1) then
         ntype = LEAF
      else if (maxNodetype_gcWanted_loc == 2) then
         ntype = ACTIVE_BLKS
      else if (maxNodetype_gcWanted_loc .GE. 3) then
         ntype = ALL_BLKS
      else if (advance_all_levels) then
         ntype = ALL_BLKS
      else
         ntype = ACTIVE_BLKS
      end if
      call Timers_start("gr_mpiAmrComm")
      Call gr_mpiAmrComm(mype,nprocs,                             &
                              lguard,lprolong,lflux,ledge,lrestrict,   & 
                              lfulltree,                               & 
                              iopt,lcc,lfc,lec,lnc,tag_offset,         &
                              getter,                                  &
                              ntype, level=lev,                        &
                              nlayersx=nlayersx,nlayersy=nlayersy,nlayersz=nlayersz)
      call Timers_stop("gr_mpiAmrComm")


      Endif ! If (no_permanent_guardcells)

      itContext % mype = mype
      itContext % iopt = iopt

      Return
    End Subroutine gr_amrGuardcellHead
end Module gr_amrGuardcellHead_mod
