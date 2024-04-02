!----------------------------------------------------------------------
!!  This file is from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright 2023 UChicago Argonne, LLC and contributors
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!  2023-03-16 K. Weide  Pass ioff,joff,koff to amr_restrict_unk_fun
!!  2023-03-17 K. Weide  Special handling for interp_mask_unk_res 40

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

#ifdef DEBUG_ALL
#define DEBUG
#endif

! This routine should be moved to 'utilities', not part of the PARAMESH kernal

      subroutine mpi_amr_restrict_fulltree(mype,iopt,lcc,lfc,lec,lnc)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use paramesh_comm_data

      use paramesh_interfaces, only : amr_1blk_copy_soln, & 
     &                                amr_1blk_guardcell_reset, & 
     &                                amr_perm_to_1blk, & 
     &                                amr_1blk_guardcell, & 
     &                                amr_restrict_unk_fun, & 
     &                                amr_restrict_nc_fun, & 
     &                                amr_restrict_fc_fun, & 
     &                                amr_restrict_ec_fun, & 
     &                                amr_restrict_work_fun, & 
     &                                amr_restrict_work_fun_recip, & 
     &                                comm_int_max_to_all, & 
     &                                comm_int_min_to_all & 
     &                                ,amr_block_geometry

      use paramesh_mpi_interfaces, only :  & 
     &                                mpi_amr_comm_setup
  use gr_flashHook_interfaces

!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data. 
! The parent get data for each child and them applies
! the restriction operator to it. Guardcell data may be needed for the
! child blocks, depending on the particular restriction operator being used.
! Thus amr_1blk_guardcell is called below.
! This routine calls a user provided routine called restrict_fun
! which defines the pattern of restriction which the user wishes to
! apply.
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------
#include "Flashx_mpi_implicitNone.fh"
      integer, intent(in)  :: mype,iopt
      logical, intent(in)  :: lcc,lfc,lec,lnc

#ifdef REAL8
      real,parameter :: tiny = 1.d-40
#else
      real,parameter :: tiny = 1.e-25
#endif

!-----Local Variables and Arrays
      integer nguard0,nguard_work0

!------------------------------------
! local arrays

      real temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
      real send(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      real recvn0(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d, & 
     &                                      kl_bnd:ku_bnd+k3d)
    Real tempn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       &
                                          kl_bnd1:ku_bnd1+k3d)
    Real sendn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,       &
                                          kl_bnd1:ku_bnd1+k3d)

      integer ::  maxbnd
  Real,Allocatable :: tempf(:,:,:,:)
  Real,Allocatable :: sendf(:,:,:,:)

      logical l_srl_only,ldiag
      logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      logical :: lfound

      integer remote_pe0,remote_block0
      integer remote_pe,remote_block,icoord,nprocs, ierror
      integer cnodetype,cempty

      integer llrefine_min,llrefine_max

      integer :: lb,level,ich,jchild,ioff,joff,koff,nlayers
      integer :: idest,i,j,k,ii,jj,kk,ivar,iopt0,jface,ng0
      integer :: ia,ja,ka,ib,jb,kb,isa,isb,jsa,jsb,ksa,ksb

      integer :: tag_offset
      integer :: iblk, iblock, ilays, jlays, klays
      integer :: id, jd, kd, is, js, ks
      integer :: ip1, jp1, kp1, ip3, jp3, kp3
      integer :: ng1, ndel
      integer :: i1,i2,j1,j2,k1,k2
  logical,allocatable :: curvilinear_cons_mask_unk(:)

!------------------------------------

!-----Begin Executable Code
#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'entered mpi_amr_restrict_fulltree: pe ',mype, & 
     &           ' iopt ',iopt
#endif /* DEBUG_FLOW_TRACE */


  nguard0 = nguard*npgs
  nguard_work0 = nguard_work*npgs

      maxbnd = max(nbndvare,nbndvarc,nbndvar)
      allocate(  & 
     & tempf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &       kl_bnd1:ku_bnd1+k3d))
      allocate( & 
     & sendf(maxbnd,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &       kl_bnd1:ku_bnd1+k3d))

  if (iopt == 1) then
     allocate(curvilinear_cons_mask_unk(size(interp_mask_unk_res,1)))
     if (lcc .AND. curvilinear_conserve) then
        curvilinear_cons_mask_unk(:) = (interp_mask_unk_res(:) .NE. 40)
     else
        curvilinear_cons_mask_unk(:) = .FALSE.
     end if
  end if

      if (.not.diagonals) then
         write(*,*) 'amr_1blk_restrict:  diagonals off'
      end if

! For cell-corner data, during a restriction operation, the
! data on a block boundary shared with a neighbor at the same
! refinement level, needs to be acquired during the operation
! of amr_1blk_guardcell_srl for the parent during the call to
! amr_1blk_guardcell_c_to_f. This flag tells amr_1blk_guardcell_srl
! to get this data. Ordinarily amr_1blk_guardcell_srl does not 
! get this data.
      lrestrict_in_progress = .true.

  call amr_1blk_guardcell_reset()

  call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierror)

!
! Make sure the gt_unk, gt_facevarx, etc copy of the solution exists
! if NO_PERMANENT_GUARDCELL is defined. This may be needed to fill
! guardcells if the restriction operator needs guardcell data.
  level = -1
      call amr_1blk_copy_soln(level)




!-----Cycle through parents in decreasing order of refinement
      llrefine_max = maxval(lrefine)
      llrefine_min = llrefine_max
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      llrefine_min = min(lrefine(lb),llrefine_min)
      enddo
      endif
      call comm_int_max_to_all1 (llrefine_max)
      call comm_int_min_to_all1 (llrefine_min)


      if(llrefine_max.gt.llrefine_min) then
     Do level = llrefine_max-1,llrefine_min,-1


! Now parents of leaf nodes get data
! from their children and then perform restriction on it.

! Setup the communications required to perform the necessary guardcell
! filling and restriction operations.
! Note: as coded, the call to mpi_morton_bnd_restrict(in the calling
! routine mpi_amr_restrict) will select
! all communications needed for restriction at all levels. We
! could try to modify mpi_morton_bnd_restrict so it acted on
! a specified level, which would eliminate some unnecessary
! communications.
! 
        tag_offset = 100

      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .true.
      lfulltree = .false.  ! Must be set to .false. in order to use
                           ! the correct communication pattern!
        Call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,           &
     &                        lflux,ledge,lrestrict,lfulltree, & 
     &                        iopt,lcc,lfc,lec,lnc,tag_offset)

        If (lnblocks > 0) Then
           Do lb = 1,lnblocks


!-----Is this a parent block of at least one leaf node?
      if(nodetype(lb).ge.2.and.lrefine(lb).eq.level) then


!-----If yes then cycle through its children.
                 Do ich=1,nchild

        jchild = ich


!-------Is this child a leaf block? If it is then fetch its data.
                    remote_pe     = child(2,ich,lb)
                    remote_block  = child(1,ich,lb)

#ifdef DEBUG
         write(*,*) 'pe ',mype,' blk ',lb,' searching for child ', & 
     &     jchild ,' at input address ',remote_block,remote_pe, & 
     &     ' lnblocks ',lnblocks
#endif /* DEBUG */


!-------if (remote_block,remote_pe) is not a local block then it must have a 
!-------local copy available in the buffer space at the end of the local
!-------block list.
        if(remote_pe.ne.mype) then
        lfound = .false.
        iblk = strt_buffer
!          do iblk = strt_buffer,last_buffer
          do while(.not.lfound.and.iblk.le.last_buffer)
#ifdef DEBUG
             write(*,*) 'mpi_amr_restrict_fulltree: pe ', & 
     &       mype,' blk ',lb, & 
     &     ' searching buffer for ',ich,' child ', & 
     &     child(:,ich,lb),' current buffer entry ', & 
     &          ' iblk ',iblk,' laddress ',laddress(:,iblk) & 
     &          ,' child(:,:,lb) ',child(:,:,lb)
#endif /* DEBUG */
            if(remote_block.eq.laddress(1,iblk).and. & 
     &           remote_pe .eq.laddress(2,iblk) ) then
              remote_block = iblk
              remote_pe    = mype
              lfound = .true.
#ifdef DEBUG
             write(*,*) 'mpi_amr_restrict_fulltree: pe ',mype, & 
     &         ' child ', & 
     &          child(:,ich,lb),' located in buffer slot ', & 
     &          iblk
#endif /* DEBUG */
            endif
            iblk = iblk+1
          enddo 
        endif

                    cnodetype = nodetype(remote_block)
                    cempty = empty(remote_block)


!--------compute the offset in the parent block appropriate for this child
         ioff = mod(jchild-1,2)*nxb/2
         joff = mod((jchild-1)/2,2)*nyb/2
         koff = mod((jchild-1)/4,2)*nzb/2



! Get the child blocks data and fill its guardcells, putting the result
! into the current working block
           if(iopt.eq.1) then
                nlayers = nguard
           elseif(iopt.ge.2) then
                nlayers = nguard_work
           endif


           if((.not.lnc).and.(.not.lec)) then
            if( lcc.or.lfc.or.(iopt.ge.2)) then
!
! Put child blocks data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding. Note, for even grid sizes the guardcells
! do not need to be filled.
             idest = 1
             call amr_perm_to_1blk(lcc,lfc,lec,lnc, & 
     &                             remote_block,remote_pe, & 
     &                             iopt,idest)


            endif
           elseif(lnc.or.lec)  then
!
! Put child blocks data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding, and fill its guardcells with valid
! data.

             ldiag = diagonals
             icoord = 0                          ! fill in all coord directions
!             l_srl_only = .true.
             l_srl_only = .false.
             call amr_1blk_guardcell(mype,iopt,nlayers, & 
     &                               remote_block,remote_pe, & 
     &                               lcc,lfc,lec,lnc, & 
     &                               l_srl_only,icoord,ldiag)
           endif
           if (iopt.eq.1) call flash_convert_cc_hook(unk1(:,:,:,:,1), nvar, &
                il_bnd1,iu_bnd1, jl_bnd1,ju_bnd1, kl_bnd1,ku_bnd1, &
                why=gr_callReason_RESTRICT)


!-----------------------
         if (curvilinear) then
! compute geometry variables for the child block (remote_block,remote_pe)
#ifdef DEBUG
         write(*,*) 'pe ',mype,' blk ',lb,' searching for child ', & 
     &     jchild ,' at address ',remote_block,remote_pe,' for geom '
#endif /* DEBUG */
         call amr_block_geometry(remote_block,remote_pe)

         if (curvilinear_conserve) then

                             where(interp_mask_unk_res(:) .NE. 40) &
                                   interp_mask_unk_res(:) = 1
            interp_mask_work_res(:) = 1
            interp_mask_facex_res(:) = 1
            interp_mask_facey_res(:) = 1
            interp_mask_facez_res(:) = 1
            interp_mask_ec_res(:) = 1
            interp_mask_nc_res(:) = 1

         if(iopt.eq.1) then

           i1 = 1 + nguard
           i2 = nxb + nguard
           j1 = 1 + k2d*nguard
           j2 = nyb + k2d*nguard
           k1 = 1 + k3d*nguard
           k2 = nzb + k3d*nguard

! Compute volume weighted cell center data for conservative restriction
           do ivar = 1,nvar
                                   If (int_gcell_on_cc(ivar) &
                                        .AND.curvilinear_cons_mask_unk(ivar))          &
                    unk1(ivar,i1:i2,j1:j2,k1:k2,1) = & 
     &                  unk1(ivar,i1:i2,j1:j2,k1:k2,1) & 
     &                  *cell_vol(i1:i2,j1:j2,k1:k2)

           enddo
! Compute area weighted cell face-center data for conservative restriction
           do ivar = 1,nfacevar
             if(int_gcell_on_fc(1,ivar)) & 
     &           facevarx1(ivar,i1:i2+1,j1:j2,k1:k2,1) = & 
     &              facevarx1(ivar,i1:i2+1,j1:j2,k1:k2,1) & 
     &               *cell_area1(i1:i2+1,j1:j2,k1:k2)
             if(int_gcell_on_fc(2,ivar)) & 
     &           facevary1(ivar,i1:i2,j1:j2+k2d,k1:k2,1) = & 
     &              facevary1(ivar,i1:i2,j1:j2+k2d,k1:k2,1) & 
     &               *cell_area2(i1:i2,j1:j2+k2d,k1:k2)
             if(int_gcell_on_fc(2,ivar)) & 
     &           facevarz1(ivar,i1:i2,j1:j2,k1:k2+k3d,1) = & 
     &              facevarz1(ivar,i1:i2,j1:j2,k1:k2+k3d,1) & 
     &               *cell_area3(i1:i2,j1:j2,k1:k2+k3d)
           enddo
! Compute distance weighted cell edge-center data for conservative restriction
           do ivar = 1,nvaredge
             if(int_gcell_on_ec(1,ivar)) & 
     &             unk_e_x1(ivar,i1:i2,j1:j2+k2d,k1:k2+k3d,1) = & 
     &              unk_e_x1(ivar,i1:i2,j1:j2+k2d,k1:k2+k3d,1) & 
     &               *cell_leng1(i1:i2,j1:j2+k2d,k1:k2+k3d)
             if(int_gcell_on_ec(2,ivar)) & 
     &             unk_e_y1(ivar,i1:i2+1,j1:j2,k1:k2+k3d,1) = & 
     &              unk_e_y1(ivar,i1:i2+1,j1:j2,k1:k2+k3d,1) & 
     &               *cell_leng2(i1:i2+1,j1:j2,k1:k2+k3d)
             if(int_gcell_on_ec(3,ivar)) & 
     &             unk_e_z1(ivar,i1:i2+1,j1:j2+k2d,k1:k2,1) = & 
     &              unk_e_z1(ivar,i1:i2+1,j1:j2+k2d,k1:k2,1) & 
     &               *cell_leng3(i1:i2+1,j1:j2+k2d,k1:k2)
           enddo

         else

! Compute volume weighted cell center data for conservative restriction
!  of work1.
           ndel = nguard_work - nguard
           do k=kl_bnd1+nguard*k3d,ku_bnd1+nguard*k3d
           do j=jl_bnd1+nguard*k2d,ju_bnd1+nguard*k2d
           do i=il_bnd1+nguard    ,iu_bnd1+nguard
             work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1) =  & 
     &            work1(i+ndel,j+ndel*k2d,k+ndel*k3d,1)  & 
     &                                        *cell_vol(i,j,k)
           enddo
           enddo
           enddo

         endif
         endif

! now reset geometry factors to appropriate values for the current block lb
         call amr_block_geometry(lb,mype)

         endif
!-----------------------



         if(iopt.eq.1) then

! Compute restricted cell-centered data from the data in the buffer
           if(lcc) then
             call amr_restrict_unk_fun(unk1(:,:,:,:,1),temp,ioff,joff,koff)

           do k=1+nguard*k3d,nzb+nguard*k3d,2
             kk = (k-nguard*k3d)/2+1+nguard*k3d
             do j=1+nguard*k2d,nyb+nguard*k2d,2
               jj = (j-nguard*k2d)/2+1+nguard*k2d
               do i=1+nguard,nxb+nguard,2
                 ii = (i-nguard)/2+1+nguard
                 do ivar=1,nvar
                   if(int_gcell_on_cc(ivar)) then
                   send(ivar,ii,jj,kk) = temp(ivar,i,j,k)
                   endif
                 enddo
               enddo
             enddo
           enddo

! update the parent block
           do k=1,nzb+(-nzb/2)*k3d
             do j=1,nyb+(-nyb/2)*k2d
               do i=1,nxb-nxb/2
                 do ivar=1,nvar
                                         If (int_gcell_on_cc(ivar)) Then
                   if (curvilinear .AND. curvilinear_cons_mask_unk(ivar)) Then
                   unk(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                      k+nguard0*k3d+koff,lb) = & 
     &                send(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_vol(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
                                            Else
                   unk(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                      k+nguard0*k3d+koff,lb) = & 
     &                send(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                                            End If
                                         End If
                 enddo
               enddo
             enddo
           enddo

                          End If  ! End If (lcc)



! Compute restricted cell corner data from the data in the buffer
           if(lnc) then


             call amr_restrict_nc_fun( unk_n1(:,:,:,:,1), & 
     &                                 tempn )



           do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
             kk = (k-nguard*k3d)/2+1+nguard*k3d
             do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
               jj = (j-nguard*k2d)/2+1+nguard*k2d
               do i=1+nguard,nxb+nguard+1,2
                 ii = (i-nguard)/2+1+nguard
                 do ivar=1,nvarcorn
                   sendn(ivar,ii,jj,kk) = tempn(ivar,i,j,k)
                 enddo
               enddo
             enddo
           enddo


! update the parent block
           do k=1,nzb+(-nzb/2+1)*k3d
             do j=1,nyb+(-nyb/2+1)*k2d
               do i=1,nxb-nxb/2+1
                 do ivar=1,nvarcorn
                   unk_n(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                           k+nguard0*k3d+koff,lb)= & 
     &              sendn(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
                 enddo
               enddo
             enddo
           enddo



           endif                 ! end of lnc if test
           endif                ! end of iopt iftest


! Compute restricted cell-face-centered data from the data in the buffer
           if(lfc) then
! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                                   kl_bnd1:ku_bnd1) & 
     &   =  facevarx1(:nfacevar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                  kl_bnd1:ku_bnd1,1)

       call amr_restrict_fc_fun(sendf,tempf,1)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nfacevar

               if (curvilinear) then
               if (curvilinear_conserve) then
               facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_area1(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
               else
               facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               facevarx(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo


! y face next
       if(ndim.ge.2) then

! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                        kl_bnd1:ku_bnd1) & 
     &   = facevary1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                 kl_bnd1:ku_bnd1,1)
       call amr_restrict_fc_fun(sendf,tempf,2)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nfacevar
               if (curvilinear) then
               if (curvilinear_conserve) then
               facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / (cell_area2(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff) + tiny)
               else
               facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               facevary(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo
       endif                      ! end of ndim>=2 test

! z face last
       if(ndim.eq.3) then

! Compute restricted data from the data in the buffer
       sendf(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                                 kl_bnd1:ku_bnd1+k3d) & 
     &    = facevarz1(:nfacevar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, & 
     &                  kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_fc_fun(sendf,tempf,3)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nfacevar
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nfacevar
               if (curvilinear) then
               if (curvilinear_conserve) then
               facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_area3(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
               else
               facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               facevarz(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo

       endif                      ! end of ndim=3 test

           endif                 ! end of lfc if test

           if (ndim > 1) then
! Compute restricted cell-edge-centered data from the data in the buffer
           if(lec) then

! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                         kl_bnd1:ku_bnd1+k3d) & 
     &   =  unk_e_x1(:nvaredge,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, & 
     &                  kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_ec_fun(sendf,tempf,1)


       sendf = 0.
       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo


! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2
             do ivar=1,nvaredge
               if (curvilinear) then
               if (curvilinear_conserve) then
               unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_leng1(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
               else
               unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               unk_e_x(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo


! y edge next
! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d) & 
     &   = unk_e_y1(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1, & 
     &                 kl_bnd1:ku_bnd1+k3d,1)
       call amr_restrict_ec_fun(sendf,tempf,2)

       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2+1)*k3d
         do j=1,nyb+(-nyb/2)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nvaredge
               if (curvilinear) then
               if (curvilinear_conserve) then
               unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_leng2(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
               else
               unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               unk_e_y(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo

       if (ndim == 3) then
! z edge last
! Compute restricted data from the data in the buffer
       sendf(:nvaredge,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d, & 
     &                                   kl_bnd1:ku_bnd1) & 
     &    = unk_e_z1(:nvaredge,il_bnd1:iu_bnd1+1, & 
     &       jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1,1)
       call amr_restrict_ec_fun(sendf,tempf,3)

       do k=1+nguard*k3d,nzb+nguard*k3d,2
         kk = (k-nguard*k3d)/2+1+nguard*k3d
         do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
           jj = (j-nguard*k2d)/2+1+nguard*k2d
           do i=1+nguard,nxb+nguard+1,2
             ii = (i-nguard)/2+1+nguard
             do ivar=1,nvaredge
               sendf(ivar,ii,jj,kk) = tempf(ivar,i,j,k)
             enddo
           enddo
         enddo
       enddo

! update the parent block
       do k=1,nzb+(-nzb/2)*k3d
         do j=1,nyb+(-nyb/2+1)*k2d
           do i=1,nxb-nxb/2+1
             do ivar=1,nvaredge
               if (curvilinear) then
               if (curvilinear_conserve) then
               unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d) & 
     &           / cell_leng3(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
               else
               unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
               else
               unk_e_z(ivar,i+nguard0+ioff,j+nguard0*k2d+joff, & 
     &                       k+nguard0*k3d+koff,lb)= & 
     &          sendf(ivar,i+nguard,j+nguard*k2d,k+nguard*k3d)
               endif
             enddo
           enddo
         enddo
       enddo
       end if ! if (ndim == 3)

       endif                 ! end of lec if test
       end if ! if (ndim > 1)


         if(iopt.ge.2) then
           iopt0 = iopt-1


! Compute restricted cell-centered workspace data from the data in the buffer


           if(mod(nxb,2).eq.0) then
             call amr_restrict_work_fun(work1(:,:,:,1),tempw1,iopt)
           else
             call amr_restrict_work_fun_recip(work1(:,:,:,1),tempw1)
           endif


           do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
             kk = (k-nguard_work*k3d)/2+1+nguard_work*k3d
             do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
               jj = (j-nguard_work*k2d)/2+1+nguard_work*k2d
               do i=1+nguard_work,nxb+nguard_work,2
                 ii = (i-nguard_work)/2+1+nguard_work
                 recvw1(ii,jj,kk,1) = tempw1(i,j,k)
               enddo
             enddo
           enddo


! update the parent block
           do k=1,nzb+(-nzb/2)*k3d
             do j=1,nyb+(-nyb/2)*k2d
               do i=1,nxb-nxb/2
                 if (curvilinear) then
                 if (curvilinear_conserve) then
                 work(i+nguard_work0+ioff,j+nguard_work0*k2d+joff, & 
     &                k+nguard_work0*k3d+koff,lb,iopt0) =   & 
     &             recvw1(i+nguard_work,j+nguard_work*k2d, & 
     &                    k+nguard_work*k3d,1) & 
     &           / cell_vol(i+nguard+ioff,j+nguard*k2d+joff, & 
     &                      k+nguard*k3d+koff)
                 else
                 work(i+nguard_work0+ioff,j+nguard_work0*k2d+joff, & 
     &                k+nguard_work0*k3d+koff,lb,iopt0) =   & 
     &             recvw1(i+nguard_work,j+nguard_work*k2d, & 
     &                    k+nguard_work*k3d,1)
                 endif
                 else
                 work(i+nguard_work0+ioff,j+nguard_work0*k2d+joff, & 
     &                k+nguard_work0*k3d+koff,lb,iopt0) =   & 
     &             recvw1(i+nguard_work,j+nguard_work*k2d, & 
     &                    k+nguard_work*k3d,1)
                 endif
               enddo
             enddo
           enddo

       endif                      ! end of iopt if test


      enddo                       ! end of loop over children

      if (iopt.eq.1) call flash_unconvert_cc_hook(unk(:,:,:,:,lb), nvar, &
                il_bnd,iu_bnd, jl_bnd,ju_bnd, kl_bnd,ku_bnd, &
                where=gr_cells_INTERIOR, why=gr_callReason_RESTRICT, &
                nlayers_in_data=nguard0)

! If using odd sized grid blocks then parent copies any face bounding
! a leaf block
      if(iopt.eq.1) then
      if(lnc) then


! cycle through parents neighbors
        do jface = 1,nfaces

! get this neighbors nodetype
        remote_pe     = neigh(2,jface,lb)
        remote_block  = neigh(1,jface,lb)

        if(remote_block.gt.0) then

          cnodetype = -1
          if(remote_pe.eq.mype) cnodetype = nodetype(remote_block)

! if (remote_block,remote_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
        remote_pe0     = remote_pe
        remote_block0  = remote_block

        if(remote_pe0.ne.mype) then
        lfound = .false.
        iblk = strt_buffer
!          do iblk = strt_buffer,last_buffer
          do while(.not.lfound.and.iblk.le.last_buffer)
#ifdef DEBUG
             write(*,*) 'mpi_amr_1blk_restrict: pe ',mype,' blk ',lb, & 
     &     ' searching buffer for neigh ', & 
     &     neigh(:,jface,lb),' current buffer entry ', & 
     &          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
            if(remote_block0.eq.laddress(1,iblk).and. & 
     &           remote_pe0 .eq.laddress(2,iblk) ) then
              remote_block0 = iblk
              remote_pe0    = mype
              lfound = .true.
#ifdef DEBUG
             write(*,*) 'mpi_amr_1blk_restrict: pe ',mype,' neigh ', & 
     &          neigh(:,jface,lb),' located in buffer slot ', & 
     &          iblk,' unk of neigh ',unk(1,:,:,:,iblk)
#endif /* DEBUG */
            endif
            iblk = iblk+1
          enddo

        if(lfound) then
         cnodetype = nodetype(remote_block0)
        endif

        endif

        if(cnodetype.eq.1) then

        if(iopt.eq.1) then
             ng0 = nguard0
        elseif(iopt.ge.2) then
             ng0 = nguard_work0
        endif
        ia = 1+ng0
        ib = nxb+ng0+1
        ja = 1+ng0*k2d
        jb = nyb+(ng0+1)*k2d
        ka = 1+ng0*k3d
        kb = nzb+(ng0+1)*k3d
        isa = 1+ng0
        isb = nxb+ng0+1
        jsa = 1+ng0*k2d
        jsb = nyb+(ng0+1)*k2d
        ksa = 1+ng0*k3d
        ksb = nzb+(ng0+1)*k3d

         if(jface.eq.1) then
           ib  = ia
           isa = isb
         elseif(jface.eq.2) then
           ia  = ib
           isb = isa
         elseif(jface.eq.3) then
           jb  = ja
           jsa = jsb
         elseif(jface.eq.4) then
           ja  = jb
           jsb = jsa
         elseif(jface.eq.5) then
           kb  = ka
           ksa = ksb
         elseif(jface.eq.6) then
           ka  = kb
           ksb = ksa
         endif


! Copy neighbor face into this block
         if(iopt.eq.1) then
           if(remote_block.le.lnblocks & 
     &                .and.remote_pe.eq.mype) then
             unk_n(:,ia:ib,ja:jb,ka:kb,lb) =  & 
     &                      unk_n(:,isa:isb,jsa:jsb,ksa:ksb,remote_block)

           else


! The next section is largely borrowed from amr_1blk_guardcell_srl.
! It sets the index ranges for copying the common face unk_n data.
! Only 1 layer is required here
!-----------------------------
             iblock = 1

! Range - source indeces are initially computed as though there
! are no permanent guardcells.
        ilays = nxb
        jlays = nyb*k2d
        klays = nzb*k3d
! Starting indeces on destination working block
        id = 1 + nguard
        jd = 1 + nguard*k2d
        kd = 1 + nguard*k3d
! Starting indeces on source block
        is = 1 + nguard0
        js = 1 + nguard0*k2d
        ks = 1 + nguard0*k3d

        ip1 = 0
        jp1 = 0
        kp1 = 0
        ip3 = 0
        jp3 = 0
        kp3 = 0

        if(jface.eq.1) then
          is = 1 + nxb + ng0
          ilays = 0
        elseif(jface.eq.2) then
          id = 1 + nxb + nguard
          is = 1 + ng0
          ilays = 0
        elseif(jface.eq.3) then
          js = 1 + (nyb+ng0)*k2d
          jlays = 0
        elseif(jface.eq.4) then
          jd = 1 + (nyb + nguard)*k2d
          js = 1 + ng0*k2d
          jlays = 0
        elseif(jface.eq.5) then
          ks = 1 + (nzb+ng0)*k3d
          klays = 0
        elseif(jface.eq.6) then
          kd = 1 + (nzb + nguard)*k3d
          ks = 1 + ng0*k3d
          klays = 0
        endif

!-----------------------------
! now we use this index info to extract remote data directly from the
! receive buffer. amr_1blk_nc_cp_remote write to unk_n1 so after it exits
! we need to copy the result to unk_n.
             call amr_1blk_nc_cp_remote( & 
     &                 mype,remote_pe,remote_block,iblock, & 
     &                 id,jd,kd,is,js,ks, & 
     &                 ilays,jlays,klays, & 
     &                 ip1,jp1,kp1,ip3,jp3,kp3,0)

             ng1 = nguard*(1-npgs)
             unk_n(:,id-ng1:id+ilays-ng1, & 
     &               jd-ng1*k2d:jd+(jlays-ng1)*k2d, & 
     &               kd-ng1*k3d:kd+(klays-ng1)*k3d,lb) = & 
     &        unk_n1(:,id:id+ilays, & 
     &                 jd:jd+jlays*k2d, & 
     &                 kd:kd+klays*k3d,iblock)


         endif

           endif

        endif                     ! end of cnodetype if test

      endif                       ! end of remote_block if test
      enddo                       ! end of jface loop
      endif                       ! end of lnc if test
       endif                      ! end of iopt if test


      endif
      enddo                       ! end of loop over blocks

      endif



! Make sure that the global copy of the newly restricted data is
! up to date.
      call amr_1blk_copy_soln(level)

      enddo                      ! end of loop over refinement levels
      endif

  lrestrict_in_progress = .false.

  if (allocated(curvilinear_cons_mask_unk)) then
     deallocate(curvilinear_cons_mask_unk)
  end if

      deallocate(tempf)
      deallocate(sendf)

#ifdef DEBUG_FLOW_TRACE
      write(*,*) 'exiting mpi_amr_restrict_fulltree: pe ',mype, & 
     &           ' iopt ',iopt
#endif /* DEBUG_FLOW_TRACE */

      return
      end subroutine mpi_amr_restrict_fulltree

