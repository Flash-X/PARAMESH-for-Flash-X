!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
! Copyright (C) 2012, 2013, 2016 The University of Chicago
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!
!     Petros Tzeferacos and Dongwook Lee, July - August 2012
!        Customizations for cylindrical geometry.
!        New prolongation method (modified Li&Li, 2004) and
!        corrected amr_minmod (was not minmod at all)
!        for balsara prolongation.
!
!     Klaus Weide, April 2013
!        Indexing fixes for cases in which i,j,k are <= nguard
!
!     Klaus Weide, May 2016
!        moved amr_minmod inside subroutine balpro


       subroutine prol_fc_dbz_init(n,i_divf_fc_vars)
        use prolong_arrays, only :  & 
     &   prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n
        use physicaldata, only : interp_mask_facex,  & 
     &   interp_mask_facey,interp_mask_facez

        implicit none

        integer, intent(in) :: n, i_divf_fc_vars(:,:)
        integer i,iface
        prol_fc_dbz_n = n ! n should equal nbndvar (or nfacevar?)
        If (.Not.allocated(prol_fc_dbz_ivar))                          &
            allocate(prol_fc_dbz_ivar(3,prol_fc_dbz_n))
        do i = 1,prol_fc_dbz_n
        do iface=1,3
           prol_fc_dbz_ivar(iface,i) = i_divf_fc_vars(iface,i)
        end do
        interp_mask_facex(i_divf_fc_vars(1,i)) = -200 ! Corresponds to fc_dbz.
        interp_mask_facey(i_divf_fc_vars(2,i)) = -200 ! Corresponds to fc_dbz.
        interp_mask_facez(i_divf_fc_vars(3,i)) = -200 ! Corresponds to fc_dbz.
        end do
        prol_fc_dbz = .true.
       end subroutine prol_fc_dbz_init

       function prol_fc_dbz_varp(ivar, iface) result(ldbz)
        use prolong_arrays, only :  & 
     &    prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n

        implicit none

        integer, intent(in) :: ivar, iface ! iface from {1 2 3} id'd {x y z}
        logical :: ldbz ! true if this var is prolonged with dbz routines
        integer :: i
        ldbz = .false.
        do i = 1, prol_fc_dbz_n
         ldbz = (prol_fc_dbz_ivar(iface,i) == ivar).or.ldbz
        end do
        return
       end function prol_fc_dbz_varp

       
     !**************************************************************
     ! Main Driver for the Balsara prolongation routine
     ! Sets indices and calls balpro which does the interp.
     !
       subroutine amr_1blk_fc_prol_dbz(  & 
     &        recvfx, recvfy, recvfz,        & 
     &        nfacevar_in, iv1, iv2, iv3,  & 
     &        ia,ib,ja,jb,ka,kb,    & 
     &        idest,ioff,joff,koff,          & 
     &        mype,lb,parent_pe,parent_blk   & 
     & )
      !**************************************************************
        use paramesh_dimensions 
        use physicaldata
        use tree
        ! use prolong_arrays 

        implicit none

        real, intent(inout), dimension(:,:,:,:) :: recvfx,recvfy,recvfz
        integer, intent(in) :: nfacevar_in
        integer, intent(in) :: iv1,iv2,iv3
        integer, intent(in) :: ia,ib,ja,jb,ka,kb
        integer, intent(in) :: idest, ioff, joff, koff
        integer, intent(in) :: mype, lb, parent_pe, parent_blk

        integer :: nv,n1,n2,n3
        integer :: icl, icu, jcl, jcu, kcl, kcu
        integer :: i, j, k
        real :: x0, y0, z0, x1, y1, z1, x2, y2, z2
        real :: dx, dy, dz, ddx, ddy, ddz
        real :: fx_0, fy_0, fz_0, fx_1, fy_1, fz_1
        real :: tx2, ty2, tz2
        real ::  scrh1,scrh2,scrh3
        logical :: left, right
        real :: Br_c(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, kl_bnd1:ku_bnd1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mlrdbg    
        integer :: ii,jj,kk
        integer :: i1, j1, k1
        logical :: mlrdbg,mlrdbg1,mlrdbg2,mlrdbg3
        logical :: mlrdbg100

        integer :: plb,plp
        real  :: tx1,ty1,tz1,tdx,tdy,tdz
        real  :: x_coll,xL,xR,r0,r1,divb, r3
        real  :: cell_face_coord1_parent(il_bnd1:iu_bnd1+1),&
                 cell_face_coord2_parent(jl_bnd1:ju_bnd1+1),&
                 cell_face_coord3_parent(kl_bnd1:ku_bnd1+1)

        real  :: cell_face_coord1_child (il_bnd1:iu_bnd1+1),&
                 cell_face_coord2_child (jl_bnd1:ju_bnd1+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        mlrdbg =.false. ! .true.
        mlrdbg1=.false. ! .true.
        mlrdbg2=.false. ! .true.
        mlrdbg3=.false. ! .true.
        mlrdbg100=.false. ! .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        nv= nfacevar
        n1=iu_bnd1-il_bnd1+1
        n2=ju_bnd1-jl_bnd1+k2d
        n3=ku_bnd1-kl_bnd1+k3d

        plb = parent_blk
        plp = parent_pe

        if(nfacevar.ne.nfacevar_in)then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'a1fpd:nfacevar inconsistent!'
          call amr_abort
        end if

        !------------------------------------

        ! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb

        ! DL added this call for cylindrical
        if (cylindrical_pm) then
           Call amr_block_geometry(plb,plp)
           cell_face_coord1_parent = cell_face_coord1
           cell_face_coord2_parent = cell_face_coord2
           if (ndim == 3) cell_face_coord3_parent = cell_face_coord3
        
           Call amr_block_geometry(lb,mype)
           cell_face_coord1_child = cell_face_coord1
           cell_face_coord2_child = cell_face_coord2
        endif

        
        
        
        
        if (cylindrical_pm .eqv. .FALSE.) then
        
    ! Interpolation loop for all geometries but cylindrical.
    ! The latter is treated seperately after this.
    !
    ! Note that the range of indeces used in the facevar plane differs
    ! depending on the value of iface_off. This assumes that the face values
    ! corresponding to index 0 (ie nguard faces to the left of the block
    ! boundary) are never needed, when iface_off=-1. 

    ! Iterate up to upper_bound-1 and flag prolong operator to calculate values
    ! for right hands of intervals.
    !

        ! NOTE: i,j,k ARE THE INDICES FOR A NEWLY CREATED CHILD BLOCK
        kloop: do k=kcl,kcu+k3d
         jloop: do j=jcl,jcu+k2d
          iloop: do i=icl,icu+1
             ! 
             ! Endpoints for the box.  Points to be prolonged to are, e.g. [x0,(y0+y1)/2,(z0+z1)/2]
             ! which are on the faces of the box. Must be reconsidered for non-Cartesian grids.
             !
             ! Choose convenient (local) coordinates.
             !

             ! NOTE: i1,j1,k1 ARE INDICES FOR THE PARENT BLOCK
             !       i, j, k  run from 5 to 13 (9 of them)
             !       ii,jj,kk run from 10 to 18 (9 of them)
             !       ioff is either 0 or 4 (amr_prolong_face_fun_init.F90)
             !       using ioff, i1 runs from either 5 to 9 or 9 to 13
             ii = i + nguard + 1
             i1 = ii/2 + ioff
             jj = j + nguard + 1
             j1 = jj/2 + joff
             if (ndim == 3) then
                kk = k + nguard + 1
                k1 = kk/2 + koff
             else
                kk = k
                k1 = k
             end if

             dx = 2.0*bsize(1,lb)/real(nxb) !Note: bsize is a newly created child block (see gr_createBlock.F90)
             dy = 2.0*bsize(2,lb)/real(nyb)
             dz = 2.0*bsize(3,lb)/real(nzb)

             x0 = -0.5*2.0*bsize(1,lb)
             y0 = -0.5*2.0*bsize(2,lb)
             z0 = -0.5*2.0*bsize(3,lb)

             ! Left hand side of intervals, *off included for pedantism
             x1 = x0 + dx*real(i1-1-nguard)
             y1 = y0 + dy*real(j1-1-nguard)
             z1 = z0 + dz*real(k1-1-nguard)

             ! Right hand side of intervals, *off included for pedantism
             x2 = x1 + dx
             y2 = y1 + dy
             z2 = z1 + dz

             ! Distance between coarse points
             ddx = dx
             ddy = dy
             if (ndim == 3) then
                ddz = dz
             else
                ddz = 1.
             end if

             tdx = bsize(1,lb)/real(nxb)
             if (ndim >= 2) then
                tdy = bsize(2,lb)/real(nyb)
             else
                tdy = 0.
             end if
             if (ndim == 3) then
                tdz = bsize(3,lb)/real(nzb)
             else
                tdz = 0.
             end if


             !! dx is a grid delta for parent block; tdx is a grid delta for children blocks
             ! Require blocks be cut in half.
             tx1 = -tdx*mod(real(i),2.0)
             ty1 = -tdy*mod(real(j),2.0)
             tz1 = -tdz*mod(real(k),2.0)

             tx2 = tx1 + tdx
             ty2 = ty1 + tdy
             tz2 = tz1 + tdz

             r0=1.
             r1=1.

             left=.true. ! always calculate left sides of box
!!!             right=((i.eq.icu).or.(j.eq.jcu).or.(k.eq.kcu)) ! on the last box, calculate the other side
             ! compute interpolated values at points referred to above.


             call balpro( r0,r1,& 
     &                   iv1,iv2,iv3,           & ! Which variables to use in each array.
     &                   i1,j1,k1,              & ! Coords of the cell
     &                   nv, n1, n2, n3,        & ! 
     &                   tx1,ty1,tz1,           & ! 
     &                   tx2,ty2,tz2,           & ! 
     &                   recvfx, recvfy, recvfz,& ! The field arrays.
     &                   fx_0, fy_0, fz_0,      & ! The output field at p0 points (0-faces).
     &                   fx_1, fy_1, fz_1,      & ! The output field at p1 points (1-faces).
     &                   ddx, ddy, ddz,         & ! Size of cell in each of these directions.
     &                   left, .true.,          & ! Flags to return f*_0 and f*_1 respectively.
     &                   k2d, k3d, mype )



             if (j <= jcu .And. k <= kcu) then
                facevarx1(iv1,i,j,k,idest) = fx_0
             end if
             if (i <= icu .And. k <= kcu) then
              facevary1(iv2,i,j,k,idest) = fy_0
             end if
             if (i <= icu .And. j <= jcu) then
              facevarz1(iv3,i,j,k,idest) = fz_0
             end if

!             if(left)then
!              facevarx1(iv1,i,j,k,idest) = fx_0
!              facevary1(iv2,i,j,k,idest) = fy_0
!              facevarz1(iv3,i,j,k,idest) = fz_0
!             end if
       
!             if(i == icu)then
!              facevarx1(iv1,i+1,j,k,idest)   = fx_1
!             end if
!             if (j == jcu) then
!              facevary1(iv2,i,j+k2d,k,idest) = fy_1 
!             end if
!             if (k == kcu) then
!              facevarz1(iv3,i,j,k+k3d,idest) = fz_1
!             end if
       
            enddo iloop
           enddo jloop
          enddo kloop
    endif           
          
          
          
          
    if (cylindrical_pm) then
          
    ! Interpolation loop for cylindrical geometry.
    !
        ! NOTE: i,j,k ARE THE INDICES FOR A NEWLY CREATED CHILD BLOCK
        kloopC1: do k=kcl,kcu+k3d
         jloopC1: do j=jcl,jcu+k2d
          iloopC1: do i=icl,icu+1
             !
             ! NOTE: i1,j1,k1 ARE INDICES FOR THE PARENT BLOCK
             !       i, j, k  run from 5 to 13 (9 of them)
             !       ii,jj,kk run from 10 to 18 (9 of them)
             !       ioff is either 0 or 4 (amr_prolong_face_fun_init.F90)
             !       using ioff, i1 runs from either 5 to 9 or 9 to 13
             ii = i + nguard + 1
             i1 = ii/2 + ioff
             jj = j + nguard + 1
             j1 = jj/2 + joff
             if (ndim == 3) then
                kk = k + nguard + 1
                k1 = kk/2 + koff
             else
                kk = k
                k1 = k
             end if

             r0 = cell_face_coord1_parent(i1)
             r1 = cell_face_coord1_parent(i1+1)
             dx = 0.5*(r1**2-r0**2)

             if (ndim >= 2) dy = cell_face_coord2_parent(j1+k2d)-cell_face_coord2_parent(j1)
             if (ndim == 3) dz = cell_face_coord3_parent(k1+ 1 )-cell_face_coord3_parent(k1)
 
             ddx = dx
             ddy = dy

             if (ndim == 3) then
               ddz = dz
             else
               ddz = 1.
             end if
             
             xL     = 0.5*r0**2
             xR     = 0.5*r1**2

             x_coll = 0.5*(r0+r1)
             x_coll = x_coll*x_coll*0.5

             tx1    = (xL-x_coll)*mod(real(i),2.0)
             tx2    = (xR-x_coll)*(1.0-mod(real(i),2.0))

             if (ndim >= 2) then
                tdy = bsize(2,lb)/real(nyb)

             else
                tdy = 0.
             end if
             if (ndim == 3) then
                tdz = bsize(3,lb)/real(nzb)
             else
                tdz = 0.
             end if
             ty1 = -tdy*mod(real(j),2.0)
 
             tz1 = -tdz*mod(real(k),2.0)

             ty2 = ty1 + tdy

             tz2 = tz1 + tdz


             left=.true. ! always calculate left sides of box
!!!             right=((i.eq.icu).or.(j.eq.jcu).or.(k.eq.kcu)) ! on the last box, calculate the other side
             ! compute interpolated values at points referred to above.


             call balpro( r0,r1,& 
     &                   iv1,iv2,iv3,           & ! Which variables to use in each array.
     &                   i1,j1,k1,              & ! Coords of the cell
     &                   nv, n1, n2, n3,        & ! 
     &                   tx1,ty1,tz1,           & ! 
     &                   tx2,ty2,tz2,           & ! 
     &                   recvfx, recvfy, recvfz,& ! The field arrays.
     &                   fx_0, fy_0, fz_0,      & ! The output field at p0 points (0-faces).
     &                   fx_1, fy_1, fz_1,      & ! The output field at p1 points (1-faces).
     &                   ddx, ddy, ddz,         & ! Size of cell in each of these directions.
     &                   left, .true.,          & ! Flags to return f*_0 and f*_1 respectively.
     &                   k2d, k3d, mype )


             r0 = cell_face_coord1_child( i )
             r1 = cell_face_coord1_child(i+1)

             if (j <= jcu .And. k <= kcu) then
!                   fx_0=fx_0/r0
!                   facevarx1(iv1,i,j,k,idest) = fx_0
             end if
             if (i <= icu .And. k <= kcu) then
              facevary1(iv2,i,j,k,idest) = fy_0
!              facevary1(iv2,i,j+k2d,k,idest) = fy_1 
             end if
             if (i <= icu .And. j <= jcu) then
              facevarz1(iv3,i,j,k,idest) = fz_0
             end if
       
            enddo iloopC1
           enddo jloopC1
          enddo kloopC1

             kloopC2: do k=kcl,kcu
                jloopC2: do j=jcl,jcu
                   iloopC2: do i=icl,icu

                       ii = i + nguard + 1
                       i1 = ii/2 + ioff
                       jj = j + nguard + 1
                       j1 = jj/2 + joff
                       if (ndim == 3) then
                         kk = k + nguard + 1
                         k1 = kk/2 + koff
                       else
                         kk = k
                         k1 = k
                       end if

                      r0=cell_face_coord1_child(i)
                      r1=cell_face_coord1_child(i+1)
                      dx = (r1**2 - r0**2)*0.5
                      dy = (cell_face_coord2_child(j+1)-cell_face_coord2_child(j))


!close divB
if (mod(i,2)==1 .and. mod(j,2)==1) then

Br_c(iv1,i1,j1,k) = (r0/r1)*recvfx(iv1,i1,j1,k) &
     - (facevary1(iv2,i,j+2,k,idest) - facevary1(iv2,i,j,k,idest))*dx/(2.*dy*r1)
endif


                   enddo iloopC2
                enddo jloopC2
             enddo kloopC2
          
             kloopC3: do k=kcl,kcu
                jloopC3: do j=jcl,jcu
                   iloopC3: do i=icl,icu

                       ii = i + nguard + 1
                       i1 = ii/2 + ioff
                       jj = j + nguard + 1
                       j1 = jj/2 + joff
                       if (ndim == 3) then
                         kk = k + nguard + 1
                         k1 = kk/2 + koff
                       else
                         kk = k
                         k1 = k
                       end if

                      r0=cell_face_coord1_child(i)
                      r1=cell_face_coord1_child(i+1)
                      dx = (r1**2 - r0**2)*0.5
                      dy = (cell_face_coord2_child(j+1)-cell_face_coord2_child(j))
          
!injections          
          
                       if (mod(i,2) == 1 .and. mod(j,2)==1) then
                         facevarx1(iv1,i,j,k,idest) = recvfx(iv1,i1,j1,k)
                         facevarx1(iv1,i,j+1,k,idest) = recvfx(iv1,i1,j1,k)
                         facevarx1(iv1,i+1,j,k,idest) = Br_c(iv1,i1,j1,k)
                         facevarx1(iv1,i+1,j+1,k,idest) = Br_c(iv1,i1,j1,k)
                         facevarx1(iv1,i+2,j,k,idest) = recvfx(iv1,i1+1,j1,k)
                         facevarx1(iv1,i+2,j+1,k,idest) = recvfx(iv1,i1+1,j1,k)
                       endif
                   enddo iloopC3
                enddo jloopC3
             enddo kloopC3
          
          
             kloopC4: do k=kcl,kcu+k3d
                jloopC4: do j=jcl,jcu+k2d
                   iloopC4: do i=icl,icu+1

                       ii = i + nguard + 1
                       i1 = ii/2 + ioff
                       jj = j + nguard + 1
                       j1 = jj/2 + joff
                       if (ndim == 3) then
                         kk = k + nguard + 1
                         k1 = kk/2 + koff
                       else
                         kk = k
                         k1 = k
                       end if
                   
                      r0=cell_face_coord1_child(i)
                      r1=cell_face_coord1_child(i+1)
                      dx = (r1**2 - r0**2)*0.5
                      dy = cell_face_coord2_child(j+1)-cell_face_coord2_child(j)
!close divb

if (mod(j,2)==1) then
        facevary1(iv2,i,j+1,k,idest) = facevary1(iv2,i,j,k,idest) &
        -(r1*facevarx1(iv1,i+1,j,k,idest) &
        - r0*facevarx1(iv1,i,j,k,idest))*dy/dx 
endif

                   enddo iloopC4
                enddo jloopC4
             enddo kloopC4
     endif


         end subroutine amr_1blk_fc_prol_dbz

         subroutine balpro(  r0,r1,& 
     & iv1,iv2,iv3,      & ! Which variables to use in each array.
     & i,j,k,            & ! Coords of the cell
     & nv, n1, n2, n3,   & ! 
     & x0,y0,z0,         & ! 
     & x1,y1,z1,         & ! 
     & fx, fy, fz,       & ! The field arrays.
     & fx_0, fy_0, fz_0, & 
     & fx_1, fy_1, fz_1, & 
     & ddx, ddy, ddz,    & ! Size of cell in each of these directions.
     & left, right,      & ! Side(s) of cell to calculate
     & k2d, k3d, mype )

           use paramesh_dimensions, only: il_bnd1,iu_bnd1,jl_bnd1,ju_bnd1,kl_bnd1,ku_bnd1,ndim
           use physicaldata

         implicit none


         real :: x0,y0,z0,x1,y1,z1, scrh1,scrh2,scrh3,xpos,ypos
         
         logical :: left, right
         integer, intent(in) :: k2d,k3d
         integer :: hp
         real    :: r0,r1
         integer :: iv1,iv2,iv3,i,j,k
         integer :: nv,n1,n2,n3
         real ::   fx(nv,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,    kl_bnd1:ku_bnd1    ),&
                   fy(nv,il_bnd1:iu_bnd1,  jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1    ),&
                   fz(nv,il_bnd1:iu_bnd1,  jl_bnd1:ju_bnd1,    kl_bnd1:ku_bnd1+k3d)
         real :: fx_0, fy_0, fz_0
         real :: fx_1, fy_1, fz_1

         real :: ddx, ddy, ddz
         real :: ddx2, ddy2, ddz2

         real :: a0,b0,c0
         real :: ax, ay, az, bx, by, bz, cx, cy, cz

         real :: axx, axy, axz, ayz
         real :: bxy, bxz, byy, byz
         real :: cxy, cxz, cyz, czz

         real :: axyz, bxyz, cxyz
         real :: axxy, axxz, byyx, byyz, cxzz, cyzz

         real :: fxp, fyp, fzp
         real :: fxm, fym, fzm

         real ::  & 
     & dy_fxp, dz_fxp, dx_fyp, dz_fyp, dx_fzp, dy_fzp,  & 
     & dy_fxm, dz_fxm, dx_fym, dz_fym, dx_fzm, dy_fzm 

         real :: dxy_fzp1, dxy_fzp2, dxy_fzp3

         real :: dxy_fzp, dxy_fzm, dyz_fxp,  & 
     &                             dyz_fxm, dxz_fyp, dxz_fym

         integer :: ihm, jhm, khm, ihp, jhp, khp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!mlrdbg    
         integer :: mype
         logical :: mlrdbg,baldbg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         mlrdbg=.false. ! .true.
         baldbg=.false. ! .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The p-variables can refer to the greater end of the current cell, and the m-variables
    ! refer to the current cell, or the p-variables can refer to the current cell while the
    ! m-variables refer to the previous cell.
         !!! legacy hm=0; hp=1

         ihm= 0;    jhm= 0;      khm= 0
         !! ihp = hp;    jhp = hp;      khp = hp
         ihp= 1;    jhp= k2d;    khp= k3d  ! previously, these were all set to hp=1.

         fxp = fx(iv1,i+ihp,j,    k    )*r1
         fyp = fy(iv2,i,    j+jhp,k    )
         fzp = fz(iv3,i,    j,    k+khp)

         fxm = fx(iv1,i+ihm,j,    k    )*r0
         fym = fy(iv2,i,    j+jhm,k    )
         fzm = fz(iv3,i,    j,    k+khm)


         dy_fxp = amr_minmod(r1*fx(iv1,i+ihp,j+1,k)-r1*fx(iv1,i+ihp,j,  k), & 
     &                       r1*fx(iv1,i+ihp,j,  k)-r1*fx(iv1,i+ihp,j-1,k)) ! dy(i+ihp,[j,j-1],k)
         if (ndim>2) then
         dz_fxp = amr_minmod(r1*fx(iv1,i+ihp,j,k+1)-r1*fx(iv1,i+ihp,j,k  ), & 
     &                       r1*fx(iv1,i+ihp,j,k  )-r1*fx(iv1,i+ihp,j,k-1)) ! dz(i+ihp,j,[k,k-1])
         else
         dz_fxp = 0.
         endif

         dx_fyp = amr_minmod(fy(iv2,i+1,j+jhp,k)-fy(iv2,i,  j+jhp,k), & 
     &                       fy(iv2,i,  j+jhp,k)-fy(iv2,i-1,j+jhp,k)) ! dx([i,i-1],j+jhp,k)
         if (ndim>2) then
         dz_fyp = amr_minmod(fy(iv2,i,j+jhp,k+1)-fy(iv2,i,j+jhp,k  ), & 
     &                       fy(iv2,i,j+jhp,k  )-fy(iv2,i,j+jhp,k-1)) ! dz(i,j+jhp,[k,k-1])
         else
         dz_fyp = 0.
         endif

         if (ndim>2) then
         dx_fzp = amr_minmod(fz(iv3,i+1,j,k+khp)-fz(iv3,i,  j,k+khp), & 
     &                       fz(iv3,i,j,  k+khp)-fz(iv3,i-1,j,k+khp)) ! dx([i,i-1],j,k+khp)
         dy_fzp = amr_minmod(fz(iv3,i,j+1,k+khp)-fz(iv3,i,j,  k+khp), & 
     &                       fz(iv3,i,j,  k+khp)-fz(iv3,i,j-1,k+khp)) ! dy(i,[j,j-1],k+khp)
         else
         dx_fzp = 0.
         dy_fzp = 0.
         endif

         dy_fxm = amr_minmod(r0*fx(iv1,i+ihm,j+1,k)-r0*fx(iv1,i+ihm,j,  k), & 
     &                       r0*fx(iv1,i+ihm,j,  k)-r0*fx(iv1,i+ihm,j-1,k)) ! dy(i+ihm,[j,j-1],k)
         if (ndim>2) then
         dz_fxm = amr_minmod(r0*fx(iv1,i+ihm,j,k+1)-r0*fx(iv1,i+ihm,j,k  ), & 
     &                       r0*fx(iv1,i+ihm,j,k  )-r0*fx(iv1,i+ihm,j,k-1)) ! dz(i+ihm,j,[k,k-1])
         else
         dz_fxm = 0.
         endif

         dx_fym = amr_minmod(fy(iv2,i+1,j+jhm,k)-fy(iv2,i,  j+jhm,k), & 
     &                       fy(iv2,i,  j+jhm,k)-fy(iv2,i-1,j+jhm,k)) ! dx([i,i-1],j+jhm,k)
         if (ndim>2) then
         dz_fym = amr_minmod(fy(iv2,i,j+jhm,k+1)-fy(iv2,i,j+jhm,k  ), & 
     &                       fy(iv2,i,j+jhm,k  )-fy(iv2,i,j+jhm,k-1)) ! dz(i,j+jhm,[k,k-1])

         dx_fzm = amr_minmod(fz(iv3,i+1,j,k+khm)-fz(iv3,i,  j,k+khm), & 
     &                       fz(iv3,i,j,  k+khm)-fz(iv3,i-1,j,k+khm)) ! dx([i,i-1],j,k+khm)
         dy_fzm = amr_minmod(fz(iv3,i,j+1,k+khm)-fz(iv3,i,j,  k+khm), & 
     &                       fz(iv3,i,j,  k+khm)-fz(iv3,i,j-1,k+khm)) ! dy(i,[j,j-1],k+khm)
         else
         dz_fym = 0.
         dx_fzm = 0.
         dy_fzm = 0.
         endif
         ! Is the following correct?
         if (ndim>2) then
         dxy_fzp =  & 
     & amr_minmod(  & 
     &  fz(iv3,i+1,j+1,k+khp) + fz(iv3,i,  j,  k+khp) -  & 
     &  fz(iv3,i+1,j,  k+khp) - fz(iv3,i,  j+1,k+khp),  & 
     &  fz(iv3,i,  j,  k+khp) + fz(iv3,i-1,j-1,k+khp) -  & 
     &  fz(iv3,i,  j-1,k+khp) - fz(iv3,i-1,j,  k+khp)  & 
     &  )

         dxy_fzm =  & 
     & amr_minmod(  & 
     &  fz(iv3,i+1,j+1,k+khm) + fz(iv3,i,  j,  k+khm) -  & 
     &  fz(iv3,i+1,j,  k+khm) - fz(iv3,i,  j+1,k+khm),  & 
     &  fz(iv3,i,  j,  k+khm) + fz(iv3,i-1,j-1,k+khm) -  & 
     &  fz(iv3,i,  j-1,k+khm) - fz(iv3,i-1,j,  k+khm)  & 
     &  )

         dyz_fxp =  & 
     & amr_minmod(  & 
     &  r1*fx(iv1,i+ihp,j+1,k+1) + r1*fx(iv1,i+ihp,j,  k  ) -  & 
     &  r1*fx(iv1,i+ihp,j+1,k  ) - r1*fx(iv1,i+ihp,j,  k+1),  & 
     &  r1*fx(iv1,i+ihp,j,  k  ) + r1*fx(iv1,i+ihp,j-1,k-1) -  & 
     &  r1*fx(iv1,i+ihp,j,  k-1) - r1*fx(iv1,i+ihp,j-1,k  )  & 
     & )

         dyz_fxm =  & 
     & amr_minmod(  & 
     &  r0*fx(iv1,i+ihm,j+1,k+1) + r0*fx(iv1,i+ihm,j,  k  ) -  & 
     &  r0*fx(iv1,i+ihm,j+1,k  ) - r0*fx(iv1,i+ihm,j,  k+1),  & 
     &  r0*fx(iv1,i+ihm,j,  k  ) + r0*fx(iv1,i+ihm,j-1,k-1) -  & 
     &  r0*fx(iv1,i+ihm,j,  k-1) - r0*fx(iv1,i+ihm,j-1,k  ) & 
     & )

         dxz_fyp =  & 
     & amr_minmod(  & 
     &  fy(iv2,i+1,j+jhp,k+1) + fy(iv2,i,  j+jhp,k  ) -  & 
     &  fy(iv2,i+1,j+jhp,k  ) - fy(iv2,i,  j+jhp,k+1), & 
     &  fy(iv2,i,  j+jhp,k  ) + fy(iv2,i-1,j+jhp,k-1) -  & 
     &  fy(iv2,i,  j+jhp,k-1) - fy(iv2,i-1,j+jhp,k  )  & 
     & )

         dxz_fym =  & 
     & amr_minmod(  & 
     &  fy(iv2,i+1,j+jhm,k+1) + fy(iv2,i,  j+jhm,k  ) -  & 
     &  fy(iv2,i+1,j+jhm,k  ) - fy(iv2,i,  j+jhm,k+1),  & 
     &  fy(iv2,i,  j+jhm,k  ) + fy(iv2,i-1,j+jhm,k-1) -  & 
     &  fy(iv2,i,  j+jhm,k-1) - fy(iv2,i-1,j+jhm,k  ) & 
     & )
         else
         dxy_fzp = 0.
         dxy_fzm = 0.
         dyz_fxp = 0.
         dyz_fxm = 0.
         dxz_fyp = 0.
         dxz_fym = 0.
         endif
        ! ddx := ddx(xp-xm)
        ! For our case, hm == 0, so ddx(...) := ddx(i1,j1,k1)

        ddx2=ddx*ddx
        ddy2=ddy*ddy
        ddz2=ddz*ddz

        axyz =       ( ( dyz_fxp - dyz_fxm ) / (ddy * ddz) ) / ddx !
        bxyz =       ( ( dxz_fyp - dxz_fym ) / (ddz * ddx) ) / ddy !
        cxyz =       ( ( dxy_fzp - dxy_fzm ) / (ddx * ddy) ) / ddz !

        ayz  = 0.5 * ( dyz_fxp + dyz_fxm ) / (ddy * ddz) !
        bxz  = 0.5 * ( dxz_fyp + dxz_fym ) / (ddz * ddx) !
        cxy  = 0.5 * ( dxy_fzp + dxy_fzm ) / (ddx * ddy) !

        ax   =       ( fxp - fxm ) / ddx   !
        by   =       ( fyp - fym ) / ddy   !
        cz   =       ( fzp - fzm ) / ddz   !

        axy  =       ( ( dy_fxp - dy_fxm ) / ddy ) / ddx  !
        byz  =       ( ( dz_fyp - dz_fym ) / ddz ) / ddy  !
        cxz  =       ( ( dx_fzp - dx_fzm ) / ddx ) / ddz  !

        axz  =       ( ( dz_fxp - dz_fxm ) / ddz ) / ddx  !
        bxy  =       ( ( dx_fyp - dx_fym ) / ddx ) / ddy  !
        cyz  =       ( ( dy_fzp - dy_fzm ) / ddy ) / ddz  !

        ay   = 0.5 * ( dy_fxp + dy_fxm ) / ddy + cxyz * ddx2 / 16.0 !
        bz   = 0.5 * ( dz_fyp + dz_fym ) / ddz + axyz * ddy2 / 16.0 !
        cx   = 0.5 * ( dx_fzp + dx_fzm ) / ddx + bxyz * ddz2 / 16.0 !

        az   = 0.5 * ( dz_fxp + dz_fxm ) / ddz + bxyz * ddx2 / 16.0 !
        bx   = 0.5 * ( dx_fyp + dx_fym ) / ddx + cxyz * ddy2 / 16.0 !
        cy   = 0.5 * ( dy_fzp + dy_fzm ) / ddy + axyz * ddz2 / 16.0 !

        axx  = - 0.5 * ( bxy + cxz )  !
        byy  = - 0.5 * ( cyz + axy )  !
        czz  = - 0.5 * ( axz + byz )  !

        a0   = 0.5 * ( fxp + fxm ) - 0.25 * axx * ddx2 !
        b0   = 0.5 * ( fyp + fym ) - 0.25 * byy * ddy2 !
        c0   = 0.5 * ( fzp + fzm ) - 0.25 * czz * ddz2 !

        cyzz = -axyz * 0.25 !
        byyz = cyzz !

        axxz = -bxyz * 0.25 !
        cxzz = axxz !

        byyx = -cxyz * 0.25 
        axxy = byyx

        if(left)then
          fx_0 = ( a0 + ax * x0 + axx * x0 * x0 )  & 
     &     + ( ay + axy * x0 + axxy * x0 * x0 ) * 0.5 * (y1+y0)  & 
     &     + ( az + axz * x0 + axxz * x0 * x0 ) * 0.5 * (z1+z0)  & 
     &     + (ayz + axyz * x0) * 0.25 * (y1 + y0) * (z1 + z0)

          fy_0 = ( b0 + by * y0 + byy * y0 * y0 )  & 
     &     + ( bz + byz * y0 + byyz * y0 * y0) * 0.5 * (z1+z0)  & 
     &     + ( bx + bxy * y0 + byyx * y0 * y0) * 0.5 * (x1+x0) &
     &     + (bxz + bxyz * y0) * 0.25 * (z1 + z0) * (x1 + x0)

          fz_0 = ( c0 + cz * z0 + czz * z0 * z0 )  & 
     &     + ( cx + cxz * z0 + cxzz * z0 * z0) * 0.5 * (x1+x0)  & 
     &     + ( cy + cyz * z0 + cyzz * z0 * z0) * 0.5 * (y1+y0)  & 
     &     + (cxy + cxyz * z0) * 0.25 * (x1 + x0) * (y1 + y0)
         end if

    ! 
    ! These may have to be reconsidered for non-orthogonal cells.
    ! 
         if(right)then
           fx_1 = ( a0 + ax * x1 + axx * x1 * x1 )  & 
     &      + ( ay + axy * x1 + axxy * x1 * x1) * 0.5 * (y1+y0)  & 
     &      + ( az + axz * x1 + axxz * x1 * x1) * 0.5 * (z1+z0)  & 
     &      + (ayz + axyz * x1) * 0.25 * (y1 + y0) * (z1 + z0)

           fy_1 = ( b0 + by * y1 + byy * y1 * y1 )  &
     &      + ( bz + byz * y1 + byyz * y1 * y1) * 0.5 * (z1+z0)  & 
     &      + ( bx + bxy * y1 + byyx * y1 * y1) * 0.5 * (x1+x0)  &
     &      + (bxz + bxyz * y1) * 0.25 * (z1 + z0) * (x1 + x0)
       
           fz_1 = ( c0 + cz * z1 + czz * z1 * z1 )  &
     &      + ( cx + cxz * z1 + cxzz * z1 * z1) * 0.5 * (x1+x0)  & 
     &      + ( cy + cyz * z1 + cyzz * z1 * z1) * 0.5 * (y1+y0)  & 
     &      + (cxy + cxyz * z1) * 0.25 * (x1 + x0) * (y1 + y0)
          end if

        contains
          function amr_minmod(a,b) result(mm)
            use physicaldata, ONLY: cylindrical_pm
            implicit none

            real, intent(in) :: a,b
            real :: mm
            if (a*b > 0.0) then 
               if(abs(a)<abs(b))then
                  mm=a
               else
                  mm=b
               endif
            else
               mm = 0.0
            endif
            if (cylindrical_pm) mm =0.0 !first order !! DEV: Horrible Hack (but avoids code repetition in balpro ;-) P.T.)!!! 
          end function amr_minmod

         end subroutine balpro
