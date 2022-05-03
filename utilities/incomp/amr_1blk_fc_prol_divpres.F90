!----------------------------------------------------------------------
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2009 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!!REORDER(5): facevar[xyz]1
!!REORDER(4): recv[uvw], recvf[xyz]

#include "paramesh_preprocessor.fh"

Module gr_pmDivpres_mod

  ! DIVPRES variables
  public :: prol_fc_divpres, prol_fc_divpres_ivar, prol_fc_divpres_n
  logical, save :: prol_fc_divpres = .false.
  integer, allocatable, save :: prol_fc_divpres_ivar(:,:)
  integer, save :: prol_fc_divpres_n = 0

  interface
     subroutine gr_pmDivpresApply(recvfx,recvfy,recvfz, &
          ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
          lb)
       use paramesh_dimensions, ONLY: il_bnd1,jl_bnd1,kl_bnd1
       implicit none
       Real,intent(inout) :: recvfx(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
       Real,intent(inout) :: recvfy(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
       Real,intent(inout) :: recvfz(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
       integer, intent(in) :: ia,ib,ja,jb,ka,kb
       integer, intent(in) :: idest,ioff,joff,koff
       integer, intent(in) :: lb
     end subroutine gr_pmDivpresApply
  end interface

Contains

      subroutine amr_1blk_fc_prol_divpres     &
       (recvu,recvv,recvw,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
        iv1,iv2,iv3,lb)


!
!------------------------------------------------------------------------
!
! This is a stub routine and is meant to be a place holder to allow
! a user to write their own prolongation routine for cell centered 
! data.  
!
! NOTE: To use this feature you must define interp_mask_face's to be >= 20.
!
! NOTE2: Use one of the other routines which are provided for doing this
! operation as an example.
!
! Marcos Vanella, July 2009
!------------------------------------------------------------------------
!

      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

      implicit none

      include 'mpif.h'

      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff
      integer, intent(in) :: iv1,iv2,iv3,lb
      real,    intent(inout) :: recvu(:,:,:,:), &
               recvv(:,:,:,:),  recvw(:,:,:,:)


      integer :: i,j,k
      integer :: offi,offj,offk
      integer :: ii,jj,kk
      integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
      integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax

      real :: dx,dy,dz,Del,InvAx(3,3),Ui(3,1),a(3,1)
      real :: b(5,1),b2(4,1),usol(4,1)
      real :: U1,U1p,V1,V1p,u_1,u_1p,v_1,v_1p, &
              u_2,u_2p,v_2,v_2p,u_10,u_12,DivU
      real, save :: InvAxi(3,3)
      real, save :: A2T(4,5),InvA2(4,4),f,fin  

      real, save :: InvAxi3d(9,9),fyz,fxz,A3T(8,9),InvA3(8,8)
      real :: InvAx3d(9,9),Ui3(9,1),a3(9,1),A_1(1,9)
      real :: b3(9,1),b3i(8,1),usol3(8,1)
      integer :: i_f,j_f,k_f
      real :: Del1,Del2,W1,W1p,xi,yi
      real :: u_1v(1,1),u_2v(1,1),u_3v(1,1),u_4v(1,1)

      integer,parameter :: largei = 100

      logical,save :: first_call = .true.
      !---------------------------------------------------------------

      ! Parents Cell sizes:
      dx = 2.*bsize(1,lb)/real(nxb)
      dy = 2.*bsize(2,lb)/real(nyb)
      dz = 2.*bsize(3,lb)/real(nzb)

 
      if (first_call) then
         first_call = .false.


      if (ndim .eq. 2) then
      ! InvAxi: 2D
      InvAxi(1,1:3)= (/      23./12.,       -4./3.,       5./12. /)
      InvAxi(2,1:3)= (/     -21./10.,       16./5.,     -11./10. /)
      InvAxi(3,1:3)= (/       8./15.,     -16./15.,       8./15. /)           

      f   = dy/dx ! Ratio dy/dx is always the same.
      fin = 1./f
 
      ! Divergence equations + u1m:
      A2T(1,1:5) = (/  1., -1., 0.,  0., 1. /)
      A2T(2,1:5) = (/  0.,  0., -1., 1., 0. /)
      A2T(3,1:5) = (/  fin,  0., 0., -fin, 0. /)
      A2T(4,1:5) = (/  0.,  fin, -fin, 0., 0. /)

      InvA2(1,1:4) = (/  1.,    -1.,          -f,           f /)
      InvA2(2,1:4) = (/ -1.,     2.,       1.5*f,      -1.5*f /)
      InvA2(3,1:4) = (/  -f,  1.5*f,  1.75*f**2., -1.25*f**2. /)      
      InvA2(4,1:4) = (/   f, -1.5*f, -1.25*f**2.,  1.75*f**2. /)

      elseif (ndim .eq. 3)  then

      ! InvAxi: 3D
      InvAxi3d(1,1:9) = (/   253./72.,     -115./48.,     101./144., &   
                            -115./48.,        16./9.,      -25./48., &
                            101./144.,      -25./48.,       5./36. /)

      InvAxi3d(2,1:9) = (/ -899./240.,        23./4.,    -449./240., &   
                               21./8.,      -64./15.,        11./8., &
                           -179./240.,         5./4.,    -89./240. /)

      InvAxi3d(3,1:9) = (/ -899./240.,        21./8.,    -179./240., &
                               23./4.,      -64./15.,         5./4., &
                           -449./240.,        11./8.,    -89./240. /)

      InvAxi3d(4,1:9) = (/  399./100.,      -63./10.,     199./100., &
                             -63./10.,      256./25.,      -33./10., &
                            199./100.,      -33./10.,     99./100. /)

      InvAxi3d(5,1:9) = (/  337./360.,     -23./12.,      337./360., &
                               -2./3.,      64./45.,         -2./3., &
                             67./360.,      -5./12.,      67./360. /)

      InvAxi3d(6,1:9) = (/  337./360.,       -2./3.,       67./360., &
                             -23./12.,      64./45.,        -5./12., &
                            337./360.,       -2./3.,      67./360. /)

      InvAxi3d(7,1:9) = (/ -299./300.,      21./10.,     -299./300., &
                                8./5.,    -256./75.,          8./5., &
                           -149./300.,      11./10.,    -149./300. /)

      InvAxi3d(8,1:9) = (/ -299./300.,        8./5.,     -149./300., &
                              21./10.,    -256./75.,        11./10., &
                           -299./300.,        8./5.,    -149./300. /)

      InvAxi3d(9,1:9) = (/   56./225.,      -8./15.,       56./225., &
                              -8./15.,    256./225.,        -8./15., &
                             56./225.,      -8./15.,      56./225. /)

      fyz = dy/dz
      fxz = dx/dz

      ! Divergence equations + Wi,j+1,k+! interpolation
      A3T(1,1:9)=(/ 0.,  0.,  0.,  0., 1./fxz, -1./fxz,  0.,  0.,  0. /)
      A3T(2,1:9)=(/ 0.,  0., 1./fxz, -1./fxz,  0.,  0.,  0.,  0.,  0. /)
      A3T(3,1:9)=(/ 0., 1./fyz,  0., -1./fyz,  0.,  0.,  0.,  0.,  0. /)
      A3T(4,1:9)=(/ 0.,  0.,  0.,  0., 1./fyz, 0., -1./fyz,  0.,   0. /)
      A3T(5,1:9)=(/ 1.,  0.,  0.,  0.,   -1.,  0.,  0.,  0.,  0. /)
      A3T(6,1:9)=(/ 0.,  1.,  0.,  0.,    0., -1.,  0.,  0.,  0. /)
      A3T(7,1:9)=(/ 0.,  0.,  1.,  0.,    0.,  0., -1.,  0.,  1. /)
      A3T(8,1:9)=(/ 0.,  0.,  0.,  1.,    0.,  0.,  0., -1.,  0. /)

      InvA3(1,1:8)=(/ 23./8.*fxz**2., -11./8.*fxz**2., 17./8.*fxz*fyz,  & 
               -13./8.*fxz*fyz, 5./8.*fxz, -5./2.*fxz, fxz, 3./8.*fxz /)
      InvA3(2,1:8)=(/ -11./8.*fxz**2., 15./8.*fxz**2., -13./8.*fxz*fyz, &
               9./8.*fxz*fyz, -1./8.*fxz, 3./2.*fxz, -fxz, 1./8.*fxz  /)
      InvA3(3,1:8)=(/ 17./8.*fxz*fyz, -13./8.*fxz*fyz, 23./8.*fyz**2.,  &
               -11./8.*fyz**2., 3./8.*fyz, -5./2.*fyz, fyz, 5./8.*fyz /)
      InvA3(4,1:8)=(/ -13./8.*fxz*fyz, 9./8.*fxz*fyz, -11./8.*fyz**2.,  &
               15./8.*fyz**2., 1./8.*fyz, 3./2.*fyz, -fyz, -1./8.*fyz /)
      InvA3(5,1:8)=(/  5./8.*fxz, -1./8.*fxz, 3./8.*fyz, 1./8.*fyz,     &
               7./8.,       -1./2.,            0.,          1./8.     /)
      InvA3(6,1:8)=(/ -5./2.*fxz, 3./2.*fxz, -5./2.*fyz, 3./2.*fyz,     &
              -1./2.,           3.,           -1.,         -1./2.     /)
      InvA3(7,1:8)=(/        fxz,      -fxz,        fyz,      -fyz,     &
                  0.,          -1.,            1.,             0.     /)
      InvA3(8,1:8)=(/  3./8.*fxz, 1./8.*fxz,  5./8.*fyz,  -1./8.*fyz,   &
               1./8.,       -1./2.,            0.,          7./8.     /)


      endif

      endif

      ifmin=ia
      ifmax=ib
      jfmin=ja
      jfmax=jb
      kfmin=ka
      kfmax=kb

      offi = 0
      offj = 0
      offk = 0
      if(ioff.gt.0) offi = nxb/2
      if(joff.gt.0) offj = nyb*k2d/2
      if(koff.gt.0) offk = nzb*k3d/2

      kcmin = ((kfmin-nguard-1+largei)/2 +       &
                      nguard - largei/2 )*k3d +  &
                      1 + offk
      kcmax = ((kfmax-nguard-1+largei)/2 +       &
                      nguard - largei/2 )*k3d +  &
                      1 + offk
      jcmin = ((jfmin-nguard-1+largei)/2 +       &
                      nguard - largei/2 )*k2d +  &
                      1 + offj
      jcmax = ((jfmax-nguard-1+largei)/2 +       &
                      nguard - largei/2 )*k2d +  &
                      1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +       &
                      nguard - largei/2 ) +      &
                      1 + offi
      icmax = ((ifmax-nguard-1+largei)/2 +       &
                      nguard - largei/2 ) +      &
                      1 + offi




      if (ndim .eq. 2) then       ! 2D CASE
 

      ! Interpolate faces in X:
      Del = dy    
      InvAx(1,1:3) =              InvAxi(1,1:3)
      InvAx(2,1:3) = 1./Del     * InvAxi(2,1:3)
      InvAx(3,1:3) = 1./Del**2. * InvAxi(3,1:3)

!!$c$$$      InvAx(1,1:3) = (/      23./12.,       -4./3.,       5./12. /)
!!$c$$$      InvAx(2,1:3) = (/ -21./10./Del,   16./5./Del, -11./10./Del /)
!!$c$$$      InvAx(3,1:3) = 
!!$c$$$&            (/ 8./15./Del**2., -16./15./Del**2., 8./15./Del**2. /)

       
      do k = kcmin,kcmax
         do j = jcmin,jcmax
            do i = icmin,icmax+1

            Ui(1,1) = recvu(iv1,i,j-1,k); !U(i,j-1);
            Ui(2,1) = recvu(iv1,i,j,k);   !U(i,j);
            Ui(3,1) = recvu(iv1,i,j+1,k); !U(i,j+1);
                  
            a = MATMUL(InvAx,Ui);
 
           facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest) &
            = a(1,1) + a(2,1)*(5./4.)*Del + a(3,1)*(5./4.*Del)**2.


           facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+1,k,idest) &
            = a(1,1) + a(2,1)*(7./4.)*Del + a(3,1)*(7./4.*Del)**2.
        
           enddo
        enddo
      enddo


      ! Y direction
      Del = dx
      InvAx(1,1:3) =              InvAxi(1,1:3)
      InvAx(2,1:3) = 1./Del     * InvAxi(2,1:3)
      InvAx(3,1:3) = 1./Del**2. * InvAxi(3,1:3)

!!$c$$$      InvAx(1,1:3) = (/      23./12.,       -4./3.,       5./12. /)
!!$c$$$      InvAx(2,1:3) = (/ -21./10./Del,   16./5./Del, -11./10./Del /)
!!$c$$$      InvAx(3,1:3) = 
!!$c$$$&            (/ 8./15./Del**2., -16./15./Del**2., 8./15./Del**2. /)


      do k = kcmin,kcmax
         do j = jcmin,jcmax+1
            do i = icmin,icmax

            Ui(1,1) = recvv(iv2,i-1,j,k); !V(i-1,j);
            Ui(2,1) = recvv(iv2,i,j,k);   !V(i,j);
            Ui(3,1) = recvv(iv2,i+1,j,k); !V(i+1,j);
                  
            a = MATMUL(InvAx,Ui);
 
           facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest) &
            = a(1,1) + a(2,1)*(5./4.)*Del + a(3,1)*(5./4.*Del)**2.


           facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin),k,idest) &
            = a(1,1) + a(2,1)*(7./4.)*Del + a(3,1)*(7./4.*Del)**2.

           enddo
        enddo
      enddo


!!$c$$$      ! Divergence on Parent Cells
!!$c$$$      do k = kcmin,kcmax
!!$c$$$         do j = jcmin,jcmax
!!$c$$$            do i = icmin,icmax
!!$c$$$            ! Fine grid values of velocities on faces of parent Cell
!!$c$$$            u_1p = 
!!$c$$$&           facevarx1(iv1,ifmin+2*(i-icmin)+2,jfmin+2*(j-jcmin),k,idest)
!!$c$$$
!!$c$$$            u_2p =
!!$c$$$&         facevarx1(iv1,ifmin+2*(i-icmin)+2,jfmin+2*(j-jcmin)+1,k,idest)
!!$c$$$
!!$c$$$            u_1 = 
!!$c$$$&           facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest)
!!$c$$$
!!$c$$$            u_2 =
!!$c$$$&         facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+1,k,idest)
!!$c$$$
!!$c$$$            v_1p =
!!$c$$$&           facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+2,k,idest)
!!$c$$$
!!$c$$$            v_2p =
!!$c$$$&         facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin)+2,k,idest)
!!$c$$$
!!$c$$$            v_1 =
!!$c$$$&           facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest)
!!$c$$$
!!$c$$$            v_2 =
!!$c$$$&         facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin),k,idest)
!!$c$$$
!!$c$$$            DivU = 0.5*((u_1p+u_2p)-(u_1+u_2))/dx + 
!!$c$$$&                  0.5*((v_1p+v_2p)-(v_1+v_2))/dy
!!$c$$$
!!$c$$$c$$$          write(*,*) 'DivU=',DivU
!!$c$$$
!!$c$$$           enddo
!!$c$$$        enddo
!!$c$$$      enddo


      ! Fill internal variables:
      ! Divergence equations + u1m:
!!$c$$$      f   = dy/dx
!!$c$$$      fin = 1./f
!!$c$$$      A2T(1,1:5) = (/  1., -1., 0.,  0., 1. /)
!!$c$$$      A2T(2,1:5) = (/  0.,  0., -1., 1., 0. /)
!!$c$$$      A2T(3,1:5) = (/  fin,  0., 0., -fin, 0. /)
!!$c$$$      A2T(4,1:5) = (/  0.,  fin, -fin, 0., 0. /)
!!$c$$$
!!$c$$$      InvA2(1,1:4) = (/  1.,    -1.,          -f,           f /)
!!$c$$$      InvA2(2,1:4) = (/ -1.,     2.,       1.5*f,      -1.5*f /)
!!$c$$$      InvA2(3,1:4) = (/  -f,  1.5*f,  1.75*f**2., -1.25*f**2. /)      
!!$c$$$      InvA2(4,1:4) = (/   f, -1.5*f, -1.25*f**2.,  1.75*f**2. /)


      do k = kcmin,kcmax
         do j = jcmin,jcmax
            do i = icmin,icmax

            ! Values of U and V on parent Cell
            U1  = recvu(iv1,i,j,k);  
            U1p = recvu(iv1,i+1,j,k);
        
            V1  = recvv(iv2,i,j,k);
            V1p = recvv(iv2,i,j+1,k);

            DivU = (U1p-U1)/dx + (V1p-V1)/dy

       
            ! Fine grid values of velocities on faces of parent Cell
            u_1p = &
            facevarx1(iv1,ifmin+2*(i-icmin)+2,jfmin+2*(j-jcmin),k,idest)

            u_2p = &
          facevarx1(iv1,ifmin+2*(i-icmin)+2,jfmin+2*(j-jcmin)+1,k,idest)

            u_1 =  &
            facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest)

            u_2 =  &
          facevarx1(iv1,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+1,k,idest)

            v_1p = &
            facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+2,k,idest)

            v_2p = &
          facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin)+2,k,idest)

            v_1 =  &
            facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin),k,idest)

            v_2 =  &
          facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin),k,idest)

            u_10=  &
            facevarx1(iv1,ifmin+2*(i-icmin)-2,jfmin+2*(j-jcmin),k,idest)

            u_12=  &
            facevarx1(iv1,ifmin+2*(i-icmin)+4,jfmin+2*(j-jcmin),k,idest)            

            ! Build Right Hand Side:
            ! Using Divergence equations + u1m linear expression
            ! Solve for pseudo-inverse (least-squares)        
            b(1,1) = dx*(DivU/2. + (  u_1/dx +  v_1/dy));
            b(2,1) = dx*(DivU/2. + (-u_1p/dx +  v_2/dy));
            b(3,1) = dx*(DivU/2. - ( u_2p/dx + v_2p/dy));
            b(4,1) = dx*(DivU/2. + (  u_2/dx - v_1p/dy));
!!$c$$$            b(5,1) =(u_1+u_1p)/2.

            if (i .eq. icmin) then
             b(5,1) = 0.375*u_1+0.75*u_1p-0.125*u_12
            else
             b(5,1) =-0.125*u_10+0.75*u_1+0.375*u_1p
            endif
        
            b2   =  MATMUL(A2T,b)

            usol =  MATMUL(InvA2,b2)


            ! Dump result on internal velocities of parent cell
            facevarx1(iv1,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin),k,idest) &
            = usol(1,1)

          facevarx1(iv1,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin)+1,k,idest) &
            = usol(2,1)

            facevary1(iv2,ifmin+2*(i-icmin),jfmin+2*(j-jcmin)+1,k,idest) &
            = usol(3,1)

          facevary1(iv2,ifmin+2*(i-icmin)+1,jfmin+2*(j-jcmin)+1,k,idest) &
            = usol(4,1)


           enddo
        enddo
      enddo
  

!!$      ! Check final divergence
!!$c$$$      write(*,*) 'Block Number=',lb
!!$c$$$      do k = kfmin,kfmax
!!$c$$$         do j = jfmin,jfmax
!!$c$$$            do i = ifmin,ifmax
!!$c$$$
!!$c$$$
!!$c$$$            ! Fine grid values of velocities on faces of parent Cell
!!$c$$$            u_1p = 
!!$c$$$&           facevarx1(iv1,i+1,j,k,idest)
!!$c$$$
!!$c$$$            u_1 = 
!!$c$$$&           facevarx1(iv1,i,j,k,idest)
!!$c$$$
!!$c$$$            v_1p =
!!$c$$$&           facevary1(iv2,i,j+1,k,idest)
!!$c$$$
!!$c$$$            v_1 =
!!$c$$$&           facevary1(iv2,i,j,k,idest)
!!$c$$$
!!$c$$$            DivU = (u_1p-u_1)/dx + 
!!$c$$$&                  (v_1p-v_1)/dy
!!$c$$$
!!$c$$$!          write(*,*) 'DivU=',DivU
!!$c$$$
!!$c$$$           enddo
!!$c$$$        enddo
!!$c$$$      enddo

      elseif(ndim .eq. 3) then    ! 3D CASE


         ! X interpolate, 2nd order flux preserving:
         Del1 = dy
         Del2 = dz

         InvAx3d(1,1:9) =                          InvAxi3d(1,1:9)
         InvAx3d(2,1:9) = 1./Del1             *    InvAxi3d(2,1:9)
         InvAx3d(3,1:9) = 1./Del2             *    InvAxi3d(3,1:9)
         InvAx3d(4,1:9) = 1./(Del1*Del2)      *    InvAxi3d(4,1:9)
         InvAx3d(5,1:9) = 1./(Del1**2.)       *    InvAxi3d(5,1:9)
         InvAx3d(6,1:9) = 1./(Del2**2.)       *    InvAxi3d(6,1:9)
         InvAx3d(7,1:9) = 1./(Del1**2.*Del2)  *    InvAxi3d(7,1:9)
         InvAx3d(8,1:9) = 1./(Del1*Del2**2.)  *    InvAxi3d(8,1:9)
         InvAx3d(9,1:9) = 1./(Del1**2. * Del2**2.)*InvAxi3d(9,1:9)

         do k = kcmin,kcmax
           do j = jcmin,jcmax
              do i = icmin,icmax+1

            i_f = ifmin+2*(i-icmin);
            j_f = jfmin+2*(j-jcmin);
            k_f = kfmin+2*(k-kcmin);

            Ui3(1,1) = recvu(iv1,i,j-1,k-1); !U(i,j-1,k-1);
            Ui3(2,1) = recvu(iv1,i,j,k-1);   !U(i,j,k-1);
            Ui3(3,1) = recvu(iv1,i,j+1,k-1); !U(i,j+1,k-1);

            Ui3(4,1) = recvu(iv1,i,j-1,k);   !U(i,j-1,k);
            Ui3(5,1) = recvu(iv1,i,j,k);     !U(i,j,k);
            Ui3(6,1) = recvu(iv1,i,j+1,k);   !U(i,j+1,k);

            Ui3(7,1) = recvu(iv1,i,j-1,k+1); !U(i,j-1,k+1);
            Ui3(8,1) = recvu(iv1,i,j,k+1);   !U(i,j,k+1);
            Ui3(9,1) = recvu(iv1,i,j+1,k+1); !U(i,j+1,k+1);

            a3 = MATMUL(InvAx3d,Ui3);

            ! Point u1, x=5/4*Delx; y=5/4*Dely
            xi=1.25*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_1v = MATMUL(A_1,a3)
            facevarx1(iv1,i_f,j_f,k_f,idest) = u_1v(1,1) 

            ! Point u2, x=7/4*Delx; y=5/4*Dely
            xi=1.75*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_2v = MATMUL(A_1,a3)
            facevarx1(iv1,i_f,j_f+1,k_f,idest) = u_2v(1,1) 

            ! Point u3, x=5/4*Delx; y=7/4*Dely
            xi=1.25*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_3v = MATMUL(A_1,a3)
            facevarx1(iv1,i_f,j_f,k_f+1,idest) = u_3v(1,1) 

            ! Point u4, x=7/4*Delx; y=7/4*Dely
            xi=1.75*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_4v = MATMUL(A_1,a3)
            facevarx1(iv1,i_f,j_f+1,k_f+1,idest) = u_4v(1,1) 
        
            ! Store in Uf:
            !Uf(x_i+2*(i-ng-1),y_i+2*(j-ng-1),z_i+2*(k-ng-1))     = u1;
            !Uf(x_i+2*(i-ng-1),y_i+2*(j-ng-1)+1,z_i+2*(k-ng-1))   = u2;      
            !Uf(x_i+2*(i-ng-1),y_i+2*(j-ng-1),z_i+2*(k-ng-1)+1)   = u3;
            !Uf(x_i+2*(i-ng-1),y_i+2*(j-ng-1)+1,z_i+2*(k-ng-1)+1) = u4;             

             enddo
           enddo
         enddo


         ! Y interpolate, 2nd order flux preserving:
         Del1 = dx
         Del2 = dz

         InvAx3d(1,1:9) =                          InvAxi3d(1,1:9)
         InvAx3d(2,1:9) = 1./Del1             *    InvAxi3d(2,1:9)
         InvAx3d(3,1:9) = 1./Del2             *    InvAxi3d(3,1:9)
         InvAx3d(4,1:9) = 1./(Del1*Del2)      *    InvAxi3d(4,1:9)
         InvAx3d(5,1:9) = 1./(Del1**2.)       *    InvAxi3d(5,1:9)
         InvAx3d(6,1:9) = 1./(Del2**2.)       *    InvAxi3d(6,1:9)
         InvAx3d(7,1:9) = 1./(Del1**2.*Del2)  *    InvAxi3d(7,1:9)
         InvAx3d(8,1:9) = 1./(Del1*Del2**2.)  *    InvAxi3d(8,1:9)
         InvAx3d(9,1:9) = 1./(Del1**2. * Del2**2.)*InvAxi3d(9,1:9)

         do k = kcmin,kcmax
           do j = jcmin,jcmax+1
              do i = icmin,icmax

            i_f = ifmin+2*(i-icmin);
            j_f = jfmin+2*(j-jcmin);
            k_f = kfmin+2*(k-kcmin);

            Ui3(1,1) = recvv(iv2,i-1,j,k-1); !V(i-1,j,k-1);
            Ui3(2,1) = recvv(iv2,i,j,k-1);   !V(i,j,k-1);
            Ui3(3,1) = recvv(iv2,i+1,j,k-1); !V(i+1,j,k-1);

            Ui3(4,1) = recvv(iv2,i-1,j,k);   !V(i-1,j,k);
            Ui3(5,1) = recvv(iv2,i,j,k);     !V(i,j,k);
            Ui3(6,1) = recvv(iv2,i+1,j,k);   !V(i+1,j,k);

            Ui3(7,1) = recvv(iv2,i-1,j,k+1); !V(i-1,j,k+1);
            Ui3(8,1) = recvv(iv2,i,j,k+1);   !V(i,j,k+1);
            Ui3(9,1) = recvv(iv2,i+1,j,k+1); !V(i+1,j,k+1);

            a3 = MATMUL(InvAx3d,Ui3);

            ! Point u1, x=5/4*Delx; y=5/4*Dely
            xi=1.25*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_1v = MATMUL(A_1,a3)
            facevary1(iv2,i_f,j_f,k_f,idest) = u_1v(1,1) 

            ! Point u2, x=7/4*Delx; y=5/4*Dely
            xi=1.75*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_2v = MATMUL(A_1,a3)
            facevary1(iv2,i_f+1,j_f,k_f,idest) = u_2v(1,1) 

            ! Point u3, x=5/4*Delx; y=7/4*Dely
            xi=1.25*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_3v = MATMUL(A_1,a3)
            facevary1(iv2,i_f,j_f,k_f+1,idest) = u_3v(1,1) 

            ! Point u4, x=7/4*Delx; y=7/4*Dely
            xi=1.75*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_4v = MATMUL(A_1,a3)
            facevary1(iv2,i_f+1,j_f,k_f+1,idest) = u_4v(1,1) 
        
            ! Store in Vf:
            !Vf(x_i+2*(i-ng-1)  ,y_i+2*(j-ng-1),z_i+2*(k-ng-1))   = u1;
            !Vf(x_i+2*(i-ng-1)+1,y_i+2*(j-ng-1),z_i+2*(k-ng-1))   = u2;      
            !Vf(x_i+2*(i-ng-1)  ,y_i+2*(j-ng-1),z_i+2*(k-ng-1)+1) = u3;
            !Vf(x_i+2*(i-ng-1)+1,y_i+2*(j-ng-1),z_i+2*(k-ng-1)+1) = u4;                    

             enddo
           enddo
         enddo


         ! Z interpolate, 2nd order flux preserving:
         Del1 = dx
         Del2 = dy

         InvAx3d(1,1:9) =                          InvAxi3d(1,1:9)
         InvAx3d(2,1:9) = 1./Del1             *    InvAxi3d(2,1:9)
         InvAx3d(3,1:9) = 1./Del2             *    InvAxi3d(3,1:9)
         InvAx3d(4,1:9) = 1./(Del1*Del2)      *    InvAxi3d(4,1:9)
         InvAx3d(5,1:9) = 1./(Del1**2.)       *    InvAxi3d(5,1:9)
         InvAx3d(6,1:9) = 1./(Del2**2.)       *    InvAxi3d(6,1:9)
         InvAx3d(7,1:9) = 1./(Del1**2.*Del2)  *    InvAxi3d(7,1:9)
         InvAx3d(8,1:9) = 1./(Del1*Del2**2.)  *    InvAxi3d(8,1:9)
         InvAx3d(9,1:9) = 1./(Del1**2. * Del2**2.)*InvAxi3d(9,1:9)

         do k = kcmin,kcmax+1
           do j = jcmin,jcmax
              do i = icmin,icmax

            i_f = ifmin+2*(i-icmin);
            j_f = jfmin+2*(j-jcmin);
            k_f = kfmin+2*(k-kcmin);

            Ui3(1,1) = recvw(iv3,i-1,j-1,k); !W(i-1,j-1,k);
            Ui3(2,1) = recvw(iv3,i,j-1,k);   !W(i,j-1,k);
            Ui3(3,1) = recvw(iv3,i+1,j-1,k); !W(i+1,j-1,k);

            Ui3(4,1) = recvw(iv3,i-1,j,k);   !W(i-1,j,k);
            Ui3(5,1) = recvw(iv3,i,j,k);     !W(i,j,k);
            Ui3(6,1) = recvw(iv3,i+1,j,k);   !W(i+1,j,k);

            Ui3(7,1) = recvw(iv3,i-1,j+1,k); !W(i-1,j+1,k);
            Ui3(8,1) = recvw(iv3,i,j+1,k);   !W(i,j+1,k);
            Ui3(9,1) = recvw(iv3,i+1,j+1,k); !W(i+1,j+1,k);

            a3 = MATMUL(InvAx3d,Ui3);

            ! Point u1, x=5/4*Delx; y=5/4*Dely
            xi=1.25*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., & 
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_1v = MATMUL(A_1,a3)
            facevarz1(iv3,i_f,j_f,k_f,idest) = u_1v(1,1) 

            ! Point u2, x=7/4*Delx; y=5/4*Dely
            xi=1.75*Del1; yi=1.25*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_2v = MATMUL(A_1,a3)
            facevarz1(iv3,i_f+1,j_f,k_f,idest) = u_2v(1,1)

            ! Point u3, x=5/4*Delx; y=7/4*Dely
            xi=1.25*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_3v = MATMUL(A_1,a3)
            facevarz1(iv3,i_f,j_f+1,k_f,idest) = u_3v(1,1)

            ! Point u4, x=7/4*Delx; y=7/4*Dely
            xi=1.75*Del1; yi=1.75*Del2;
            A_1(1,1:9)  = (/ 1., xi, yi, xi*yi, xi**2., yi**2., &
                   xi**2.*yi, xi*yi**2., xi**2*yi**2 /)
            u_4v = MATMUL(A_1,a3)
            facevarz1(iv3,i_f+1,j_f+1,k_f,idest) = u_4v(1,1)
        
            ! Store in Wf:
            !Wf(x_i+2*(i-ng-1)  ,y_i+2*(j-ng-1),z_i+2*(k-ng-1))   = u1;
            !Wf(x_i+2*(i-ng-1)+1,y_i+2*(j-ng-1),z_i+2*(k-ng-1))   = u2;      
            !Wf(x_i+2*(i-ng-1)  ,y_i+2*(j-ng-1)+1,z_i+2*(k-ng-1)) = u3;
            !Wf(x_i+2*(i-ng-1)+1,y_i+2*(j-ng-1)+1,z_i+2*(k-ng-1)) = u4;  

             enddo
           enddo
         enddo

!!$c$$$         DivU = 0.
!!$c$$$         do k = kcmin,kcmax
!!$c$$$           do j = jcmin,jcmax
!!$c$$$              do i = icmin,icmax
!!$c$$$
!!$c$$$            i_f = ifmin+2*(i-icmin);
!!$c$$$            j_f = jfmin+2*(j-jcmin);
!!$c$$$            k_f = kfmin+2*(k-kcmin);
!!$c$$$
!!$c$$$              U1 =0.25*((facevarx1(1,i_f+2,j_f,k_f,idest)   +
!!$c$$$&                         facevarx1(1,i_f+2,j_f+1,k_f,idest) +
!!$c$$$&                         facevarx1(1,i_f+2,j_f,k_f+1,idest) +
!!$c$$$&                         facevarx1(1,i_f+2,j_f+1,k_f+1,idest))-
!!$c$$$&                        (facevarx1(1,i_f,j_f,k_f,idest)   +
!!$c$$$&                         facevarx1(1,i_f,j_f+1,k_f,idest) +
!!$c$$$&                         facevarx1(1,i_f,j_f,k_f+1,idest) +
!!$c$$$&                         facevarx1(1,i_f,j_f+1,k_f+1,idest)))/dx +
!!$c$$$&                  0.25*((facevary1(1,i_f,j_f+2,k_f,idest)   +
!!$c$$$&                         facevary1(1,i_f+1,j_f+2,k_f,idest) +
!!$c$$$&                         facevary1(1,i_f,j_f+2,k_f+1,idest) +
!!$c$$$&                         facevary1(1,i_f+1,j_f+2,k_f+1,idest))-
!!$c$$$&                        (facevary1(1,i_f,j_f,k_f,idest)   +
!!$c$$$&                         facevary1(1,i_f+1,j_f,k_f,idest) +
!!$c$$$&                         facevary1(1,i_f,j_f,k_f+1,idest) +
!!$c$$$&                         facevary1(1,i_f+1,j_f,k_f+1,idest)))/dy +
!!$c$$$&                  0.25*((facevarz1(1,i_f,j_f,k_f+2,idest)   +
!!$c$$$&                         facevarz1(1,i_f+1,j_f,k_f+2,idest) +
!!$c$$$&                         facevarz1(1,i_f,j_f+1,k_f+2,idest) +
!!$c$$$&                         facevarz1(1,i_f+1,j_f+1,k_f+2,idest))-
!!$c$$$&                        (facevarz1(1,i_f,j_f,k_f,idest)   +
!!$c$$$&                         facevarz1(1,i_f+1,j_f,k_f,idest) +
!!$c$$$&                         facevarz1(1,i_f,j_f+1,k_f,idest) +
!!$c$$$&                         facevarz1(1,i_f+1,j_f+1,k_f,idest)))/dz                  
!!$c$$$
!!$c$$$               DivU = max(DivU,abs(U1))
!!$c$$$
!!$c$$$              enddo
!!$c$$$           enddo
!!$c$$$        enddo
!!$c$$$        write(*,*) 'lb=',lb,', DivU parent-Interpolated=',DivU



         ! Fill MetaCell internal velocities, divergence preserving:
         do k = kcmin,kcmax
           do j = jcmin,jcmax
              do i = icmin,icmax

            i_f = ifmin+2*(i-icmin);
            j_f = jfmin+2*(j-jcmin);
            k_f = kfmin+2*(k-kcmin);         


            ! Divergence of U:
            U1  = recvu(iv1,i,j,k);  
            U1p = recvu(iv1,i+1,j,k);
        
            V1  = recvv(iv2,i,j,k);
            V1p = recvv(iv2,i,j+1,k);

            W1  = recvw(iv3,i,j,k);
            W1p = recvw(iv3,i,j,k+1);

            DivU = (U1p-U1)/dx + (V1p-V1)/dy + (W1p-W1)/dz

            
            if (i .eq. icmin) then
!!$c$$$            ! Ui+1,j,k
!!$c$$$            facevarx1(iv1,i_f+1,j_f,k_f,idest) = 
!!$c$$$&           0.5*(facevarx1(iv1,i_f+2,j_f,k_f,idest)+
!!$c$$$&                facevarx1(iv1,i_f,j_f,k_f,idest))
!!$c$$$        
!!$c$$$            ! Ui+1,j+1,k+1
!!$c$$$            facevarx1(iv1,i_f+1,j_f+1,k_f+1,idest) = 
!!$c$$$&           0.5*(facevarx1(iv1,i_f+2,j_f+1,k_f+1,idest)+
!!$c$$$&                facevarx1(iv1,i_f,j_f+1,k_f+1,idest))

            ! Ui+1,j,k
            facevarx1(iv1,i_f+1,j_f,k_f,idest) =        &
            0.375*facevarx1(iv1,i_f,j_f,k_f,idest)   +  &
            0.75 *facevarx1(iv1,i_f+2,j_f,k_f,idest) -  &
            0.125*facevarx1(iv1,i_f+4,j_f,k_f,idest)
        
            ! Ui+1,j+1,k+1
            facevarx1(iv1,i_f+1,j_f+1,k_f+1,idest) =    &
            0.375*facevarx1(iv1,i_f,j_f+1,k_f+1,idest)   +  &
            0.75 *facevarx1(iv1,i_f+2,j_f+1,k_f+1,idest) -  &
            0.125*facevarx1(iv1,i_f+4,j_f+1,k_f+1,idest)
       
            else

            ! Ui+1,j,k
            facevarx1(iv1,i_f+1,j_f,k_f,idest) =        &
            -0.125*facevarx1(iv1,i_f+2,j_f,k_f,idest) + &
             0.75*facevarx1(iv1,i_f,j_f,k_f,idest)    + &
             0.375 *facevarx1(iv1,i_f+2,j_f,k_f,idest) 

        
            ! Ui+1,j+1,k+1
            facevarx1(iv1,i_f+1,j_f+1,k_f+1,idest) =        &
            -0.125*facevarx1(iv1,i_f-2,j_f+1,k_f+1,idest) + &
             0.75*facevarx1(iv1,i_f,j_f+1,k_f+1,idest)    + &
             0.375 *facevarx1(iv1,i_f+2,j_f+1,k_f+1,idest) 

            endif


            if (j .eq. jcmin) then
!!$c$$$            ! Vi,j+1,k
!!$c$$$            facevary1(iv2,i_f,j_f+1,k_f,idest) = 
!!$c$$$&           0.5*(facevary1(iv2,i_f,j_f+2,k_f,idest)+
!!$c$$$&                facevary1(iv2,i_f,j_f,k_f,idest))
!!$c$$$        
!!$c$$$            ! Vi+1,j+1,k+1
!!$c$$$            facevary1(iv2,i_f+1,j_f+1,k_f+1,idest) = 
!!$c$$$&           0.5*(facevary1(iv2,i_f+1,j_f+2,k_f+1,idest)+
!!$c$$$&                facevary1(iv2,i_f+1,j_f,k_f+1,idest));


            ! Vi,j+1,k
            facevary1(iv2,i_f,j_f+1,k_f,idest) = &
            0.375*facevary1(iv2,i_f,j_f,k_f,idest)   + &
            0.75 *facevary1(iv2,i_f,j_f+2,k_f,idest) - &
            0.125*facevary1(iv2,i_f,j_f+4,k_f,idest)
        
            ! Vi+1,j+1,k+1
            facevary1(iv2,i_f+1,j_f+1,k_f+1,idest) =  &
            0.375*facevary1(iv2,i_f+1,j_f,k_f+1,idest)   + &
            0.75 *facevary1(iv2,i_f+1,j_f+2,k_f+1,idest) - &
            0.125*facevary1(iv2,i_f+1,j_f+4,k_f+1,idest)

            else

            ! Vi,j+1,k
            facevary1(iv2,i_f,j_f+1,k_f,idest) =  &
            -0.125*facevary1(iv2,i_f,j_f-2,k_f,idest)  + &
             0.75 *facevary1(iv2,i_f,j_f,k_f,idest)    + &
             0.375 *facevary1(iv2,i_f,j_f+2,k_f,idest) 

        
            ! Vi+1,j+1,k+1
            facevary1(iv2,i_f+1,j_f+1,k_f+1,idest) =  &
            -0.125*facevary1(iv2,i_f+1,j_f-2,k_f+1,idest) + &
             0.75*facevary1(iv2,i_f+1,j_f,k_f+1,idest)    + &
             0.375 *facevary1(iv2,i_f+1,j_f+2,k_f+1,idest) 

            endif


            ! b values
            b3(1,1) =  DivU - (facevarx1(iv1,i_f+1,j_f,k_f,idest)-  &
                               facevarx1(iv1,i_f,j_f,k_f,idest))/dx &
                            - (facevary1(iv2,i_f,j_f+1,k_f,idest)-  &
                               facevary1(iv2,i_f,j_f,k_f,idest))/dy &
                            +  facevarz1(iv3,i_f,j_f,k_f,idest)/dz
                              
            b3(2,1) =  DivU - (facevarx1(iv1,i_f+2,j_f,k_f,idest)-    &
                               facevarx1(iv1,i_f+1,j_f,k_f,idest))/dx &
                            +  facevarz1(iv3,i_f+1,j_f,k_f,idest)/dz  &
                            +  facevary1(iv2,i_f+1,j_f,k_f,idest)/dy          
        
            b3(3,1) =  DivU - (facevary1(iv2,i_f,j_f+2,k_f,idest)- &
                               facevary1(iv2,i_f,j_f+1,k_f,idest))/dy &  
                            +  facevarz1(iv3,i_f,j_f+1,k_f,idest)/dz  &
                            +  facevarx1(iv1,i_f,j_f+1,k_f,idest)/dx
                   
            b3(4,1) =  DivU -  facevary1(iv2,i_f+1,j_f+2,k_f,idest)/dy &
                            +  facevarz1(iv3,i_f+1,j_f+1,k_f,idest)/dz &
                            -  facevarx1(iv1,i_f+2,j_f+1,k_f,idest)/dx 
        
            b3(5,1) =  DivU +  facevary1(iv2,i_f,j_f,k_f+1,idest)/dy &
                            -  facevarz1(iv3,i_f,j_f,k_f+2,idest)/dz &
                            +  facevarx1(iv1,i_f,j_f,k_f+1,idest)/dx
                   
            b3(6,1) =  DivU - (facevary1(iv2,i_f+1,j_f+1,k_f+1,idest)- &
                               facevary1(iv2,i_f+1,j_f,k_f+1,idest))/dy & 
                            -  facevarz1(iv3,i_f+1,j_f,k_f+2,idest)/dz  &
                            -  facevarx1(iv1,i_f+2,j_f,k_f+1,idest)/dx
        
            b3(7,1) =  DivU - (facevarx1(iv1,i_f+1,j_f+1,k_f+1,idest)- &
                               facevarx1(iv1,i_f,j_f+1,k_f+1,idest))/dx &
                            -  facevarz1(iv3,i_f,j_f+1,k_f+2,idest)/dz  &
                            -  facevary1(iv2,i_f,j_f+2,k_f+1,idest)/dy
                   
            b3(8,1) =  DivU - (facevarx1(iv1,i_f+2,j_f+1,k_f+1,idest)- &
                              facevarx1(iv1,i_f+1,j_f+1,k_f+1,idest))/dx &
                            - (facevary1(iv2,i_f+1,j_f+2,k_f+1,idest)-   &
                              facevary1(iv2,i_f+1,j_f+1,k_f+1,idest))/dy &
                            -  facevarz1(iv3,i_f+1,j_f+1,k_f+2,idest)/dz
            b3(1:8,1) = dz*b3(1:8,1)


            if (k .eq. kcmin) then
!!$c$$$            ! Wi,j+1,k+1
!!$c$$$            b3(9,1) = 0.5*(facevarz1(iv3,i_f,j_f+1,k_f+2,idest)+
!!$c$$$&                          facevarz1(iv3,i_f,j_f+1,k_f,idest))

            ! Wi,j+1,k+1
            b3(9,1) = &
             0.375*facevarz1(iv3,i_f,j_f+1,k_f,idest)   + &
             0.75 *facevarz1(iv3,i_f,j_f+1,k_f+2,idest) - &
             0.125*facevarz1(iv3,i_f,j_f+1,k_f+4,idest)

            else

            ! Wi,j+1,k+1
            b3(9,1) = &
            -0.125*facevarz1(iv3,i_f,j_f+1,k_f-2,idest) + &
             0.75 *facevarz1(iv3,i_f,j_f+1,k_f,idest)   + &
             0.375*facevarz1(iv3,i_f,j_f+1,k_f+2,idest) 

            endif

       
            ! Solve for the rest of internal velocities using divergence equations:
            b3i= MATMUL(A3T,b3)
            usol3 = MATMUL(InvA3,b3i)
        
            ! U values
            facevarx1(iv1,i_f+1,j_f,k_f+1,idest) =   usol3(1,1)
            facevarx1(iv1,i_f+1,j_f+1,k_f,idest) =   usol3(2,1)
            ! V values
            facevary1(iv2,i_f+1,j_f+1,k_f,idest) =   usol3(3,1)
            facevary1(iv2,i_f,j_f+1,k_f+1,idest) =   usol3(4,1)
            ! W values
            facevarz1(iv3,i_f,j_f,k_f+1,idest)   =   usol3(5,1)
            facevarz1(iv3,i_f+1,j_f,k_f+1,idest) =   usol3(6,1)
            facevarz1(iv3,i_f,j_f+1,k_f+1,idest) =   usol3(7,1)
            facevarz1(iv3,i_f+1,j_f+1,k_f+1,idest) = usol3(8,1)

              enddo
           enddo
         enddo

!!$c$$$         DivU = 0.
!!$c$$$         do k = kcmin,kcmax
!!$c$$$           do j = jcmin,jcmax
!!$c$$$              do i = icmin,icmax
!!$c$$$
!!$c$$$                 U1   = (recvu(iv1,i+1,j,k)-recvu(iv1,i,j,k))/dx
!!$c$$$&                     + (recvv(iv2,i,j+1,k)-recvv(iv2,i,j,k))/dy
!!$c$$$&                     + (recvw(iv3,i,j,k+1)-recvw(iv3,i,j,k))/dz  
!!$c$$$
!!$c$$$                 DivU = max(DivU,abs(U1))
!!$c$$$
!!$c$$$
!!$c$$$              enddo
!!$c$$$           enddo
!!$c$$$        enddo
!!$c$$$
!!$c$$$        write(*,*) 'lb=',lb,', Div U parent=',DivU
!!$c$$$
!!$c$$$
!!$c$$$
!!$c$$$      ! Check final divergence
!!$c$$$      DivU = 0.
!!$c$$$      do k = kfmin,kfmax
!!$c$$$         do j = jfmin,jfmax
!!$c$$$            do i = ifmin,ifmax
!!$c$$$
!!$c$$$
!!$c$$$            ! Fine grid values of velocities on faces of parent Cell
!!$c$$$            u_1p = 
!!$c$$$&           facevarx1(iv1,i+1,j,k,idest)
!!$c$$$
!!$c$$$            u_1 = 
!!$c$$$&           facevarx1(iv1,i,j,k,idest)
!!$c$$$
!!$c$$$            v_1p =
!!$c$$$&           facevary1(iv2,i,j+1,k,idest)
!!$c$$$
!!$c$$$            v_1 =
!!$c$$$&           facevary1(iv2,i,j,k,idest)
!!$c$$$
!!$c$$$            W1p = facevarz1(iv3,i,j,k+1,idest)
!!$c$$$            W1  = facevarz1(iv3,i,j,k,idest)
!!$c$$$
!!$c$$$            U1  = (u_1p-u_1)/dx + 
!!$c$$$&                 (v_1p-v_1)/dy +
!!$c$$$&                 ( W1p- W1)/dz
!!$c$$$
!!$c$$$            DivU = max(DivU,abs(2.*U1))
!!$c$$$
!!$c$$$           enddo
!!$c$$$        enddo
!!$c$$$      enddo
!!$c$$$
!!$c$$$        write(*,*) 'lb=',lb,', Div U Child=',DivU
!!$c$$$
!!$c$$$!        call mpi_finalize()

      endif

      return
      end subroutine amr_1blk_fc_prol_divpres


      subroutine prol_fc_divpres_init(n,tf,i_divf_fc_vars)
      
        use paramesh_dimensions, only: ndim
        use physicaldata, only : interp_mask_facex, &
        interp_mask_facey,interp_mask_facez

      implicit none

      integer, intent(in) :: n, tf, i_divf_fc_vars(tf,n)
      integer i,iface

      prol_fc_divpres_n = n 

      if (.not. allocated(prol_fc_divpres_ivar)) &
      allocate(prol_fc_divpres_ivar(tf,prol_fc_divpres_n))

      do i = 1,prol_fc_divpres_n
        do iface=1,ndim
         prol_fc_divpres_ivar(iface,i) = i_divf_fc_vars(iface,i)
        end do

        interp_mask_facex(i_divf_fc_vars(1,i)) = -200 ! Corresponds to fc_divpres.
        if(ndim .gt. 1) interp_mask_facey(i_divf_fc_vars(2,i)) = -200 ! Corresponds to fc_divpres.
        if(ndim .eq. 3) interp_mask_facez(i_divf_fc_vars(3,i)) = -200 ! Corresponds to fc_divpres.
      end do

      prol_fc_divpres = .true.

      end subroutine prol_fc_divpres_init


end Module gr_pmDivpres_mod

  subroutine gr_pmDivpresApply(recvfx,recvfy,recvfz, &
       ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
       lb)
    use paramesh_dimensions, ONLY: il_bnd1,jl_bnd1,kl_bnd1
    use gr_pmDivpres_mod, ONLY: &
         prol_fc_divpres_n, prol_fc_divpres_ivar
    use gr_pmDivpres_mod, ONLY: amr_1blk_fc_prol_divpres
    use paramesh_dimensions, ONLY : &
         nfacevar, &
         iu_bnd1, ju_bnd1, ku_bnd1, &
         ndim, k2d, k3d
    use physicaldata, ONLY: facevarx1, facevary1, facevarz1

    implicit none
    Real,intent(inout) :: recvfx(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfy(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    Real,intent(inout) :: recvfz(:, il_bnd1:, jl_bnd1:,         kl_bnd1:)
    integer, intent(in) :: ia,ib,ja,jb,ka,kb
    integer, intent(in) :: idest,ioff,joff,koff
    integer, intent(in) :: lb

    Integer :: iprol, iv1, iv2, iv3

    recvfx(1:nfacevar, il_bnd1:iu_bnd1+1,          &
         jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)          &
         = facevarx1(1:nfacevar, il_bnd1:iu_bnd1+1,  &
         jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1, 1)


    Divpres_ndim:if (ndim .eq. 2) then

       recvfy(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, &
            kl_bnd1:ku_bnd1)                                 &
            = facevary1(1:nfacevar, il_bnd1:iu_bnd1,           &
            jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)

       recvfz = 0.

       do iprol = 1, prol_fc_divpres_n
          iv1 = prol_fc_divpres_ivar(1,iprol)
          iv2 = prol_fc_divpres_ivar(2,iprol)
          iv3 = 0

!!$               write(*,*) 'iv1,iv2,divp=',iv1,iv2

          call amr_1blk_fc_prol_divpres                   &
               (recvfx,recvfy,recvfz,ia,ib,ja,jb,ka,kb,idest,  &
               ioff,joff,koff,iv1,iv2,iv3,lb)

       enddo

    elseif(ndim .eq. 3) then

       recvfy(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d, &
            kl_bnd1:ku_bnd1)                                 &
            = facevary1(1:nfacevar, il_bnd1:iu_bnd1,             &
            jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)

       recvfz(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,           &
            kl_bnd1:ku_bnd1+k3d)                                   &
            = facevarz1(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1, &
            kl_bnd1:ku_bnd1+k3d, 1)

       do iprol = 1, prol_fc_divpres_n
          iv1 = prol_fc_divpres_ivar(1,iprol)
          iv2 = prol_fc_divpres_ivar(2,iprol)
          iv3 = prol_fc_divpres_ivar(3,iprol)

          call amr_1blk_fc_prol_divpres                    &
               (recvfx,recvfy,recvfz,ia,ib,ja,jb,ka,kb,idest,   &
               ioff,joff,koff,iv1,iv2,iv3,lb)

       enddo


    endif Divpres_ndim

  end subroutine gr_pmDivpresApply
