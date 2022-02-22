!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
! Copyright (C) 2016 The University of Chicago
! Copyright (C) 2022 UChicago Argonne, LLC and contributors
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "Simulation.h"


      subroutine amr_1blk_cc_prol_work_mg                    &
       (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,         &
        lb,order)


!
!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted
! from the solution array work, and performs a prolongation operation
! on it, between the bounds ranges ia to ib, ja to jb, and ka to kb.
! The data in recv is from a parent block and the
! result of the prolongation operation is written directly into one
! layer of the working block array work1(...,idest).
! The position of the child within the parent block is specified by
! the ioff, joff and koff arguments.
!
! This particular prolongation uses a more general interpolation proceedure
! than some of the other routines provided with PARAMESH.  Any 'order' of
! interpolation can be selected for any variable (as described below).
! This routine works by using a general polynomial fitting algorithm.
! The interpolations are performed first in the `x'
! direction.  `Y' interapolations follow, but use the interpolated
! data from the `x' sweep.  The 'Z' sweep is similarly performed.
!
! To select the `order' (we use the term order here loosely) of interpolation
! the array interp_mask_work should have data in it that is >= 30.
! Since the interpolation scheme is general, one can select
! different orders of interpolation for different variables as.  For instance,
! if
!   interp_mask_work(1) = 30
!   interp_mask_work(2) = 31
!   interp_mask_work(3) = 32
! then this routine will be called by amr_1blk_cc_prol_gen_work_fun with
! 'order' set to 0, 1, and 2, respectively, with the effect that WORK
! variable 1 will be prolongated used simple direct injection, work variable 2
! will be prolongated using linear interpolation and work variable 3 will be
! prolongated using quadratic interpolation.
!
! Finally, the `order' of interpolation must be equal or less than nguard_work.
! This ensures that enough guardcells space is available to compute
! the interpolation weights for the polynomial fits.
!
! It is applied to all WORK variables whose corresponding element
! of interp_mask_work is in the range 30..32.  Usually there is at most one
! work variable, so what matters is usually the value of interp_mask_work(1).
!
! NOTE: This routine may not be as effcient as some of the other, similar
!       routine provided for prolongation. So, if you don't need the
!       flexibility of this routine, you might want to consider using another
!       or writing another yourself.
!
! NOTE2:  This routine does NOT guarantee conservative prologation at
!         refinement jumps.  This is described in the documentation.
!
! Written :     Kevin Olson,  March 2002 and based on similar routines
!               by Peter MacNeice.
!------------------------------------------------------------------------
! HISTORY   2010    modified to do direct injection for first cell in from a face,
!                   and named amr_1blk_cc_prol_work_mg            - Marcos Vanella
!           2016    moved into Grid/GridMain                      - Klaus Weide
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,lb
      integer, intent(in) :: order
      real,    intent(inout) :: recv(:,:,:)

!------------------------------------
! local arrays

      real :: weight_right
      real :: weight_left
      real :: tempx, tempy

      integer :: i,j,k
      integer :: offi,offj,offk
      integer,parameter :: largei = 100
      integer,parameter :: maxorder = 3
      integer :: ii,iorder
      integer :: icmin,jcmin,kcmin
      integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      integer :: ipar, jpar, kpar

      integer,save,allocatable :: imina(:,:), &
                                  imaxa(:,:)
      real,save,allocatable    :: weightx(:,:,:)

      integer,save,allocatable :: jmina(:,:), &
                                  jmaxa(:,:)
      real,save,allocatable    :: weighty(:,:,:)

      integer,save,allocatable :: kmina(:,:), &
                                  kmaxa(:,:)
      real,save,allocatable    :: weightz(:,:,:)
      
      logical,save :: first_call = .true.

!------------------------------------

      if (first_call) then
         !$omp critical (CritRegAmr_1blk_cc_prol_work_mg)
         first_call = .false.

         if (.NOT. allocated(imina)) then
            allocate (imina(iuw1,0:maxorder))
            allocate (imaxa(iuw1,0:maxorder))
            allocate (weightx(iuw1,iuw1,0:maxorder))

            allocate (jmina(juw1,0:maxorder))
            allocate (jmaxa(juw1,0:maxorder))
            allocate (weighty(juw1,juw1,0:maxorder))

            allocate (kmina(kuw1,0:maxorder))
            allocate (kmaxa(kuw1,0:maxorder))
            allocate (weightz(kuw1,kuw1,0:maxorder))
         endif

         do iorder = 0,maxorder

!!! XXXXX !!!

         i = ((1-nguard_work-1+largei)/2 + &
               nguard_work - largei/2 ) + 1

         !write(*,*) 'iuw1',iuw1,'i=',i

         do ii = 1,iuw1

            if ((mod(ii,2) .ne. 0 .and. mod(nguard_work,2) .ne. 0) .or. &
                (mod(ii,2)  ==  0 .and. mod(nguard_work,2)  ==  0)) then

                                ! right point

               if (ii < nguard + nxb/2) then
                  imina(ii,iorder) = i
                  imaxa(ii,iorder) = i + iorder
               else
                  imina(ii,iorder) = i + iorder/2 - iorder
                  imaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = &
                         weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightx(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1

            else
                                ! left point

               if (ii < nguard + nxb/2) then
                  imina(ii,iorder) = i - iorder/2
                  imaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  imina(ii,iorder) = i - iorder
                  imaxa(ii,iorder) = i
               end if

               do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = &
                         weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  !write(*,*)'ipar=',ipar,'iorder',iorder,'ii=',ii
                  !write(*,*) 'imina=',imina(ii,iorder),'imaxa=',imaxa(ii,iorder)
                  weightx(ipar,ii,iorder) = weight_left
               end do

            end if

         end do                 ! end loop over ii

!!! YYYYY !!!

         if (ndim >= 2) then

         i = ((1-nguard_work-1+largei)/2 + &
               nguard_work - largei/2 ) + 1
         do ii = 1,juw1

            if ((mod(ii,2) .ne. 0 .and. mod(nguard_work,2) .ne. 0) .or. &
                (mod(ii,2)  ==  0 .and. mod(nguard_work,2)  ==  0)) then

                                ! right point

               if (ii < nguard + nyb/2) then
                  jmina(ii,iorder) = i
                  jmaxa(ii,iorder) = i + iorder
               else
                  jmina(ii,iorder) = i + iorder/2 - iorder
                  jmaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = &
                        weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weighty(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1

            else
                                ! left point

               if (ii < nguard + nyb/2) then
                  jmina(ii,iorder) = i - iorder/2
                  jmaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  jmina(ii,iorder) = i - iorder
                  jmaxa(ii,iorder) = i
               end if

               do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = &
                        weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weighty(ipar,ii,iorder) = weight_left
               end do

            end if

         end do                 ! end loop over ii

         end if                 ! end if (ndim

!!! ZZZZZ !!!

         if (ndim == 3) then

         i = ((1-nguard_work-1+largei)/2 + &
               nguard_work - largei/2 ) + 1
         do ii = 1,kuw1

            if ((mod(ii,2) .ne. 0 .and. mod(nguard_work,2) .ne. 0) .or. &
                (mod(ii,2)  ==  0 .and. mod(nguard_work,2)  ==  0)) then

                                ! right point

               if (ii < nguard + nzb/2) then
                  kmina(ii,iorder) = i
                  kmaxa(ii,iorder) = i + iorder
               else
                  kmina(ii,iorder) = i + iorder/2 - iorder
                  kmaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = &
                        weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightz(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1

            else
                                ! left point

               if (ii < nguard + nzb/2) then
                  kmina(ii,iorder) = i - iorder/2
                  kmaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  kmina(ii,iorder) = i - iorder
                  kmaxa(ii,iorder) = i
               end if

               do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = &
                        weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightz(ipar,ii,iorder) = weight_left
               end do

            end if

         end do                 ! end loop over ii

         end if                 ! end if (ndim

         end do                 ! end loop over iorder
         !$omp end critical (CritRegAmr_1blk_cc_prol_work_mg)

      end if                    ! end if (first_call



! Set the bounds on the loop controlling the interpolation.
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

      kcmin = ((kfmin-nguard_work-1+largei)/2 + &
                     nguard_work - largei/2 )*k3d + &
                     1 + offk
      jcmin = ((jfmin-nguard_work-1+largei)/2 + &
                      nguard_work - largei/2 )*k2d + &
                      1 + offj
      icmin = ((ifmin-nguard_work-1+largei)/2 + &
                      nguard_work - largei/2 ) + &
                      1 + offi





! Main Interpolation loop.






! Interpolate !




      do k = kfmin,kfmax
      do j = jfmin,jfmax
      do i = ifmin,ifmax
         
         work1(i,j,k,idest) = 0.

         if (ndim == 3) then

! DO direct injection for first cell in from a face
               
            if ((neigh(1,1,lb) <= -20 .and. i == nguard_work+1) .or. &
                (neigh(1,3,lb) <= -20 .and. j == nguard_work+1) .or. &
                (neigh(1,5,lb) <= -20 .and. k == nguard_work+1) .or. &
                (neigh(1,2,lb) <= -20 .and. i == nguard_work+nxb) .or. &
                (neigh(1,4,lb) <= -20 .and. j == nguard_work+nyb) .or. &
                (neigh(1,6,lb) <= -20 .and. k == nguard_work+nzb))   &
                then

               kpar = ((k-nguard_work-1+largei)/2 + &
                        nguard_work - largei/2 )*k3d + & 
                        1 + offk
               jpar = ((j-nguard_work-1+largei)/2 + &
                        nguard_work - largei/2 )*k2d + & 
                        1 + offj
               ipar = ((i-nguard_work-1+largei)/2 + &
                        nguard_work - largei/2 ) +  &
                        1 + offi

               work1(i,j,k,idest) = recv(ipar,jpar,kpar)

            else

            do kpar = kmina(k,order),kmaxa(k,order)
            tempy = 0.
            do jpar = jmina(j,order),jmaxa(j,order)
            tempx = 0.
            do ipar = imina(i,order),imaxa(i,order)
               tempx = tempx + &
                    weightx(ipar,i,order)* &
                    recv(ipar+offi,jpar+offj,kpar+offk)
            end do
               tempy = tempy + &
                    weighty(jpar,j,order)*tempx
            end do
               work1(i,j,k,idest) = work1(i,j,k,idest) + &
                   weightz(kpar,k,order)*tempy
            end do

            end if

         elseif (ndim == 2) then

            kpar = 1
            do jpar = jmina(j,order),jmaxa(j,order)
            do ipar = imina(i,order),imaxa(i,order)
               work1(i,j,k,idest) = work1(i,j,k,idest) + &
                   weightx(ipar,i,order)* &
                   weighty(jpar,j,order)* &
                   recv(ipar+offi,jpar+offj,kpar+offk)
            end do
            end do

         elseif (ndim == 1) then

            kpar = 1
            jpar = 1
            do ipar = imina(i,order),imaxa(i,order)
               work1(i,j,k,idest) = work1(i,j,k,idest) + &
                   weightx(ipar,i,order)* &
                   recv(ipar+offi,jpar+offj,kpar+offk)
            end do

         end if

      end do                    ! end loop over i
      end do                    ! end loop over j
      end do                    ! end loop over k


 2    return
      end subroutine amr_1blk_cc_prol_work_mg

