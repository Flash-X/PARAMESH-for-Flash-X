!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003, 2004 United States Government as represented by the
! National Aeronautics and Space Administration, Goddard Space Flight
! Center.  All Rights Reserved.
! Copyright (C) 2019 The University of Chicago
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------


!     Klaus Weide,  September 2019, *clean_divb*
!        Dummy version for this unused interface, for FLASH5
!
       subroutine prol_fc_clean_divb_init(n,i_divf_fc_vars)
       use prolong_arrays, only : prol_fc_clean_divb
       implicit none
       integer, intent(in) :: n, i_divf_fc_vars(3,n)
       prol_fc_clean_divb = .FALSE.
       end subroutine prol_fc_clean_divb_init

       subroutine amr_1blk_fc_clean_divb(  & 
     &        nfacevar_in,  & 
     &        ia,ib,ja,jb,ka,kb,    & 
     &        ionea,ioneb,  & 
     &        jonea,joneb,  & 
     &        konea,koneb,  & 
     &        idest,ioff,joff,koff,          & 
     &        mype,lb,parent_pe,parent_blk   & 
     & )

       implicit none

       integer, intent(in) :: nfacevar_in
       integer, intent(in) :: ia,ib,ja,jb,ka,kb
       integer, intent(in) :: ionea,ioneb
       integer, intent(in) :: jonea,joneb
       integer, intent(in) :: konea,koneb
       integer, intent(in) :: idest, ioff, joff, koff
       integer, intent(in) :: mype, lb, parent_pe, parent_blk

       end subroutine amr_1blk_fc_clean_divb
