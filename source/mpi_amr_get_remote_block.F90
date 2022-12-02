!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/mpi_amr_get_remote_block
!! NAME
!!
!!   mpi_amr_get_remote_block
!!
!! SYNOPSIS
!!
!!   Call mpi_amr_get_remote_block (mype, remote_pe, remote_block, idest,
!!                                   iopt, lcc, lfc, lec, lnc,
!!                                   nlayersx,nlayersy,nlayersz )
!!
!!   Call mpi_amr_get_remote_block (integer, integer, integer, integer,
!!                                   integer, logical, logical, logical, logical,
!!                                   optional integer, optional integer, 
!!                                   optional integer)
!!
!! ARGUMENTS
!!
!!   integer :: mype             
!!     The local processor
!!
!!   integer :: remote_pe        
!!     The remote processor.
!!
!!   integer :: remote_block     
!!     The local block id of the block to be copied from
!!     the remote processor.
!!    
!!   integer :: idest            
!!     Selects the storage space in the 1blk data structures which is to
!!     be used in this Call. If the leaf node is having its
!!     guardcells filled then set this to 1, if its parent
!!     is being filled set it to 2.
!!
!!   integer :: iopt             
!!     A switch to control which data source is to be used
!!     iopt=1 will use 'unk', 'facevarx,y,z', 'unk_e_x,y,z', and 'unk_n'
!!     iopt>=2 will use 'work'
!!
!!   logical :: lcc              
!!     A logical switch which requests cell centered data.
!!
!!   logical :: lfc              
!!     A logical switch which requests cell-face centered data.
!!
!!   logical :: lec              
!!     A logical switch which requests cell-edge centered data.
!!
!!   logical :: lnc              
!!     A logical switch which requests cell-corner data.
!!
!!   optional, interger :: nlayersx, nlayersy, nlayersz
!!     Optional argumentes which indicate the number of guardcells to use
!!     in the x, y, and z directions.
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!   tree
!!   workspace
!!   mpi_morton
!!   paramesh_interfaces
!!   paramesh_mpi_interfaces
!!
!! CALLS
!!
!!   amr_mpi_find_blk_in_buffer
!!   mpi_set_message_limits
!!
!! RETURNS
!!
!!   Nothing returned. 
!!
!! DESCRIPTION
!! 
!!   This routine copies guard cell information to face iface in layer
!!   idest of the working block, from the appropriate face of the neighboring 
!!   block, assuming that the neighboring block is on a different processor.
!! 
!! AUTHORS
!!
!!  Written by Peter MacNeice (July 1997).  Modified by Kevin Olson for
!!  directional guardcell filling.
!!
!! MODIFICATIONS
!!  2022-11-08 K. Weide  Subroutine now takes argument 'pdg';
!!                       added ONLY to some USE statements;
!!                       include "Flashx_mpi_implicitNone.fh";
!!                       tweaked / added some debugging output.
!!  2022-12-01 K. Weide  Removed TAB characters
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine mpi_amr_get_remote_block(mype,remote_pe,remote_block, & 
          idest,iopt,lcc,lfc,lec,lnc, pdg,ig,                          &
          nlayersx,nlayersy,nlayersz)

!-----Use statements
      use gr_pmPdgDecl, ONLY : pdg_t
      Use paramesh_dimensions, ONLY: gr_thePdgDimens
      Use paramesh_dimensions, ONLY: ndim, k2d, k3d, nbndvar, nbndvare, nvarcorn, &
                                     npgs, nguard_work
      Use physicaldata, only: gt_unk
      Use physicaldata, only: facevarx,   facevary,   facevarz, &
                              gt_facevarx,gt_facevary,gt_facevarz, &
                              facevarx1,  facevary1,  facevarz1
      Use physicaldata, only: unk_e_x,    unk_e_y,    unk_e_z, &
                              gt_unk_e_x, gt_unk_e_y, gt_unk_e_z, &
                              unk_e_x1,  unk_e_y1,  unk_e_z1
      Use physicaldata, only: unk_n, gt_unk_n, unk_n1
      Use physicaldata, only: int_gcell_on_cc,int_gcell_on_fc,int_gcell_on_ec,int_gcell_on_nc
      Use physicaldata, only: ngcell_on_cc,       ngcell_on_fc,   ngcell_on_ec,   ngcell_on_nc
      Use physicaldata, only: gcell_on_cc_pointer,gcell_on_fc_pointer,gcell_on_ec_pointer,gcell_on_nc_pointer
      Use tree
      Use workspace
      Use mpi_morton
      Use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      Use paramesh_mpi_interfaces, only : mpiSet_message_limits
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

!-----Include statements.
#include "Flashx_mpi_implicitNone.fh"

!-----Input/Output arguments.
      Integer, intent(in) :: mype,remote_pe,remote_block
      Integer, intent(in) :: idest,iopt
      Logical, intent(in) :: lcc,lfc,lec,lnc
      type(pdg_t), intent(INOUT) :: pdg
      integer, intent(in) :: ig
      Integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!-----Local arrays and variables.
      Integer :: iopt0
      Integer :: nguard0 
      Integer :: nguard_work0
      Integer :: ierrorcode,ierr
      Integer :: dtype
      Integer :: vtype
      Integer :: index,index0
      Integer :: ia, ib, ja, jb, ka, kb, i, j, k
      Integer :: ivar, ivar_next
      Logical :: lfound

#ifndef FL_NON_PERMANENT_GUARDCELLS
# define no_permanent_guardcells .FALSE.      
#endif

!-----Begin executable code.
  ASSOCIATE(nguard       => gr_thePdgDimens(ig) % nguard, &
            nxb          => gr_thePdgDimens(ig) % nxb,    &
            nyb          => gr_thePdgDimens(ig) % nyb,    &
            nzb          => gr_thePdgDimens(ig) % nzb,    &
            nvar         => gr_thePdgDimens(ig) % nvar,   &
            unk          => pdg % unk,                    &
            unk1         => pdg % unk1)

      nguard0 = nguard*npgs
      nguard_work0 = nguard_work*npgs

      If (remote_block <= lnblocks.And.remote_pe == mype) Then

      If (iopt == 1) Then

      If (lcc) Then

        Do ivar=1,nvar
          If (int_gcell_on_cc(ivar)) Then
!---------Copy complete remote block into a buffer block Called recv.
          If (no_permanent_guardcells) Then
          unk1(ivar,                                                   & 
               1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,        & 
               1+nguard*k3d:nzb+nguard*k3d,idest) =                    & 
            gt_unk(ivar,                                               & 
               1+nguard0:nxb+nguard0,                                  & 
               1+nguard0*k2d:nyb+nguard0*k2d,                          & 
               1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
          Else ! no_permanent_guardcells
          unk1(ivar,nguard+1:nguard+nxb,                               & 
                    nguard*k2d+1:nguard*k2d+nyb,                       & 
                    nguard*k3d+1:nguard*k3d+nzb,idest) =               & 
           unk(ivar,nguard0+1:nguard0+nxb,                             & 
                    nguard0*k2d+1:nguard0*k2d+nyb,                     & 
                    nguard0*k3d+1:nguard0*k3d+nzb,remote_block)
          End If
          End If
        End Do  ! End Do ivar=1,nvar

      End If  ! End If (lcc)

      If (lfc) Then

       Do ivar=1,nbndvar
       If (no_permanent_guardcells) Then
       If (int_gcell_on_fc(1,ivar)) Then
       facevarx1(ivar, & 
                 1+nguard:nxb+nguard+1,1+nguard*k2d:nyb+nguard*k2d,    & 
                 1+nguard*k3d:nzb+nguard*k3d,idest) =                  & 
        gt_facevarx(ivar,                                              & 
                    1+nguard0:nxb+nguard0+1,                           & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       If (ndim >= 2) Then
       If (int_gcell_on_fc(2,ivar)) Then
       facevary1(ivar,                                                 & 
                 1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d+k2d,  & 
                 1+nguard*k3d:nzb+nguard*k3d,idest) =                  & 
        gt_facevary(ivar,                                              & 
                    1+nguard0:nxb+nguard0,                             & 
                    1+nguard0*k2d:nyb+nguard0*k2d+k2d,                 & 
                    1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       End If
       If (ndim == 3) Then
       If (int_gcell_on_fc(3,ivar)) Then
       facevarz1(ivar,                                                 & 
                 1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,      & 
                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =              & 
        gt_facevarz(ivar,                                              & 
                    1+nguard0:nxb+nguard0,                             & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       End If

       Else ! no_permanent_guardcells

       If (int_gcell_on_fc(1,ivar)) Then
       facevarx1(ivar,                                                 & 
                 1+nguard:nxb+nguard+1,                                & 
                 1+nguard*k2d:nyb+nguard*k2d,                          & 
                 1+nguard*k3d:nzb+nguard*k3d,idest) =                  & 
           facevarx(ivar,                                              & 
                    1+nguard0:nxb+nguard0+1,                           & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       If (ndim >= 2) Then
       If (int_gcell_on_fc(2,ivar)) Then
       facevary1(ivar,                                                 & 
                 1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d+k2d,  & 
                 1+nguard*k3d:nzb+nguard*k3d,idest) =                  & 
           facevary(ivar,                                              & 
                    1+nguard0:nxb+nguard0,                             & 
                    1+nguard0*k2d:nyb+nguard0*k2d+k2d,                 & 
                    1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       End If
       If (ndim == 3) Then
       If (int_gcell_on_fc(3,ivar)) Then
       facevarz1(ivar, & 
                 1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,      & 
                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =              & 
           facevarz(ivar,                                              & 
                    1+nguard0:nxb+nguard0,                             & 
                    1+nguard0*k2d:nyb+nguard0*k2d,                     & 
                    1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       End If

       End If  ! End If (no_permanent_guardcells)
       End Do  ! End Do ivar=1,nbndvar

      End If  ! End If (lfc)


      If (ndim > 1) Then
      If (lec) Then
!-----Copy complete remote block into a buffer block Called recv.
       Do ivar=1,nbndvare
       If (no_permanent_guardcells) Then
       If (int_gcell_on_ec(1,ivar)) Then
       unk_e_x1(ivar,                                                  & 
                1+nguard:nxb+nguard,                                   & 
                1+nguard*k2d:nyb+nguard*k2d+k2d,                       & 
                1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =               & 
        gt_unk_e_x(ivar,                                               & 
                   1+nguard0:nxb+nguard0,                              & 
                   1+nguard0*k2d:nyb+nguard0*k2d+k2d,                  & 
                   1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       If (int_gcell_on_ec(2,ivar)) Then
       unk_e_y1(ivar,                                                  & 
                1+nguard:nxb+nguard+1,                                 & 
                1+nguard*k2d:nyb+nguard*k2d,                           & 
                1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =               & 
        gt_unk_e_y(ivar,                                               & 
                   1+nguard0:nxb+nguard0+1,                            & 
                   1+nguard0*k2d:nyb+nguard0*k2d,                      & 
                   1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       If (ndim == 3) Then
       If (int_gcell_on_ec(3,ivar)) Then
       unk_e_z1(ivar,                                                  & 
                1+nguard:nxb+nguard+1,                                 & 
                1+nguard*k2d:nyb+nguard*k2d+k2d,                       & 
                1+nguard*k3d:nzb+nguard*k3d,idest) =                   & 
        gt_unk_e_z(ivar,                                               & 
                   1+nguard0:nxb+nguard0+1,                            & 
                   1+nguard0*k2d:nyb+nguard0*k2d+k2d,                  & 
                   1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       End If ! If (ndim == 3)

       Else ! no_permanent_guardcells

       If (int_gcell_on_ec(1,ivar)) Then
       unk_e_x1(ivar,                                                  & 
                1+nguard:nxb+nguard,                                   & 
                1+nguard*k2d:nyb+nguard*k2d+k2d,                       & 
                1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =               & 
           unk_e_x(ivar, & 
                   1+nguard0:nxb+nguard0, & 
                   1+nguard0*k2d:nyb+nguard0*k2d+k2d, & 
                   1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       If (int_gcell_on_ec(2,ivar)) Then
       unk_e_y1(ivar,                                                  & 
                1+nguard:nxb+nguard+1,                                 & 
                1+nguard*k2d:nyb+nguard*k2d,                           & 
                1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =               & 
           unk_e_y(ivar,                                               & 
                   1+nguard0:nxb+nguard0+1,                            & 
                   1+nguard0*k2d:nyb+nguard0*k2d,                      & 
                   1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If
       If (ndim == 3) Then
       If (int_gcell_on_ec(3,ivar)) Then
       unk_e_z1(ivar,                                                  & 
                1+nguard:nxb+nguard+1,                                 & 
                1+nguard*k2d:nyb+nguard*k2d+k2d,                       & 
                1+nguard*k3d:nzb+nguard*k3d,idest) =                   & 
           unk_e_z(ivar,                                               & 
                   1+nguard0:nxb+nguard0+1,                            & 
                   1+nguard0*k2d:nyb+nguard0*k2d+k2d,                  & 
                   1+nguard0*k3d:nzb+nguard0*k3d,remote_block)
       End If
       End If  ! If (ndim == 3)

       End If  ! End If (no_permanent_guardcells)
       End Do  ! End Do ivar=1,nbndvare

      End If  ! End If (lec)                             
      End If  ! End If (ndim > 1)

      If (lnc) Then
!------Copy complete remote block into a buffer block Called recv.
       Do ivar=1,nvarcorn
       If (int_gcell_on_nc(ivar)) Then
       If (no_permanent_guardcells) Then
       unk_n1(ivar,                                                    & 
              1+nguard:nxb+nguard+1,                                   & 
              1+nguard*k2d:nyb+nguard*k2d+k2d,                         & 
              1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =                 & 
        gt_unk_n(ivar,                                                 & 
                 1+nguard0:nxb+nguard0+1,                              & 
                 1+nguard0*k2d:nyb+nguard0*k2d+k2d,                    & 
                 1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)

       Else ! no_permanent_guardcells

       unk_n1(ivar,                                                    & 
              1+nguard:nxb+nguard+1,                                   & 
              1+nguard*k2d:nyb+nguard*k2d+k2d,                         & 
              1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =                 & 
           unk_n(ivar,                                                 & 
                 1+nguard0:nxb+nguard0+1,                              & 
                 1+nguard0*k2d:nyb+nguard0*k2d+k2d,                    & 
                 1+nguard0*k3d:nzb+nguard0*k3d+k3d,remote_block)
       End If  ! End If (no_permanent_guardcells)
       End If  ! End If (int_gcell_on_nc(ivar))
       End Do  ! End Do ivar=1,nvarcorn
      End If  ! End If (lnc)

      ElseIf (iopt >= 2) Then

      If (lcc) Then

        iopt0 = iopt-1
!-------Copy complete remote block into a buffer block Called recvw.
        work1(1+nguard_work:nxb+nguard_work,                           & 
              1+nguard_work*k2d:nyb+nguard_work*k2d,                   & 
              1+nguard_work*k3d:nzb+nguard_work*k3d,idest) =           & 
          work(                                                        & 
               1+nguard_work0:nxb+nguard_work0,                        & 
               1+nguard_work0*k2d:nyb+nguard_work0*k2d,                & 
               1+nguard_work0*k3d:nzb+nguard_work0*k3d,                & 
               remote_block,iopt0)

      End If  ! End If (lcc)

      End If  ! End If (iopt == 1)

#ifdef DEBUG
399   format(1x,'@ ',I0,'  LOCAL: ',I3,'@',I0,' idest ',I0,' -- lcc ',&
                 L1,' copied to unk1 per var: ',I0,' Reals.')
      print 399,mype,remote_block,remote_pe,idest,lcc,&
           size(unk1(ivar,nguard+1:nguard+nxb,                               &
                    nguard*k2d+1:nguard*k2d+nyb,                       &
                    nguard*k3d+1:nguard*k3d+nzb,idest))
#endif

      Else  

      Call amr_mpi_find_blk_in_buffer(mype,remote_block,               & 
                                      remote_pe,idest,dtype,index0,    &
                                      lfound)

#ifdef DEBUG
      400 format('>@ ',I0,' remote: ',I3,'@',I0,' idest ',I0,' -- found ',&
                 L1,' dtype ',I0,' index0 ',I0)
      print 400,mype,remote_block,remote_pe,idest,lfound,dtype,index0
#endif
      If ((.Not.lfound).Or.(dtype.ne.14.And.dtype.ne.14+27)) Then
          write(*,*) 'Paramesh error : pe ',mype,' needed full blk ',  & 
            remote_block,remote_pe,' but could not find it or only ',  & 
            ' found part of it in the message buffer.',                & 
            '  Found',lfound,', dtype=',dtype,'.'
          Call mpi_abort(amr_mpi_meshComm,ierrorcode,ierr)
      End If

      If (lcc) Then

        If (iopt == 1) Then
          vtype = 1
          index = index0
          Call mpiSet_message_limits(                                 & 
                       dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                  & 
                       nlayersx,nlayersy,nlayersz)

          If (no_permanent_guardcells) Then
             ia = ia + nguard
             ib = ib + nguard
             ja = ja + nguard*k2d
             jb = jb + nguard*k2d
             ka = ka + nguard*k3d
             kb = kb + nguard*k3d
          End If

          Do k = ka,kb
          Do j = ja,jb
          Do i = ia,ib
            Do ivar=1,ngcell_on_cc
              ivar_next = gcell_on_cc_pointer(ivar)
              unk1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_cc
          End Do
          End Do
          End Do

        ElseIf (iopt > 1) Then

          vtype = 0
          index = index0
          Call mpiSet_message_limits(                                 & 
                       dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                  & 
                       nlayersx,nlayersy,nlayersz)

          If (no_permanent_guardcells) Then
             ia = ia + nguard_work
             ib = ib + nguard_work
             ja = ja + nguard_work*k2d
             jb = jb + nguard_work*k2d
             ka = ka + nguard_work*k3d
             kb = kb + nguard_work*k3d
          End If

          Do k = ka,kb
          Do j = ja,jb
          Do i = ia,ib
            work1(i,j,k,idest) =  temprecv_buf(index+1)
            index = index+1
          End Do
          End Do
          End Do

        End If  ! End If (iopt == 1)
      
      End If  ! End If (lcc)             

      If (lfc) Then

        if (iopt == 1) then

!-------starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                             index + ngcell_on_cc*message_size_cc(dtype)

        vtype = 2
        Call mpiSet_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)

        If (no_permanent_guardcells) Then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
        End If

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
            Do ivar=1,ngcell_on_fc(1)
             ivar_next = gcell_on_fc_pointer(1,ivar)
             facevarx1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_fc(1)
        End Do
        End Do
        End Do

        If (ndim >= 2) Then
         vtype = 3
         Call mpiSet_message_limits(                                  & 
                    dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                     & 
                    nlayersx,nlayersy,nlayersz)
 
         If (no_permanent_guardcells) Then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
         End If

         Do k = ka,kb
         Do j = ja,jb
         Do i = ia,ib
            Do ivar=1,ngcell_on_fc(2)
             ivar_next = gcell_on_fc_pointer(2,ivar)
             facevary1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_fc(2)
         End Do
         End Do
         End Do
        End If ! End If (ndim >= 2)

        If (ndim == 3) Then
         vtype = 4
         Call mpiSet_message_limits(                                  & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)

         If (no_permanent_guardcells) Then
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
         End If

         Do k = ka,kb
         Do j = ja,jb
         Do i = ia,ib
            Do ivar=1,ngcell_on_fc(3)
             ivar_next = gcell_on_fc_pointer(3,ivar)
             facevarz1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_fc(3)
         End Do
         End Do
         End Do
        End If  ! End If (ndim == 3)

        End If  ! End If (iopt == 1)

      End If  ! End If (lfc)

      If (ndim > 1) Then
      If (lec) Then

        If (iopt == 1) Then

!-------starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                             index + ngcell_on_cc*message_size_cc(dtype)
        If(l_datapacked(3)) index =                                    &
                            index + ngcell_on_fc(1) *                  &
                                    message_size_fcx(dtype)            &
                                  + ngcell_on_fc(2) *                  &
                                    message_size_fcy(dtype)            &
                                  + ngcell_on_fc(3) *                  &
                                    message_size_fcz(dtype)   
        vtype = 5
        Call mpiSet_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)

        If (no_permanent_guardcells) Then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        End If

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
            Do ivar=1,ngcell_on_ec(1)
              ivar_next = gcell_on_ec_pointer(1,ivar)
              unk_e_x1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_ec(1)
        End Do
        End Do
        End Do

        If (ndim >= 2) Then
        vtype = 6
        Call mpiSet_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)

        If (no_permanent_guardcells) Then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
          kb = kb + nguard*k3d
        End If

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
            Do ivar=1,ngcell_on_ec(2)
              ivar_next = gcell_on_ec_pointer(2,ivar)
              unk_e_y1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_ec(2)
        End Do
        End Do
        End Do
        End If  ! End If (ndim >= 2)

        If (ndim == 3) Then
        vtype =7
        Call mpiSet_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)

        If (no_permanent_guardcells) Then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        End If

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
            Do ivar=1,ngcell_on_ec(3)
              ivar_next = gcell_on_ec_pointer(3,ivar)
              unk_e_z1(ivar_next,i,j,k,idest) =                        & 
                  temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_ec(3)
        End Do
        End Do
        End Do
       End If  ! End If (ndim == 3)

       End If  ! End If (iopt ==  1)

      End If  ! End If (lec)
      End If  ! End If (ndim > 1)

      If (lnc) Then

        If (iopt == 1) Then

! starting index if cell-centered data is also included in recv_buf
        index = index0
        If (l_datapacked(2)) index =                                   & 
                            index + ngcell_on_cc*message_size_cc(dtype)
        If(l_datapacked(3)) index =                                    &
                            index + ngcell_on_fc(1) *                  &
                                    message_size_fcx(dtype)            &
                                  + ngcell_on_fc(2) *                  &
                                    message_size_fcy(dtype)            &
                                  + ngcell_on_fc(3) *                  &
                                    message_size_fcz(dtype)
           If (l_datapacked(4)) index =                                & 
                            index + maxval(ngcell_on_ec(1:ndim))       & 
                                     *message_size_ec(dtype)

        vtype = 8
        Call mpiSet_message_limits(                                   & 
                     dtype,ia,ib,ja,jb,ka,kb,vtype,ig,                    & 
                     nlayersx,nlayersy,nlayersz)


        If (no_permanent_guardcells) Then
           ia = ia + nguard
           ib = ib + nguard
           ja = ja + nguard*k2d
           jb = jb + nguard*k2d
           ka = ka + nguard*k3d
           kb = kb + nguard*k3d
        End If

        Do k = ka,kb
        Do j = ja,jb
        Do i = ia,ib
            Do ivar=1,ngcell_on_nc
              ivar_next = gcell_on_nc_pointer(ivar)
              unk_n1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            End Do
            index = index+ngcell_on_nc
        End Do
        End Do
        End Do

        End If  ! End If (iopt == 1)

      End If  ! End If (lnc)

      End If  ! End If (int_gcell_on_nc(ivar))
  end ASSOCIATE
      Return
      End Subroutine mpi_amr_get_remote_block
