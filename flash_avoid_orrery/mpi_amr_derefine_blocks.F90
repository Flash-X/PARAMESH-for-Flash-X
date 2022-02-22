!!****if* source/mpi_amr_derefine_blocks
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
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
!!   mpi_amr_derefine_blocks
!! 
!! SYNOPSIS
!!
!!   Call mpi_amr_derefine_blocks (lnblocks_old, mype)
!!
!!   Call mpi_amr_derefine_blocks (integer, integer)
!! 
!! ARGUMENTS      
!!
!!    integer, intent(in) :: lnblocks_old, mype
!!      lnbllocks_old -> The number of local blocks before derefinement.
!!      mype   -> The local processor id.
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
!!   constants
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine handles derefinement to the mesh.  It removes blocks 
!!   whic are marked for derinement.
!!
!! AUTHORS
!!
!!   Kevin Olson, 1997
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      Subroutine amr_derefine_blocks(lnblocks_old, mype)

!-----Use statements
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use constants
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

      Implicit None

!-----Input/Output arguments.
      Integer, intent(inout) :: lnblocks_old
      Integer, intent(in)    :: mype

!-----Include statements.
      Include 'mpif.h'

!-----Local arrays and variables
      Real    :: eps,accuracy
      Integer :: new_loc(maxblocks_tr)
      Integer :: i,j,k,jsend
      Integer :: lnblocks2
      Integer :: neight(2,mfaces,maxblocks_tr)
      Integer :: childt(2,mchild,maxblocks_tr)
      Integer :: parentt(2,maxblocks_tr)
      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      Integer :: reqr(maxblocks_tr)
      Integer :: ierr,nsend,nrecv
      Integer :: nodetype_chi(nchild,maxblocks_tr)
      Integer :: jr0
      Logical :: lsend
      Logical :: useFlashCustomVersion

!-----Begin Executable code.
      useFlashCustomVersion = .false.
      if (use_flash_surr_blks_fill) then
         !The custom subroutine may only be used if there is valid
         !information in surr_blks data structure.  The standard 
         !PARAMESH subroutine works with (and updates) neigh data
         !structure, so we use this until the "orrery" has 
         !transferred neigh data into surr_blks.  (The surr_blks
         !flag is actually flipped by FLASH because FLASH
         !initialization is not well integrated with PARAMESH.)
         useFlashCustomVersion = surr_blks_valid
      end if

      If (useFlashCustomVersion) then

         call amr_derefine_blocks_flash (lnblocks_old, mype)

      Else

      accuracy = 10./10.**precision(accuracy)
      eps = accuracy

!-----remove blocks marked for derefinement by packing the data
      Do i = 1,maxblocks_tr
         new_loc(i) = -1
      End Do 

!-----Compute new_loc, new_loc marks where each block will end up after the
!-----derefinement is Do ne
      k = 1
      Do i = 1,lnblocks
         If (.not.derefine(i)) Then
            new_loc(i) = k
            k = k + 1
          End If
      End Do 

!-----reconnect all pointers
      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nrecv = nrecv + 1
             Call MPI_IRECV(parentt(1,i),1,MPI_INTEGER,                & 
                  parent(2,i),i,amr_mpi_meshComm,                        & 
                  reqr(nrecv),ierr)
           Else
             parentt(1,i) = new_loc(parent(1,i))
           End If
         End If
       End Do 
       
       nsend = 0
       Do i = 1,lnblocks
         Do j = 1,nchild
           If (child(1,j,i) > 0) Then
             If (child(2,j,i).ne.mype) Then
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,               & 
                    child(2,j,i),child(1,j,i),amr_mpi_meshComm,          & 
                    ierr)
             End If
           End If
         End Do 
       End Do 

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      nrecv = 0
      Do i = 1,lnblocks
        Do j = 1,nchild
          If (child(1,j,i) > 0) Then
            If (child(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(childt(1,j,i),1,MPI_INTEGER,              & 
                   child(2,j,i),child(1,j,i),amr_mpi_meshComm,           & 
                   reqr(nrecv),ierr)
            Else
              childt(1,j,i) = new_loc(child(1,j,i))
            End If
          End If
        End Do 
       End Do 
       
       nsend = 0
       Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nsend = nsend + 1
             Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,                 & 
                  parent(2,i),i,amr_mpi_meshComm,                        & 
                  ierr)
           End If
         End If
       End Do 

       If (nrecv > 0) Then
         Call MPI_WAITALL(nrecv,reqr,statr,ierr)
       End If

       Do j = 1,nfaces

         If (mod(j,2) == 0) Then
            jsend = j - 1
         Else
            jsend = j + 1
         End If
            
         nrecv = 0
         Do i = 1,lnblocks
            If (neigh(1,j,i) > 0) Then
               If (neigh(2,j,i).ne.mype) Then
                  nrecv = nrecv + 1
                  Call MPI_IRECV(neight(1,j,i),1,MPI_INTEGER,          & 
                       neigh(2,j,i),neigh(1,j,i),amr_mpi_meshComm,       & 
                       reqr(nrecv),ierr)
               Else
                  neight(1,j,i) = new_loc(neigh(1,j,i))
               End If
            End If
         End Do 
      
         nsend = 0
         Do i = 1,lnblocks

         lsend=.True.

         If (spherical_pm) Then
          If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lsend=.False.
          If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lsend=.False.
         End If

         If (lsend) Then
           If (neigh(1,jsend,i) > 0) Then
             If (neigh(2,jsend,i).ne.mype) Then
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,               & 
                    neigh(2,jsend,i),i,amr_mpi_meshComm,                 & 
                    ierr)
             End If
           End If
         End If

         If (spherical_pm) Then
           lsend = .True.
           jr0 = jsend
           If (j == 3.and.abs(bnd_box(1,2,i)) < eps) jr0 = 3
           If (j == 4.and.abs(bnd_box(1,2,i)) < eps) lsend = .False.
           If (j == 4.and.abs(bnd_box(2,2,i)-pi) < eps) jr0 = 4
           If (j == 3.and.abs(bnd_box(2,2,i)-pi) < eps) lsend = .False.
           If (abs(bnd_box(1,2,i)) < eps.and.                          & 
               abs(bnd_box(2,2,i)-pi) < eps) Then
             write(*,*) 'both poles in blk ',i,abs(bnd_box(1,2,i)),    &  
               abs(bnd_box(2,2,i)-pi)
             lsend = .True.
           End If

           If (lsend.and.jr0 == j) Then
             If (neigh(1,jr0,i) > 0) Then
               If (neigh(2,jr0,i).ne.mype) Then
                  nsend = nsend + 1
                  Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,            & 
                       neigh(2,jr0,i),i,amr_mpi_meshComm,                & 
                       ierr)
               End If
             End If
            End If
         End If  ! End If (spherical_pm)

         End Do  ! End Do i = 1,lnblocks

         If (nrecv > 0) Then
            Call MPI_WAITALL(nrecv,reqr,statr,ierr)
         End If

      End Do  ! End Do j = 1,nfaces

      Do i = 1,lnblocks_old
        If (new_loc(i).ne.i.and.new_loc(i) > 0) Then
          If (nvar > 0) unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
          If (nfacevar > 0) Then
             facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
             facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
             facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
          End If
          If (nvaredge > 0) Then
             unk_e_x(:,:,:,:,new_loc(i)) = unk_e_x(:,:,:,:,i)
             unk_e_y(:,:,:,:,new_loc(i)) = unk_e_y(:,:,:,:,i)
             unk_e_z(:,:,:,:,new_loc(i)) = unk_e_z(:,:,:,:,i)
          End If
          If (nvarcorn > 0) unk_n(:,:,:,:,new_loc(i)) =                & 
                                             unk_n(:,:,:,:,i)
        End If
      End Do 

      parent(1,1:lnblocks) = parentt(1,1:lnblocks)
      child(1,:,1:lnblocks) = childt(1,:,1:lnblocks)
      neigh(1,:,1:lnblocks) = neight(1,:,1:lnblocks)

      k = 1
      lnblocks2 = lnblocks
      Do i = 1,lnblocks
         
         If (.not.derefine(i)) Then
            
            If (k.ne.i) Then
               Do j = 1,nchild
                  child(1,j,k) = child(1,j,i)
                  child(2,j,k) = child(2,j,i)
               End Do 
               parent(1,k) = parent(1,i)
               parent(2,k) = parent(2,i)
               Do j = 1,nfaces
                  neigh(1,j,k) = neigh(1,j,i)
                  neigh(2,j,k) = neigh(2,j,i)
               End Do 
               Do j = 1,ndim
                  coord(j,k) = coord(j,i)
                  bnd_box(1,j,k) = bnd_box(1,j,i)
                  bnd_box(2,j,k) = bnd_box(2,j,i)
               End Do 
               bsize(:,k) = bsize(:,i)
               newchild(k) = newchild(i)
               which_child(k) = which_child(i)
               lrefine(k) = lrefine(i)
               bflags(:,k) = bflags(:,i)
               work_block(k) = work_block(i)
               If (empty_cells) Then
                  empty(k) = empty(i)
               End If
               
            End If

            k = k + 1
            
         Else
            
            lnblocks2 = lnblocks2 - 1
            lnblocks_old = lnblocks_old - 1
            
         End If  ! End If (.not.derefine(i))
         
      End Do  ! End Do i = 1,lnblocks

!-----overwrite old locations
      Do i = lnblocks2+1,lnblocks
         
         derefine(i) = .FALSE.
         Do j = 1,nchild
            child(1,j,i) = -1
            child(2,j,i) = -1
         End Do 
         parent(1,i) = -1
         parent(2,i) = -1
         Do j = 1,nfaces
            neigh(1,j,i) = -1
            neigh(2,j,i) = -1
         End Do 
         Do j = 1,ndim
            coord(j,i) = -1.
            bnd_box(1,j,i) = -1.
            bnd_box(2,j,i) = -1.
         End Do 
         bsize(:,i) = -1.
         nodetype(i) = -1
         which_child(i) = -1
         newchild(i) = .FALSE.
         lrefine(i) = -1
         bflags(:,i) = -1
         work_block(i) = 0.
         If (empty_cells) Then
            empty(i) = 0
         End If
         
      End Do  ! End Do i = lnblocks2+1,lnblocks
      
      lnblocks = lnblocks2

!-----reset node types
      Do i = 1,lnblocks
         nodetype(i) = 3
         If (child(1,1,i) <= -1) Then
            nodetype(i) = 1
         End If
      End Do 
      nrecv = 0
      Do i = 1,lnblocks
         Do j = 1,nchild
            nodetype_chi(j,i) = -1 
            If (child(1,j,i) > -1) Then
            If (child(2,j,i).ne.mype) Then
               nrecv = nrecv + 1
               Call MPI_IRECV(nodetype_chi(j,i),                       & 
                              1,                                       & 
                              MPI_INTEGER,                             & 
                              child(2,j,i),                            & 
                              child(1,j,i),                            & 
                              amr_mpi_meshComm,                          & 
                              reqr(nrecv),                             & 
                              ierr)
            Else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            End If
            End If
         End Do 
      End Do  ! End Do i = 1, lnblocks

      nsend = 0
      Do i = 1,lnblocks
!--------send nodetype to your parent
         If (parent(1,i) >= 1) Then
         If (parent(2,i).ne.mype) Then
            nsend = nsend + 1
            Call MPI_SSEND(nodetype(i),                                & 
                           1,                                          & 
                           MPI_INTEGER,                                & 
                           parent(2,i),                                & 
                           i,                                          & 
                           amr_mpi_meshComm,                             & 
                           ierr)
         End If
         End If
      End Do   ! End Do i = 1, lnblocks

      If (nrecv > 0) Then
         Call MPI_WAITALL (nrecv, reqr, statr, ierr)
      End If

      Do i = 1,lnblocks
         Do j = 1,nchild
            If (nodetype_chi(j,i) == 1) nodetype(i) = 2
         End Do 
      End Do 

!-----reset derefine flags
      Do i = 1,maxblocks_tr
         derefine(i) = .FALSE.
      End Do

      End If

      Return
      End Subroutine amr_derefine_blocks


!FLASH's custom version of the same subroutine.  There are sufficient
!difference between this and the original that it is cleaner just
!to create a new subroutine.
Subroutine amr_derefine_blocks_flash (lnblocks_old, mype)

  !-----Use statements
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use constants
      Use Paramesh_comm_data, ONLY : amr_mpi_meshComm

  Implicit None

  !-----Input/Output arguments.
  Integer, intent(inout) :: lnblocks_old
  Integer, intent(in)    :: mype

  !-----Include statements.
  Include 'mpif.h'

  !-----Local arrays and variables
  Integer :: new_loc(maxblocks_tr)
  Integer :: i,j,k,jsend
  Integer :: lnblocks2
  Integer :: surr_blkst(2,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
  Integer :: childt(2,mchild,maxblocks_tr)
  Integer :: parentt(2,maxblocks_tr)
  Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
  Integer :: reqr(maxblocks_tr)
  Integer :: ierr,nsend,nrecv
  Integer :: nodetype_chi(nchild,maxblocks_tr)     
  Integer :: errorcode
  Integer :: li, lj, lk, faceAxis, faceSide, irf
  Integer, dimension(mdim) :: gCell, gs, gr

  !-----Begin Executable code.

  !Do some sanity checks about the geometry of the mesh.
  if ( lsingular_line .and. ndim > 1 .and. &
       (spherical_pm .or. polar_pm) ) then
     print *, "FLASH's amr_derefine_blocks not appropriate here."
     call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
  end if


  !By the time this subroutine is called we need correct data in 
  !surr_blks for all the top level blocks.
  If (.not.surr_blks_valid) then
     print *, "surr_blks_valid flag indicates top level not set"
     call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
  end if



  !-----remove blocks marked for derefinement by packing the data
  Do i = 1,maxblocks_tr
     new_loc(i) = -1
  End Do

  !-----Compute new_loc, new_loc marks where each block will end up after the
  !-----derefinement is done
  k = 1
  Do i = 1,lnblocks
     If (.not.derefine(i)) Then
        new_loc(i) = k
        k = k + 1
     End If
  End Do


  !-----reconnect all pointers
  parentt(:,1:lnblocks) = parent(:,1:lnblocks)
  childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
  surr_blks(3,:,:,:,:) = -999 !To make it clear that data in the 
  !nodetype index is nonsense until we call mpi_amr_exchange_nodetype.
  surr_blkst(1:2,:,:,:,1:lnblocks) = surr_blks(1:2,:,:,:,1:lnblocks)

  nrecv = 0
  Do i = 1,lnblocks
     If (parent(1,i) > 0) Then
        If (parent(2,i).ne.mype) Then
           nrecv = nrecv + 1
           Call MPI_IRECV(parentt(1,i),1,MPI_INTEGER,                & 
                parent(2,i),i,amr_mpi_meshComm,                        & 
                reqr(nrecv),ierr)
        Else
           parentt(1,i) = new_loc(parent(1,i))
        End If
     End If
  End Do

  nsend = 0
  Do i = 1,lnblocks
     Do j = 1,nchild
        If (child(1,j,i) > 0) Then
           If (child(2,j,i).ne.mype) Then
              nsend = nsend + 1
              Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,               & 
                   child(2,j,i),child(1,j,i),amr_mpi_meshComm,          & 
                   ierr)
           End If
        End If
     End Do
  End Do

  If (nrecv > 0) Then
     Call MPI_WAITALL(nrecv,reqr,statr,ierr)
  End If

  nrecv = 0
  Do i = 1,lnblocks
     Do j = 1,nchild
        If (child(1,j,i) > 0) Then
           If (child(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(childt(1,j,i),1,MPI_INTEGER,              & 
                   child(2,j,i),child(1,j,i),amr_mpi_meshComm,           & 
                   reqr(nrecv),ierr)
           Else
              childt(1,j,i) = new_loc(child(1,j,i))
           End If
        End If
     End Do
  End Do

  nsend = 0
  Do i = 1,lnblocks
     If (parent(1,i) > 0) Then
        If (parent(2,i).ne.mype) Then
           nsend = nsend + 1
           Call MPI_SSEND (new_loc(i),1,MPI_INTEGER,                 & 
                parent(2,i),i,amr_mpi_meshComm,                        & 
                ierr)
        End If
     End If
  End Do

  If (nrecv > 0) Then
     Call MPI_WAITALL(nrecv,reqr,statr,ierr)
  End If


  do lk = -k3d, k3d
     do lj = -k2d, k2d
        do li = -k1d, k1d
           gCell = (/li,lj,lk/)
           validRegion: if ( (0 /= sum(abs(gCell))) ) then
              gr = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)
              gs = (/(-li+k1d+1),(-lj+k2d+1),(-lk+k3d+1)/)

              nrecv = 0
              Do i = 1,lnblocks
                 If (surr_blks(1,gr(1),gr(2),gr(3),i) > 0) Then
                    If (surr_blks(2,gr(1),gr(2),gr(3),i).ne.mype) Then
                       nrecv = nrecv + 1
                       Call MPI_IRECV(surr_blkst(1,gr(1),gr(2),gr(3),i), 1, &
                            MPI_INTEGER, surr_blks(2,gr(1),gr(2),gr(3),i), &
                            surr_blks(1,gr(1),gr(2),gr(3),i), &
                            amr_mpi_meshComm, reqr(nrecv), ierr)
                    Else
                       surr_blkst(1,gr(1),gr(2),gr(3),i) = &
                            new_loc(surr_blks(1,gr(1),gr(2),gr(3),i))
                    End If
                 End If
              End Do

              nsend = 0
              Do i = 1,lnblocks
                 If (surr_blks(1,gs(1),gs(2),gs(3),i) > 0) Then
                    If (surr_blks(2,gs(1),gs(2),gs(3),i).ne.mype) Then
                       nsend = nsend + 1
                       Call MPI_SSEND(new_loc(i), 1, MPI_INTEGER, &
                            surr_blks(2,gs(1),gs(2),gs(3),i), i, & 
                            amr_mpi_meshComm, ierr)
                    End If
                 End If
              End Do

              If (nrecv > 0) Then
                 Call MPI_WAITALL(nrecv,reqr,statr,ierr)
              End If

           end if validRegion
        end do
     end do
  end do


  Do i = 1,lnblocks_old
     If (new_loc(i).ne.i.and.new_loc(i) > 0) Then
        If (nvar > 0) unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
        If (nfacevar > 0) Then
           facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
           facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
           facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
        End If
        If (nvaredge > 0) Then
           unk_e_x(:,:,:,:,new_loc(i)) = unk_e_x(:,:,:,:,i)
           unk_e_y(:,:,:,:,new_loc(i)) = unk_e_y(:,:,:,:,i)
           unk_e_z(:,:,:,:,new_loc(i)) = unk_e_z(:,:,:,:,i)
        End If
        If (nvarcorn > 0) unk_n(:,:,:,:,new_loc(i)) =                & 
             unk_n(:,:,:,:,i)
     End If
  End Do

  parent(1,1:lnblocks) = parentt(1,1:lnblocks)
  child(1,:,1:lnblocks) = childt(1,:,1:lnblocks)
  surr_blks(1,:,:,:,1:lnblocks) = surr_blkst(1,:,:,:,1:lnblocks)


  k = 1
  lnblocks2 = lnblocks
  Do i = 1,lnblocks

     If (.not.derefine(i)) Then

        If (k.ne.i) Then
           Do j = 1,nchild
              child(1,j,k) = child(1,j,i)
              child(2,j,k) = child(2,j,i)
           End Do
           parent(1,k) = parent(1,i)
           parent(2,k) = parent(2,i)
           do lk = -k3d, k3d
              do lj = -k2d, k2d
                 do li = -k1d, k1d
                    gr = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)
                    surr_blks(1:2,gr(1),gr(2),gr(3),k) = &
                         surr_blks(1:2,gr(1),gr(2),gr(3),i)
                 end do
              end do
           end do
           Do j = 1,ndim
              coord(j,k) = coord(j,i)
              bnd_box(1,j,k) = bnd_box(1,j,i)
              bnd_box(2,j,k) = bnd_box(2,j,i)
           End Do
           bsize(:,k) = bsize(:,i)
           newchild(k) = newchild(i)
           which_child(k) = which_child(i)
           lrefine(k) = lrefine(i)
           bflags(:,k) = bflags(:,i)
           work_block(k) = work_block(i)
           If (empty_cells) Then
              empty(k) = empty(i)
           End If

        End If

        k = k + 1

     Else

        lnblocks2 = lnblocks2 - 1
        lnblocks_old = lnblocks_old - 1

     End If  ! End If (.not.derefine(i))

  End Do  ! End Do i = 1,lnblocks

  !-----overwrite old locations
  Do i = lnblocks2+1,lnblocks

     derefine(i) = .FALSE.
     Do j = 1,nchild
        child(1,j,i) = -1
        child(2,j,i) = -1
     End Do
     parent(1,i) = -1
     parent(2,i) = -1
     do lk = -k3d, k3d
        do lj = -k2d, k2d
           do li = -k1d, k1d
              gr = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)
              surr_blks(1:2,gr(1),gr(2),gr(3),i) = -1
           end do
        end do
     end do
     Do j = 1,ndim
        coord(j,i) = -1.
        bnd_box(1,j,i) = -1.
        bnd_box(2,j,i) = -1.
     End Do
     bsize(:,i) = -1.
     nodetype(i) = -1
     which_child(i) = -1
     newchild(i) = .FALSE.
     lrefine(i) = -1
     bflags(:,i) = -1
     work_block(i) = 0.
     If (empty_cells) Then
        empty(i) = 0
     End If

  End Do  ! End Do i = lnblocks2+1,lnblocks


  lnblocks = lnblocks2

  !-----reset node types
  Do i = 1,lnblocks
     nodetype(i) = 3
     If (child(1,1,i) <= -1) Then
        nodetype(i) = 1
     End If
  End Do
  nrecv = 0
  Do i = 1,lnblocks
     Do j = 1,nchild
        nodetype_chi(j,i) = -1 
        If (child(1,j,i) > -1) Then
           If (child(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(nodetype_chi(j,i),             & 
                   1,                                       & 
                   MPI_INTEGER,                             & 
                   child(2,j,i),                            & 
                   child(1,j,i),                            & 
                   amr_mpi_meshComm,                          & 
                   reqr(nrecv),                             & 
                   ierr)
           Else
              nodetype_chi(j,i) = nodetype(child(1,j,i))
           End If
        End If
     End Do
  End Do  ! End Do i = 1, lnblocks

  nsend = 0
  Do i = 1,lnblocks
     !--------send nodetype to your parent
     If (parent(1,i) >= 1) Then
        If (parent(2,i).ne.mype) Then
           nsend = nsend + 1
           Call MPI_SSEND(nodetype(i),                      & 
                1,                                          & 
                MPI_INTEGER,                                & 
                parent(2,i),                                & 
                i,                                          & 
                amr_mpi_meshComm,                             & 
                ierr)
        End If
     End If
  End Do   ! End Do i = 1, lnblocks

  If (nrecv > 0) Then
     Call MPI_WAITALL (nrecv, reqr, statr, ierr)
  End If

  Do i = 1,lnblocks
     Do j = 1,nchild
        If (nodetype_chi(j,i) == 1) nodetype(i) = 2
     End Do
  End Do


  !-----reset derefine flags
  Do i = 1,maxblocks_tr
     derefine(i) = .FALSE.
  End Do

  Return
End Subroutine amr_derefine_blocks_flash
