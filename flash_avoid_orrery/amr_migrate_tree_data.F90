!!****fi source/amr_migrate_tree_data
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2010 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!   amr_migrate_tree_data
!!
!! SYNOPSIS
!!
!!   call amr_migrate_tree_data (new_loc,nprocs,mype)
!!   call amr_migrate_tree_data (integer array, integer, integer)
!!
!! ARGUMENTS
!!   
!!   integer, intent(inout) :: new_loc(:,:)
!!     new locations (processor and location in morton list) to migrate data to
!!   integer, intent(in) :: nprocs, mype
!!     number of procs. and proc. id.
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
!!   io
!!   paramesh_comm_data
!!
!! CALLS
!!
!!   fill_old_loc
!! 
!! RETURNS
!!
!!   Nothing returned
!!
!! DESCRIPTION
!!
!!   This subroutine move tree data and reconnect all pointers according to the 
!!   new locations which are passed in in the argument 'new_loc'.  
!!   'new_loc' is constructed during the morton sort and ordering by work 
!!   before it is passed to this routine.
!!
!! AUTHORS
!!
!!   Kevin Olson (1996-2001).
!!   Chris Daley (2010) -    modifications, variant amr_migrate_tree_data_flash
!!
!!***

#include "paramesh_preprocessor.fh"

      Subroutine amr_migrate_tree_data (new_loc,nprocs,mype)

!-----Use statements.
      Use paramesh_dimensions
      Use physicaldata
      Use tree
      Use io
      Use paramesh_comm_data
      Use paramesh_interfaces, only : fill_old_loc
      Use paramesh_mpi_interfaces, only : MPI_int_SSEND

#include "Flashx_mpi_implicitNone.fh"

!-----Input/Output variables.
      Integer, Intent(inout) :: new_loc(2,maxblocks_tr)
      Integer, Intent(in)    :: nprocs,mype

!-----Local variables and arrays
      Integer, Parameter :: buf_size = mdim+mdim+2*mdim
      Integer, Parameter :: ibuf_size = 2*mfaces+2*mchild+2+5+10
      Real    :: buffer(buf_size)
      Real    :: buffert(buf_size,maxblocks)
      Integer :: ibuffer(ibuf_size)
      Integer :: ibuffert(ibuf_size,maxblocks)
      Integer :: neight(2,mfaces,maxblocks_tr)
      Integer :: childt(2,mchild,maxblocks_tr)
      Integer :: parentt(2,maxblocks_tr)
      Integer :: i,j,k,jj
      Integer :: old_loc(2,maxblocks_tr)
      Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      Integer :: reqr(maxblocks_tr)
      Integer :: ierr,nsend,nrecv
      Logical :: newchildt(maxblocks_tr)
      Logical :: useFlashCustomVersion

      !It does not hurt to have an explicit interface here, but keep it consistent!
      interface
         Subroutine amr_migrate_tree_data_flash (new_loc,nprocs,mype)
           use tree, ONLY: maxblocks_tr
           implicit none
           Integer, Intent(inout) :: new_loc(2,maxblocks_tr)
           Integer, Intent(in)    :: nprocs,mype
         End Subroutine amr_migrate_tree_data_flash
      end interface


!-----Begin executable code.
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

         call amr_migrate_tree_data_flash (new_loc,nprocs,mype)

      Else

      Call fill_old_loc(new_loc,old_loc,nprocs,mype)

!-----count no. of new blocks
      new_lnblocks = 0
      Do i = 1,maxblocks
        If (old_loc(1,i) > -1) Then
          new_lnblocks = new_lnblocks + 1
        End If
      End Do

!-----update pointers to parents, children and neighbors
      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
             nrecv = nrecv + 1
             Call MPI_IRECV(parentt(1,i),2,MPI_INTEGER,                & 
                  parent(2,i),i,amr_mpi_meshComm,                        & 
                  reqr(nrecv),ierr)
           Else
             parentt(:,i) = new_loc(:,parent(1,i))
           End If
         End If
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         Do j = 1,nchild
           If (child(1,j,i) > 0) Then
             If (child(2,j,i).ne.mype) Then
               ! parent is sending to all its children
               nsend = nsend + 1
               Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,             & 
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
              Call MPI_IRECV(childt(1,j,i),2,MPI_INTEGER,              & 
                   child(2,j,i),child(1,j,i),amr_mpi_meshComm,           & 
                   reqr(nrecv),ierr)
            Else
              childt(:,j,i) = new_loc(:,child(1,j,i))
            End If
          End If
        End Do
       End Do
       
       nsend = 0
       Do i = 1,lnblocks
         If (parent(1,i) > 0) Then
           If (parent(2,i).ne.mype) Then
!------------child is sending to its parent
             nsend = nsend + 1
             Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,               & 
                  parent(2,i),i,amr_mpi_meshComm,                        & 
                  ierr)
           End If
         End If
       End Do

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      nrecv = 0
      Do i = 1,lnblocks
        Do j = 1,nfaces
          If (neigh(1,j,i) > 0) Then
            If (neigh(2,j,i).ne.mype) Then
              nrecv = nrecv + 1
              Call MPI_IRECV(neight(1,j,i),2,MPI_INTEGER,              & 
                   neigh(2,j,i),neigh(1,j,i),amr_mpi_meshComm,           & 
                   reqr(nrecv),ierr)
           Else
             neight(:,j,i) = new_loc(:,neigh(1,j,i))
            End If
          End If
        End Do
      End Do

      nsend = 0
      Do i = 1,lnblocks
        Do j = 1,nfaces
          If (neigh(1,j,i) > 0) Then
            If (neigh(2,j,i).ne.mype) Then
              nsend = nsend + 1
              Call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,          & 
                   neigh(2,j,i),i,amr_mpi_meshComm,                      & 
                   ierr)
            End If
          End If
        End Do
      End Do

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

      parent(:,1:lnblocks) = parentt(:,1:lnblocks)
      child(:,:,1:lnblocks) = childt(:,:,1:lnblocks)
      neigh(:,:,1:lnblocks) = neight(:,:,1:lnblocks)

!-----initialize temp buffer array
      ibuffer(:) = -1
      Do i = 1,maxblocks
         buffert(:,i) = -1.
         ibuffert(:,i) = -1
         newchildt(i) = .FALSE.
      End Do

      nrecv = 0
      Do i = 1,new_lnblocks
        If (old_loc(2,i).ne.mype) Then
          nrecv = nrecv + 1
          Call MPI_IRECV(buffert(1,i),buf_size,                        & 
               amr_mpi_real,                                           & 
               old_loc(2,i),i,amr_mpi_meshComm,                          & 
               reqr(nrecv),ierr)
          nrecv = nrecv + 1
          Call MPI_IRECV(ibuffert(1,i),ibuf_size,MPI_INTEGER,          & 
               old_loc(2,i),i+maxblocks_tr,amr_mpi_meshComm,             & 
               reqr(nrecv),ierr)
          nrecv = nrecv + 1
          Call MPI_IRECV(newchildt(i),1,MPI_LOGICAL,                   & 
               old_loc(2,i),i+2*maxblocks_tr,amr_mpi_meshComm,           & 
               reqr(nrecv),ierr)
        End If
      End Do

      nsend = 0

      Do i = 1,lnblocks

!-------pack buffer for sending
        k = 0
        Do j = 1,mdim
          k = k + 1
          buffer(k) = coord(j,i)
        End Do
        Do j = 1,mdim
          Do jj = 1,2
            k = k + 1
            buffer(k) = bnd_box(jj,j,i)
          End Do
        End Do
        Do j = 1,mdim
          k = k + 1
          buffer(k) = bsize(j,i)
        End Do

        k = 0
        Do j = 1,mchild
          Do jj = 1,2
            k = k + 1
            ibuffer(k) = child(jj,j,i)
          End Do
        End Do
        Do j = 1,mfaces
          Do jj = 1,2
            k = k + 1
            ibuffer(k) = neigh(jj,j,i)
          End Do
        End Do
        Do j = 1,2
          k = k + 1
          ibuffer(k) = parent(j,i)
        End Do
        k = k + 1
        ibuffer(k) = lrefine(i)
        k = k + 1
        ibuffer(k) = nodetype(i)
        k = k + 1
        ibuffer(k) = which_child(i)
        k = k + 1
        ibuffer(k) = empty(i)
        k = k + 1
        ibuffer(k) = work_block(i)
        Do j=1,mflags
        k = k + 1
        ibuffer(k) = bflags(j,i)
        End Do

        If (new_loc(2,i).ne.mype) Then
          nsend = nsend + 1
          Call MPI_SSEND(buffer,buf_size,amr_mpi_real,              &
               new_loc(2,i),new_loc(1,i),                              & 
               amr_mpi_meshComm,ierr)
          nsend = nsend + 1
          Call MPI_SSEND(ibuffer,ibuf_size,MPI_INTEGER,             &
               new_loc(2,i),new_loc(1,i)+maxblocks_tr,                 & 
               amr_mpi_meshComm,ierr)
          nsend = nsend + 1
          Call MPI_SSEND(newchild(i),1,MPI_LOGICAL,                    & 
               new_loc(2,i),new_loc(1,i)+2*maxblocks_tr,               & 
               amr_mpi_meshComm,ierr)
        Else
          buffert(1:buf_size,new_loc(1,i)) = buffer(1:buf_size)
          ibuffert(1:ibuf_size,new_loc(1,i)) = ibuffer(1:ibuf_size)
          newchildt(new_loc(1,i)) = newchild(i)
        End If

      End Do  ! End Do i = 1,lnblocks

      If (nrecv > 0) Then
        Call MPI_WAITALL(nrecv,reqr,statr,ierr)
      End If

!-----now unpack the buffer
      Do i = 1,new_lnblocks

        k = 0
        Do j = 1,mdim
          k = k + 1
          coord(j,i) = buffert(k,i)
        End Do
        Do j = 1,mdim
          Do jj = 1,2
            k = k + 1
            bnd_box(jj,j,i) = buffert(k,i)
          End Do
        End Do
        Do j = 1,mdim
          k = k + 1
          bsize(j,i) = buffert(k,i)
        End Do
        k = 0
        Do j = 1,mchild
          Do jj = 1,2
            k = k + 1
            child(jj,j,i) = ibuffert(k,i)
          End Do
        End Do
        Do j = 1,mfaces
          Do jj = 1,2
            k = k + 1
            neigh(jj,j,i) = ibuffert(k,i)
          End Do
        End Do
        Do j = 1,2
          k = k + 1
          parent(j,i) = ibuffert(k,i)
        End Do
        k = k + 1
        lrefine(i) = ibuffert(k,i)
        k = k + 1
        nodetype(i) = ibuffert(k,i)
        k = k + 1
        which_child(i) = ibuffert(k,i)
        k = k + 1
        empty(i) = max(ibuffert(k,i),0)     ! empty must be either 1 or 0.
                                            ! It cannot be -1.
        k = k + 1
        work_block(i) = ibuffert(k,i)
        Do j=1,mflags
        k = k + 1
        bflags(j,i) = ibuffert(k,i)
        End Do

        newchild(i) = newchildt(i)

      End Do  ! End Do i = 1,new_lnblocks

      End If

      Return
      End Subroutine amr_migrate_tree_data



!FLASH's custom version of the same subroutine.  There are sufficient
!difference between this and the original that it is cleaner just
!to create a new subroutine.
Subroutine amr_migrate_tree_data_flash (new_loc,nprocs,mype)

  !-----Use statements.
  Use paramesh_dimensions
  Use physicaldata
  Use tree
  Use io
  Use paramesh_comm_data
  Use paramesh_interfaces, only : fill_old_loc
  Use paramesh_mpi_interfaces, only : MPI_int_SSEND


#include "Flashx_mpi_implicitNone.fh"

  !-----Input/Output variables.
  Integer, Intent(inout) :: new_loc(2,maxblocks_tr)
  Integer, Intent(in)    :: nprocs,mype

  !-----Local variables and arrays
  Integer, Parameter :: buf_size = mdim+mdim+2*mdim
  Integer, Parameter :: ibuf_size = 2*mfaces+2*mchild+2+5+10+(3*3**mdim)
  Real    :: buffer(buf_size)
  Real    :: buffert(buf_size,maxblocks)
  Integer :: ibuffer(ibuf_size)
  Integer :: ibuffert(ibuf_size,maxblocks)
  Integer :: childt(2,mchild,maxblocks_tr)
  Integer :: parentt(2,maxblocks_tr)
  Integer :: i,j,k,jj
  Integer :: old_loc(2,maxblocks_tr)
  Integer :: statr(MPI_STATUS_SIZE,maxblocks_tr * (3**ndim - 1))
  Integer :: reqr(maxblocks_tr * (3**ndim - 1))
  Integer :: ierr,nsend,nrecv
  Logical :: newchildt(maxblocks_tr)

  !Declarations for FLASH's custom implementation.
  integer, dimension(mdim) :: gCell, gs, gr
  integer :: li, lj, lk, f, faceAxis, faceSide, irf

  !Note 1st dim is 2 not 3. We fill nodetype much later.
  Integer :: surr_blkst(2,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
  Integer :: errorcode

  !-----Begin executable code.

  !Do some sanity checks about the geometry of the mesh.
  if ( lsingular_line .and. ndim > 1 .and. &
       (spherical_pm .or. polar_pm) ) then
     print *, "FLASH's amr_migrate_tree_data not appropriate here."
     call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
  end if


  !By the time this subroutine is called we need correct data in 
  !surr_blks for all the top level blocks.
  If (.not.surr_blks_valid) then
     print *, "surr_blks_valid flag indicates top level not set"
     call mpi_abort(amr_mpi_meshComm,errorcode,ierr)
  end if


  Call fill_old_loc(new_loc,old_loc,nprocs,mype)

  !-----count no. of new blocks
  new_lnblocks = 0
  Do i = 1,maxblocks
     If (old_loc(1,i) > -1) Then
        new_lnblocks = new_lnblocks + 1
     End If
  End Do

  !-----update pointers to parents, children and neighbors
  parentt(:,1:lnblocks) = parent(:,1:lnblocks)
  childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
  surr_blkst(1:2,:,:,:,1:lnblocks) = surr_blks(1:2,:,:,:,1:lnblocks)

  nrecv = 0
  Do i = 1,lnblocks
     If (parent(1,i) > 0) Then
        If (parent(2,i).ne.mype) Then
           nrecv = nrecv + 1
           Call MPI_IRECV(parentt(1,i),2,MPI_INTEGER,                & 
                parent(2,i),i,amr_mpi_meshComm,                        & 
                reqr(nrecv),ierr)
        Else
           parentt(:,i) = new_loc(:,parent(1,i))
        End If
     End If
  End Do

  nsend = 0
  Do i = 1,lnblocks
     Do j = 1,nchild
        If (child(1,j,i) > 0) Then
           If (child(2,j,i).ne.mype) Then
              ! parent is sending to all its children
              nsend = nsend + 1
              Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,             & 
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
              Call MPI_IRECV(childt(1,j,i),2,MPI_INTEGER,              & 
                   child(2,j,i),child(1,j,i),amr_mpi_meshComm,           & 
                   reqr(nrecv),ierr)
           Else
              childt(:,j,i) = new_loc(:,child(1,j,i))
           End If
        End If
     End Do
  End Do

  nsend = 0
  Do i = 1,lnblocks
     If (parent(1,i) > 0) Then
        If (parent(2,i).ne.mype) Then
           !------------child is sending to its parent
           nsend = nsend + 1
           Call MPI_SSEND (new_loc(1,i),2,MPI_INTEGER,               & 
                parent(2,i),i,amr_mpi_meshComm,                        & 
                ierr)
        End If
     End If
  End Do

  If (nrecv > 0) Then
     Call MPI_WAITALL(nrecv,reqr,statr,ierr)
  End If



  nrecv = 0
  Do i = 1,lnblocks
     do lk = -k3d, k3d
        do lj = -k2d, k2d
           do li = -k1d, k1d
              gCell = (/li,lj,lk/)
              if ( (0 /= sum(abs(gCell))) ) then
                 gr = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)

                 if (surr_blks(1,gr(1),gr(2),gr(3),i) > 0) then
                    if (surr_blks(2,gr(1),gr(2),gr(3),i).ne.myPE) then
                       nrecv = nrecv + 1
                       Call MPI_IRECV(surr_blkst(1,gr(1),gr(2),gr(3),i), &
                            2, MPI_INTEGER, &
                            surr_blks(2,gr(1),gr(2),gr(3),i), &
                            surr_blks(1,gr(1),gr(2),gr(3),i), &
                            amr_mpi_meshComm, reqr(nrecv),ierr)
                    else
                       surr_blkst(1:2,gr(1),gr(2),gr(3),i) = &
                            new_loc(1:2,surr_blks(1,gr(1),gr(2),gr(3),i))
                    end if
                 end if
              end if
           end do
        end do
     end do
  end do


  nsend = 0
  Do i = 1,lnblocks
     do lk = -k3d, k3d
        do lj = -k2d, k2d
           do li = -k1d, k1d
              gCell = (/li,lj,lk/)
              if ( (0 /= sum(abs(gCell))) ) then
                 gs = (/(-li+k1d+1),(-lj+k2d+1),(-lk+k3d+1)/)

                 if (surr_blks(1,gs(1),gs(2),gs(3),i) > 0) then
                    if (surr_blks(2,gs(1),gs(2),gs(3),i).ne.myPE) then
                       nsend = nsend + 1
                       Call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,          & 
                            surr_blks(2,gs(1),gs(2),gs(3),i),&
                            i,amr_mpi_meshComm,ierr)
                    end if
                 end if
              end if
           end do
        end do
     end do
  end do

  If (nrecv > 0) Then
     Call MPI_WAITALL(nrecv,reqr,statr,ierr)
  End If


  parent(:,1:lnblocks) = parentt(:,1:lnblocks)
  child(:,:,1:lnblocks) = childt(:,:,1:lnblocks)
  surr_blks(1:2,:,:,:,1:lnblocks) = surr_blkst(1:2,:,:,:,1:lnblocks)


  !-----initialize temp buffer array
  Do i = 1,maxblocks
     buffert(:,i) = -1.
     ibuffert(:,i) = -1
     newchildt(i) = .FALSE.
  End Do

  nrecv = 0
  Do i = 1,new_lnblocks
     If (old_loc(2,i).ne.mype) Then
        nrecv = nrecv + 1
        Call MPI_IRECV(buffert(1,i),buf_size,                        & 
             amr_mpi_real,                                           & 
             old_loc(2,i),i,amr_mpi_meshComm,                          & 
             reqr(nrecv),ierr)
        nrecv = nrecv + 1
        Call MPI_IRECV(ibuffert(1,i),ibuf_size,MPI_INTEGER,          & 
             old_loc(2,i),i+maxblocks_tr,amr_mpi_meshComm,             & 
             reqr(nrecv),ierr)
        nrecv = nrecv + 1
        Call MPI_IRECV(newchildt(i),1,MPI_LOGICAL,                   & 
             old_loc(2,i),i+2*maxblocks_tr,amr_mpi_meshComm,           & 
             reqr(nrecv),ierr)
     End If
  End Do

  nsend = 0

  Do i = 1,lnblocks

     !-------pack buffer for sending
     k = 0
     Do j = 1,mdim
        k = k + 1
        buffer(k) = coord(j,i)
     End Do
     Do j = 1,mdim
        Do jj = 1,2
           k = k + 1
           buffer(k) = bnd_box(jj,j,i)
        End Do
     End Do
     Do j = 1,mdim
        k = k + 1
        buffer(k) = bsize(j,i)
     End Do

     k = 0
     Do j = 1,mchild
        Do jj = 1,2
           k = k + 1
           ibuffer(k) = child(jj,j,i)
        End Do
     End Do
     Do j = 1,2
        k = k + 1
        ibuffer(k) = parent(j,i)
     End Do
     k = k + 1
     ibuffer(k) = lrefine(i)
     k = k + 1
     ibuffer(k) = nodetype(i)
     k = k + 1
     ibuffer(k) = which_child(i)
     k = k + 1
     ibuffer(k) = empty(i)
     k = k + 1
     ibuffer(k) = work_block(i) !CD: This should be in real array not integer array.
     Do j=1,mflags
        k = k + 1
        ibuffer(k) = bflags(j,i)
     End Do

     !Pack surr_blks.
     do lk = 1, 1+2*k3d
        do lj = 1, 1+2*k2d
           do li = 1, 1+2*k1d
              do f = 1, 2
                 k = k + 1
                 ibuffer(k) = surr_blks(f,li,lj,lk,i)
              end do
           end do
        end do
     end do

     !CD: So that we do not send any uninitialized data.  It is never
     !looked at by the receiver, but it causes memory checkers to
     !complain.  Need to make ibuf_size smaller eventually.
     if (k < ibuf_size) then
        ibuffer(k+1:ibuf_size) = -999
     end if


     If (new_loc(2,i).ne.mype) Then
        nsend = nsend + 1
        Call MPI_SSEND(buffer,buf_size,amr_mpi_real,              &
             new_loc(2,i),new_loc(1,i),                              & 
             amr_mpi_meshComm,ierr)
        nsend = nsend + 1
        Call MPI_SSEND(ibuffer,ibuf_size,MPI_INTEGER,             &
             new_loc(2,i),new_loc(1,i)+maxblocks_tr,                 & 
             amr_mpi_meshComm,ierr)
        nsend = nsend + 1
        Call MPI_SSEND(newchild(i),1,MPI_LOGICAL,                    & 
             new_loc(2,i),new_loc(1,i)+2*maxblocks_tr,               & 
             amr_mpi_meshComm,ierr)
     Else
        buffert(1:buf_size,new_loc(1,i)) = buffer(1:buf_size)
        ibuffert(1:ibuf_size,new_loc(1,i)) = ibuffer(1:ibuf_size)
        newchildt(new_loc(1,i)) = newchild(i)
     End If

  End Do  ! End Do i = 1,lnblocks

  If (nrecv > 0) Then
     Call MPI_WAITALL(nrecv,reqr,statr,ierr)
  End If

  !-----now unpack the buffer
  Do i = 1,new_lnblocks

     k = 0
     Do j = 1,mdim
        k = k + 1
        coord(j,i) = buffert(k,i)
     End Do
     Do j = 1,mdim
        Do jj = 1,2
           k = k + 1
           bnd_box(jj,j,i) = buffert(k,i)
        End Do
     End Do
     Do j = 1,mdim
        k = k + 1
        bsize(j,i) = buffert(k,i)
     End Do
     k = 0
     Do j = 1,mchild
        Do jj = 1,2
           k = k + 1
           child(jj,j,i) = ibuffert(k,i)
        End Do
     End Do
     Do j = 1,2
        k = k + 1
        parent(j,i) = ibuffert(k,i)
     End Do
     k = k + 1
     lrefine(i) = ibuffert(k,i)
     k = k + 1
     nodetype(i) = ibuffert(k,i)
     k = k + 1
     which_child(i) = ibuffert(k,i)
     k = k + 1
     empty(i) = max(ibuffert(k,i),0)     ! empty must be either 1 or 0.
     ! It cannot be -1.
     k = k + 1
     work_block(i) = ibuffert(k,i)
     Do j=1,mflags
        k = k + 1
        bflags(j,i) = ibuffert(k,i)
     End Do

     !Unpack surr_blks.
     do lk = 1, 1+2*k3d
        do lj = 1, 1+2*k2d
           do li = 1, 1+2*k1d
              do f = 1, 2
                 k = k + 1
                 surr_blks(f,li,lj,lk,i) = ibuffert(k,i)
              end do
           end do
        end do
     end do
     !Update central surr_blks index to myself.
     surr_blks(1,1+k1d,1+k2d,1+k3d,i) = i
     surr_blks(2,1+k1d,1+k2d,1+k3d,i) = myPE
     surr_blks(3,1+k1d,1+k2d,1+k3d,i) = nodetype(i)

     newchild(i) = newchildt(i)

  End Do  ! End Do i = 1,new_lnblocks


  call mpi_amr_exchange_nodetype (mype, new_lnblocks)


  !CD: At this point we copy into neigh (to preserve original behavior)
  do i = 1, new_lnblocks
     do lk = -k3d, k3d
        do lj = -k2d, k2d
           do li = -k1d, k1d
              gCell = (/li,lj,lk/)
              if ( (0 /= sum(abs(gCell))) ) then  !valid region
                 gr = (/(li+k1d+1),(lj+k2d+1),(lk+k3d+1)/)

                 !The orrery calculation sets block and processor slot
                 !to the negative number, so do this here too.
                 If (surr_blks(1,gr(1),gr(2),gr(3),i) <= -1) Then
                    surr_blks(2,gr(1),gr(2),gr(3),i) = &
                         surr_blks(1,gr(1),gr(2),gr(3),i)
                 End If

                 if ( (1 == sum(abs(gCell))) ) then
                    faceAxis = sum((/1,2,3/) * abs(gCell)) !x:1, y:2, z:3
                    faceSide = sum(gCell) !low side:-1, high side:+1
                    irf = faceAxis*2 - dim(0,faceSide)
                    neigh(1:2,irf,i) = surr_blks(1:2,gr(1),gr(2),gr(3),i)
                 end if
              end if
           end do
        end do
     end do
  end do

  Return
End Subroutine amr_migrate_tree_data_flash
