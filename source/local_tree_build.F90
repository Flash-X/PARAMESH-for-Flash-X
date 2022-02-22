!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/local_tree_build
!! NAME
!!
!!   local_tree_build
!!
!! SYNOPSIS
!!
!!   Call local_tree_build()
!!
!! ARGUMENTS
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!   local_tree_common
!!   physicaldata
!!   tree
!!   paramesh_dimensions
!!   paramesh_comm_data
!!   mpi_morton
!!
!! CALLS
!!
!!   add_block_to_tree
!!
!! RETURNS
!!
!!   Does not return anything.  Upon exit a local tree has been constructed.
!!
!! DESCRIPTION
!!
!!   This routine constucts a unique tree of blocks local to each processor.  
!!   Each tree represents the entire physical domain, but each is pruned at 
!!   different points depending on the local list of blocks stored on any 
!!   particular processor.  The pruning is done such that each local tree
!!   contains just enough information
!!   so that is can be searched by the local blocks to find their surrounding
!!   blocks.  A digital orrey is used to pass the local block lists from
!!   processor to processor.  As new data arrive they are added (or not) to
!!   the tree by calling the routine 'add_block_to_tree'.
!!
!! AUTHORS
!!
!!   Kevin Olson
!!
!!***

#include "paramesh_preprocessor.fh"

Subroutine local_tree_build ()
  
  !-----Use Statements.
  Use local_tree_common
  Use physicaldata
  Use tree
  Use paramesh_dimensions
  Use paramesh_comm_data
  Use mpi_morton, Only : lperiodicx, lperiodicy, lperiodicz
  
  Implicit None
  
  !-----Include Statements.
  Include 'mpif.h'
  
  !-----Local Variables
  
  Type(node), Pointer :: tpar
  Real :: lcoord(3)
  Real, Allocatable :: coord1(:,:), bsize1(:,:)
  Integer :: lnblocks_max
  Integer :: lnblocks_old
  Integer :: proc, ipass
  Integer :: status(MPI_STATUS_SIZE) 
  Integer :: idest, lb, isrc, iproc, ierr, nprocs, mype
  Integer :: no_of_level1_blocks, no_of_level1_blocks_loc
  Integer :: i,j,k,imin,imax,jmin,jmax,kmin,kmax
  Integer :: n_children_at_0
  Integer :: nxface, nyface, nzface
  Integer :: lb1, lnblocks1, lnblocks_max1
  Integer, Allocatable :: nodetype1(:), neigh1(:,:,:), block_id(:)
  Logical :: off_proc, searched
  
  !-----Begin Executable Code
  
  Call MPI_COMM_SIZE (amr_mpi_meshComm,nprocs,ierr)
  Call MPI_COMM_RANK (amr_mpi_meshComm,mype,ierr)
  
  Nullify(local_tree)
  Nullify(tpar)
  
  Call MPI_ALLREDUCE(lnblocks, lnblocks_max, 1,                    &
       MPI_INTEGER, MPI_MAX,                         &
       amr_mpi_meshComm, ierr)
  
  !-----Find the total number of level 1 blocks
  no_of_level1_blocks_loc = 0
  no_of_level1_blocks     = 0
  Do lb = 1, lnblocks
     If (lrefine(lb) == 1) Then
        no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
        !-----------ADD extra level 1 'ghost' blocks if boundaries are periodic
        nxface = 0
        If (coord(1,lb) - bsize(1,lb) < grid_xmin .and.            &
             neigh(1,1,lb) > 0) Then ! periodic
           no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
           nxface = nxface + 1
        End If
        If (coord(1,lb) + bsize(1,lb) > grid_xmax .and.            &
             neigh(1,2,lb) > 0) Then ! periodic
           no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
           nxface = nxface + 1
        End If
        
        nyface = 0
        If (ndim >= 2) Then 
           If (coord(2,lb) - bsize(2,lb) < grid_ymin .and.            &
                neigh(1,3,lb) > 0) Then ! periodic
              no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
              nyface = nyface + 1
           End If
           If (coord(2,lb) + bsize(2,lb) > grid_ymax .and.            &
                neigh(1,4,lb) > 0) Then ! periodic
              no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
              nyface = nyface + 1
           End If
        End If  ! End If (ndim >= 2)
        
        nzface = 0
        If (ndim == 3) Then
           If (coord(3,lb) - bsize(3,lb) < grid_zmin .and.            &
                neigh(1,5,lb) > 0) Then ! periodic
              no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
              nzface = nzface + 1
           End If
           If (coord(3,lb) + bsize(3,lb) > grid_zmax .and.   &
                neigh(1,6,lb) > 0) Then ! periodic
              no_of_level1_blocks_loc = no_of_level1_blocks_loc + 1
              nzface = nzface + 1
           End If
        End If  ! End If (ndim == 3)
        
        !-----------Edges
        no_of_level1_blocks_loc = no_of_level1_blocks_loc +        &
             nxface*nyface
        no_of_level1_blocks_loc = no_of_level1_blocks_loc +        &
             nxface*nzface
        no_of_level1_blocks_loc = no_of_level1_blocks_loc +        &
             nyface*nzface
        !-----------Corners
        no_of_level1_blocks_loc = no_of_level1_blocks_loc +        &
             nxface*nyface*nzface
        
     end if
  end do
  Call MPI_ALLREDUCE(no_of_level1_blocks_loc, no_of_level1_blocks, &
       1,                                            &
       MPI_INTEGER, MPI_SUM,                         &
       amr_mpi_meshComm, ierr)
  
  !-----Store data in temp arrays about level 1 blocks so that all the data is not
  !-----passed around the digital orrey for ipass = 1
  lnblocks1 = 0
  Do lb = 1, lnblocks
     If (lrefine(lb) == 1) Then
        lnblocks1 = lnblocks1 + 1
     End if
  End Do
  Call MPI_ALLREDUCE(lnblocks1, lnblocks_max1, 1,                  &
       MPI_INTEGER, MPI_MAX,                         &
       amr_mpi_meshComm, ierr)
  
  !-----Allocate temp space for arrays needed.
  Allocate(nodetype1(lnblocks_max1))
  Allocate(coord1(mdim,lnblocks_max1))
  Allocate(bsize1(mdim,lnblocks_max1))
  Allocate(neigh1(2,mfaces,lnblocks_max1))
  Allocate(block_id(lnblocks_max1))
  nodetype1 = -9999
  coord1 = -9999.9
  bsize1 = -9999.9
  neigh1 = -9999
  block_id = -9999


  
  !-----Store the data from level 1 blocks in temp arrays
  lb1 = 1
  Do lb = 1,lnblocks
     If (lrefine(lb) == 1) Then
        nodetype1(lb1)  = nodetype(lb)
        coord1(:,lb1)   = coord(:,lb)
        bsize1(:,lb1)   = bsize(:,lb)
        neigh1(:,:,lb1) = neigh(:,:,lb)
        block_id(lb1) = lb
        lb1 = lb1 + 1
     End If
  End do
  
  n_children_at_0 = 0
  !-----First pass only add lrefine = 1 blocks to tree(s)
  !-----Second pass add the rest of the blocks.
  Do ipass = 1,2  
     
     lnblocks_old = lnblocks
     proc = mype
     !-----Loop through all processors
     Do iproc = 0, nprocs-1
        
        If (iproc == 0) Then
           off_proc = .False.
        Else
           off_proc = .True.
        End If
        
        If (ipass == 1) Then
           Do lb1 = 1, lnblocks1
              off_proc = .False.
              searched = .False.
              
              kmin = 0
              kmax = 0
              jmin = 0
              jmax = 0
              imin = 0
              imax = 0
              
              !-----------Tests for periodic boundary conditions and reset of 
              !-----------of loop indeces
              If (coord1(1,lb1) - bsize1(1,lb1) < grid_xmin .and.        &
                   neigh1(1,1,lb1) > 0) Then ! periodic
                 imax = 1
              End If
              If (coord1(1,lb1) + bsize1(1,lb1) > grid_xmax .and.        &
                   neigh1(1,2,lb1) > 0) Then ! periodic
                 imin = -1
              End If
              
              If (ndim >= 2) Then 
                 If (coord1(2,lb1) - bsize1(2,lb1) < grid_ymin .and.     &
                      neigh1(1,3,lb1) > 0) Then ! periodic
                    jmax = 1
                 End If
                 If (coord1(2,lb1) + bsize1(2,lb1) > grid_ymax .and.     &
                      neigh1(1,4,lb1) > 0) Then ! periodic
                    jmin = -1
                 End If
              End If  ! End If (ndim >= 2)
              
              If (ndim == 3) Then
                 If (coord1(3,lb1) - bsize1(3,lb1) < grid_zmin .and.     &
                      neigh1(1,5,lb1) > 0) Then ! periodic
                    kmax = 1
                 End If
                 If (coord1(3,lb1) + bsize1(3,lb1) > grid_zmax .and.     &
                      neigh1(1,6,lb1) > 0) Then ! periodic
                    kmin = -1
                 End If
              End If
              
              !-----------Loop through possible neighbors and add block to local tree
              Do k = kmin, kmax
                 lcoord(3) = coord1(3,lb1) + k*(grid_zmax-grid_zmin)
                 Do j = jmin, jmax
                    lcoord(2) = coord1(2,lb1) + j*(grid_ymax-grid_ymin)
                    Do i = imin, imax
                       lcoord(1) = coord1(1,lb1) + i*(grid_xmax-grid_xmin)
                       Call add_block_to_tree (local_tree, tpar, local_tree,   &
                            lcoord(:),bsize1(:,lb1),1,                         & 
                            block_id(lb1), proc, nodetype1(lb1),               & 
                            1, n_children_at_0, no_of_level1_blocks,           &
                            off_proc, searched)
                    End Do
                 End Do
              End Do
           End Do  ! End Do lb1 = 1, lnblocks1
           
        Else
           
           Do lb = 1, lnblocks
              If (lrefine(lb) .ne. 1) Then
                 searched = .False.
                 
                 kmin = 0
                 kmax = 0
                 jmin = 0
                 jmax = 0
                 imin = 0
                 imax = 0
                 
                 !--------------Tests for periodic boundary conditions and reset of 
                 !--------------of loop indeces
                 !CD: Changed neigh(1,1:6,lb) > 0 to neigh(1,1:6,lb) >= -1
                 !(4 Feb 2010).  This fixed problem in 2D periodic simulation
                 !where lower-left block of computational domain did not
                 !locate upper-right block of computational domain.  Face 1
                 !and face 3 of lower-left block was -1 indicating that there
                 !was no neighbor at same refinement along face directions.
                 !However, there was a same refinement level neighbor along
                 !the corner direction (i.e. upper-right block).  Making the
                 !change means that we search for the corner neighbor.
                 If (coord(1,lb) - bsize(1,lb) < grid_xmin .and.         &
                      neigh(1,1,lb) >= -1) Then ! periodic
                    imax = 1
                 End If
                 If (coord(1,lb) + bsize(1,lb) > grid_xmax .and.         &
                      neigh(1,2,lb) >= -1) Then ! periodic
                    imin = -1
                 End If
                 
                 If (ndim >= 2) Then 
                    If (coord(2,lb) - bsize(2,lb) < grid_ymin .and.         &
                         neigh(1,3,lb) >= -1) Then ! periodic
                       jmax = 1
                    End If
                    If (coord(2,lb) + bsize(2,lb) > grid_ymax .and.         &
                         neigh(1,4,lb) >= -1) Then ! periodic
                       jmin = -1
                    End If
                 End If  ! End If (ndim >= 2)
                 
                 If (ndim == 3) Then
                    If (coord(3,lb) - bsize(3,lb) < grid_zmin .and.         &
                         neigh(1,5,lb) >= -1) Then ! periodic
                       kmax = 1
                    End If
                    If (coord(3,lb) + bsize(3,lb) > grid_zmax .and.         &
                         neigh(1,6,lb) >= -1) Then ! periodic
                       kmin = -1
                    End If
                 End If  ! End If (ndim == 3)
                 
                 !--------------Loop through possible neighbors and add block to local tree
                 Do k = kmin, kmax
                    lcoord(3) = coord(3,lb) + k*(grid_zmax-grid_zmin)
                    Do j = jmin, jmax
                       lcoord(2) = coord(2,lb) + j*(grid_ymax-grid_ymin)
                       Do i = imin, imax
                          lcoord(1) = coord(1,lb) + i*(grid_xmax-grid_xmin)
                          Call add_block_to_tree (local_tree, tpar, local_tree,& 
                               lcoord(:),bsize(:,lb),lrefine(lb),& 
                               lb, proc, nodetype(lb),           & 
                               1, n_children_at_0, no_of_level1_blocks,&
                               off_proc, searched)
                       End Do
                    End Do
                 End Do
                 
              End If  ! End If (lrefine(lb) .ne. 1)
           End Do   ! End Do lb = 1, lnblocks
        End If   ! End If (ipass == 1) 
        
        !-----Now perform digital orrey communication.
        
        !-----Destination processor to send data to.
        idest = mype + 1
        If (idest > nprocs-1) idest = 0
        !-----Source processor to receive data from.
        isrc = mype - 1
        If (isrc < 0) isrc = nprocs-1
        
        proc = proc - 1
        If (proc < 0) proc = nprocs-1
        
        If (ipass == 1) Then
           
           Call MPI_SENDRECV_REPLACE (lnblocks1, 1, MPI_INTEGER,            &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (block_id, lnblocks_max1, MPI_INTEGER, &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (nodetype1, lnblocks_max1, MPI_INTEGER,&
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (neigh1, 2*mfaces*lnblocks_max1,       &
                MPI_INTEGER,                          &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (coord1, mdim*lnblocks_max1,           &
                amr_mpi_real,                         &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (bsize1, mdim*lnblocks_max1,           &
                amr_mpi_real,                         &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
        Else
           
           Call MPI_SENDRECV_REPLACE (lnblocks, 1, MPI_INTEGER,             &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (lrefine, lnblocks_max, MPI_INTEGER,   &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (nodetype, lnblocks_max, MPI_INTEGER,  &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (neigh, 2*mfaces*lnblocks_max,         &
                MPI_INTEGER,                          &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (coord, mdim*lnblocks_max,             &
                amr_mpi_real,                         &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
           
           Call MPI_SENDRECV_REPLACE (bsize, mdim*lnblocks_max,             &
                amr_mpi_real,                         &
                idest, mype,                          &
                isrc,  isrc,                          &
                amr_mpi_meshComm, status, ierr)
        End if  !  end if (ipass == 1)
        
     End Do  ! End Do iproc = 0, nprocs-1
     
  End Do  ! End Do ipass = 1,2
  
  Deallocate(nodetype1)
  Deallocate(coord1)
  Deallocate(bsize1)
  Deallocate(neigh1)
  Deallocate(block_id)
  
  If (mype == 0) Then
     If (lnblocks .ne. lnblocks_old) Then
        Print *,' ERROR after digital orrey !!!!!!!!!!!!!!!!!!!!! '
        Print *,' lnblocks = ',lnblocks,lnblocks_old
     End If
  End If
  
End Subroutine local_tree_build

