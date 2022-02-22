!!****if* source/Grid/GridMain/AMR/Paramesh4/bittree/amr_verify_bittree
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
!!  amr_verify_bittree
!!
!! SYNOPSIS
!!
!!  call amr_verify_bittree()
!!
!! DESCRIPTION
!!
!!  Have bittree do a sanity check by making sure the information in tree agrees
!!  with bittree's calculations.
!!
!! ARGUMENTS
!!  none
!!
!! HISTORY
!!
!!  2021 - 2022  Tom Klosterman
!!***

subroutine amr_verify_bittree()
  
  use bittree, only : gr_btIdentify, bittree_block_count, &
                      gr_getIntCoords
  use Paramesh_comm_data, only: amr_mpi_meshComm
  use paramesh_dimensions, only: ndim
  use tree, only: lnblocks, bsize, coord, lrefine, &
                  grid_xmin, grid_ymin, grid_zmin, &
                  grid_xmax, grid_ymax, grid_zmax
  use iso_c_binding, only: c_bool, c_int
  use Driver_interface, only: Driver_abort
  implicit none
  
#include "Flashx_mpi.h"

  integer :: b, locb, lev, proc0, proc1, nprocs, ierr
  integer(c_int) :: nb
  integer :: i
  integer :: lcoord(3)
  logical :: invalid
  logical(c_bool) :: cfalse = .FALSE.
  
  integer, allocatable :: all_nblock(:), all_recv(:), all_disp(:)
  integer, allocatable :: coord_list(:,:), all_coords(:,:)
  
  call MPI_COMM_RANK(amr_mpi_meshComm, proc0, ierr)
  call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)

!-Make sure bittree returns correct proc/locblk for all local blocks
  invalid = .false.
  do b=1, lnblocks
    call gr_getIntCoords(b,lcoord)
    lev = lrefine(b)
    call gr_btIdentify(nprocs,lev,lcoord,proc1,locb)

    if(lrefine(b) /= lev .or. proc0 /= proc1 .or. b /= locb) then
      if(proc0 .eq. 0) then
        if(.not. invalid) print *,'BITTREE IS WRONG!!! proc=', proc0
        print *, ' actual/bittree: proc:',proc0,proc1,'locb:',b,locb,'lev:',lrefine(b),lev
      end if
      invalid = .true.
    end if

  end do
  
  call MPI_Allreduce(MPI_IN_PLACE, invalid, 1, FLASH_LOGICAL, MPI_LOR, amr_mpi_meshComm, ierr)
 
!-If any blocks were invalid on any processor, output errors in bittree.misery.log and abort
  if(invalid) then
    if(proc0 == 0) allocate(all_nblock(nprocs))
    call MPI_Gather(lnblocks, 1, FLASH_INTEGER, all_nblock, 1, FLASH_INTEGER, 0, amr_mpi_meshComm, ierr)
    if(proc0 == 0) print *, 'gather 1'
    
    allocate(coord_list(1+ndim,lnblocks))
    do b=1, lnblocks
      coord_list(1,b) = lrefine(b)
      call gr_getIntCoords(b,lcoord)
      coord_list(2:1+ndim,b) = lcoord(1:ndim)
    end do
    
    if(proc0 == 0) then
      allocate(all_recv(nprocs))
      allocate(all_disp(nprocs))
      allocate(all_coords(1+ndim,sum(all_nblock)))
      all_recv = all_nblock*(1+ndim)
      all_disp(1) = 0
      do i=2, nprocs
        all_disp(i) = all_disp(i-1) + all_recv(i-1)
      end do
    end if
    
    call MPI_Gatherv(coord_list, (1+ndim)*lnblocks, FLASH_INTEGER, all_coords, &
      all_recv, all_disp, FLASH_INTEGER, 0, amr_mpi_meshComm, ierr)
    if(proc0 == 0) print *, 'gather 2'
    
    if(proc0 == 0) then
      all_disp = all_disp/(1+ndim) ! change all_disp to starting 0-based block
      
      open(1349,file='bittree.misery.log',action='write')
      
      call bittree_block_count(cfalse,nb)
      write(1349,*) 'procs=', nprocs, ' actual blocks=', ubound(all_coords,2), ' bittree blocks=', nb
      write(1349,*) ''
      
      do b=1, ubound(all_coords,2)
        if(proc0 < nprocs-1) then
          if(b-1 >= all_disp(proc0+2)) then
            proc0 = proc0 + 1
          end if
        end if
        lev = all_coords(1,b)
        call gr_btIdentify(nprocs, lev, all_coords(2:,b), proc1, locb)
        if(lev /= all_coords(1,b) .or. proc0 /= proc1 .or. b-all_disp(proc0+1) /= locb) then
          write(1349,*) 'lev=', all_coords(1,b), ' ijk=', all_coords(2:,b), &
            ' proc=',proc0, ' lblk=', b-all_disp(proc0+1), &
            ' BITTREE lev=',lev,' proc=',proc1,' lblk=', locb
        else
          write(1349,*) 'lev=', all_coords(1,b), ' ijk=', all_coords(2:,b), &
            ' proc=',proc0, ' lblk=', b-all_disp(proc0+1)
        end if
      end do
      close(1349)
      proc0 = 0
    end if
    
    deallocate(coord_list)
    if(proc0 == 0) then
      deallocate(all_nblock)
      deallocate(all_recv)
      deallocate(all_disp)
      deallocate(all_coords)
      call Driver_abort('Bittree suicide, see bittree.misery.log for discrepancies.')
    end if
  end if
end subroutine
