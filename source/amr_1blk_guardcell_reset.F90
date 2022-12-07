!!****if* source/amr_1blk_guardcell_reset
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
!!   amr_1blk_guardcell_reset
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_guardcell_reset()
!!
!! ARGUMENTS
!!
!!  None
!!
!! INCLUDES
!!
!! USES
!!
!!   paramesh_dimensions
!!   physicaldata
!!
!! CALLS
!!
!! RETURNS
!!
!!   Nothing returned.
!!
!! DESCRIPTION
!!
!! This routine resets some variables which manage the guardcell
!! filling operation when operating in 1-block mode. It should be called
!! from amr_initialize, and at a synchronization point which separates
!! guardcell filling of data with a given time stamp, and guardcell
!! filling at the next time stamp.
!!
!! AUTHORS
!!
!!   Peter MacNeice, February 1999.
!!
!! MODIFICATIONS
!!  2022-05-13 K. Weide  Small cleanup
!!  2022-05-20 K. Weide  Added "pcache_gcregions" variables
!!  2022-11-30 Klaus Weide  USE statements with updated ONLY lists.
!!***

      Subroutine amr_1blk_guardcell_reset

      Use physicaldata, ONLY: lnew_parent, &
           pcache_blk_u, pcache_pe_u, &
           pcache_blk_w, pcache_pe_w
      Use physicaldata, ONLY: pcache_gcregions_u, pcache_gcregions_w

      Implicit None

!-----Begin Executable Code

! reset id recording parent block which is currently in cache
          lnew_parent = .True.
          pcache_blk_u = -1
          pcache_pe_u  =  -1
          pcache_gcregions_u  =  -1
          pcache_blk_w = -1
          pcache_pe_w  =  -1
          pcache_gcregions_w  =  -1


      Return
      End Subroutine amr_1blk_guardcell_reset
