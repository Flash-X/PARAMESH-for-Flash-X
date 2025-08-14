!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!     Klaus Weide, May 2025, *io.F90*
!        Simplify
!

      Module io

        Integer, Save :: iu_log = 6

        Public :: output_dir, amr_log_file
        Character (Len=80) :: output_dir
        Character (Len=80) :: amr_log_file

      End Module io
