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
!

      module constants

      implicit none
  
      real(kind=kind(1.0d0)), parameter ::  & 
     &        half=0.5d0,  & 
     &        zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0, & 
     &        pi=3.141592653589793d0

      end module constants
