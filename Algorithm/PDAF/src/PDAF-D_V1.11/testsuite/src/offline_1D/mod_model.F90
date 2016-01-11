!$Id: mod_model.F90 786 2009-12-07 10:33:52Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_model

! !DESCRIPTION:
! This module provides shared variables for the model.
! This variant is for the offline mode of PDAF.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE

! !PUBLIC DATA MEMBERS:  
!    ! Control model run - available as command line options
  INTEGER :: dim_state           ! Model state dimension

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
  INTEGER, ALLOCATABLE :: local_dims(:)  ! Array for local state dimensions
!EOP

END MODULE mod_model
