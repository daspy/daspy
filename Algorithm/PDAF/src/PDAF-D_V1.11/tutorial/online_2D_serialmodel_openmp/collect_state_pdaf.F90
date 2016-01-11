!$Id: collect_state_pdaf.F90 1409 2013-09-25 11:47:03Z lnerger $
!BOP
!
! !ROUTINE: collect_state_pdaf --- Initialize state vector from model fields
!
! !INTERFACE:
SUBROUTINE collect_state_pdaf(dim_p, state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! For the offline mode of PDAF this routine only needs
! to exist for linking. It is never called.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code based on online implementation
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: nx, ny, field

  IMPLICIT NONE
  
! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! local state vector

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_X    (as U_coll_state)
! Called by: PDAF_assimilate_X   (as U_coll_state)
!EOP

! *** local variables ***
  INTEGER :: i, j         ! Counters
  

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

  DO j = 1, nx
     state_p(1 + (j-1)*ny : j*ny) = field(1:ny, j)
  END DO
  
END SUBROUTINE collect_state_pdaf
