!$Id: init_dim_obs_local.F90 1160 2011-09-14 09:32:08Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_local --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_local(domain, step, dim_obs, dim_obs_l)

! !DESCRIPTION:
! User-supplied routine for PDAF (LSEIK):
!
! The routine is called during the loop over
! all local analysis domains. It has to set 
! the dimension of the local observation vector 
! for the current local analsis domain.
!
! This variant is for the Lorenz96 model without
! parallelization. A local observation 
! domain is used that is defined by the cut-off 
! distance lseik\_range around the current grid
! point that is updated.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_model, &
       ONLY: dim_state
  USE mod_parallel, &
       ONLY: mype_filter
  USE mod_assimilation, &
       ONLY: local_range, local_range2

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: domain     ! Current local analysis domain
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(in)  :: dim_obs    ! Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
!EOP


! *** local variables ***
  INTEGER :: i                       ! Counter
  INTEGER :: domain_g                ! Global domain index


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  dim_obs_l = 1 + local_range + local_range2

  ! local dimension is larger than state dimension:
  ! reset to state_dimension
  IF (dim_obs_l > dim_state) dim_obs_l = dim_state

END SUBROUTINE init_dim_obs_local

