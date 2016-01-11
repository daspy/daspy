!$Id: obs_op_f_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update
! before the loop over all local analysis domains
! is entered.  The routine has to perform the
! operation of the observation operator acting on
! a state vector.  The full vector of all
! observations required for the localized analysis
! on the PE-local domain has to be initialized.
! This is usually data on the PE-local domain plus
! some region surrounding the PE-local domain.
! This data is gathered by MPI operations. The
! gathering has to be done here, since in the loop
! through all local analysis domains, no global
! MPI operations can be performed, because the
! number of local analysis domains can vary from
! PE to PE.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index, &
       STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, &
       XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, R_Local, H_Local, &
       FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, &
       XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
! Called by: PDAF_lestkf_update  (as U_obs_op)
! Called by: PDAF_letkf_update   (as U_obs_op)
!EOP

! *** local variables ***
  INTEGER :: i       ! Counter
  INTEGER :: member

! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  !print*,"dim_obs_f",dim_obs_f
  DO i = 1, dim_obs_f
     CALL PDAF_get_obsmemberid(member)
     IF (member > 0) THEN
         m_state_f(i) = DOT_PRODUCT(H_Local(i,:),HXF_NC(:,member))
         !print*,"m_state_f",i,member,H_Local(i,:),HXF_NC(:,member)
     END IF
  END DO

END SUBROUTINE obs_op_f_pdaf
