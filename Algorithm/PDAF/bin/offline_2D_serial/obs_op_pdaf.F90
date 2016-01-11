!$Id: obs_op_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to
! provide the observed sub-state for the PE-local
! domain.
!
! Implementation for the 2D offline example
! with or without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs, obs_index, &
       STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, &
       XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, R_Local, H_Local, &
       FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, &
       XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***
  INTEGER :: i,j       ! Counter
  INTEGER :: member
  REAL    :: invdimens

! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  !print*,"dim_obs_p",dim_obs_p
  DO i = 1, dim_obs_p
     m_state_p(i) = 0
     CALL PDAF_get_obsmemberid(member)
     IF (member > 0) THEN
         m_state_p(i) = DOT_PRODUCT(H_Local(i,:),HXF_NC(:,member))
         !print*,"m_state_p",i,member,m_state_p(i),H_Local(i,:),HXF_NC(:,member)
     ELSE IF (member == 0) THEN ! For ETKF analysis
        invdimens = 1.0 / REAL(ENSEMBLE_NUMBER)
        DO j=1,ENSEMBLE_NUMBER
            !print*,member,i,j,m_state_p(i),DOT_PRODUCT(H_Local(i,:),HXF_NC(:,j))
            m_state_p(i) = m_state_p(i) + invdimens * DOT_PRODUCT(H_Local(i,:),HXF_NC(:,j))
        END DO
        !IF(m_state_p(i) == 0) THEN
        !    print*,"H_Local(i,:)",maxval(H_Local(i,:))
        !    print*,"HXF_NC(:,j)",maxval(HXF_NC(:,j))
        !END IF
     END IF
  END DO

END SUBROUTINE obs_op_pdaf
