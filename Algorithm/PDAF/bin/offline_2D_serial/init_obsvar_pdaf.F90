!$Id: init_obsvar_pdaf.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: init_obsvar_pdaf --- Get mean observation error variance
!
! !INTERFACE:
SUBROUTINE init_obsvar_pdaf(step, dim_obs_p, obs_p, meanvar)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF

! This routine will only be called, if the adaptive
! forgetting factor feature is used. Please note that
! this is an experimental feature.
!
! The routine is called in global filters (like SEIK)
! during the analysis or in local filters (e.g. LSEIK)
! before the loop over local analysis domains 
! by the routine PDAF\_set\_forget that estimates an 
! adaptive forgetting factor.  The routine has to 
! initialize the mean observation error variance.  
! For global filters this should be the global mean,
! while for local filters it should be the mean for the
! PE-local  sub-domain.  (See init\_obsvar\_l_pdaf()
! for a localized variant for local filters.)
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
       ONLY: rms_obs, obs_index, &
       STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, &
       XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, R_Local, H_Local, &
       FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, &
       XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step          ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p     ! PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) ! PE-local observation vector
  REAL, INTENT(out)   :: meanvar       ! Mean observation error variance

! !CALLING SEQUENCE:
! Called by: PDAF_set_forget    (as U_init_init_obs_covar)
!EOP

! *** local variables ***
  INTEGER :: i          ! counter

! *****************************
! *** Compute mean variance ***
! *****************************

  ! We assume that all observations have the same error.
  ! Thus, the mean variance is the error variance of each single observation.

  meanvar = 0.0
  DO i = 1, dim_obs_p
    meanvar = meanvar + R_Local(i,i) * (1.0 / REAL(dim_obs_p))
  END DO

END SUBROUTINE init_obsvar_pdaf
