!$Id: localB_seq.F90 1450 2013-12-12 09:30:00Z lnerger $
!BOP
!
! !ROUTINE: localB --- apply localization matrix B
!
! !INTERFACE:
SUBROUTINE localB(dim, dim_obs, HP, HPH)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF (EnSKF)
    !
    ! This routine applies a localization matrix B
    ! to the matrices HP and HPH^T of the EnSKF.
    !
    ! !REVISION HISTORY:
    ! 2010-03 - Lars Nerger
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_assimilation, &
        ONLY: B, B_step, local_range, locweight, use_obs_mask, obs_index, coords_obs, coords_l, &
        STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, Normal_Score_Trans, &
        XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, &
        FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, Normal_Score_Trans_NAME, &
        XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME, &
        State_DIM_Single_Layer, Parameter_DIM, Par_Uniform_STD, Alpha_Inflation, &
        State_DIM_Single_Layer_NAME, Parameter_DIM_NAME, Par_Uniform_STD_NAME, Alpha_Inflation_NAME, &
        Parameter_Optimization_Flag, Parameter_Optimization_Flag_NAME, &
        State_DIM_Single_Column, State_DIM_Single_Column_NAME, &
        X_Left, X_Right, Y_Lower, Y_Upper, &
        lmbda_DIM, minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint, &
        minimize_lbfgsb_factr, minimize_lbfgsb_pgtol, &
        lmbda_state, lmbda_parameter, lmbda_bias, hlmbda, minimize_lbfgsb_epsilon_in

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: dim     ! State dimension
    INTEGER, INTENT(in) :: dim_obs ! number of observations
    REAL, INTENT(inout) :: HP(dim_obs, dim) ! Matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs)     ! Matrix HPH

    ! *** local variables ***
    INTEGER :: i, j          ! Index of observation component
    INTEGER :: ilow, iup     ! Lower and upper bounds of observation domain
    REAL    :: ivariance_obs ! Inverse of variance of the observations
    REAL    :: distance      ! Distance between points in the domain
    REAL    :: weight        ! Localization weight
    REAL    :: cfaci         ! Half-weight distance
	REAL    :: halfrange     ! Half range
	REAL    :: spp_r, HP_Temp
	
    ! **********************
    ! *** INITIALIZATION ***
    ! **********************
	spp_r = REAL(local_range)
	use_obs_mask = .TRUE.
	
    IF (locweight == 2) THEN

        WRITE (*,'(8x,a,f6.1)') &
            '---   Apply localization, support radus ', local_range

        IF (.NOT. use_obs_mask) THEN
            ! This is for global obervations
            DO j = 1, dim_obs
                DO i = 1, dim
                    HP(j,i) = HP(j,i) * weight
                END DO
            END DO
            DO j = 1, dim_obs
                DO i = 1, dim_obs
                    HPH(j,i) = HPH(j,i) * weight
                END DO
            END DO
        ELSE
            ! This is for gappy obervations
			DO j = 1, dim_obs
				DO i = 1, dim
                    !print*,"obs_index(j)",obs_index(j)
                    !print*,"coords_obs(:, j)",coords_obs(:, j)
			
                    ! distance between analysis point and current observation
                    coords_l(1) = XF_COORD_NC(i,1)
                    coords_l(2) = XF_COORD_NC(i,2)
                    !print*,"coords_l",coords_l
                    distance = SQRT(REAL((coords_l(1) - coords_obs(1, j))**2 + &
                        (coords_l(2) - coords_obs(2, j))**2))
                    !print*,"i",i,"j",j,"distance",distance

					!IF(distance < abs(GridSize_Sys)) THEN
					!	continue
					!	weight = 1.0
					!END IF
					
                    halfrange = spp_r / 2.0
                    !print*,"halfrange",halfrange
			
                    ! Compute weight
                    !         cutoff: IF (distance <= srange) THEN
                    IF (distance <= halfrange) THEN
                        weight = -0.25 * (distance / halfrange)**5 &
                            + 0.5 * (distance / halfrange)**4 &
                            + 5.0 / 8.0 * (distance / halfrange)**3 &
                            - 5.0 / 3.0 * (distance / halfrange)**2 + 1.0

                    ELSEIF (distance > halfrange .AND. distance < spp_r) THEN
                        weight = 1.0 / 12.0 * (distance / halfrange)**5 &
                            - 0.5 * (distance / halfrange)**4 &
                            + 5.0 / 8.0 * (distance / halfrange)**3 &
                            + 5.0 / 3.0 * (distance / halfrange)**2 &
                            - 5.0 * (distance / halfrange) &
                            + 4.0 - 2.0 / 3.0 * halfrange / distance

                        ! Check if weight >0 (Could be <0 due to numerical precision)
                        IF (weight < 0.0) THEN
                            weight = 0.0
                        END IF
                    ELSE
                        weight=0.0
                    ENDIF
					HP(j,i) = HP(j,i) * weight
                END DO
            END DO
			
			! Only assimilate 1 obs for the model grid cells with observation, re-assign HP
			!DO j = 1, dim_obs
			!	HP_Temp = HP(j,obs_index(j))
			!	HP(:,obs_index(j)) = 0.0
			!	HP(j,obs_index(j)) = HP_Temp
			!END DO
			! Only assimilate 1 obs for the model grid cells with observation
			
			
            DO j = 1, dim_obs
                DO i = 1, dim_obs
                    distance = SQRT(REAL((coords_obs(1, j) - coords_obs(1, i))**2 + &
                        (coords_obs(2, j) - coords_obs(2, i))**2))
                    !print*,"j",j,"i",i,"distance",distance
			
                    halfrange = spp_r / 2.0
                    !print*,"halfrange",halfrange
			
                    ! Compute weight
                    !         cutoff: IF (distance <= srange) THEN
                    IF (distance <= halfrange) THEN
                        weight = -0.25 * (distance / halfrange)**5 &
                            + 0.5 * (distance / halfrange)**4 &
                            + 5.0 / 8.0 * (distance / halfrange)**3 &
                            - 5.0 / 3.0 * (distance / halfrange)**2 + 1.0

                    ELSEIF (distance > halfrange .AND. distance < spp_r) THEN
                        weight = 1.0 / 12.0 * (distance / halfrange)**5 &
                            - 0.5 * (distance / halfrange)**4 &
                            + 5.0 / 8.0 * (distance / halfrange)**3 &
                            + 5.0 / 3.0 * (distance / halfrange)**2 &
                            - 5.0 * (distance / halfrange) &
                            + 4.0 - 2.0 / 3.0 * halfrange / distance

                        ! Check if weight >0 (Could be <0 due to numerical precision)
                        IF (weight < 0.0) THEN
                            weight = 0.0
                        END IF
                    ELSE
                        weight=0.0
                    ENDIF
                    HPH(j,i) = HPH(j,i) * weight
                END DO
            END DO
        END IF
    END IF

END SUBROUTINE localB
