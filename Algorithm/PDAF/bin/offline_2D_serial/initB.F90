!$Id: initB.F90 1450 2013-12-12 09:30:00Z lnerger $
!BOP
!
! !ROUTINE: initB - Initialize localization matrix B
!
! !INTERFACE:
SUBROUTINE initB(locweight)

! !DESCRIPTION:
! User-supplied routine for PDAF (constant forecast error covariance) 
!
! This implementation is still suboptimal. We initialize a matrix of
! size dim_state x dim_state. While this works for the Lorenz96 model,,
! it should not be done for large-scale models.
!
! !REVISION HISTORY:
! 2010 - Tijana Janjic
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: local_range, B, B_step, srange, obs_index, coords_obs, coords_l, &
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
		
  USE mod_parallel, &
       ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: locweight  ! Type of localization
!EOP

! *** local variables ***
  INTEGER :: i, j          ! Index of observation component
  REAL    :: distance      ! Distance between points in the domain 
  REAL    :: halfrange     ! Half range
  INTEGER :: dim_half      ! Size of first half of state vector
  REAL    :: spp_r

	
! **********************
! *** INITIALIZATION ***
! **********************
  ALLOCATE(B(STATE_DIM,STATE_DIM))
  ALLOCATE(B_step(STATE_DIM,STATE_DIM))

  spp_r = REAL(local_range)

  IF (locweight == 1) THEN
     WRITE(*,'(9x,a,f10.2)') &
          '--- LOCALIZATION: Initialize B with support radius: ', spp_r
     WRITE(*,'(9x,a)') &
          '--- LOCALIZATION: Use 5th-order polynomial'
!      IF (srange < spp_r) THEN
!         WRITE(*,'(9x,a,f10.2)') &
!              '--- LOCALIZATION: Use cut-off radius: ', srange
!      END IF
  ELSE IF (locweight == 2) THEN
     WRITE(*,'(9x,a,f10.2)') &
          '--- LOCALIZATION: Initialize B with support radius: ', spp_r
     WRITE(*,'(9x,a)') &
          '--- LOCALIZATION: Use unit weight up to support radius'
  ELSE
     WRITE(*,'(9x,a)') &
          '--- LOCALIZATION: Initialize B with unit weight'
!      WRITE(*,'(9x,a,f10.2)') &
!           '--- LOCALIZATION: Use cut-off radius: ', srange
  END IF


  IF (locweight==2) THEN
     ! Size of first half of state vector
     dim_half = FLOOR(REAL(STATE_DIM)/2.0)+1

     ! Initialize fist column of B
     B = 0.0
     B_step = 0.0
     firsthalf: DO i = 1, dim_half
        ! *** Compute localizing weighting term according to
        ! *** distance of the observation from current grid point
		
		
		print*,"obs_index(i)",obs_index(i)
		print*,"coords_obs(1, obs_index(i))",coords_obs(:, obs_index(i))
		print*,"distance",distance

        ! distance between analysis point and current observation
		coords_l(1) = XF_COORD_NC(i,1)
		coords_l(2) = XF_COORD_NC(i,2)
		print*,"coords_l",coords_l
		distance = SQRT(REAL((coords_l(1) - coords_obs(1, obs_index(i)))**2 + &
            (coords_l(2) - coords_obs(2, obs_index(i)))**2))
		
		
        halfrange = spp_r / 2.0
		print*,"halfrange",halfrange
        
        ! Compute weight
!         cutoff: IF (distance <= srange) THEN
           distances: IF (distance <= halfrange) THEN
              B(i,1) = -0.25 * (distance / halfrange)**5 &
                   + 0.5 * (distance / halfrange)**4 &
                   + 5.0 / 8.0 * (distance / halfrange)**3 &
                   - 5.0 / 3.0 * (distance / halfrange)**2 + 1.0

           ELSEIF (distance > halfrange .AND. distance < spp_r) THEN
              B(i,1) = 1.0 / 12.0 * (distance / halfrange)**5 &
                   - 0.5 * (distance / halfrange)**4 &
                   + 5.0 / 8.0 * (distance / halfrange)**3 &
                   + 5.0 / 3.0 * (distance / halfrange)**2 &
                   - 5.0 * (distance / halfrange) &
                   + 4.0 - 2.0 / 3.0 * halfrange / distance

              ! Check if weight >0 (Could be <0 due to numerical precision)
              IF (B(i, 1) < 0.0) THEN
                 B(i, 1) = 0.0
              END IF
           ELSE
              B(i, 1)=0.0
           ENDIF distances
           IF (distance < spp_r) THEN
              B_step(i,1) = 1.0
           END IF
!         ELSE
!            B(i, 1)=0.0
!         END IF cutoff
     END DO firsthalf

     ! Initialize second half by periodic BC
     secondhalf: DO i = dim_half+1, STATE_DIM
        B(i, 1) = B(STATE_DIM-i+2, 1) 
        B_step(i, 1) = B_step(STATE_DIM-i+2, 1) 
     END DO secondhalf


     ! *** Fill columns 2 to STATE_DIM

     ! lower triangle and diagonal of matrix
     DO j = 2, STATE_DIM
        DO i = j, STATE_DIM
           B(i, j) = B(i-j+1, 1)
           B_step(i, j) = B_step(i-j+1, 1)
        END DO
     END DO

     ! upper triangle of matrix
     DO j = 2, STATE_DIM
        DO i = 1, j - 1
           B(i,j) = B(STATE_DIM-j+i+1, 1)
           B_step(i,j) = B_step(STATE_DIM-j+i+1, 1)
        END DO
     END DO

  ELSE IF (locweight == 1) THEN
     ! Unitialize with unit weight up to support radius

     ! Size of first half of state vector
     dim_half = FLOOR(REAL(STATE_DIM)/2.0)+1

     ! Initialize fist column of B
     B = 0.0
     firsthalf2: DO i = 1, dim_half
        ! *** Compute localizing weighting term according to
        ! *** distance of the observation from current grid point
        
        ! Compute weight
        distances2: IF (distance < spp_r) THEN
           B(i, 1) = 1.0
        ELSE
           B(i, 1) = 0.0
        ENDIF distances2
     END DO firsthalf2

     ! Initialize second half by periodic BC
     secondhalf2: DO i = dim_half+1, STATE_DIM
        B(i, 1) = B(STATE_DIM-i+2, 1) 
     END DO secondhalf2


     ! *** Fill columns 2 to STATE_DIM

     ! lower triangle and diagonal of matrix
     DO j = 2, STATE_DIM
        DO i = j, STATE_DIM
           B(i, j) = B(i-j+1, 1)
        END DO
     END DO

     ! upper triangle of matrix
     DO j = 2, STATE_DIM
        DO i = 1, j - 1
           B(i,j) = B(STATE_DIM-j+i+1, 1)
        END DO
     END DO

     B_step = B

  ELSE
     B = 1.0
     B_step = B
  END IF



END SUBROUTINE initB 
