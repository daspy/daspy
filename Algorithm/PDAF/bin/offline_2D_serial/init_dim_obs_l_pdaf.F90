!$Id: init_dim_obs_l_pdaf.F90 1399 2013-05-06 09:21:15Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_l_pdaf --- Set dimension of local observation vector
!
! !INTERFACE:
SUBROUTINE init_dim_obs_l_pdaf(domain_p, step, dim_obs_f, dim_obs_l)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: LSEIK/LETKF/LESTKF
    !
    ! The routine is called during the loop over
    ! all local analysis domains. It has to set
    ! the dimension of the local observation vector
    ! for the current local analysis domain.
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
        ONLY: nx, ny, obs, local_range, coords_obs, coords_l, obs_index_l, &
        STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, &
        XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, R_Local, R_Local_l, coords_obs_l, &
        FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, &
        XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME, &
        X_Left, X_Right, Y_Lower, Y_Upper, GridSize_Sys, GridSize_Obs

    use netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in)  :: domain_p   ! Current local analysis domain
    INTEGER, INTENT(in)  :: step       ! Current time step
    INTEGER, INTENT(in)  :: dim_obs_f  ! Full dimension of observation vector
    INTEGER, INTENT(out) :: dim_obs_l  ! Local dimension of observation vector

    ! !CALLING SEQUENCE:
    ! Called by: PDAF_lseik_update   (as U_init_dim_obs_l)
    ! Called by: PDAF_lestkf_update  (as U_init_dim_obs_l)
    ! Called by: PDAF_letkf_update   (as U_init_dim_obs_l)
    !EOP


    ! *** local variables ***
    INTEGER :: i, cnt                   ! Counters
    INTEGER :: limits_x(2), limits_y(2) ! Coordinate limites for observation domain
    REAL :: distance                    ! Distance between observation and analysis domain

    INTEGER :: num_local_obs_temp
    REAL :: dist(dim_obs_f) !Distance between all observations and analysis domain
    INTEGER, ALLOCATABLE :: localobs(:)
    REAL(8), ALLOCATABLE :: dist_copy(:)
    INTEGER ::coeffs_whole_index
    INTEGER :: seed(12), NINT_Correlation

    NINT_Correlation = NINT(local_range/REAL(abs(GridSize_Sys),8))
    NINT_Correlation = max(NINT_Correlation,1)

    ! Comment this line, because it does not work under parallelpython
    !print*,"NINT_Correlation",NINT_Correlation,local_range,REAL(abs(GridSize_Sys),8)


    ! **********************************************
    ! *** Initialize local observation dimension ***
    ! **********************************************

    ! Coordinates of local analysis domain
    coords_l(1) = XF_COORD_NC(domain_p,1)
    coords_l(2) = XF_COORD_NC(domain_p,2)

    !Determine coordinate limits for observation domain
    limits_x(1) = coords_l(1) - local_range
    if (limits_x(1) < X_Left) limits_x(1) = X_Left
    limits_x(2) = coords_l(1) + local_range
    if (limits_x(2) > X_Right) limits_x(2) = X_Right

    limits_y(1) = coords_l(2) - local_range
    if (limits_y(1) < Y_Lower) limits_y(1) = Y_Lower
    limits_y(2) = coords_l(2) + local_range
    if (limits_y(2) > Y_Upper) limits_y(2) = Y_Upper

    ! Count observations within local_range
    dim_obs_l = 0
    DO i = 1, dim_obs_f
        IF (obs(i) .ne. -9999.0) THEN
            IF (coords_obs(1, i) >= limits_x(1) .AND. coords_obs(1, i) <= limits_x(2) .AND. &
                coords_obs(2, i) >= limits_y(1) .AND. coords_obs(2, i) <= limits_y(2)) THEN
                distance = SQRT(REAL((coords_l(1) - coords_obs(1,i))**2 + &
                    (coords_l(2) - coords_obs(2,i))**2))
                dist(i) = distance
                IF (distance <= REAL(local_range)) THEN
                    dim_obs_l = dim_obs_l + 1
                END IF
            END IF
        END IF
    END DO

    ! Find the Closest Obs for Model Grids with Observations
    ALLOCATE(dist_copy(dim_obs_f))
    ALLOCATE(localobs(dim_obs_f))

    localobs(:) = 0
    dist_copy(:) = dist

    IF (size(pack(dist,dist < abs(GridSize_Sys))) > 0) THEN
        num_local_obs_temp = 1
    ELSE
        num_local_obs_temp = dim_obs_l
    END IF

    IF (num_local_obs_temp > 0) THEN
        DO WHILE (count(localobs /= 0) < num_local_obs_temp)
            coeffs_whole_index = minloc(dist_copy(:), DIM=1, MASK = dist_copy < NINT_Correlation*abs(GridSize_Sys))
            localobs(coeffs_whole_index) = 1
            dist_copy(coeffs_whole_index) = NINT_Correlation*abs(GridSize_Sys)
        END DO
    END IF

    dim_obs_l = num_local_obs_temp
    !print*,"num_local_obs_temp",num_local_obs_temp
    !print*,"localobs",localobs

    DEALLOCATE(dist_copy)
    DEALLOCATE(localobs)
    ! Stop Find the Closest Obs for Model Grids with Observations

    ! Initialize index array for local observations in full observed vector
    IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
    ALLOCATE(obs_index_l(dim_obs_l))

    IF (ALLOCATED(R_Local_l)) DEALLOCATE(R_Local_l)
    IF (ALLOCATED(coords_obs_l)) DEALLOCATE(coords_obs_l)
    ALLOCATE(R_Local_l(dim_obs_l,dim_obs_l))
    ALLOCATE(coords_obs_l(2,dim_obs_l))

    cnt = 0
    !!!!!!!!!!!!!!!!!!!!!!! For Model Grids with Observation, we only use 1 local obs
    IF (dim_obs_l == 1) THEN
        cnt = cnt + 1
        obs_index_l(cnt) = coeffs_whole_index
        R_Local_l(cnt,cnt) = R_NC(coeffs_whole_index,coeffs_whole_index)
        coords_obs_l(:,cnt) = coords_obs(:, coeffs_whole_index)
    ELSE
        DO i = 1, dim_obs_f
            IF (obs(i) .ne. -9999.0) THEN
                IF (coords_obs(1, i) >= limits_x(1) .AND. coords_obs(1, i) <= limits_x(2) .AND. &
                    coords_obs(2, i) >= limits_y(1) .AND. coords_obs(2, i) <= limits_y(2)) THEN
                    distance = SQRT(REAL((coords_l(1) - coords_obs(1,i))**2 + &
                        (coords_l(2) - coords_obs(2,i))**2))

                    IF (distance <= REAL(local_range)) THEN
                        cnt = cnt + 1
                        obs_index_l(cnt) = i
                        R_Local_l(cnt,cnt) = R_NC(i,i)
                        coords_obs_l(:,cnt) = coords_obs(:, i)
                    END IF
                END IF
            END IF
        END DO
    END IF
    !print*,"obs_index_l",obs_index_l,"dim_obs_f",dim_obs_f,"local_range",local_range,"limits_x",limits_x,"limits_y",limits_y
    !print*,"coords_l",coords_l,"coords_obs",coords_obs

contains
    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

END SUBROUTINE init_dim_obs_l_pdaf

