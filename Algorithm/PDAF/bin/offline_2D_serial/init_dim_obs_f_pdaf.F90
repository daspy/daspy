!$Id: init_dim_obs_f_pdaf.F90 1421 2013-09-25 15:10:05Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_f_pdaf --- Set full dimension of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_f_pdaf(step, dim_obs_f)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: LSEIK/LETKF/LESTKF
    !
    ! The routine is called in PDAF\_lseik\_update
    ! at the beginning of the analysis step before
    ! the loop through all local analysis domains.
    ! It has to determine the dimension of the
    ! observation vector according to the current
    ! time step for all observations required for
    ! the analyses in the loop over all local
    ! analysis domains on the PE-local state domain.
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
        ONLY : nx, ny, obs, obs_index, coords_obs, &
        STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, Normal_Score_Trans, &
        XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, R_Local, H_Local, &
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

    use netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in)  :: step      ! Current time step
    INTEGER, INTENT(out) :: dim_obs_f ! Dimension of full observation vector

    ! !CALLING SEQUENCE:
    ! Called by: PDAF_lseik_update   (as U_init_dim_obs)
    ! Called by: PDAF_lestkf_update  (as U_init_dim_obs)
    ! Called by: PDAF_letkf_update   (as U_init_dim_obs)
    !EOP

    ! *** Local variables
    INTEGER :: i, j                     ! Counters
    INTEGER :: cnt, cnt0                ! Counters

    REAL, ALLOCATABLE :: HXF_NNscore(:), HXF_Trans(:)
    REAL, ALLOCATABLE :: Observation_NNscore(:), Observation_Trans(:)
	REAL, ALLOCATABLE :: Observation_NNscore_Perb(:), Observation_NNscore_Ens(:), Observation_NNscore_Ens_Array(:), Observation_Trans_Ens(:), Observation_Noise(:)
	REAL :: Inflation_Fac, hxf_var, HXF_Trans_Var, HXF_NNscore_Var, Observation_Trans_Var
	
    ! *********************************************
    ! *** Initialize full observation dimension ***
    ! *********************************************

    ! Count observations
    cnt = 0
    DO i = 1, STATE_DIM
        IF (OBS_NC(i) .ne. -9999.0) cnt = cnt + 1
    END DO

    ! Set number of observations
    dim_obs_f = cnt

    ! Initialize vector of observations and index array
    IF (ALLOCATED(obs_index)) DEALLOCATE(obs_index)
    IF (ALLOCATED(obs)) DEALLOCATE(obs)
    IF (ALLOCATED(coords_obs)) DEALLOCATE(coords_obs)
    ALLOCATE(obs_index(dim_obs_f))
    ALLOCATE(obs(dim_obs_f))
    ALLOCATE(coords_obs(2, dim_obs_f))

    IF (ALLOCATED(R_Local)) DEALLOCATE(R_Local)
    IF (ALLOCATED(H_Local)) DEALLOCATE(H_Local)
    ALLOCATE(R_Local(dim_obs_f,dim_obs_f))
    ALLOCATE(H_Local(dim_obs_f,STATE_DIM))
	
    ALLOCATE(hlmbda(minimize_lbfgsb_n))
	
	ALLOCATE(HXF_NNscore(cnt*ENSEMBLE_NUMBER))
    ALLOCATE(HXF_Trans(cnt*ENSEMBLE_NUMBER))
	
	ALLOCATE(Observation_NNscore(cnt))
	ALLOCATE(Observation_Trans(cnt))
	
	!ALLOCATE(Observation_NNscore_Perb(cnt))
	!ALLOCATE(Observation_NNscore_Ens_Array(cnt*ENSEMBLE_NUMBER))
	!ALLOCATE(Observation_Noise(ENSEMBLE_NUMBER))
	
    cnt = 0
    cnt0 = 0
    DO i = 1, STATE_DIM
        cnt0 = cnt0 + 1
        IF (OBS_NC(i) .ne. -9999.0) THEN
            cnt = cnt + 1
            obs_index(cnt) = cnt0
            IF (Normal_Score_Trans == 1) THEN
                HXF_NNscore(((cnt-1)*ENSEMBLE_NUMBER+1):cnt*ENSEMBLE_NUMBER) = HXF_NC(i,:)
                Observation_NNscore(cnt) = OBS_NC(i)
				!CALL normal(ENSEMBLE_NUMBER,0.0,sqrt(R_NC(i,i)),Observation_Noise,OBS_NC(i))
				!Observation_NNscore_Ens_Array(((cnt-1)*ENSEMBLE_NUMBER+1):cnt*ENSEMBLE_NUMBER) = OBS_NC(i) + Observation_Noise
				!print*,"Observation_Noise",Observation_Noise
				!Observation_NNscore_Perb(cnt) = min(Observation_NNscore(cnt), minval(Observation_NNscore_Ens_Array(((cnt-1)*ENSEMBLE_NUMBER+1):cnt*ENSEMBLE_NUMBER)))
			END IF
        END IF
    END DO
	
	cnt = 0
    cnt0 = 0
    DO i = 1, STATE_DIM
        cnt0 = cnt0 + 1

        IF (OBS_NC(i) .ne. -9999.0) THEN
            cnt = cnt + 1
            obs_index(cnt) = cnt0

            

            !modify the R according to the model variance
            CALL com_covar(ENSEMBLE_NUMBER,HXF_NC(i,:),HXF_NC(i,:),hxf_var)
            IF (R_NC(cnt,cnt) > hxf_var) THEN
                Inflation_Fac = R_NC(cnt,cnt) * max(0.1, min(1.0,(Alpha_Inflation * (hxf_var - R_NC(cnt,cnt)) / R_NC(cnt,cnt) + 1.0)))  !RTPS (multiplicative)
                R_NC(cnt,cnt) = Inflation_Fac
            END IF

            obs(cnt) = OBS_NC(i)
            R_Local(cnt,cnt) = R_NC(cnt,cnt)
            H_Local(cnt,:) = H_NC(cnt,:)
            coords_obs(1, cnt) = OBS_COORD_NC(cnt,1)
            coords_obs(2, cnt) = OBS_COORD_NC(cnt,2)
        END IF
    END DO
	DEALLOCATE(HXF_NNscore)
    DEALLOCATE(HXF_Trans)
	DEALLOCATE(Observation_NNscore)
	DEALLOCATE(Observation_Trans)
	!IF (Normal_Score_Trans == 1) THEN
	!	DEALLOCATE(Observation_NNscore_Perb)
	!	DEALLOCATE(Observation_NNscore_Ens_Array)
	!	DEALLOCATE(Observation_NNscore_Ens)
	!	DEALLOCATE(Observation_Trans_Ens)
	!END IF
    ! *** Clean up ***

    DEALLOCATE(H_NC)
    DEALLOCATE(OBS_NC)
    !DEALLOCATE(R_NC)
    DEALLOCATE(OBS_COORD_NC)

contains
    !-----------------------------------------------------------------------
    ! Mean
    !-----------------------------------------------------------------------
    SUBROUTINE com_mean(ndim,var,amean)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: ndim
        REAL,INTENT(IN) :: var(ndim)
        REAL,INTENT(OUT) :: amean

        INTEGER :: i

        amean = 0.0d0
        DO i=1,ndim
            amean = amean + var(i)
        END DO
        amean = amean / REAL(ndim,8)

        RETURN
    END SUBROUTINE com_mean

    !-----------------------------------------------------------------------
    ! Standard deviation
    !-----------------------------------------------------------------------
    SUBROUTINE com_stdev(ndim,var,aout)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: ndim
        REAL,INTENT(IN) :: var(ndim)
        REAL,INTENT(OUT) :: aout

        REAL :: amean
        REAL :: dev(ndim)

        CALL com_mean(ndim,var,amean)

        dev(:) = var(:) - amean

        aout = SQRT( SUM(dev*dev) / REAL(ndim-1,8) )

        RETURN
    END SUBROUTINE com_stdev
    !-----------------------------------------------------------------------
    ! Covariance
    !-----------------------------------------------------------------------
    SUBROUTINE com_covar(ndim,var1,var2,cov)
        IMPLICIT NONE

        INTEGER,INTENT(IN) :: ndim
        REAL(8),INTENT(IN) :: var1(ndim)
        REAL(8),INTENT(IN) :: var2(ndim)
        REAL(8),INTENT(OUT) :: cov

        REAL(8) :: amean1,amean2
        REAL(8) :: dev1(ndim),dev2(ndim)

        CALL com_mean(ndim,var1,amean1)
        CALL com_mean(ndim,var2,amean2)

        dev1(:) = var1(:) - amean1
        dev2(:) = var2(:) - amean2

        cov = SUM( dev1*dev2 ) / REAL(ndim-1,8)

        RETURN
    END SUBROUTINE com_covar

    SUBROUTINE normal(nbv,mean,sigma,array,seed_in) !returns a normal distribution
    implicit none
    INTEGER, intent(in) :: nbv
	real(8), intent(in) :: seed_in,mean,sigma
	  INTEGER :: i
	real(8), intent(out) :: array(nbv)
	 INTEGER :: seed(12)
	 real(8) :: pi, temp
	 
	  pi = 4.0*ATAN(1.0)
	  seed(:) = seed_in
      CALL RANDOM_SEED(PUT = seed)
	  CALL RANDOM_NUMBER(array) ! Uniform distribution
	 
	! Now convert to normal distribution
	  DO i = 1, nbv-1, 2
		temp = sigma * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
		array(i+1) = sigma * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
		array(i) = temp
	  END DO
    return
	end SUBROUTINE normal

    subroutine init_random_seed()
        implicit none
        integer, allocatable :: seed(:)
        integer :: i, n, un, istat, dt(8), pid, t(2), s
        integer(8) :: count, tms

        call random_seed(size = n)
        allocate(seed(n))
        ! First try if the OS provides a random number generator
        open(newunit=un, file="/dev/urandom", access="stream", &
            form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(un) seed
            close(un)
        else
            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            call system_clock(count)
            if (count /= 0) then
                t = transfer(count, t)
            else
                call date_and_time(values=dt)
                tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
                t = transfer(tms, t)
            end if
            s = ieor(t(1), t(2))
            pid = getpid() + 1099279 ! Add a prime
            s = ieor(s, pid)
            if (n >= 3) then
                seed(1) = t(1) + 36269
                seed(2) = t(2) + 72551
                seed(3) = pid
                if (n > 3) then
                    seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
                end if
            else
                seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
            end if
        end if
        call random_seed(put=seed)
    end subroutine init_random_seed


END SUBROUTINE init_dim_obs_f_pdaf

