!$Id: prepoststep_ens_offline.F90 1399 2013-05-06 09:21:15Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_offline --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
    state_p, Uinv, ens_p, flag)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! The routine is called for global filters (e.g. SEIK)
    ! before the analysis and after the ensemble transformation.
    ! For local filters (e.g. LSEIK) the routine is called
    ! before and after the loop over all local analysis
    ! domains.
    ! The routine provides full access to the state
    ! estimate and the state ensemble to the user.
    ! Thus, user-controlled pre- and poststep
    ! operations can be performed here. For example
    ! the forecast and the analysis states and ensemble
    ! covariance matrix can be analized, e.g. by
    ! computing the estimated variances.
    ! For the offline mode, this routine is the place
    ! in which the writing of the analysis ensemble
    ! can be performed.
    !
    ! If a user considers to perform adjustments to the
    ! estimates (e.g. for balances), this routine is
    ! the right place for it.
    !
    ! Implementation for the 2D offline example
    ! without parallelization.
    !
    ! !REVISION HISTORY:
    ! 2013-02 - Lars Nerger - Initial code based on offline_1D
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_assimilation, &
        ONLY: nx, ny, incremental, filename, subtype, covartype, &
        STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, Normal_Score_Trans, &
        XF_NC, XF_NC_Copy, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, &
        FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, Normal_Score_Trans_NAME, &
        XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME, &
        State_DIM_Single_Layer, Parameter_DIM, Par_Sens_Dim, Par_Uniform_STD, Alpha_Inflation, &
        State_DIM_Single_Layer_NAME, Parameter_DIM_NAME, Par_Uniform_STD_NAME, Alpha_Inflation_NAME, &
        Parameter_Optimization_Flag, Parameter_Optimization_Flag_NAME, &
        State_DIM_Single_Column, State_DIM_Single_Column_NAME, &
        X_Left, X_Right, Y_Lower, Y_Upper, &
        Bias_Model_Dim, Bias_Obs_Dim, Bias_Model_Uniform_STD,Bias_Model_Uniform_STD_NAME, &
        Bias_Obs_Uniform_STD, Bias_Obs_Uniform_STD_NAME, &
        Bias_Forecast_Model_Option, Bias_Forecast_Model_Option_NAME, &
        Bias_Observation_Model_Option, Bias_Observation_Model_Option_NAME, &
        Model_Inflation_Uniform_STD, Model_Inflation_Uniform_STD_NAME, &
        lmbda_DIM, minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint, &
        minimize_lbfgsb_factr, minimize_lbfgsb_pgtol, &
        lmbda_state, lmbda_parameter, lmbda_bias, hlmbda, minimize_lbfgsb_epsilon_in, &
        Par_Sens_Dim, State_DIM_Single_Column, Bias_Model_Dim, Bias_Obs_Dim

    use netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: step        ! Current time step (not relevant for offline mode)
    INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
    INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
    INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
    REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
    ! The array 'state_p' is not generally not initialized in the case of SEIK.
    ! It can be used freely here.
    REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
    REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
    INTEGER, INTENT(in) :: flag        ! PDAF status flag

    ! !CALLING SEQUENCE:
    ! Called by: PDAF_get_state      (as U_prepoststep)
    ! Called by: PDAF_seik_update    (as U_prepoststep)
    ! Called by: PDAF_lseik_update    (as U_prepoststep)
    ! Calls: PDAF_add_increment
    ! Calls: PDAF_seik_TtimesA
    ! Calls: memcount
    ! Calls: dgemm (BLAS)
    ! Calls: dgesv (LAPACK)
    ! Calls: MPI_send
    ! Calls: MPI_recv
    !EOP

    ! *** local variables ***
    INTEGER :: i, j, member         ! counters
    INTEGER, SAVE :: allocflag = 0       ! Flag for memory counting
    LOGICAL, SAVE :: firstio = .TRUE.    ! File output is peformed for first time?
    LOGICAL, SAVE :: firsttime = .TRUE.    ! Routine is called for first time?
    REAL :: invdim_ens                   ! Inverse ensemble size
    REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
    REAL :: rmserror_est                 ! estimated RMS error
    REAL, ALLOCATABLE :: variance(:)     ! model state variances
    REAL, ALLOCATABLE :: field(:,:)     ! global model field
    CHARACTER(len=2) :: ensstr          ! String for ensemble member

    ! For Normal Score
    REAL :: XF_NNscore(ENSEMBLE_NUMBER), XF_Trans(ENSEMBLE_NUMBER)
    REAL :: Analysis_NNscore(ENSEMBLE_NUMBER), Analysis_Trans(ENSEMBLE_NUMBER)

    ! For Inflation
    REAL :: dxf_dev(dim_p), dxa_dev(dim_p), xf_xm, xa_xm, dxa_min, dxa_max, xa_median
    REAL :: dxf(dim_p,ENSEMBLE_NUMBER), dxa(dim_p,ENSEMBLE_NUMBER), U(ENSEMBLE_NUMBER)
    INTEGER :: Par_STD_Index, Bias_Model_Dim_LETKF, Bias_Obs_Dim_LETKF
    INTEGER :: seed(12)

    ! because we defined the default values of Bias_Model_Dim and Bias_Obs_Dim are non-zeros, so...
    IF (Bias_Forecast_Model_Option > 0) THEN
        Bias_Model_Dim_LETKF = Bias_Model_Dim
    ELSE
        Bias_Model_Dim_LETKF = 0
    END IF

    IF (Bias_Observation_Model_Option > 0) THEN
        Bias_Obs_Dim_LETKF = Bias_Obs_Dim
    ELSE
        Bias_Obs_Dim_LETKF = 0
    END IF

    seed(:) = 123456
    CALL RANDOM_SEED(PUT = seed)
    CALL RANDOM_NUMBER(U)
    U = 0.8 * (1.0 - U) + 1.2 * U
    !print*,"U",U

    ! **********************
    ! *** INITIALIZATION ***
    ! **********************

    IF (firsttime) THEN
        WRITE (*, '(8x, a)') 'Analize forecasted state ensemble'
    ELSE
        WRITE (*, '(8x, a)') 'Analize and write assimilated state ensemble'
    END IF

    ! Allocate fields
    ALLOCATE(variance(dim_p))

    ! Initialize numbers
    rmserror_est  = 0.0
    invdim_ens    = 1.0 / REAL(dim_ens)
    invdim_ensm1  = 1.0 / REAL(dim_ens - 1)


    ! **************************************************************
    ! *** Perform prepoststep for SEIK with re-inititialization. ***
    ! *** The state and error information is completely in the   ***
    ! *** ensemble.                                              ***
    ! *** Also performed for SEIK without re-init at the initial ***
    ! *** time.                                                  ***
    ! **************************************************************

    IF (.not. firsttime) THEN

        DO i = 1, dim_p

            !print*,"ens_p",i,"Row",ens_p(i,:)
            !print*,"XF_NC_Copy",i,"Row",XF_NC_Copy(i,:)

            ! Inflation
            CALL com_mean(ENSEMBLE_NUMBER,XF_NC_Copy(i,:),xf_xm)
            CALL com_mean(ENSEMBLE_NUMBER,ens_p(i,:),xa_xm)
            dxa(i,:) = ens_p(i,:) - xa_xm
            dxf(i,:) = XF_NC_Copy(i,:) - xf_xm

            !print*,"minval(dxf)",minval(dxf),"maxval(dxf)",maxval(dxf)
            !print*,"minval(dxa)",minval(dxa),"maxval(dxa)",maxval(dxa)

            CALL com_stdev(ENSEMBLE_NUMBER,dxf(i,:),dxf_dev(i))
            CALL com_stdev(ENSEMBLE_NUMBER,dxa(i,:),dxa_dev(i))

            !print*,"Parameter_Optimization_Flag",Parameter_Optimization_Flag,Alpha_Inflation,dxf_dev(i),dxa_dev(i)
            !print*,i,nx,Bias_Forecast_Model_Option,Bias_Observation_Model_Option,Bias_Model_Dim,Bias_Obs_Dim
            ! Parameter Inflation
            !IF((Parameter_Optimization_Flag>0) .AND. (Alpha_Inflation .NE. 0.0) .AND. (i>State_DIM_Single_Layer) .AND. (dxa_dev(i) .GT. 0.0) .AND. (dxa_dev(i) < dxf_dev(i))) THEN
            IF((Parameter_Optimization_Flag>0) .AND. (Alpha_Inflation .NE. 0.0) .AND. &
                (Par_Sens_Dim >= 1) .AND. (i>State_DIM_Single_Layer)) THEN
                !print*,"inflation",i,(i-State_DIM_Single_Column*State_DIM_Single_Layer)/State_DIM_Single_Layer+1,State_DIM_Single_Column

                IF ((State_DIM_Single_Layer) > 1 .AND. (i>State_DIM_Single_Column*State_DIM_Single_Layer)) THEN
                    !Par_STD_Index = (i-State_DIM_Single_Column*State_DIM_Single_Layer-1)/State_DIM_Single_Layer+1
                    Par_STD_Index = (i-1-State_DIM_Single_Column*State_DIM_Single_Layer)/State_DIM_Single_Layer+1
                ELSE
                    Par_STD_Index = i-State_DIM_Single_Column
                END iF

                !Inflation
                !print*,i,"Par_STD_Index",Par_STD_Index,"dxa_dev(i)",dxa_dev(i),"Par_Uniform_STD(Par_STD_Index)",Par_Uniform_STD(Par_STD_Index)
                !print*,ens_p(i,:)
                IF (dxa_dev(i) < Par_Uniform_STD(Par_STD_Index)) THEN
                    IF (dxa_dev(i) .GT. 0.0) THEN
                        ens_p(i,:) = xa_xm + dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * &
                            (Par_Uniform_STD(Par_STD_Index) - dxa_dev(i)) / dxa_dev(i) + 1.0)))  !RTPS (multiplicative)
                    ELSE
                        seed(:) = 123456
                        CALL RANDOM_SEED(PUT = seed)
                        CALL RANDOM_NUMBER(U)
                        U = 0.8 * (1.0 - U) + 1.2 * U
                        ens_p(i,:) = ens_p(i,:) * U
                    END IF
                    !print*,xa_xm,i,State_DIM_Single_Layer,Par_Uniform_STD(mod(i,State_DIM_Single_Layer))
                    !ens_p(i,:) = ens_p(i,:) + dxa(i,:) * (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(i)) / dxa_dev(i))  !RTPS (multiplicative)
                    !ens_p(i,:) = ens_p(i,:) * (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(i)) / dxa_dev(i) + 1.0)  !RTPS (multiplicative)
                   !ens_p(i,:) = xa_xm + dxa(i,:) * (1.0 - Alpha_Inflation) + Alpha_Inflation * dxf(i,:)    !RTPP (multiplicative+additave)
                   !print*,ens_p(i,:), xa_xm
                   !print*,dxf(i,:)
                   !print*,"**************************************"
                END IF
                !print*,ens_p(i,:)
                !print*,"Par_STD_Index,xa_xm,dxa_dev(i)",Par_STD_Index,xa_xm,dxa_dev(i),dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(i)) / dxa_dev(i) + 1.0)))

                !linspace of parameter ensembls to keep uniform spread
                !CALL MDIAN(ens_p(i,:),Ensemble_Number,xa_median)

                !dxa_max = xa_median + Par_Uniform_STD(Par_STD_Index)/(sqrt(1/12.0)*2.0)
                !dxa_min = 2.0 * xa_median - dxa_max
                !print*,"i,xa_median,Par_Uniform_STD(Par_STD_Index),dxa_min,dxa_max",i,xa_median,Par_Uniform_STD(Par_STD_Index),dxa_min,dxa_max

                !seed(:) = 123456
                !CALL RANDOM_SEED(PUT = seed)
                !CALL RANDOM_NUMBER(U)

                !ens_p(i,:) = (U*(dxa_max-dxa_min))+dxa_min
                !print*,"ens_p(i,:)",ens_p(i,:)

            ! State Inflation
            ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Forecast_Model_Option == 0) .AND. &
                (Bias_Observation_Model_Option == 0) .AND. (Alpha_Inflation .NE. 0.0) .AND. (dxa_dev(i) < dxf_dev(i))) THEN

                IF (dxa_dev(i) .GT. 0.0) THEN
                    ens_p(i,:) = xa_xm + dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * &
                        (dxf_dev(i) - dxa_dev(i)) / dxa_dev(i) + 1.0)))  !RTPS (multiplicative)
                ELSE
                    seed(:) = 123456
                    CALL RANDOM_SEED(PUT = seed)
                    CALL RANDOM_NUMBER(U)
                    U = 0.8 * (1.0 - U) + 1.2 * U
                    ens_p(i,:) = ens_p(i,:) * U
                END IF
                !future inflation
                !CALL com_mean(Ensemble_Number,ens_p(i,:),xa_xm)
                !dxa(i,:) = ens_p(i,:) - xa_xm
                !CALL com_stdev(Ensemble_Number,dxa(i,:),dxa_dev(i))
                !!print*,"dxa_dev(i),Model_Inflation_Uniform_STD(i)",dxa_dev(i),Model_Inflation_Uniform_STD(i)
                !IF ((dxa_dev(i) < Model_Inflation_Uniform_STD(i)) .AND. (dxa_dev(i) .GT. 0.0)) THEN
                !    ens_p(i,:) = xa_xm + dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * (Model_Inflation_Uniform_STD(i) - dxa_dev(i)) / dxa_dev(i) + 1.0)))  !RTPS (multiplicative)
                !END IF
                !CALL com_mean(Ensemble_Number,ens_p(i,:),xa_xm)
                !dxa(i,:) = ens_p(i,:) - xa_xm
                !CALL com_stdev(Ensemble_Number,dxa(i,:),dxa_dev(i))
                !print*,"dxa_dev(i)",dxa_dev(i)

            ! Model Bias Inflation
            ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Forecast_Model_Option > 0) .AND. &
                (Alpha_Inflation .NE. 0.0) .AND. (i>State_DIM_Single_Layer) .AND. (i<=(State_DIM_Single_Layer+Bias_Model_Dim_LETKF)) .AND. (Bias_Model_Dim_LETKF >= 1)) THEN
                IF (dxa_dev(i) < Bias_Model_Uniform_STD(i-State_DIM_Single_Layer)) THEN
                    IF (dxa_dev(i) .GT. 0.0) THEN
                        ens_p(i,:) = xa_xm + dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * &
                            (Bias_Model_Uniform_STD(i-State_DIM_Single_Layer) - dxa_dev(i)) / dxa_dev(i) + 1.0)))  !RTPS (multiplicative)
                        !ens_p(i,:) = ens_p(i,:) * (Alpha_Inflation * (Bias_Model_Uniform_STD(i-State_DIM_Single_Layer) - dxa_dev(i)) / dxa_dev(i) + 1.0)  !RTPS (multiplicative)
                    ELSE
                        seed(:) = 123456
                        CALL RANDOM_SEED(PUT = seed)
                        CALL RANDOM_NUMBER(U)
                        U = 0.8 * (1.0 - U) + 1.2 * U
                        ens_p(i,:) = ens_p(i,:) * U
                    END IF
                END IF
                !linspace of parameter ensembls to keep uniform spread
                CALL MDIAN(ens_p(i,:),Ensemble_Number,xa_median)

                dxa_max = xa_median + Bias_Model_Uniform_STD(i-State_DIM_Single_Layer)/(sqrt(1/12.0)*2.0)
                dxa_min = 2.0 * xa_median - dxa_max
                !print*,"i,xa_median,Bias_Model_Uniform_STD(i-State_DIM_Single_Layer),dxa_min,dxa_max",i,xa_median,Bias_Model_Uniform_STD(i-State_DIM_Single_Layer),dxa_min,dxa_max

                seed(:) = 123456
                CALL RANDOM_SEED(PUT = seed)
                CALL RANDOM_NUMBER(U)

                ens_p(i,:) = (U*(dxa_max-dxa_min))+dxa_min
                !print*,"ens_p(i,:)",ens_p(i,:)

            ! Observation Bias Inflation
            ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Observation_Model_Option > 0) .AND. &
                (Alpha_Inflation .NE. 0.0) .AND. (i>(State_DIM_Single_Layer+Bias_Model_Dim_LETKF)) .AND. (Bias_Obs_Dim_LETKF >= 1)) THEN
                IF (dxa_dev(i) < Bias_Obs_Uniform_STD(i-(State_DIM_Single_Layer+Bias_Model_Dim_LETKF))) THEN
                    IF (dxa_dev(i) .GT. 0.0) THEN
                        ens_p(i,:) = xa_xm + dxa(i,:) * max(1.0, min(2.0, (Alpha_Inflation * &
                            (Bias_Obs_Uniform_STD(i-(State_DIM_Single_Layer+Bias_Model_Dim_LETKF)) - dxa_dev(i)) / dxa_dev(i) + 1.0)))  !RTPS (multiplicative)
                    ELSE
                        seed(:) = 123456
                        CALL RANDOM_SEED(PUT = seed)
                        CALL RANDOM_NUMBER(U)
                        U = 0.8 * (1.0 - U) + 1.2 * U
                        ens_p(i,:) = ens_p(i,:) * U
                    END IF
                END IF
                !linspace of parameter ensembls to keep uniform spread
                CALL MDIAN(ens_p(i,:),Ensemble_Number,xa_median)

                dxa_max = xa_median + Bias_Obs_Uniform_STD(i-(State_DIM_Single_Layer+Bias_Model_Dim_LETKF))/(sqrt(1/12.0)*2.0)
                dxa_min = 2.0 * xa_median - dxa_max
                !print*,"i,xa_median,Bias_Obs_Uniform_STD(i-(State_DIM_Single_Layer+Bias_Model_Dim_LETKF),dxa_min,dxa_max",i,xa_median,Bias_Obs_Uniform_STD(i-(State_DIM_Single_Layer+Bias_Model_Dim_LETKF),dxa_min,dxa_max

                seed(:) = 123456
                CALL RANDOM_SEED(PUT = seed)
                CALL RANDOM_NUMBER(U)

                ens_p(i,:) = (U*(dxa_max-dxa_min))+dxa_min
                !print*,"ens_p(i,:)",ens_p(i,:)
            END IF

        END DO

        DEALLOCATE(XF_NC_Copy)

    END IF

    ! *** Compute mean state
    WRITE (*, '(8x, a)') '--- compute ensemble mean'

    state_p = 0.0
    DO member = 1, dim_ens
        DO i = 1, dim_p
            state_p(i) = state_p(i) + ens_p(i, member)
        END DO
    END DO
    state_p(:) = invdim_ens * state_p(:)

    ! *** Compute sampled variances ***
    variance(:) = 0.0
    DO member = 1, dim_ens
        DO j = 1, dim_p
            variance(j) = variance(j) &
                + (ens_p(j, member) - state_p(j)) &
                * (ens_p(j, member) - state_p(j))
        END DO
    END DO
    variance(:) = invdim_ensm1 * variance(:)

    ! ************************************************************
    ! *** Compute RMS errors according to sampled covar matrix ***
    ! ************************************************************

    ! total estimated RMS error
    DO i = 1, dim_p
        rmserror_est = rmserror_est + variance(i)
    ENDDO
    rmserror_est = SQRT(rmserror_est / dim_p)


    ! *****************
    ! *** Screen IO ***
    ! *****************

    ! Output RMS errors given by sampled covar matrix
    WRITE (*, '(12x, a, es12.4)') &
        'RMS error according to sampled variance: ', rmserror_est

  
    ! *******************
    ! *** File output ***
    ! *******************

    IF (.not. firsttime) THEN

        print *,"*** Writing file ", FILE_NAME, "! "
        ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
        ! the file.
        call check( nf90_open(FILE_NAME, NF90_WRITE, ncid) )

        ! Get the varid of the data variable, based on its name.
        call check( nf90_inq_varid(ncid, XM_NAME, varid) )
        ! Write the pretend data to the file. Although netCDF supports
        ! reading and writing subsets of data, in this case we write all the
        ! data in one operation.
        call check( nf90_put_var(ncid, varid, state_p, start = (/ 1 /), count = (/ STATE_DIM /)) )

        ! Get the varid of the data variable, based on its name.
        call check( nf90_inq_varid(ncid, XA_NAME, varid) )
        ! Write the pretend data to the file. Although netCDF supports
        ! reading and writing subsets of data, in this case we write all the
        ! data in one operation.
        call check( nf90_put_var(ncid, varid, ens_p, start = (/ 1, 1 /), count = (/ STATE_DIM, ENSEMBLE_NUMBER /)) )

        !print*,"state_p",state_p
        !print*,"ens_p",ens_p(1,:)

        ! Close the file. This frees up any internal netCDF resources
        ! associated with the file.
        call check( nf90_close(ncid) )

        print *,"*** SUCCESS write file ", FILE_NAME, "! "
        print *,""
    END IF

    ! ********************
    ! *** finishing up ***
    ! ********************

    DEALLOCATE(variance)

    firsttime = .FALSE.

contains
    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

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

        !*******************************************************
    !* Given an array X of N numbers, returns their median *
    !* value XMED. The array X is not modified, and is     *
    !* accessed sequentially in each consecutive pass.     *
    !*******************************************************
    SUBROUTINE MDIAN(X,N,XMED)
        integer, INTENT(IN) :: N
        real(8), INTENT(IN) :: X(N)
        real(8), INTENT(INOUT) :: XMED
        real(8), parameter :: BIG = 1.e30, AFAC=1.5, AMP=1.5
        integer N2

        call hpsort(N,X)
        N2=N/2
        if (2*N2.eq.N) then
            XMED = 0.5*(X(N2)+X(N2+1))
        else
            XMED = X(N2+1)
        endif
        return
    END SUBROUTINE MDIAN

    !*****************************************************
    !*  Sorts an array RA of length N in ascending order *
    !*                by the Heapsort method             *
    !* ------------------------------------------------- *
    !* INPUTS:                                           *
    !*      N     size of table RA                       *
    !*      RA    table to be sorted                     *
    !* OUTPUT:                                           *
    !*      RA    table sorted in ascending order        *
    !*                                                   *
    !* NOTE: The Heapsort method is a N Log2 N routine,  *
    !*       and can be used for very large arrays.      *
    !*****************************************************
    SUBROUTINE HPSORT(N,RA)
        integer, INTENT(IN) :: N
        real(8) RA(N),RRA
        integer L, IR, I, J
        L=N/2+1
        IR=N
    !The index L will be decremented from its initial value during the
    !"hiring" (heap creation) phase. Once it reaches 1, the index IR
    !will be decremented from its initial value down to 1 during the
    !"retirement-and-promotion" (heap selection) phase.
10  continue
    if(L > 1)then
        L=L-1
        RRA=RA(L)
    else
        RRA=RA(IR)
        RA(IR)=RA(1)
        IR=IR-1
        if(IR.eq.1)then
            RA(1)=RRA
            return
        end if
    end if
    I=L
    J=L+L
20  if(J.le.IR)then
        if(J < IR)then
            if(RA(J) < RA(J+1))  J=J+1
        end if
        if(RRA < RA(J))then
            RA(I)=RA(J)
            I=J; J=J+J
        else
            J=IR+1
        end if
        goto 20
    end if
    RA(I)=RRA
    goto 10
END

SUBROUTINE linspace(d1,d2,n,grid)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: d1, d2
    REAL(8), DIMENSION(n), INTENT(OUT) :: grid

    INTEGER :: indxi


    grid(1) = d1
    DO indxi= 0,n-2
        grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
    END DO
    grid(n) = d2

    !MATLAB
    !grid = [d1+(0:n-2)*(d2-d1)/(floor(n)-1) d2];


END SUBROUTINE
END SUBROUTINE prepoststep_ens_offline
