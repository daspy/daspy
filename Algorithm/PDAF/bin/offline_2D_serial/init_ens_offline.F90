!$Id: init_ens_offline.F90 1421 2013-09-25 15:10:05Z lnerger $
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
    ens_p, flag)

    ! !DESCRIPTION:
    ! User-supplied routine for PDAF.
    ! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
    !
    ! The routine is called when the filter is
    ! initialized in PDAF\_filter\_init.  It has
    ! to initialize an ensemble of dim\_ens states.
    ! For the offline mode, the ensemble will be
    ! typically read-in from files.
    !
    ! The routine is called by all filter processes and
    ! initializes the ensemble for the PE-local domain.
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
        ONLY: nx, ny, locweight, use_obs_mask, &
        STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, Normal_Score_Trans, &
        XF_NC, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, &
        FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, Normal_Score_Trans_NAME, &
        XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME, &
        State_DIM_Single_Layer, Parameter_DIM, Par_Uniform_STD, Alpha_Inflation, &
        State_DIM_Single_Layer_NAME, Parameter_DIM_NAME, Par_Uniform_STD_NAME, Alpha_Inflation_NAME, &
        Parameter_Optimization_Flag, Parameter_Optimization_Flag_NAME, &
        Bias_Forecast_Model_Option, Bias_Observation_Model_Option, &
        State_DIM_Single_Column, State_DIM_Single_Column_NAME, &
        X_Left, X_Right, Y_Lower, Y_Upper, &
		lmbda_DIM, minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint, &
		minimize_lbfgsb_factr, minimize_lbfgsb_pgtol, &
		lmbda_state, lmbda_parameter, lmbda_bias, hlmbda, minimize_lbfgsb_epsilon_in, &
		Par_Sens_Dim, State_DIM_Single_Column, Bias_Model_Dim, Bias_Obs_Dim

    use netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
    ! It is not necessary to initialize the array 'state_p' for SEIK.
    ! It is available here only for convenience and can be used freely.
    REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
    REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
    INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

    ! !CALLING SEQUENCE:
    ! Called by: PDAF_filter_init    (as U_ens_init)
    !EOP

    ! *** local variables ***
    INTEGER :: i, j, member  ! Counters
    INTEGER, SAVE :: allocflag = 0      ! Flag for memory counting
    REAL, ALLOCATABLE :: field(:,:)     ! global model field
    CHARACTER(len=2) :: ensstr          ! String for ensemble member

    REAL, ALLOCATABLE :: XF_NNscore_Whole(:), XF_NNscore(:), XF_Trans(:)
	REAL :: Observation_NNscore(1), Observation_Trans(1)
	
    ! **********************
    ! *** INITIALIZATION ***
    ! **********************

    ! *** Generate full ensemble on filter-PE 0 ***
    WRITE (*, '(/9x, a)') 'Initialize state ensemble'
    WRITE (*, '(9x, a)') '--- read ensemble from files'
    WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens

    ! ********************************
    ! *** Read ensemble from files ***
    ! ********************************


    

    ens_p = XF_NC

    ! ****************
    ! *** clean up ***
    ! ****************

    DEALLOCATE(XF_NC)
    

END SUBROUTINE init_ens_offline
