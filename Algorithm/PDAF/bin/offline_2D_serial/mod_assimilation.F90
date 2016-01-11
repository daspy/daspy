!$Id: mod_assimilation.F90 1399 2013-05-06 09:21:15Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assimilation

! !DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
! Implementation for the 2D offline example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE
  SAVE
!EOP


! *** Variables specific for offline tutorial example ***

  INTEGER :: nx, ny              ! Size of 2D grid
  INTEGER :: dim_state_p         ! Model state dimension

  REAL, ALLOCATABLE    :: obs(:)          ! Vector holding all observations
  INTEGER, ALLOCATABLE :: obs_index(:)    ! Vector holding state-vector indices of observations
  INTEGER, ALLOCATABLE :: obs_index_l(:)  ! Vector holding local state-vector indices of observations
  INTEGER, ALLOCATABLE :: coords_obs(:,:) ! Array for observation coordinates
  INTEGER :: coords_l(2)                  ! Coordinates of local analysis domain

! NC Varialbes added by Xujun ***************************************************************
  INTEGER :: STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, State_DIM_Single_Layer, Correlation_Range
  INTEGER :: Parameter_DIM, Par_Sens_Dim, State_DIM_Single_Column, Bias_Model_Dim, Bias_Obs_Dim
  INTEGER :: ncid, varid, Normal_Score_Trans, Parameter_Optimization_Flag, Bias_Forecast_Model_Option, Bias_Observation_Model_Option
  REAL    :: Alpha_Inflation, X_Left, X_Right, Y_Lower, Y_Upper, GridSize_Sys, GridSize_Obs
  ! Read the data into these arrays.
  LOGICAL			:: use_obs_mask
  REAL, ALLOCATABLE :: B(:,:), B_step(:,:)
  REAL, ALLOCATABLE :: XF_NC(:,:), XF_NC_Copy(:,:), HXF_NC(:,:), Par_Uniform_STD(:)
  REAL, ALLOCATABLE :: OBS_NC(:), XF_COORD_NC(:,:), OBS_COORD_NC(:,:), R_NC(:,:)
  REAL, ALLOCATABLE :: R_Local(:,:), R_Local_l(:,:), coords_obs_l(:,:)
  REAL, ALLOCATABLE :: Bias_Model_Uniform_STD(:),Bias_Obs_Uniform_STD(:),Model_Inflation_Uniform_STD(:)
  INTEGER, ALLOCATABLE :: H_NC(:,:)
  INTEGER, ALLOCATABLE :: H_Local(:,:)
  integer			   :: lmbda_DIM, minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint
  REAL	    		   :: minimize_lbfgsb_factr, minimize_lbfgsb_pgtol
  REAL, ALLOCATABLE    :: minimize_lbfgsb_epsilon_in(:), lmbda_state(:), lmbda_parameter(:,:), lmbda_bias(:,:), hlmbda(:)

! This is the name of the data file we will read.
  character (len = *), parameter :: FILE_NAME = "NC_to_PDAF.nc"
  character (len = *), parameter :: STATE_DIM_NAME = "STATE_DIM"
  character (len = *), parameter :: OBS_DIM_NAME = "OBS_DIM"
  character (len = *), parameter :: ENSEMBLE_NUMBER_NAME = "ENSEMBLE_NUMBER"
  character (len = *), parameter :: Normal_Score_Trans_NAME = "Normal_Score_Trans"
  character (len = *), parameter :: XF_NAME = "XF_NC"
  character (len = *), parameter :: HXF_NAME = "HXF_NC"
  character (len = *), parameter :: H_NAME = "H_NC"
  character (len = *), parameter :: OBS_NAME = "OBS_NC"
  character (len = *), parameter :: XF_COORD_NAME = "XF_COORD_NC"
  character (len = *), parameter :: OBS_COORD_NAME = "OBS_COORD_NC"
  character (len = *), parameter :: GridSize_Sys_NAME = "GridSize_Sys"
  character (len = *), parameter :: GridSize_Obs_NAME = "GridSize_Obs"
  character (len = *), parameter :: R_NAME = "R_NC"
  character (len = *), parameter :: XA_NAME = "XA_NC"
  character (len = *), parameter :: XM_NAME = "XM_NC"
  character (len = *), parameter :: State_DIM_Single_Layer_NAME = "State_DIM_Single_Layer"
  character (len = *), parameter :: State_DIM_Single_Column_NAME = "State_DIM_Single_Column"
  character (len = *), parameter :: Parameter_DIM_NAME = "Parameter_DIM"
  character (len = *), parameter :: Par_Sens_Dim_NAME = "Par_Sens_Dim"
  character (len = *), parameter :: Alpha_Inflation_NAME = "Alpha_Inflation"
  character (len = *), parameter :: Par_Uniform_STD_NAME = "Par_Uniform_STD"
  character (len = *), parameter :: Bias_Model_Dim_NAME = "Bias_Model_Dim"
  character (len = *), parameter :: Bias_Obs_Dim_NAME = "Bias_Obs_Dim"
  character (len = *), parameter :: Parameter_Optimization_Flag_NAME = "Parameter_Optimization_Flag"
  character (len = *), parameter :: Bias_Model_Uniform_STD_NAME = "Bias_Model_Uniform_STD"
  character (len = *), parameter :: Bias_Obs_Uniform_STD_NAME = "Bias_Obs_Uniform_STD"
  character (len = *), parameter :: Model_Inflation_Uniform_STD_NAME = "Model_Inflation_Uniform_STD"
  character (len = *), parameter :: Bias_Forecast_Model_Option_NAME = "Bias_Forecast_Model_Option"
  character (len = *), parameter :: Bias_Observation_Model_Option_NAME = "Bias_Observation_Model_Option"
  character (len = *), parameter :: Correlation_Range_NAME = "Correlation_Range"
  character (len = *), parameter :: minimize_lbfgsb_n_NAME = "minimize_lbfgsb_n"
  character (len = *), parameter :: minimize_lbfgsb_m_NAME = "minimize_lbfgsb_m"
  character (len = *), parameter :: minimize_lbfgsb_factr_NAME = "minimize_lbfgsb_factr"
  character (len = *), parameter :: minimize_lbfgsb_pgtol_NAME = "minimize_lbfgsb_pgtol"
  character (len = *), parameter :: minimize_lbfgsb_iprint_NAME = "minimize_lbfgsb_iprint"
  character (len = *), parameter :: minimize_lbfgsb_epsilon_in_NAME = "minimize_lbfgsb_epsilon_in"
  character (len = *), parameter :: lmbda_DIM_NAME = "lmbda_DIM"
  
! NC Varialbes added by Xujun ***************************************************************


! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF_offline                 ***

! !PUBLIC MEMBER FUNCTIONS:
! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL    :: model_err_amp ! Amplitude for model error

! ! Settings for observations - available as command line options
  INTEGER :: delt_obs      ! time step interval between assimilation steps
  REAL    :: rms_obs       ! RMS error size for observation generation
  INTEGER :: dim_obs       ! Number of observations

! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
                          ! (0) no outputs, (1) progess info, (2) add timings
                          ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                          ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
                          ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
  INTEGER :: subtype      ! Subtype of filter algorithm
                          !   SEEK: 
                          !     (0) evolve normalized modes
                          !     (1) evolve scaled modes with unit U
                          !     (2) fixed basis (V); variable U matrix
                          !     (3) fixed covar matrix (V,U kept static)
                          !   SEIK:
                          !     (0) ensemble forecast; new formulation
                          !     (1) ensemble forecast; old formulation
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) SEIK with ensemble transformation
                          !   EnKF:
                          !     (0) analysis for large observation dimension
                          !     (1) analysis for small observation dimension
                          !   LSEIK:
                          !     (0) ensemble forecast;
                          !     (2) fixed error space basis
                          !     (3) fixed state covariance matrix
                          !     (4) LSEIK with ensemble transformation
                          !   ETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) ETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to SEIK subtypes 2/3
                          !   LETKF:
                          !     (0) ETKF using T-matrix like SEIK
                          !     (1) LETKF following Hunt et al. (2007)
                          !       There are no fixed basis/covariance cases, as
                          !       these are equivalent to LSEIK subtypes 2/3
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

! ! Filter settings - available as command line options
!    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
!    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
!    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
!    ! SEIK/ETKF/LSEIK/ETKFS
  INTEGER :: type_trans    ! Type of ensemble transformation
                           ! SEIK/LSEIK:
                           ! (0) use deterministic omega
                           ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
                           ! ETKF/LETKF with subtype=4:
                           ! (0) use deterministic symmetric transformation
                           ! (2) use product of (0) with random orthonormal matrix with
                           !     eigenvector (1,...,1)^T
!    ! LSEIK/LETKF
  REAL    :: local_range   ! Range for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  REAL    :: srange        ! Support range for 5th order polynomial
                           !   or radius for 1/e for exponential weighting
!    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                    !   (0) symmetric square root, (1) Cholesky decomposition

!    ! File output - available as a command line option
  CHARACTER(len=110) :: filename  ! file name for assimilation output

!    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     ! For SEIK: Definition of ensemble covar matrix
                           ! (0): Factor (r+1)^-1 (or N^-1)
                           ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
                           ! This setting is only for the model part; The definition
                           ! of P has also to be specified in PDAF_filter_init.
                           ! Only for upward-compatibility of PDAF!
  REAL    :: time          ! model time

  !$OMP THREADPRIVATE(obs_index_l, coords_l)

END MODULE mod_assimilation
