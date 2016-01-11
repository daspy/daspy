!$Id: init_pdaf.F90 1346 2013-04-10 08:50:00Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf - Interface routine to call initialization of PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables
! for PDAF as well as the call to the initialization
! routine PDAF_init.
!
! This variant is for the Lorenz96 model.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE parser, &
       ONLY: parse
  USE mod_model, &
       ONLY: step_null, dim_state, dt
  USE mod_modeltime, &
       ONLY: total_steps
  USE mod_parallel, &
       ONLY: mype_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       rms_obs, model_error, model_err_amp, incremental, covartype, &
       type_forget, forget, dim_bias, epsilon, rank_analysis_enkf, &
       locweight, local_range, local_range2, srange, int_rediag, &
       file_ini, file_obs, type_ensinit, seedset, type_trans, &
       type_sqrt, stepnull_means, dim_lag
  USE output_netcdf_asml, &
       ONLY: init_netcdf_asml, file_asml, delt_write_asml, write_states, &
       write_stats, write_ens

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_init
! Calls: parse
!EOP

! Local variables
  CHARACTER(len=32) :: handle  ! handle for command line parser
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag

  ! External subroutines
  EXTERNAL :: init_seik  ! SEIK ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  IF (mype_world == 0) THEN
     WRITE (*,'(/1x,a)') 'INITIALIZE PDAF'
  END IF


! **********************************************************
! ***               CONTROL OF PDAF                      ***
! **********************************************************

! *** IO options ***
  screen      = 2   ! Write screen output (1) for output, (2) add timings
  delt_write_asml = 1    ! Output interval for state information in assimilation intervals
  write_states = .TRUE.  ! Whether to write estimates states into the output file
  write_stats  = .FALSE. ! Whether to write time dependent ensemble statistics (skewness, kurtosis)
  write_ens  = .FALSE.   ! Whether to write full time dependent ensemble stats
  stepnull_means = 3001  ! Step at which the second computation of time mean error is started
                         ! (first computation of mean sis always starting at initial step)

! *** specifications for observations ***
  ! avg. observation error (used for assimilation)
  rms_obs = 1.0     ! This error is the standard deviation 
                    ! for the Gaussian distribution 
  delt_obs = 1      ! Time step interval between analysis/assimilation steps

! *** Filter specific variables
  filtertype = 1    ! Type of filter
                    !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
                    !   ESTKF (6), LESTKF (7)
  dim_ens = 30      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                    ! Number of EOFs to be used for SEEK
  dim_lag = 0       ! Size of lag in smoother
  subtype = 0       ! subtype of filter: 
                    !   SEEK: 
                    !     (0) evolve normalized modes
                    !     (1) evolve scaled modes with unit U
                    !     (2) fixed basis (V); variable U matrix
                    !     (3) fixed covar matrix (V,U kept static)
                    !   SEIK:
                    !     (0) mean forecast; new formulation
                    !     (1) mean forecast; old formulation
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) SEIK with ensemble transformation
                    !   EnKF:
                    !     (0) analysis for large observation dimension
                    !     (1) analysis for small observation dimension
                    !   LSEIK:
                    !     (0) mean forecast;
                    !     (2) fixed error space basis
                    !     (3) fixed state covariance matrix
                    !     (4) LSEIK with ensemble transformation
                    !   ETKF:
                    !     (0) ETKF using T-matrix like SEIK
                    !     (1) ETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to SEIK subtypes 2/3
                    !   LETKF:
                    !     (0) LETKF using T-matrix like SEIK
                    !     (1) LETKF following Hunt et al. (2007)
                    !       There are no fixed basis/covariance cases, as
                    !       these are equivalent to LSEIK subtypes 2/3
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  int_rediag = 1    ! Interval of analysis steps to perform 
                    !    re-diagonalization in SEEK
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  type_ensinit = 'eof' ! 'eof' for 2nd-order exact sampling from EOFs
                    !    'rnd' for random sampling from true state trajectory
  seedset = 1       ! Index of set of seeds to be used for init (only for 'rnd')
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    ! (0) for (r+1)^-1 (old SEIK); (1): for r^-1 (real ensemble
                    ! covariance matrix) This parameter has also to be set internally
                    ! in PDAF_init
  rank_analysis_enkf = 0   ! rank to be considered for inversion of HPH
                    ! in analysis of EnKF; (0) for analysis w/o eigendecomposition
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK
                    ! (0): fixed; (1) global adaptive; (2) local adaptive for LSEIK
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  epsilon = 1.0E-4  ! epsilon for approx. TLM evolution in SEEK
  local_range = 5   ! Range in grid points for observation domain in LSEIK
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE for observed ensemble
                    !   (2) use sqrt of 5th-order polynomial for observed ensemble
                    !   (3) exponentially decreasing with SRANGE for obs. error matrix
                    !   (4) use 5th-order polynomial for obs. error matrix
                    !   (5) use 5th-order polynomial for observed ensemble
                    !   (6) regulated localization of R with mean error variance
                    !   (7) regulated localization of R with single-point error variance
  srange = local_range  ! Support range for 5th-order polynomial
                    ! range for 1/e for exponential weighting

! *** File names
  file_asml = 'assimilation.nc' ! Output file
  file_ini = 'covar.nc'         ! Initialization file
  file_obs = 'obs.nc'           ! File holding observations


! ********************************
! ***      INITIALIZATION      ***
! ********************************

! *** Parse command line options ***

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL parse(handle, model_err_amp)

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL parse(handle, delt_obs)
  handle = 'rms_obs'                 ! Assumed uniform RMS error of the observations
  CALL parse(handle, rms_obs)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)
  handle = 'incremental'             ! Set whether to use incremental updating
  CALL parse(handle, incremental)

  ! Settings for smoother
  handle = 'dim_lag'                 ! Size of lag in smoother
  CALL parse(handle, dim_lag)

  ! General settings - external of PDAF
  handle = 'type_ensinit'            ! Define type of ensemble initialization
  CALL parse(handle, type_ensinit)   ! possible are 'eof' or 'rnd'
  handle = 'seedset'                 ! Choose set of seeds for init (only for 'rnd')
  CALL parse(handle, seedset)        ! valid values are 1 to 5

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
  CALL parse(handle, type_trans)
  handle = 'epsilon'                 ! Set EPSILON for SEEK
  CALL parse(handle, epsilon)
  handle = 'int_rediag'              ! Time step interval for rediagonalization in SEEK
  CALL parse(handle, int_rediag)
  handle = 'rank_analysis_enkf'      ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_analysis_enkf)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_sqrt'               ! Set type of transform square-root
  CALL parse(handle,type_sqrt)

  ! Settings for localization in LSEIK/LETKF
  handle = 'local_range'             ! Set range in grid points for observation domain
  CALL parse(handle, local_range)
  local_range2 = local_range
  handle = 'local_range2'            ! Set right-side range in grid points for observation domain
  CALL parse(handle, local_range2)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  srange = local_range               ! By default use local_range as support range
  handle = 'srange'                  ! Set support range in grid points
             ! for 5th-order polynomial or range for 1/e in exponential weighting
  CALL parse(handle, srange)

  ! Setting for file output
  handle = 'delt_write_asml'         ! Set write intetval for output in assimilation cycles
  CALL parse(handle, delt_write_asml)
  handle = 'write_states'            ! Define, whether state information is written to file
  CALL parse(handle, write_states)
  handle = 'write_stats'             ! Define, whether to write higher-order ensemble statistics
  CALL parse(handle, write_stats)
  handle = 'write_ens'               ! Define, whether to write full ensemble
  CALL parse(handle, write_ens)
  handle = 'stepnull_means '         ! Step at which computation of time mean error is started
  CALL parse(handle, stepnull_means)
  handle = 'file_asml'               ! Set name of output file
  CALL parse(handle, file_asml)
  handle = 'file_ini'                ! Set name of initialization file
  CALL parse(handle, file_ini)
  handle = 'file_obs'                ! Set name of file holding observations
  CALL parse(handle, file_obs)


! *** Initial Screen output ***
  screen2: IF (mype_world == 0) THEN

     IF (filtertype == 1) THEN
        WRITE (*, '(21x, a)') 'Filter: SEIK'
        IF (subtype == 2) THEN
           WRITE (*, '(6x, a)') '-- fixed error-space basis'
        ELSE IF (subtype == 3) THEN
           WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
        ELSE IF (subtype == 4) THEN
           WRITE (*, '(6x, a)') '-- use ensemble transformation'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 2) THEN
        WRITE (*, '(21x, a)') 'Filter: EnKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
        IF (rank_analysis_enkf > 0) THEN
           WRITE (*, '(6x, a, i5)') &
                'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
        END IF
     ELSE IF (filtertype == 3) THEN
        WRITE (*, '(21x, a)') 'Filter: LSEIK'
        IF (subtype == 2) THEN
           WRITE (*, '(6x, a)') '-- fixed error-space basis'
        ELSE IF (subtype == 3) THEN
           WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
        ELSE IF (subtype == 4) THEN
           WRITE (*, '(6x, a)') '-- use ensemble transformation'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 4) THEN
        WRITE (*, '(21x, a)') 'Filter: ETKF'
        IF (subtype == 0) THEN
           WRITE (*, '(17x, a)') '--> Variant using T-matrix'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(17x, a)') '--> Variant following Hunt et al. (2007)'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 5) THEN
        WRITE (*, '(21x, a)') 'Filter: LETKF'
        IF (subtype == 0) THEN
           WRITE (*, '(17x, a)') '--> Variant using T-matrix'
        ELSE IF (subtype == 1) THEN
           WRITE (*, '(17x, a)') '--> Variant following Hunt et al. (2007)'
        END IF
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 6) THEN
        WRITE (*, '(21x, a)') 'Filter: ESTKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     ELSE IF (filtertype == 7) THEN
        WRITE (*, '(21x, a)') 'Filter: LESTKF'
        WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
        IF (dim_lag > 0) WRITE (*, '(15x, a, i5)') 'smoother lag:', dim_lag
        WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
        WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
        IF (model_error) THEN
           WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
        END IF
     END IF
     WRITE (*, '(8x, a)') 'File names:'
     WRITE (*, '(11x, a, a)') 'Initialization: ', TRIM(file_ini)
     WRITE (*, '(11x, a, a)') 'Observations:   ', TRIM(file_obs)
     WRITE (*, '(11x, a, a)') 'Output:         ', TRIM(file_asml)

  END IF screen2


! *****************************************************
! *** Call filter initialization routine on all PEs ***
! *****************************************************

  whichinit: IF (filtertype == 1) THEN
     ! *** SEIK with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 2) THEN
     ! *** EnKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = rank_analysis_enkf ! Maximum rank for matrix inversion
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 3, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 3) THEN
     ! *** LSEIK with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 4) THEN
     ! *** ETKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 6, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 5) THEN
     ! *** LETKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 6, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 6) THEN
     ! *** ESTKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor

     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  ELSEIF (filtertype == 7) THEN
     ! *** LESTKF with init by 2nd order exact sampling ***
     filter_param_i(1) = dim_state   ! State dimension
     filter_param_i(2) = dim_ens     ! Size of ensemble
     filter_param_i(3) = dim_lag     ! Size of lag in smoother
     filter_param_i(4) = incremental ! Whether to perform incremental analysis
     filter_param_i(5) = type_forget ! Type of forgetting factor
     filter_param_i(6) = type_trans  ! Type of ensemble transformation
     filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
     filter_param_r(1) = forget      ! Forgetting factor
     
     CALL PDAF_init(filtertype, subtype, step_null, &
          filter_param_i, 7, &
          filter_param_r, 2, &
          COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, filterpe, init_seik, &
          screen, status_pdaf)
  END IF whichinit

! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

  ! Initialize netcdf output
  CALL init_netcdf_asml(step_null, dt, dim_state, filtertype, subtype, &
       dim_ens, forget, type_ensinit, local_range, local_range2, &
       locweight, srange, rms_obs, delt_obs, total_steps, &
       seedset, stepnull_means, dim_lag)

END SUBROUTINE init_pdaf
