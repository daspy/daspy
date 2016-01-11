!$Id: compute_truermse.F90 1160 2011-09-14 09:32:08Z lnerger $
!BOP
!
! !ROUTINE: compute_truermse --- Compute true rms errors for assimilation
!
! !INTERFACE:
SUBROUTINE compute_truermse(calltype, step, time, dim, state_est, &
       trmse, dim_ens, ens, hist_true, hist_mean, &
       skewness, kurtosis)

! !DESCRIPTION:
! Helper routine for the pre/poststep routines. The routine
! reads the true state from the trajectory file and computes
! the true rms error for the currrent state estimate. In 
! addition, rank histograms about the ensemble mean and
! the true state as well as higher-order statistical moments
! (skewness, kurtosis) are computed.
!
! !REVISION HISTORY:
! 2010-01 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_memcount, &
       ONLY: memcount
  USE mod_assimilation, &
       ONLY: fileid_state, state_true, stepnull_means
  USE output_netcdf, &
       ONLY: file_state

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  CHARACTER(len=3), INTENT(in) :: calltype ! Whether routine is called at 
        !(ini) initial time, (for) after forecast, or (ana) after analysis
  INTEGER, INTENT(in) :: step           ! Current time step
        ! (When the routine is called before the analysis -step is provided.)
  REAL, INTENT(in)    :: time           ! Currrent model time
  INTEGER, INTENT(in) :: dim            ! PE-local state dimension
  REAL, INTENT(in)    :: state_est(dim) ! Forecast/analysis state
  REAL, INTENT(out)   :: trmse          ! True rms error
  INTEGER, INTENT(in) :: dim_ens        ! Ensemble size
  REAL, INTENT(in)    :: ens(dim, dim_ens)     ! State ensemble
  ! Arrays for rank histograms
  INTEGER, INTENT(inout) :: hist_true(dim_ens+1, 2) ! Histogram about true state
  INTEGER, INTENT(inout) :: hist_mean(dim_ens+1, 2) ! Histogram about ensemble mean
  REAL, INTENT(out) :: skewness         ! Skewnes of ensemble
  REAL, INTENT(out) :: kurtosis         ! Kurtosis of ensemble

! !CALLING SEQUENCE:
! Called by: prepoststep routine
!EOP

! *** local variables ***
  INTEGER :: i, j, s                    ! Counters
  INTEGER :: id_dim, id_step            ! File dimension IDs
  INTEGER :: id_state                   ! File variable ID
  INTEGER :: stat(50)                   ! Array for status flag
  INTEGER :: nsteps_file                ! Number of time steps in trajectory file
  INTEGER :: pos(2)                     ! Position index for writing
  INTEGER :: cnt(2)                     ! Count index for writing
  INTEGER :: state_step1and2(2)         ! First and second time step index in state file
  INTEGER, SAVE :: statefile_laststep   ! Last time step stored in state file
  INTEGER, SAVE :: delt_state_file      ! Interval between sucessively stored states
  INTEGER :: read_pos                   ! Which time step to read from the state file 


! **********************
! *** INITIALIZATION ***
! **********************

  initcall: IF (calltype == 'ini') THEN
     
     ! Open Netcdf file holding true trajectory
     s = 1
     stat(s) = NF_OPEN(TRIM(file_state), NF_NOWRITE, fileid_state)

     ! Get dimensions
     s = s + 1
     stat(s) = NF_INQ_DIMID(fileid_state, 'timesteps', id_dim)
     s = s + 1
     stat(s) = NF_INQ_DIMLEN(fileid_state, id_dim, nsteps_file)

     ! Read time step information
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid_state, 'step', id_step)

     pos(1) = 1
     cnt(1) = 2
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid_state, id_step, pos, cnt, state_step1and2)
  
     pos(1) = nsteps_file
     cnt(1) = 1
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid_state, id_step, pos, cnt, statefile_laststep)
     s = s + 1

     DO i = 1,  s - 1
        IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error reading trajectory file, no.', i
     END DO

     ! Initialize observation interval in file
     delt_state_file = state_step1and2(2) - state_step1and2(1)

     ! Allocate array for true state
     ALLOCATE(state_true(dim))
     ! count allocated memory
     CALL memcount(3, 'r', dim)

  END IF initcall

  ! Initialize file position corresponding to current time step
  IF (calltype == 'for') THEN
     read_pos = -step / delt_state_file + 1
  ELSE
     read_pos =  step / delt_state_file + 1
  END IF
  write (*,'(8x,a,i6)') &
       '--- Read true state at file position', read_pos


! ******************************
! *** Compute true RMS error ***
! ******************************
  comprms: IF (abs(step) < statefile_laststep) THEN

     ! Read true state
     s = 1
     stat(s) = NF_INQ_VARID(fileid_state, 'state', id_state)

     pos(2) = read_pos
     cnt(2) = 1
     pos(1) = 1
     cnt(1) = dim
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid_state, id_state, pos, cnt, state_true)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading true state from file, no.', i
     END DO

     ! Compute RMS
     trmse = 0.0
     DO j = 1, dim
        trmse = trmse + (state_true(j) - state_est(j))**2
     END DO
     trmse = SQRT(trmse / dim)

     
     ! *** Compute histogram information
     CALL compute_hist(calltype, step, stepnull_means, dim, dim_ens, &
          state_true, state_est, ens, hist_true, hist_mean)

     ! *** Compute statistics of ensemble (skewness, kurtosis)
     CALL compute_stats(calltype, step, stepnull_means, dim, dim_ens, &
          state_true, state_est, ens, skewness, kurtosis)

  ELSE
     ! We don't have true state information any more - set rmse to zero
     trmse = 0.0

  END IF comprms

END SUBROUTINE compute_truermse
!BOP
!
! !ROUTINE: close_netcdf_state  --- close netcdf file for true state
!
! !INTERFACE:
SUBROUTINE close_netcdf_state()

! !DESCRIPTION:
! This routine closes the netcdf file holding the true trajectory.
! It is called in the routine next_observation.

! !USES:
  USE mod_assimilation, &
       ONLY: fileid_state

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
!EOP

! Local variables
  INTEGER :: stat(50)                ! Array for status flag

! Close file
  stat(1) = NF_CLOSE(fileid_state)
  IF (stat(1) /= NF_NOERR) &
       WRITE(*, *) 'NetCDF error in closing file with true states, no. 1'

END SUBROUTINE close_netcdf_state
!BOP
!
! !ROUTINE: compute_hist --- Increment information for rank histogram
!
! !INTERFACE:
SUBROUTINE compute_hist(calltype, step, stepnull_means, dim, dim_ens, &
     state_true, state_est, ens, hist_true, hist_mean)

! !DESCRIPTION:
! Helper routine for the pre/poststep routines. This routines
! increments the rank historgram information based on the 
! ensemble and the true state.
! The rank histogram is computed as follows: The rank at the current
! time is given by how many ensemble members have a first element 
! than is below the first element of the true state.
!
! !REVISION HISTORY:
! 2010-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=3), INTENT(in) :: calltype ! Whether routine is called at 
        !(ini) initial time, (for) after forecast, or (ana) after analysis
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: stepnull_means    ! Step step at which computation is started
  INTEGER, INTENT(in) :: dim               ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens           ! Ensemble size
  REAL, INTENT(in)    :: state_true(dim)   ! True state
  REAL, INTENT(in)    :: state_est(dim)    ! Estimated state
  REAL, INTENT(in)    :: ens(dim, dim_ens) ! State ensemble
  ! Arrays for rank histrograms
  INTEGER, INTENT(inout) :: hist_true(dim_ens+1, 2) ! Histogram about true state
  INTEGER, INTENT(inout) :: hist_mean(dim_ens+1, 2) ! Histogram about ensemble mean

! !CALLING SEQUENCE:
! Called by: compute_truermse
!EOP

! *** local variables ***
  INTEGER :: i           ! Counters
  INTEGER :: rank_ens    ! Rank of current ensemble


! *********************************
! *** Increment rank histograms ***
! *********************************

  IF (calltype=='ana') THEN
     rank_ens = 0
     DO i = 1, dim_ens
        if (ens(1,i) < state_true(1)) rank_ens = rank_ens + 1
     END DO
     hist_true(rank_ens+1, 1) = hist_true(rank_ens+1, 1) + 1
     IF (ABS(step) >= stepnull_means) THEN
        hist_true(rank_ens+1, 2) = hist_true(rank_ens+1, 2) + 1
     END IF
  END IF

  IF (calltype=='ana') THEN
     rank_ens = 0
     DO i = 1, dim_ens
        if (ens(1,i) < state_est(1)) rank_ens = rank_ens + 1
     END DO
     hist_mean(rank_ens+1, 1) = hist_mean(rank_ens+1, 1) + 1
     IF (ABS(step) >= stepnull_means) THEN
        hist_mean(rank_ens+1, 2) = hist_mean(rank_ens+1, 2) + 1
     END IF
  END IF

END SUBROUTINE compute_hist
!BOP
!
! !ROUTINE: compute_stats --- Compute ensemble statistics
!
! !INTERFACE:
SUBROUTINE compute_stats(calltype, step, stepnull_means, dim, dim_ens, &
     state_true, state_est, ens, skewness, kurtosis)

! !DESCRIPTION:
! Helper routine for the pre/poststep routines.
! It computes higher-order ensemble statistics (skewness, kurtosis).
!
! The definition used for kurtosis follows that used by Lawson and Hansen,
! Mon. Wea. Rev. 132 (2004) 1966.
!
! !REVISION HISTORY:
! 2010-09 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:

  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=3), INTENT(in) :: calltype ! Whether routine is called at 
        !(ini) initial time, (for) after forecast, or (ana) after analysis
  INTEGER, INTENT(in) :: step              ! Current time step
  INTEGER, INTENT(in) :: stepnull_means    ! Step step at which computation is started
  INTEGER, INTENT(in) :: dim               ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens           ! Ensemble size
  REAL, INTENT(in)    :: state_true(dim)   ! True state
  REAL, INTENT(in)    :: state_est(dim)    ! Estimated state
  REAL, INTENT(in)    :: ens(dim, dim_ens) ! State ensemble
  REAL, INTENT(out)   :: skewness          ! Skewness of ensemble
  REAL, INTENT(out)   :: kurtosis          ! Kurtosis of ensemble

! !CALLING SEQUENCE:
! Called by: compute_truermse
!EOP

! *** local variables ***
  INTEGER :: i           ! Counters
  REAL :: m2, m3, m4     ! Statistical moments


! **************************
! *** Compute statistics ***
! **************************

  m2 = 0.0
  m3 = 0.0
  m4 = 0.0
  do i=1,dim_ens
     m2 = m2 + (ens(1, i) - state_est(1))**2
     m3 = m3 + (ens(1, i) - state_est(1))**3
     m4 = m4 + (ens(1, i) - state_est(1))**4 
  end do
  m2 = m2 / real(dim_ens)
  m3 = m3 / real(dim_ens)
  m4 = m4 / real(dim_ens)
  
  skewness = m3 / sqrt(m2)**3
  kurtosis = m4 / m2**2 -3.0

END SUBROUTINE compute_stats
