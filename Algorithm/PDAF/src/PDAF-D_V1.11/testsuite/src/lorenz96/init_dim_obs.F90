!$Id: init_dim_obs.F90 836 2010-01-29 17:19:57Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs(step, dim_obs)

! !DESCRIPTION:
! User-supplied routine for PDAF (SEEK/SEIK/EnKF):
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! This variant is for the Lorenz96 model without
! parallelization. The full model state is observed. 
! Thus, the number of observations equals the global 
! state dimension.
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: file_obs, delt_obs_file, obsfile_laststep, have_obs
  USE mod_model, &
       ONLY: dim_state

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs  ! Dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** Local variables ***
  INTEGER :: i, s               ! Counters
  INTEGER :: stat(50)           ! Array for status flag
  INTEGER :: fileid             ! Id of netcdf file
  INTEGER :: id_step            ! ID for step information
  INTEGER :: nsteps_file        ! Number of time steps in file
  INTEGER :: obs_step1and2(2)   ! Array holding first and second time step index
  INTEGER :: pos(2)             ! Position index for writing
  INTEGER :: cnt(2)             ! Count index for writing
  LOGICAL, save :: firstcall = .TRUE. ! Whether routine is called the first time


! ****************************************
! *** Initialize observation dimension ***
! ****************************************

  ! Retrieve information on observations from file
  IF (firstcall) THEN

     ! Read observation information from file
     s = 1
     stat(s) = NF_OPEN(TRIM(file_obs), NF_NOWRITE, fileid)

     ! Read number of time steps in file
     s = s + 1
     stat(s) = NF_INQ_DIMID(fileid, 'timesteps', id_step)
     s = s + 1
     stat(s) = NF_INQ_DIMLEN(fileid, id_step, nsteps_file)

     ! Read time step information
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'step', id_step)

     pos(1) = 1
     cnt(1) = 2
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obs_step1and2)
  
     pos(1) = nsteps_file
     cnt(1) = 1
     s = s + 1
     stat(s) = NF_GET_VARA_INT(fileid, id_step, pos, cnt, obsfile_laststep)

     s = s + 1
     stat(s) = nf_close(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading observation step information, no.', i
     END DO

     ! Initialize observation interval in file
     delt_obs_file = obs_step1and2(2) - obs_step1and2(1)

     firstcall = .FALSE.

  END IF

  ! observation dimension 
  IF (step <= obsfile_laststep) THEN
     dim_obs = dim_state
  ELSE
     dim_obs = 0
  END IF

  ! Set dim_obs to 0, if we are at the end of a forecast phase without obs.
  IF (.NOT.have_obs) dim_obs = 0

END SUBROUTINE init_dim_obs

