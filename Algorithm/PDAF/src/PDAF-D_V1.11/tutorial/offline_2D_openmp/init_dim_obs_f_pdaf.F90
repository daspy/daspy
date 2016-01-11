!$Id: init_dim_obs_f_pdaf.F90 1428 2013-09-25 15:22:04Z lnerger $
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
       ONLY : nx, ny, obs, obs_index, coords_obs

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
  REAL, ALLOCATABLE :: obs_field(:,:) ! Array for observation field read from file


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Read observation field form file
  ALLOCATE(obs_field(ny, nx))

  OPEN (12, file='../inputs_offline/obs.txt', status='old')
  DO i = 1, ny
     READ (12, *) obs_field(i, :)
  END DO
  CLOSE (12)

  ! Count observations
  cnt = 0
  DO j = 1, nx
     DO i= 1, ny
        IF (obs_field(i,j) > -999.0) THEN
           cnt = cnt + 1
        END IF
     END DO
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

  cnt = 0
  cnt0 = 0
  DO j = 1, nx
     DO i= 1, ny
        cnt0 = cnt0 + 1
        IF (obs_field(i,j) > -999.0) THEN
           cnt = cnt + 1
           obs_index(cnt) = cnt0
           obs(cnt) = obs_field(i, j)
           coords_obs(1, cnt) = j
           coords_obs(2, cnt) = i
        END IF
     END DO
  END DO


! *** Clean up ***

  DEALLOCATE(obs_field)

END SUBROUTINE init_dim_obs_f_pdaf

