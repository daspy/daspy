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
       ONLY: nx, ny, local_range, coords_obs, coords_l, obs_index_l

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


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Coordinates of local analysis domain
  coords_l(1) = ceiling(real(domain_p)/real(ny))
  coords_l(2) = domain_p - (coords_l(1)-1)*ny

  !Determine coordinate limits for observation domain
  limits_x(1) = coords_l(1) - local_range
  if (limits_x(1) < 1) limits_x(1) = 1
  limits_x(2) = coords_l(1) + local_range
  if (limits_x(2) > nx) limits_x(2) = nx

  limits_y(1) = coords_l(2) - local_range
  if (limits_y(1) < 1) limits_y(1) = 1
  limits_y(2) = coords_l(2) + local_range
  if (limits_y(2) > ny) limits_y(2) = ny

  ! Count observations within local_range
  dim_obs_l = 0
  DO i = 1, dim_obs_f
     IF (coords_obs(1, i) >= limits_x(1) .AND. coords_obs(1, i) <= limits_x(2) .AND. &
          coords_obs(2, i) >= limits_y(1) .AND. coords_obs(2, i) <= limits_y(2)) THEN
        
        distance = SQRT(REAL((coords_l(1) - coords_obs(1,i))**2 + &
             (coords_l(2) - coords_obs(2,i))**2))
        IF (distance <= REAL(local_range)) dim_obs_l = dim_obs_l + 1

     END IF
  END DO

  ! Initialize index array for local observations in full observed vector
  IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
  ALLOCATE(obs_index_l(dim_obs_l))

  cnt = 0
  DO i = 1, dim_obs_f
     IF (coords_obs(1, i) >= limits_x(1) .AND. coords_obs(1, i) <= limits_x(2) .AND. &
          coords_obs(2, i) >= limits_y(1) .AND. coords_obs(2, i) <= limits_y(2)) THEN

        distance = SQRT(REAL((coords_l(1) - coords_obs(1,i))**2 + &
             (coords_l(2) - coords_obs(2,i))**2))
        IF (distance <= REAL(local_range)) THEN
           cnt = cnt + 1
           obs_index_l(cnt) = i
        END IF
     END IF
  END DO

END SUBROUTINE init_dim_obs_l_pdaf

