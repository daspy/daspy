! Copyright (c) 2004-2014 Lars Nerger
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id: PDAF_interfaces_module.F90 1525 2014-12-17 12:08:13Z lnerger $
!BOP
!
! !ROUTINE: PDAF_interfaces_module --- Interface definitions for PDAF
!
! !INTERFACE:
MODULE PDAF_interfaces_module

! !DESCRIPTION:
! Module providing interface definition of the PDAF routines that
! are called from the model code.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2012-05 - Lars Nerger - Initial code
! Later revisions - see svn log
!EOP
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  INTERFACE 
     SUBROUTINE PDAF_init(filtertype, subtype, stepnull, param_int, dim_pint, &
          param_real, dim_preal, COMM_model, COMM_filter, COMM_couple, &
          task_id, n_modeltasks, in_filterpe, U_init_ens, in_screen, &
          flag)
       INTEGER, INTENT(in) :: filtertype     ! Type of filter
       INTEGER, INTENT(in) :: subtype        ! Sub-type of filter
       INTEGER, INTENT(in) :: stepnull       ! Initial time step of assimilation
       INTEGER, INTENT(in) :: dim_pint       ! Number of integer parameters
       INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
       INTEGER, INTENT(in) :: dim_preal      ! Number of real parameter 
       REAL, INTENT(inout) :: param_real(dim_preal) ! Real parameter array
       INTEGER, INTENT(in) :: COMM_model     ! Model communicator
       INTEGER, INTENT(in) :: COMM_couple    ! Coupling communicator
       INTEGER, INTENT(in) :: COMM_filter    ! Filter communicator
       INTEGER, INTENT(in) :: task_id        ! Id of my ensemble task
       INTEGER, INTENT(in) :: n_modeltasks   ! Number of parallel model tasks
       LOGICAL, INTENT(in) :: in_filterpe    ! Is my PE a filter-PE?
       INTEGER, INTENT(in) :: in_screen      ! Control screen output:
       INTEGER, INTENT(out):: flag           ! Status flag, 0: no error, error codes:
       EXTERNAL :: U_init_ens  ! User-supplied routine for ensemble initialization
     END SUBROUTINE PDAF_init
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_state_si(nsteps, time, doexit, flag)
       INTEGER, INTENT(inout) :: nsteps  ! Flag and number of time steps
       REAL, INTENT(out)      :: time    ! current model time
       INTEGER, INTENT(inout) :: doexit  ! Whether to exit from forecasts
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_get_state_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_state(steps, time, doexit, U_next_observation, U_distribute_state, &
          U_prepoststep, flag)
       INTEGER, INTENT(inout) :: steps   ! Flag and number of time steps
       REAL, INTENT(out)      :: time    ! current model time
       INTEGER, INTENT(inout) :: doexit  ! Whether to exit from forecasts
       INTEGER, INTENT(inout) :: flag    ! Status flag
       EXTERNAL :: U_next_observation, & ! Routine to provide time step, time and dimension
                                         !   of next observation
            U_distribute_state, &        ! Routine to distribute a state vector
            U_prepoststep                ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_get_state
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seek(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA             ! Provide product R^-1 HV
     END SUBROUTINE PDAF_put_state_seek
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seek_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_seek_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seek(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Routine to provide time step, time and dimension
                                    !   of next observation
     END SUBROUTINE PDAF_assimilate_seek
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seek_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_seek_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_seik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_seik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_seik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seik(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Routine to provide time step, time and dimension
                                    !   of next observation
     END SUBROUTINE PDAF_assimilate_seik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_seik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_seik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_enkf(U_collect_state, U_init_dim_obs, U_obs_op,  &
          U_init_obs, U_prepoststep, U_add_obs_err, U_init_obs_covar, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_add_obs_err, &        ! Add obs error covariance R to HPH in EnKF
            U_init_obs_covar        ! Initialize obs. error cov. matrix R in EnKF
     END SUBROUTINE PDAF_put_state_enkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_enkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_enkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_enkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_add_obs_error, &
          U_init_obs_covar, U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obs_covar, &    ! Initialize obs. error cov. matrix R in EnKF
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_add_obs_error, &     ! Add obs error covariance R to HPH in EnKF
            U_next_observation, &  ! Routine to provide time step, time and dimension
                                   !   of next observation
            U_distribute_state     ! Routine to distribute a state vector
     END SUBROUTINE PDAF_assimilate_enkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_enkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_enkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lseik(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lseik(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation, &   ! Routine to provide time step, time and dimension
                                    !   of next observation
            U_distribute_state      ! Routine to distribute a state vector
     END SUBROUTINE PDAF_assimilate_lseik
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lseik_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lseik_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_etkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_etkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_etkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_etkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_etkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Routine to provide time step, time and dimension
                                    !   of next observation
     END SUBROUTINE PDAF_assimilate_etkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_etkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_etkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_letkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_letkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation, &   ! Routine to provide time step, time and dimension
                                    !   of next observation
            U_distribute_state      ! Routine to distribute a state vector
     END SUBROUTINE PDAF_assimilate_letkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_letkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_letkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_estkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_prepoststep, U_prodRinvA, U_init_obsvar, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_init_dim_obs, &      ! Initialize dimension of observation vector
            U_obs_op, &            ! Observation operator
            U_init_obsvar, &       ! Initialize mean observation error variance
            U_init_obs, &          ! Initialize observation vector
            U_prepoststep, &       ! User supplied pre/poststep routine
            U_prodRinvA            ! Provide product R^-1 A
     END SUBROUTINE PDAF_put_state_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_estkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_estkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_estkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_prepoststep, U_prodRinvA, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_distribute_state, &   ! Routine to distribute a state vector
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_obs_op, &             ! Observation operator
            U_init_obs, &           ! Initialize observation vector
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_prodRinvA, &          ! Provide product R^-1 HV
            U_next_observation      ! Routine to provide time step, time and dimension
                                    !   of next observation
     END SUBROUTINE PDAF_assimilate_estkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_estkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_estkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lestkf(U_collect_state, U_init_dim_obs, U_obs_op, &
          U_init_obs, U_init_obs_l, U_prepoststep, U_prodRinvA_l, U_init_n_domains_p, &
          U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, U_g2l_obs, &
          U_init_obsvar, U_init_obsvar_l, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep           ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lestkf(U_collect_state, U_distribute_state, &
          U_init_dim_obs, U_obs_op, U_init_obs, U_init_obs_l, U_prepoststep, &
          U_prodRinvA_l, U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, &
          U_g2l_state, U_l2g_state, U_g2l_obs, U_init_obsvar, U_init_obsvar_l, &
          U_next_observation, flag)
       INTEGER, INTENT(out) :: flag    ! Status flag
       EXTERNAL :: U_collect_state, &  ! Routine to collect a state vector
            U_obs_op, &             ! Observation operator
            U_init_n_domains_p, &   ! Provide number of local analysis domains
            U_init_dim_l, &         ! Init state dimension for local ana. domain
            U_init_dim_obs, &       ! Initialize dimension of observation vector
            U_init_dim_obs_l, &     ! Initialize dim. of obs. vector for local ana. domain
            U_init_obs, &           ! Initialize PE-local observation vector
            U_init_obs_l, &         ! Init. observation vector on local analysis domain
            U_init_obsvar, &        ! Initialize mean observation error variance
            U_init_obsvar_l, &      ! Initialize local mean observation error variance
            U_g2l_state, &          ! Get state on local ana. domain from full state
            U_l2g_state, &          ! Init full state from state on local analysis domain
            U_g2l_obs, &            ! Restrict full obs. vector to local analysis domain
            U_prodRinvA_l, &        ! Provide product R^-1 A on local analysis domain
            U_prepoststep, &        ! User supplied pre/poststep routine
            U_next_observation, &   ! Routine to provide time step, time and dimension
                                    !   of next observation
            U_distribute_state      ! Routine to distribute a state vector
     END SUBROUTINE PDAF_assimilate_lestkf
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_assimilate_lestkf_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_assimilate_lestkf_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_seik_TtimesA(rank, dim_col, A, B)
       INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
       REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
       REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)
     END SUBROUTINE PDAF_seik_TtimesA
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_etkf_Tleft(dim_ens, dim, A)
       INTEGER, INTENT(in) :: dim_ens      ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim               ! Number of columns in A and B
       REAL, INTENT(inout) :: A(dim_ens, dim)   ! Input/output matrix
     END SUBROUTINE PDAF_etkf_Tleft
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_estkf_OmegaA(rank, dim_col, A, B)
       INTEGER, INTENT(in) :: rank         ! Rank of initial covariance matrix
       INTEGER, INTENT(in) :: dim_col            ! Number of columns in A and B
       REAL, INTENT(in)    :: A(rank, dim_col)   ! Input matrix
       REAL, INTENT(out)   :: B(rank+1, dim_col) ! Output matrix (TA)
     END SUBROUTINE PDAF_estkf_OmegaA
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_enkf_omega(seed, r, dim_ens, omega, norm, &
          otype, screen)
       INTEGER, INTENT(in) :: seed(4)  ! Seed for random number generation
       INTEGER, INTENT(in) :: r        ! Approximated rank of covar matrix
       INTEGER, INTENT(in) :: dim_ens  ! Ensemble size
       REAL, INTENT(inout) :: omega(dim_ens,r)  ! Random matrix
       REAL, INTENT(inout) :: norm     ! Norm for ensemble transformation
       INTEGER, INTENT(in) :: otype    ! Type of omega
       INTEGER, INTENT(in) :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_enkf_omega
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_seik_omega(rank, omega, omegatype, screen)
       INTEGER, INTENT(in) :: rank      ! Approximated rank of covar matrix
       REAL, INTENT(inout) :: omega(rank+1, rank) ! Matrix Omega
       INTEGER, INTENT(in) :: omegatype ! Select type of omega
       INTEGER, INTENT(in) :: screen    ! Verbosity flag
     END SUBROUTINE PDAF_seik_omega
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_memberid(memberid)
       INTEGER,INTENT(inout) :: memberid ! Index in the local ensemble
     END SUBROUTINE PDAF_get_memberid
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_get_obsmemberid(memberid)
       INTEGER,INTENT(inout) :: memberid ! Index in the local observed ensemble
     END SUBROUTINE PDAF_get_obsmemberid
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_local_weight(wtype, rtype, cradius, sradius, distance, &
          nrows, ncols, A, var_obs, weight, verbose)
       INTEGER, INTENT(in) :: wtype      ! Type of weight function
       INTEGER, INTENT(in) :: rtype      ! Type of regulated weighting
       REAL, INTENT(in)    :: cradius    ! Cut-off radius
       REAL, INTENT(in)    :: sradius    ! Support radius 
       REAL, INTENT(in)    :: distance   ! Distance to observation
       INTEGER, INTENT(in) :: nrows      ! Number of rows in matrix A
       INTEGER, INTENT(in) :: ncols      ! Number of columns in matrix A
       REAL, INTENT(in) :: A(nrows, ncols) ! Input matrix
       REAL, INTENT(in)    :: var_obs    ! Observation variance
       REAL, INTENT(out)   :: weight     ! Weights
       INTEGER, INTENT(in) :: verbose    ! Verbosity flag
     END SUBROUTINE PDAF_local_weight
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_incremental(steps, U_dist_stateinc)
       INTEGER, INTENT(in) :: steps ! Time steps over which increment is distributed
       EXTERNAL :: U_dist_stateinc  ! Add state increment during integration
     END SUBROUTINE PDAF_incremental
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_add_increment(dim_p, state_p)
       INTEGER,INTENT(in) :: dim_p          ! State dimension
       REAL,INTENT(inout) :: state_p(dim_p) ! State vector
     END SUBROUTINE PDAF_add_increment
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_get_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(out)       :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_get_smootherens
  END INTERFACE

  INTERFACE 
     SUBROUTINE PDAF_set_smootherens(sens_point, maxlag, status)
       REAL, POINTER, INTENT(out) :: sens_point(:,:,:)  ! Pointer to smoother array
       INTEGER, INTENT(in)        :: maxlag  ! Number of past timesteps processed in sens
       INTEGER, INTENT(out)       :: status  ! Status flag 
     END SUBROUTINE PDAF_set_smootherens
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_prepost_si(flag)
       INTEGER, INTENT(inout) :: flag    ! Status flag
     END SUBROUTINE PDAF_put_state_prepost_si
  END INTERFACE

  INTERFACE
     SUBROUTINE PDAF_put_state_prepost(U_collect_state, U_prepoststep, flag)
       INTEGER, INTENT(out) :: flag   ! Status flag
       EXTERNAL :: U_collect_state, & ! Routine to collect a state vector
            U_prepoststep             ! User supplied pre/poststep routine
     END SUBROUTINE PDAF_put_state_prepost
  END INTERFACE

END MODULE PDAF_interfaces_module
