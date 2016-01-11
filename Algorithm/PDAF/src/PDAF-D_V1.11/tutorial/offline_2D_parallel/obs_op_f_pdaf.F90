!$Id: obs_op_f_pdaf.F90 1366 2013-04-24 16:25:05Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called in PDAF\_X\_update
! before the loop over all local analysis domains
! is entered.  The routine has to perform the 
! operation of the observation operator acting on 
! a state vector.  The full vector of all 
! observations required for the localized analysis
! on the PE-local domain has to be initialized.
! This is usually data on the PE-local domain plus 
! some region surrounding the PE-local domain. 
! This data is gathered by MPI operations. The 
! gathering has to be done here, since in the loop 
! through all local analysis domains, no global
! MPI operations can be performed, because the 
! number of local analysis domains can vary from 
! PE to PE.
!
! Implementation for the 2D offline example
! with parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: obs_index, local_dims_obs
  USE mod_parallel, &
       ONLY: mype_filter, npes_filter, COMM_filter, MPIerr, MPI_DOUBLE_PRECISION
  

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_obs_op)
! Called by: PDAF_lestkf_update  (as U_obs_op)
! Called by: PDAF_letkf_update   (as U_obs_op)
!EOP

! *** local variables ***
  INTEGER :: i                         ! Counter
  REAL, ALLOCATABLE :: m_state_tmp(:)  ! Temporary process-local state vector
  INTEGER, ALLOCATABLE :: local_dis(:) ! Displacement array for gathering


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  ! Initialize process-local observed state
  ALLOCATE(m_state_tmp(local_dims_obs(mype_filter+1)))

  DO i = 1, local_dims_obs(mype_filter+1)
     m_state_tmp(i) = state_p(obs_index(i))
  END DO

  ! Gather full observed state
  ALLOCATE(local_dis(npes_filter))

  local_dis(1) = 0
  DO i = 2, npes_filter
     local_dis(i) = local_dis(i-1) + local_dims_obs(i-1)
  END DO

  CALL MPI_AllGatherV(m_state_tmp, local_dims_obs(mype_filter+1), MPI_DOUBLE_PRECISION, &
       m_state_f, local_dims_obs, local_dis, MPI_DOUBLE_PRECISION, &
       COMM_filter, MPIerr)


  ! *** Clean up ***

  DEALLOCATE(m_state_tmp, local_dis)

END SUBROUTINE obs_op_f_pdaf
