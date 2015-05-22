!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!   OUTPUT
!     trans(nbv,nbv) : transformation matrix
!=======================================================================
SUBROUTINE letkf_core(Def_Print,nobs,nobsl,nbv,hdxb,rdiag,rloc,dep,parm_infl,trans)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: Def_Print,nobs
  INTEGER,INTENT(IN) :: nobsl
  INTEGER,INTENT(IN) :: nbv	! ensemble size
  REAL(8),INTENT(IN) :: hdxb(1:nobs,1:nbv)
  REAL(8),INTENT(IN) :: rdiag(1:nobs)
  REAL(8),INTENT(IN) :: rloc(1:nobs)
  REAL(8),INTENT(IN) :: dep(1:nobs)
  REAL(8),INTENT(INOUT) :: parm_infl
  REAL(8),INTENT(OUT) :: trans(nbv,nbv)
  REAL(8) :: hdxb_rinv(nobsl,nbv)
  REAL(8) :: eivec(nbv,nbv)
  REAL(8) :: eival(nbv)
  REAL(8) :: pa(nbv,nbv)
  REAL(8) :: work1(nbv,nbv)
  REAL(8) :: work2(nbv,nobsl)
  REAL(8) :: work3(nbv)
  REAL(8) :: rho
  REAL(8) :: parm(4),sigma_o,gain
  REAL(8),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl

  REAL(8),PARAMETER :: relax_alpha = 0.0d0  ! relaxation parameter     !GYL
  REAL(8),PARAMETER :: min_infl = 0.0d0     ! minimum inlfation factor !GYL

  INTEGER :: i,j,k

  IF(nobsl == 0) THEN
    trans = 0.0d0
    DO i=1,nbv
      trans(i,i) = SQRT(parm_infl)
    END DO
    RETURN
  ELSE

!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(hdxb_rinv),maxval(hdxb_rinv)
	print*,"-----------------------------------hdxb^T Rinv hdxb"
END IF
  CALL dgemm('t','n',nbv,nbv,nobsl,1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:),&
  & nobsl,0.0d0,work1,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(work1),maxval(work1)
	print*,"-----------------------------------hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)"
END IF
  IF (min_infl /= 0.0d0 .AND. parm_infl < min_infl) THEN !GYL
    parm_infl = min_infl                                 !GYL
  END IF                                                 !GYL
  rho = 1.0d0 / parm_infl
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,8) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(work1),maxval(work1)
	print*,"-----------------------------------eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]"
END IF
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(eivec),maxval(eivec)
	print*,"-----------------------------------Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv"
END IF
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,pa,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(pa),maxval(pa)
	print*,"----------------------------------Pa hdxb_rinv^T"
END IF
  CALL dgemm('n','t',nbv,nobsl,nbv,1.0d0,pa,nbv,hdxb_rinv,&
    & nobsl,0.0d0,work2,nbv)
!  DO j=1,nobsl
!    DO i=1,nbv
!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
!      DO k=2,nbv
!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(work3),maxval(work3)
	print*,"----------------------------------T = sqrt[(m-1)Pa]"
END IF
  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,8) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO
  CALL dgemm('n','t',nbv,nbv,nbv,1.0d0,work1,nbv,eivec,&
    & nbv,0.0d0,trans,nbv)
!  DO j=1,nbv
!    DO i=1,nbv
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO
!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  IF (relax_alpha /= 0.0d0) THEN            !GYL
    trans = (1.0d0 - relax_alpha) * trans   !GYL
    DO i=1,nbv                              !GYL
      trans(i,i) = relax_alpha + trans(i,i) !GYL
    END DO                                  !GYL
  END IF                                    !GYL
  DO j=1,nbv
    DO i=1,nbv
      trans(i,j) = trans(i,j) + work3(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  Inflation estimation
!-----------------------------------------------------------------------
IF (Def_Print >= 3) THEN
	print*,minval(trans),maxval(trans)
	print*,"----------------------------------Inflation estimation"
END IF
  parm = 0.0d0
  DO i=1,nobsl
    parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
  END DO
  DO j=1,nbv
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(nbv-1,8)
  parm(3) = SUM(rloc(1:nobsl))
  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
  parm_infl = parm_infl + gain * parm(4)

  RETURN
  END IF
END SUBROUTINE letkf_core

!================ letkf with bias correction (DdSM method)============================================================
SUBROUTINE letkf_core_DdSM(nobs,nobsl,nbv,hdxb,rdiag,rloc,dep,alpha,parm_infl,Kb)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  INTEGER,INTENT(IN) :: nbv	! ensemble size
  REAL(8),INTENT(IN) :: hdxb(1:nobs,1:nbv)
  REAL(8),INTENT(IN) :: rdiag(1:nobs)
  REAL(8),INTENT(IN) :: rloc(1:nobs)
  REAL(8),INTENT(IN) :: dep(1:nobs)
  REAL(8),INTENT(IN) :: alpha
  REAL(8),INTENT(INOUT) :: parm_infl
  REAL(8),INTENT(OUT) :: Kb(nbv,nobsl)
  REAL(8),ALLOCATABLE :: hdxb_rinv(:,:)
  REAL(8) :: eivec(nbv,nbv)
  REAL(8) :: eival(nbv)
  REAL(8) :: pa(nbv,nbv)
  REAL(8) :: work1(nbv,nbv)
  REAL(8),ALLOCATABLE :: work2(:,:)
  REAL(8) :: work3(nbv)
  REAL(8) :: rho,ba
  REAL(8) :: parm(4),sigma_o,gain
  REAL(8),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl

  INTEGER :: i,j,k

  IF(nobsl == 0) THEN
    Kb = 0.0d0
    RETURN
  ELSE
  ALLOCATE(hdxb_rinv(nobsl,nbv),work2(nbv,nobsl))

!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
      DO k=2,nobsl
        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / parm_infl
  DO i=1,nbv
    work1(i,i) = (1.0+alpha)*work1(i,i) + REAL(nbv-1,8) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO
  DO j=1,nbv
    DO i=1,nbv
      pa(i,j) = work1(i,1) * eivec(j,1)
      DO k=2,nbv
        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------
  DO j=1,nobsl
    DO i=1,nbv
      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
      DO k=2,nbv
        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
        Kb(i,j) = alpha*work2(i,j)
      END DO
    END DO
  END DO
  DEALLOCATE(hdxb_rinv,work2)
  RETURN
  END IF

END SUBROUTINE letkf_core_DdSM

!=======================================================================
!  Eigenvalue decomposition using subroutine rs
!    INPUT
!      INTEGER :: imode           : mode switch (0: only eiven values)
!      INTEGER :: n               : dimension of matrix
!      REAL(r_size) :: a(n,n)     : input matrix
!    OUTPUT
!      REAL(r_size) :: eival(n)   : eiven values in decending order
!                                   i.e. eival(1) is the largest
!      REAL(r_size) :: eivec(n,n) : eiven vectors
!      INTEGER :: nrank_eff       : number of positive eivenvalues
!=======================================================================
SUBROUTINE mtx_eigen(imode,n,a,eival,eivec,nrank_eff)
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: imode ! 0: calculate only eigen values
  INTEGER,INTENT(IN) :: n
  REAL(8),INTENT(IN) :: a(1:n,1:n)
  REAL(8),INTENT(OUT) :: eival(1:n)
  REAL(8),INTENT(OUT) :: eivec(1:n,1:n)
  INTEGER,INTENT(OUT) :: nrank_eff

  REAL(8) :: a8(n,n)
  REAL(8) :: eival8(n)
  REAL(8) :: eivec8(n,n)
  REAL(8) :: wrk1(n)
  REAL(8) :: wrk2(n)
  INTEGER :: ierr,i,j

  a8 = a
  eivec8 = 0.0d0
  CALL rs(n,n,a8,eival8,imode,eivec8,wrk1,wrk2,ierr)
  IF( ierr/=0 ) THEN
    WRITE(6,'(A)') '!!! ERROR (mtx_eigen): rs error code is ',ierr
    STOP 2
  END IF

  nrank_eff = n
  IF( eival8(n) > 0 ) THEN
    DO i=1,n
      IF( eival8(i) < ABS(eival8(n))*SQRT(EPSILON(eival8)) ) THEN
        nrank_eff = nrank_eff - 1
        eival8(i) = 0.0d0
        eivec8(:,i) = 0.0d0
      END IF
    END DO
  ELSE
    WRITE(6,'(A)') '!!! ERROR (mtx_eigen): All Eigenvalues are below 0'
    STOP 2
  END IF

  IF( nrank_eff<n .AND. eival8(1)/=0 ) THEN
    j = 0
    DO i=n,1,-1
      IF( eival8(i) == 0 ) THEN
        eival8(i) = eival8(n-nrank_eff-j)
        eivec(:,i) = eivec8(:,n-nrank_eff-j)
        eival8(n-nrank_eff-j) = 0.0d0
        eivec8(:,n-nrank_eff-j) = 0.0d0
        j = j+1
      END IF
    END DO
  END IF

  DO i=1,n
    eival(i) = eival8(n+1-i)
    eivec(:,i) = eivec8(:,n+1-i)
  END DO

  RETURN
END SUBROUTINE mtx_eigen
