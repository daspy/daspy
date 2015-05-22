SUBROUTINE CALL_LETKF_2D(nx,ny,nbv,num_local_obs,msw_infl,eps,R,Correlation_Par,&
    vario_type,pos_sys,pos_obs,xf_in,hxf_whole,h,y_in,xa,innovation,&
    parm_infl,Alpha_Inflation,increments,localization_map,nthreads,&
    Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,&
    Par_Uniform_STD, Par_Sens_Dim, State_DIM_Single_Layer, Def_Localization, &
    GridSize_Sys, Normal_Score_Trans, State_DIM_Single_Column,&
    Bias_Forecast_Model_Option,Bias_Observation_Model_Option,&
    Bias_Model_Dim,Bias_Obs_Dim,Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,&
    minimize_lbfgsb_n,minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol)

    use omp_lib

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE letkf_core(Def_Print,nobs,nobsl,nbv,hdxb,rdiag,rloc,dep,parm_infl,trans)
            INTEGER,INTENT(IN) :: Def_Print,nobs
            INTEGER,INTENT(IN) :: nobsl
            INTEGER,INTENT(IN) :: nbv    ! ensemble size
            REAL(8),INTENT(IN) :: hdxb(1:nobs,1:nbv)
            REAL(8),INTENT(IN) :: rdiag(1:nobs)
            REAL(8),INTENT(IN) :: rloc(1:nobs)
            REAL(8),INTENT(IN) :: dep(1:nobs)
            REAL(8),INTENT(INOUT) :: parm_infl
            REAL(8),INTENT(OUT) :: trans(nbv,nbv)
        END SUBROUTINE
    END INTERFACE

    INTEGER,INTENT(IN) :: nx    ! Number of Model Grids
    INTEGER,INTENT(IN) :: ny    ! Number of the Observations
    INTEGER,INTENT(IN) :: nbv    ! Number of the Ensembles
    INTEGER,INTENT(IN) :: num_local_obs    ! Number of Observations used in the Local Mode
    REAL(4),INTENT(INOUT) :: parm_infl(nx) ! inflation parameter
    REAL(8),INTENT(IN) :: Alpha_Inflation ! inflation parameter for xa
    REAL(8),INTENT(IN) :: GridSize_Sys ! grid cell size
    REAL(8) :: parm
    REAL(8),INTENT(IN) :: msw_infl ! inflation mode switch
    REAL(8),INTENT(IN) :: eps ! truncation parameter for correlation
    REAL(8),INTENT(IN) :: R(ny,ny)   ! Observation Error Covariance
    REAL(8),INTENT(IN) :: Correlation_Par(5,2)    ! Correlation Model Parameter
    CHARACTER(12) :: vario_type    ! Correlation Function Name
    REAL(8),INTENT(IN) :: pos_sys(nx,2)    ! Model Grid Position
    REAL(8),INTENT(IN) :: pos_obs(ny,2)    ! Observation Grid Position
    REAL(8),INTENT(IN) :: xf_in(nx,nbv)
    REAL(8),INTENT(IN) :: hxf_whole(nx,nbv)
    REAL(8),INTENT(OUT) :: xa(nx,nbv)
    REAL(8),INTENT(OUT) :: localization_map(nx)	! record whether we can find the local observations for each grid cell
    INTEGER,INTENT(IN) :: h(ny,nx),Parameter_Optimization_Flag,Par_Sens_Dim,State_DIM_Single_Layer,State_DIM_Single_Column,Def_Localization,Normal_Score_Trans
    REAL(8),INTENT(IN) :: y_in(ny),Parameter_Regularization, Par_Uniform_STD(Par_Sens_Dim)
    INTEGER,INTENT(IN) :: Bias_Forecast_Model_Option, Bias_Observation_Model_Option ! whether to do bias estimation
    INTEGER,INTENT(IN) :: Bias_Model_Dim,Bias_Obs_Dim
    REAL(8),INTENT(IN) :: Bias_Model_Uniform_STD(Bias_Model_Dim),Bias_Obs_Uniform_STD(Bias_Obs_Dim)   ! bias inflation std for model bias and observation bias
    REAL(8),INTENT(IN) :: Model_Inflation_Uniform_STD(nx)   ! model state inflation std

    integer,  intent(in)   :: minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint
    real(8), intent(in)    :: minimize_lbfgsb_factr, minimize_lbfgsb_pgtol, minimize_lbfgsb_epsilon_in(minimize_lbfgsb_n)
    real(8)                :: lmbda(minimize_lbfgsb_n),hlmbda(minimize_lbfgsb_n)

    INTEGER  :: Bias_Model_Dim_LETKF,Bias_Obs_Dim_LETKF

    INTEGER :: h4d(ny,nx)

    REAL(8) :: xm(nx), dxf(nx,nbv), hxf(ny,nbv), hxfm(ny), d(ny), hdxf(ny,nbv), R_Diag(ny), y(ny), hxf_var(ny)
    REAL(8) :: dxf_dev(nx), dxa_dev(nx), xa_xm, dxa_min, dxa_max, U(nbv)
    REAL(8) :: dxa(nx,nbv),dxa_mean
    REAL(8) :: trans(nbv,nbv)
    REAL(8),INTENT(OUT) :: innovation(ny,nbv), increments(nx,nbv)

    INTEGER :: ix,i,j,n,localobs_index,ny_pack,num_local_obs_temp,Par_STD_Index
    INTEGER :: it, ny_loc(nx), ii, jj
    REAL(8) :: dist(ny),coeffs_whole_init,coeffs_whole_max_value
    INTEGER :: thread_id
    INTEGER,INTENT(IN) ::nthreads, Def_Print

    INTEGER :: seed(12), NINT_Correlation

    INTEGER, ALLOCATABLE :: localobs(:) !,coeffs_whole_index(:)
    REAL(8), ALLOCATABLE :: dist_pack(:),d_pack(:),coeffs(:),coeffs_whole(:),coeffs_whole_copy(:),coeffs_whole_eps(:)
    REAL(8), ALLOCATABLE :: coeffs_pack(:),rdiag_pack(:),rloc_pack(:),hdxf_pack(:,:)
    REAL(8) :: Inflation_Fac, xa_median, xa_median_ens(nbv)
    INTEGER ::coeffs_whole_index

        ! Variables for Normal Score Transform
    REAL(8) :: xm_copy(nx), xf(nx,nbv), dxf_copy(nx,nbv)
    REAL(8), ALLOCATABLE :: XF_NNscore(:), XF_Trans(:)
    REAL(8), ALLOCATABLE :: HXF_NNscore(:), HXF_Trans(:)
    REAL(8) :: Analysis_NNscore(nbv), Analysis_Trans(nbv)
	REAL(8) :: Observation_NNscore(ny), Observation_Trans(ny)
    INTEGER :: iwt, ierror
    REAL(8) :: wt(nbv), XF_NNscore_Sort(nbv), XF_Trans_Sort(nbv), HXF_NNscore_Sort(nbv), HXF_Trans_Sort(nbv), Observation_NNscore_Sort(1), Observation_Trans_Sort(1), HXF_Trans_Var, HXF_NNscore_Var
    REAL(8) :: zmin,zmax,ltail,ltpar,utail,utpar
	
    IF (Def_Print >= 3) THEN
        print*,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization
        print*,Par_Uniform_STD, Par_Sens_Dim,State_DIM_Single_Layer
    END IF
	
	
    CALL OMP_SET_NUM_THREADS(nthreads)

    ! because we defined the default values of Bias_Model_Dim and Bias_Obs_Dim are non-zeros, so...
    IF (Bias_Forecast_Model_Option > 0) THEN
        Bias_Model_Dim_LETKF = Bias_Model_Dim
    ELSE
        Bias_Model_Dim_LETKF = 0
    END IF

    IF (Bias_Observation_Model_Option > 0) THEN
        Bias_Obs_Dim_LETKF = Bias_Obs_Dim
    ELSE
        Bias_Obs_Dim_LETKF = 0
    END IF
    ! Parameters for normal score back transformation
    zmin = minval(xf_in)
    zmax = maxval(xf_in)
    ltail=1
    ltpar=1
    utail=1
    utpar=1

    xf = xf_in
    y = y_in

    !$OMP PARALLEL SHARED(R_Diag) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! Diagonal of R
    !
    DO i=1,ny
        R_Diag(i) = R(i,i)
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    IF (Def_Print >= 3) THEN
        PRINT*,"R_Diag"
        PRINT*,minval(R_Diag),maxval(R_Diag)
    END IF
	
    !$OMP PARALLEL SHARED(hxf,Observation_NNscore,HXF_NNscore) PRIVATE(i,j)
    !$OMP DO SCHEDULE(STATIC)
    !---------------
    ! analysis step
    !---------------
    !
    ! hxf = H xf
    !
	
    DO i=1,ny
        Observation_NNscore(i) = y(i)
        DO j=1,nbv
            hxf(i,j) = DOT_PRODUCT(h(i,:), hxf_whole(:,j))
        ENDDO
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
	
    IF (Def_Print >= 3) THEN
        PRINT*,"hxf"
        PRINT*,minval(hxf),maxval(hxf)
    END IF
	
	
    !$OMP PARALLEL SHARED(xm,xm_copy,XF_NNscore) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! ensemble mean -> xm
    !
	DO i=1,nx
        IF (Normal_Score_Trans == 1) THEN
            XF_NNscore(((i-1)*nbv+1):i*nbv) = xf(i,:)			
        END IF
        CALL com_mean(nbv,xf(i,:),xm(i))
        CALL com_mean(nbv,xf_in(i,:),xm_copy(i))
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    IF (Def_Print >= 3) THEN
        PRINT*,"xm"
        PRINT*,minval(xm),maxval(xm)
    END IF
	
    !$OMP PARALLEL SHARED(dxf,dxf_copy,xf,xm) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! ensemble ptb -> dxf
    !
    DO i=1,nbv
        dxf(:,i) = xf(:,i) - xm(:)
        dxf_copy(:,i) = xf_in(:,i) - xm_copy(:)
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    IF (Def_Print >= 3) THEN
        PRINT*,"dxf"
        PRINT*,minval(dxf),maxval(dxf)
    END IF
	
    !$OMP PARALLEL SHARED(hxfm,hxf,hxf_var,R_Diag) PRIVATE(i,Inflation_Fac)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! hxfm = mean(H xf)
    !
    DO i=1,ny
        CALL com_mean(nbv,hxf(i,:),hxfm(i))

        !modify the R according to the model variance
        CALL com_covar(nbv,hxf(i,:),hxf(i,:),hxf_var(i))
        IF (R_Diag(i) > hxf_var(i)) THEN
            Inflation_Fac = R_Diag(i) * max(0.1, min(1.0,(Alpha_Inflation * (hxf_var(i) - R_Diag(i)) / R_Diag(i) + 1.0)))  !RTPS (multiplicative)
            R_Diag(i) = Inflation_Fac
        END IF

    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    IF (Def_Print >= 3) THEN
        PRINT*,"hxfm"
        PRINT*,minval(hxfm),maxval(hxfm)
    END IF

    !
    ! d = y - hxfm
    !
    d = y - hxfm
    IF (Def_Print >= 3) THEN
        PRINT*,"d"
        PRINT*,minval(d),maxval(d)
    END IF
    !$OMP PARALLEL SHARED(hdxf,hxfm) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! hdxf
    !
    DO i=1,nbv
        hdxf(:,i) = hxf(:,i) - hxfm(:)
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    IF (Def_Print >= 3) THEN
        PRINT*,"hdxf"
        PRINT*,minval(hdxf),maxval(hdxf)
    END IF

    !$OMP PARALLEL SHARED(innovation) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    !
    ! innovation
    !
    DO i=1,nbv
        innovation(:,i) = y - hxf(:,i)
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    IF (Def_Print >= 3) THEN
        PRINT*,"innovation"
        PRINT*,minval(innovation),maxval(innovation)
    END IF

    NINT_Correlation = NINT(Correlation_Par(4,1)/REAL(abs(GridSize_Sys),8))
    NINT_Correlation = max(NINT_Correlation,1)

    ! Comment this line, because it does not work under parallelpython
    !print*,"NINT_Correlation",NINT_Correlation,Correlation_Par(4,1),REAL(abs(GridSize_Sys),8)
	
    !$OMP PARALLEL SHARED(xa,parm_infl,dxf,h4d,hxf,hxfm,h,hdxf,increments,localization_map,ny_loc,dxf_dev,dxa_dev,dxa) PRIVATE(ix,xa_median,parm,ii,jj,num_local_obs_temp,dist,coeffs,dist_pack,ny_pack,coeffs_pack,rdiag_pack,rloc_pack,d_pack,hdxf_pack,coeffs_whole,coeffs_whole_eps,coeffs_whole_copy,coeffs_whole_max_value,coeffs_whole_init,coeffs_whole_index,i,localobs,j,it,localobs_index,trans,xa_xm,dxa_min,dxa_max,dxa_mean,seed,U,xa_median_ens, Analysis_NNscore, Analysis_Trans, Par_STD_Index)
    !$OMP DO SCHEDULE(STATIC)

    DO ix=1,nx
        IF (Def_Print >= 3) THEN
            PRINT*, 'The', ix, 'th Grid is being processing!'
        END IF
        !thread_id = omp_get_thread_num()
        !print*,nthreads,ix, thread_id
        ny_loc(ix) = 0
        
        parm = REAL(parm_infl(ix),8)

        ii = pos_sys(ix,1)
        jj = pos_sys(ix,2)

        dist = sqrt(abs(ii - pos_obs(:, 1)) ** 2 + abs(jj - pos_obs(:, 2)) ** 2)
        !print*,"dist",dist
        dist_pack = pack(dist,dist<NINT_Correlation*abs(GridSize_Sys))
        !print*,"dist_pack",dist_pack
        ny_pack = size(dist_pack)
        !print*,"ny_pack",ny_pack
		
        IF (ny_pack > 0) THEN
            IF (Def_Print >= 3) THEN
                PRINT*,"Correlation_Par",Correlation_Par(4,1), vario_type, Correlation_Par(5,1),"GridSize_Sys",GridSize_Sys,"ny_pack",ny_pack,"dist_pack",minval(dist_pack),maxval(dist_pack)
            END IF
			
            ALLOCATE(coeffs_whole(ny_pack))
            ALLOCATE(coeffs_whole_copy(ny_pack))
            ALLOCATE(localobs(ny_pack))
            ALLOCATE(coeffs(ny_pack))
			
            CALL Calc_Loccoeffs_Fortran(Correlation_Par(4,1), vario_type, dist_pack, Correlation_Par(5,1), ny_pack, coeffs_whole)
            !PRINT*,"coeffs_whole",coeffs_whole
			
            IF (Def_Print >= 3) THEN
                PRINT*,"coeffs_whole",coeffs_whole
                PRINT*,eps, maxval(coeffs_whole, MASK = coeffs_whole >= eps)
            END IF
			
            coeffs_whole_max_value = maxval(coeffs_whole, MASK = (coeffs_whole >= eps))
	        
            ! Normalization (I think it is wrong for the non-observed grid cells)
            coeffs_whole = coeffs_whole / coeffs_whole_max_value
            !print*,minval(coeffs_whole, MASK = coeffs_whole > 0)
			
            localobs(:) = 0
            coeffs_whole_copy(:) = coeffs_whole(:)
	        
            coeffs_whole_eps = pack(coeffs_whole,coeffs_whole >= eps)
	        
            IF (size(pack(dist,dist < abs(GridSize_Sys))) > 0) THEN
                num_local_obs_temp = 1
            ELSE IF (size(coeffs_whole_eps) > num_local_obs) THEN
                num_local_obs_temp = num_local_obs
            ELSE
                num_local_obs_temp = size(coeffs_whole_eps)
            END IF
	        
            IF (num_local_obs_temp > 0) THEN
                DO WHILE (count(localobs /= 0) < num_local_obs_temp)
                    coeffs_whole_index = maxloc(coeffs_whole_copy(:), DIM=1, MASK = coeffs_whole_copy >= eps)
                    localobs(coeffs_whole_index) = 1
                    coeffs_whole_copy(coeffs_whole_index) = 0.0
                END DO
            END IF
			
            !print*,localobs
            ny_loc(ix) = num_local_obs_temp
            IF (Def_Print >= 3) THEN
                print*, 'There are', ny_loc(ix), 'Observations Used!'
                print*, 'minval(coeffs_whole)',minval(coeffs_whole),'maxval(coeffs_whole)',maxval(coeffs_whole),'min(dist_pack)',minval(dist_pack),'max(dist_pack)',maxval(dist_pack),"Correlation_Par(:,1)",Correlation_Par(:,1)
            END IF
			
            IF (ny_loc(ix) > 0) THEN
				
                ALLOCATE(coeffs_pack(num_local_obs_temp))
                ALLOCATE(hdxf_pack(num_local_obs_temp,nbv))
                ALLOCATE(rdiag_pack(num_local_obs_temp))
                ALLOCATE(rloc_pack(num_local_obs_temp))
                ALLOCATE(d_pack(num_local_obs_temp))
				
                d_pack = pack(pack(d,dist<NINT_Correlation*abs(GridSize_Sys)),localobs == 1)
                coeffs_pack = pack(coeffs_whole,localobs == 1)
                rdiag_pack = pack(pack(R_Diag,dist<NINT_Correlation*abs(GridSize_Sys)),localobs == 1)
		    	
                IF (Def_Localization == 1) THEN
                    rdiag_pack = rdiag_pack / coeffs_pack
                END IF
		    	
                rloc_pack(:) = 1.0
		        
                !print*,"ix",ix,"dist_pack",pack(pack(dist,dist<NINT_Correlation*abs(GridSize_Sys)),localobs == 1),"d_pack",d_pack,"coeffs_pack",coeffs_pack,"rdiag_pack",rdiag_pack
                !IF(IX>=1000) THEN
                !	CALL ABORT
                !END IF
		    	
                DO i=1,nbv
                    hdxf_pack(:,i) = pack(pack(hdxf(:,i),dist<NINT_Correlation*abs(GridSize_Sys)),localobs == 1)
                END DO

                IF (Def_Print >= 3) THEN
                    print*,"letkf",hdxf_pack(:,:),rdiag_pack(:),rloc_pack(:),d_pack(:),parm
                END IF
		        
                IF(parm == 0.0d0) parm = abs(msw_infl)
		        
                CALL letkf_core(Def_Print,num_local_obs_temp,ny_loc(ix),nbv,hdxf_pack,rdiag_pack,rloc_pack,d_pack,parm,trans)
				
                IF (Def_Print >= 3) THEN
                    print*,"******************Finish letkf_core**************************************************"
                END IF
		        
                IF(msw_infl > 0.0d0) parm = msw_infl
                DO j=1,nbv
                    xa(ix,j) = xm(ix)
                    DO i=1,nbv
                        !                if ((ix == 2770) .or. (ix == 2768)) then
                        !                    print*,dxf(ix,i), trans(i,j)
                        !                end if
                        IF ((Parameter_Optimization_Flag>0) .AND. (ix>State_DIM_Single_Layer)) THEN
                            increments(ix,j) = increments(ix,j) + Parameter_Regularization * dxf(ix,i) * trans(i,j)
                            xa(ix,j) = xa(ix,j) + Parameter_Regularization * dxf(ix,i) * trans(i,j)
                        ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Forecast_Model_Option > 0) .AND. (ix>(nx-Bias_Model_Dim_LETKF-Bias_Obs_Dim_LETKF)) .AND. ((Bias_Model_Dim_LETKF >= 1))) THEN
                            increments(ix,j) = increments(ix,j) + dxf(ix,i) * trans(i,j)
                            xa(ix,j) = xa(ix,j) + dxf(ix,i) * trans(i,j)
                        ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Observation_Model_Option > 0) .AND. (ix>(nx-Bias_Model_Dim_LETKF)) .AND. ((Bias_Obs_Dim_LETKF >= 1))) THEN
                            increments(ix,j) = increments(ix,j) + dxf(ix,i) * trans(i,j)
                            xa(ix,j) = xa(ix,j) + dxf(ix,i) * trans(i,j)
                        ELSE
                            increments(ix,j) = increments(ix,j) + dxf(ix,i) * trans(i,j)
                            xa(ix,j) = xa(ix,j) + dxf(ix,i) * trans(i,j)
                        END IF
                    END DO
                    IF(ISNAN(increments(ix,j))) THEN
                        xa(ix,j) = xf(ix,j)
                    END IF
                END DO
                parm_infl(ix) = min(max(parm,1.0),1.1)
                localization_map(ix) = 1

                CALL com_mean(nbv,xa(ix,:),xa_xm)
                dxa(ix,:) = xa(ix,:) - xa_xm

                IF (Def_Print >= 3) THEN
                    print*,"=====================================Analsys is",xm(ix)
                    print*,"minval(dxf_copy(ix,:))",minval(dxf_copy(ix,:)),"maxval(dxf_copy(ix,:))",maxval(dxf_copy(ix,:))
                    print*,"minval(dxa(ix,:))",minval(dxa(ix,:)),"maxval(dxa(ix,:))",maxval(dxa(ix,:))
                END IF
			    
                                !dxa_min = minval(dxa(ix,:))
                                !dxa_max = maxval(dxa(ix,:))
                                !seed(:) = 123456
                                !CALL RANDOM_SEED(PUT = seed)
                                !CALL RANDOM_NUMBER(U)
                                !dxa(ix,:) = -1.0 * (dxa_max - dxa_min) / 2.0 * (1 - U) + (dxa_max - dxa_min) / 2.0 * U

                seed(:) = 123456
                CALL RANDOM_SEED(PUT = seed)
                CALL RANDOM_NUMBER(U)
                U = 0.8 * (1.0 - U) + 1.2 * U
                !print*,"U",U

                CALL com_stdev(nbv,dxf_copy(ix,:),dxf_dev(ix))
                CALL com_stdev(nbv,dxa(ix,:),dxa_dev(ix))
				
                IF (Def_Print >= 3) THEN
                    print*,"Parameter_Optimization_Flag",Parameter_Optimization_Flag,Alpha_Inflation,ix,State_DIM_Single_Layer,ny_loc(ix),dxa_dev(ix),Par_Sens_Dim
                    IF (ix>State_DIM_Single_Layer) THEN
                        print*,"inflation",xa_xm,dxa_dev(ix)
                    ENDIF
                ENDIF

                !print*,ix,nx,Bias_Forecast_Model_Option,Bias_Observation_Model_Option,Bias_Model_Dim_LETKF,Bias_Obs_Dim_LETKF
                ! Parameter Inflation
                !IF((Parameter_Optimization_Flag>0) .AND. (Alpha_Inflation .NE. 0.0) .AND. (ix>State_DIM_Single_Layer) .AND. (dxa_dev(ix) .GT. 0.0) .AND. (dxa_dev(ix) < dxf_dev(ix))) THEN
                IF((Parameter_Optimization_Flag>0) .AND. (Alpha_Inflation .NE. 0.0) .AND. (Par_Sens_Dim >= 1) .AND. (ix>State_DIM_Single_Layer)) THEN
                    !print*,"inflation",ix,(ix-State_DIM_Single_Column*State_DIM_Single_Layer)/State_DIM_Single_Layer+1,State_DIM_Single_Column

                    IF ((State_DIM_Single_Layer) > 1 .AND. (ix>State_DIM_Single_Column*State_DIM_Single_Layer)) THEN
                        !Par_STD_Index = (ix-State_DIM_Single_Column*State_DIM_Single_Layer-1)/State_DIM_Single_Layer+1
                        Par_STD_Index = (ix-1-State_DIM_Single_Column*State_DIM_Single_Layer)/State_DIM_Single_Layer+1
                    ELSE
                        Par_STD_Index = ix-State_DIM_Single_Column
                    END iF

                    !Inflation
                    IF (dxa_dev(ix) < Par_Uniform_STD(Par_STD_Index)) THEN
                        IF (Def_Print >= 3) THEN
                            print*,"inflation",ix/State_DIM_Single_Layer - State_DIM_Single_Column + 1,xa_xm,dxa(ix,:),Par_Uniform_STD(Par_STD_Index),dxa_dev(ix)
                        ENDIF
                        IF (dxa_dev(ix) .GT. 0.0) THEN
                            xa(ix,:) = xa_xm + dxa(ix,:) * max(1.0, min(2.0, (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(ix)) / dxa_dev(ix) + 1.0)))  !RTPS (multiplicative)
                        ELSE
                            seed(:) = 123456
                            CALL RANDOM_SEED(PUT = seed)
                            CALL RANDOM_NUMBER(U)
                            U = 0.8 * (1.0 - U) + 1.2 * U
                            xa(ix,:) = xa(ix,:) * U
                        END IF
                        !print*,xa_xm,ix,State_DIM_Single_Layer,Par_Uniform_STD(mod(ix,State_DIM_Single_Layer))
                        !xa(ix,:) = xa(ix,:) + dxa(ix,:) * (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(ix)) / dxa_dev(ix))  !RTPS (multiplicative)
                        !xa(ix,:) = xa(ix,:) * (Alpha_Inflation * (Par_Uniform_STD(Par_STD_Index) - dxa_dev(ix)) / dxa_dev(ix) + 1.0)  !RTPS (multiplicative)
                       !xa(ix,:) = xa_xm + dxa(ix,:) * (1.0 - Alpha_Inflation) + Alpha_Inflation * dxf(ix,:)    !RTPP (multiplicative+additave)
                       !print*,xa(ix,:), xa_xm
                       !print*,dxf(ix,:)
                       !print*,"**************************************"
                    END IF

                ! State Inflation
                ELSE IF ((Parameter_Optimization_Flag == 0) .AND. (Bias_Forecast_Model_Option == 0) .AND. (Bias_Observation_Model_Option == 0) .AND. (Alpha_Inflation .NE. 0.0) .AND. (dxa_dev(ix) < dxf_dev(ix))) THEN
                    IF (Def_Print >= 3) THEN
                        print*,"inflation",dxf_dev(ix),dxa_dev(ix),Alpha_Inflation * (dxf_dev(ix) - dxa_dev(ix)) / dxa_dev(ix) + 1.0
                    END IF

                    IF (dxa_dev(ix) .GT. 0.0) THEN
						!print*,"before,xa_xm",xa_xm
                        xa(ix,:) = xa_xm + dxa(ix,:) * max(1.0, min(2.0, (Alpha_Inflation * (dxf_dev(ix) - dxa_dev(ix)) / dxa_dev(ix) + 1.0)))  !RTPS (multiplicative)
						!CALL com_mean(nbv,xa(ix,:),xa_xm)
						!print*,"After,xa_xm",xa_xm
					ELSE
						!CALL com_mean(nbv,xa(ix,:),dxa_mean)
						!print*,"before,dxa_mean",dxa_mean
                        seed(:) = 123456
                        CALL RANDOM_SEED(PUT = seed)
                        CALL RANDOM_NUMBER(U)
                        U = 0.8 * (1.0 - U) + 1.2 * U
                        xa(ix,:) = xa(ix,:) * U
						!CALL com_mean(nbv,xa(ix,:),dxa_mean)
						!print*,"After,dxa_mean",dxa_mean
                    END IF

                END IF

                DEALLOCATE (coeffs_pack)
                DEALLOCATE (hdxf_pack)
                DEALLOCATE (rdiag_pack)
                DEALLOCATE (rloc_pack)
                DEALLOCATE (d_pack)

            ELSE
                DO j=1,nbv
                    xa(ix,j) = xf_in(ix,j)
                END DO
                localization_map(ix) = 0

            END IF

            DEALLOCATE (coeffs_whole)
            DEALLOCATE (coeffs_whole_copy)
            DEALLOCATE (localobs)
            DEALLOCATE (coeffs)
            DEALLOCATE (coeffs_whole_eps)

        ELSE
            DO j=1,nbv
                xa(ix,j) = xf_in(ix,j)
            END DO
            localization_map(ix) = 0
        END IF

        DEALLOCATE (dist_pack)
		
    END DO    ! End of nx
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
	
    RETURN

END SUBROUTINE CALL_LETKF_2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Calc_Loccoeffs_Fortran(radius, tag, dist, kappa, length, coeffs)

    IMPLICIT NONE

    INTERFACE
        SUBROUTINE Exponential_Fortran(length,h,r,correlation)
            INTEGER,INTENT(IN) :: length
            REAL(8),INTENT(IN) :: h(length)
            REAL(8),INTENT(IN) :: r
            REAL(8),INTENT(OUT) :: correlation(length)
        END SUBROUTINE Exponential_Fortran
        SUBROUTINE Spherical_Fortran(length,h,r,correlation)
            INTEGER,INTENT(IN) :: length
            REAL(8),INTENT(IN) :: h(length)
            REAL(8),INTENT(IN) :: r
            REAL(8),INTENT(OUT) :: correlation(length)
        END SUBROUTINE Spherical_Fortran
        SUBROUTINE Gaussian_Fortran(length,h,r,correlation)
            INTEGER,INTENT(IN) :: length
            REAL(8),INTENT(IN) :: h(length)
            REAL(8),INTENT(IN) :: r
            REAL(8),INTENT(OUT) :: correlation(length)
        END SUBROUTINE Gaussian_Fortran
        SUBROUTINE SteMat_Fortran(length,h,r,kappa,correlation)
            INTEGER,INTENT(IN) :: length
            REAL(8),INTENT(IN) :: h(length)
            REAL(8),INTENT(IN) :: r,kappa
            REAL(8),INTENT(OUT) :: correlation(length)
        END SUBROUTINE SteMat_Fortran
        SUBROUTINE Matern_Fortran(length,h,r,kappa,correlation)
            INTEGER,INTENT(IN) :: length
            REAL(8),INTENT(IN) :: h(length)
            REAL(8),INTENT(IN) :: r,kappa
            REAL(8),INTENT(OUT) :: correlation(length)
        END SUBROUTINE Matern_Fortran
    END INTERFACE

    REAL(8),INTENT(IN) :: radius,kappa
    CHARACTER(12),INTENT(IN) :: tag
    INTEGER,INTENT(IN)::length
    REAL(8),INTENT(IN) :: dist(length)
    REAL(8),INTENT(OUT) :: coeffs(length)

    IF (tag == 'Gauss') THEN
        coeffs = exp(-0.5 * (dist / radius) ** 2)

    ELSE IF (tag == 'Matern') THEN
        CALL Matern_Fortran(length,dist,radius,kappa,coeffs)

    ELSE IF (tag == 'Exponential') THEN
        CALL Exponential_Fortran(length,dist,radius,coeffs)

    ELSE IF (tag == 'Gaussian') THEN
        CALL Gaussian_Fortran(length,dist,radius,coeffs)

    ELSE IF (tag == 'Spherical') THEN
        CALL Spherical_Fortran(length,dist,radius,coeffs)

    ELSE IF (tag == 'Gaspari_Cohn') THEN
        CALL Gaspari_Cohn_Fortran(length,dist,radius,coeffs)

    ELSE
        print*, 'Wrong coorelation type name!!!'
        
    END IF
    
    RETURN

END SUBROUTINE Calc_Loccoeffs_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Exponential_Fortran(length,h_in,r,correlation)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: length
    REAL(8),INTENT(IN) :: h_in(length)
    REAL(8),INTENT(IN) :: r
    REAL(8),INTENT(OUT) :: correlation(length)
    REAL(8) :: h(length)
    INTEGER::i
	
    h = h_in
	
    !$OMP PARALLEL SHARED(correlation) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,length
        IF(h(i).LE.0.0) then
            correlation(i)=1.0
        ELSE 
            correlation(i)=exp(-1.0*h(i)/r)
        END IF
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    RETURN

END SUBROUTINE Exponential_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Spherical_Fortran(length,h_in,r,correlation)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: length
    REAL(8),INTENT(IN) :: h_in(length)
    REAL(8),INTENT(IN) :: r
    REAL(8),INTENT(OUT) :: correlation(length)
    REAL(8) :: h(length)
    REAL(8) :: hphi
    INTEGER :: i
    
    h = h_in
    
    !$OMP PARALLEL SHARED(correlation) PRIVATE(i,hphi)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,length
        IF(h(i).LE.0.0) then
            correlation(i)=1.0
        ELSE IF(h(i) >=r) then
            correlation(i)=0.0
        ELSE
            hphi=h(i)/r
            correlation(i)=1.0 - (1.5 * hphi) + (0.5 * (hphi ** 3))
        END IF
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    RETURN

END SUBROUTINE Spherical_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Gaussian_Fortran(length,h_in,r,correlation)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: length
    REAL(8),INTENT(IN) :: h_in(length)
    REAL(8),INTENT(IN) :: r
    REAL(8),INTENT(OUT) :: correlation(length)
    REAL(8) :: h(length)
    REAL(8) :: hphi
    INTEGER :: i
	
    h = h_in
	
    !$OMP PARALLEL SHARED(correlation) PRIVATE(i,hphi)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,length
        IF(h(i).LE.0.0) then
            correlation(i)=1.0
        ELSE 
            hphi=h(i)/r
            correlation(i)=exp(-1.0*(hphi * hphi))
        END IF
    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    RETURN

END SUBROUTINE Gaussian_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Gaspari_Cohn_Fortran(length,h_in,r,correlation)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: length
    REAL(8),INTENT(IN) :: h_in(length)
    REAL(8),INTENT(IN) :: r
    REAL(8),INTENT(OUT) :: correlation(length)
    REAL(8) :: h(length)
    REAL(8) :: cfaci
    INTEGER :: i

    h = h_in

    cfaci = REAL(r) / 2.0

    !$OMP PARALLEL SHARED(correlation) PRIVATE(i)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,length
        IF(h(i).LE.0.0) then
            correlation(i)=1.0
        ELSEIF (h(i) <= r / 2) THEN
            correlation(i) = -0.25 * (h(i) / cfaci)**5 &
                + 0.5 * (h(i) / cfaci)**4 &
                + 5.0 / 8.0 * (h(i) / cfaci)**3 &
                - 5.0 / 3.0 * (h(i) / cfaci)**2 + 1.0
        ELSEIF (h(i) > r / 2 .AND. h(i) < r) THEN
            correlation(i) = 1.0 / 12.0 * (h(i) / cfaci)**5 &
                - 0.5 * (h(i) / cfaci)**4 &
                + 5.0 / 8.0 * (h(i) / cfaci)**3 &
                + 5.0 / 3.0 * (h(i) / cfaci)**2 &
                - 5.0 * (h(i) / cfaci) &
                + 4.0 - 2.0 / 3.0 * cfaci / h(i)
        ELSE
            correlation(i) = 0.0
        ENDIF

    END DO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    RETURN

END SUBROUTINE Gaspari_Cohn_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Matern_Fortran(length,h_in,r,kappa,correlation)

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: length
    REAL(8),INTENT(IN) :: h_in(length)
    REAL(8),INTENT(IN) :: r,kappa
    REAL(8),INTENT(OUT) :: correlation(length)
    REAL(8) :: h(length)
    REAL(8) :: hphi, alpha, cte
    INTEGER :: i, nb, ize, ncalc

    REAL(8) :: bk(10)
	
    h = h_in
	
    alpha = kappa
    IF (alpha < 0.0) THEN
        alpha = alpha * (-1.0)
    ENDIF
    nb = 1 + floor(alpha) !/* nb-1 <= |alpha| < nb */
    alpha = alpha - (nb-1)
    ize = 1
	
    !print*,"***********************h",maxval(h),minval(h)
    !print*,alpha,length,h,r,kappa
	
    !$OMP PARALLEL SHARED(correlation) PRIVATE(i,hphi,cte,bk,ncalc)
    !$OMP DO SCHEDULE(STATIC)
    DO i=1,length
        hphi = h(i) / r
        IF(hphi.eq.0.0) THEN
            correlation(i)=1.0
        ELSE IF (kappa .eq. 0.5) THEN
            correlation(i) = exp(-1.0 * hphi)
            correlation(i) = 0.0
        ELSE
            cte = (2.0 ** (-1.0 * (kappa - 1.0)))/gamma(kappa)
            !print*,"alpha",alpha
            call rkbesl(hphi, alpha, nb, ize, bk, ncalc)
            !print*,ncalc
            IF(ncalc .NE. nb) THEN
                PRINT*,"Matern_Fortran Failed!!!!!!!!!!",i,h(i),hphi,-1.0*h(i),r,exp(-1.0*h(i)/r),"***********************h",maxval(h),minval(h)
                correlation(i) = exp(-1.0*h(i)/r)
            ELSE
                correlation(i) = cte * (hphi ** kappa) * bk(nb)
            ENDIF
        END IF
        IF(ncalc .NE. nb) THEN
            print*,i,length,hphi,h(i),r,kappa,alpha,nb,cte,correlation(i),bk(nb),ncalc
        END IF
    END DO
    !print*,i,length,hphi,r,kappa,alpha,nb,cte,correlation(i),bk(nb),ncalc

    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    RETURN

END SUBROUTINE Matern_Fortran

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LOGICAL FUNCTION isfinite(a)
    REAL(8) :: a
    isfinite = (a-a).EQ.0
    RETURN
END FUNCTION isfinite

!-----------------------------------------------------------------------
! Mean
!-----------------------------------------------------------------------
SUBROUTINE com_mean(ndim,var,amean)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    REAL(8),INTENT(IN) :: var(ndim)
    REAL(8),INTENT(OUT) :: amean

    INTEGER :: i

    amean = 0.0d0
    DO i=1,ndim
        amean = amean + var(i)
    END DO
    amean = amean / REAL(ndim,8)

    RETURN
END SUBROUTINE com_mean

!-----------------------------------------------------------------------
! Standard deviation
!-----------------------------------------------------------------------
SUBROUTINE com_stdev(ndim,var,aout)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    REAL(8),INTENT(IN) :: var(ndim)
    REAL(8),INTENT(OUT) :: aout

    REAL(8) :: amean
    REAL(8) :: dev(ndim)

    CALL com_mean(ndim,var,amean)

    dev(:) = var(:) - amean

    aout = SQRT( SUM(dev*dev) / REAL(ndim-1,8) )

    RETURN
END SUBROUTINE com_stdev
!-----------------------------------------------------------------------
! Covariance
!-----------------------------------------------------------------------
SUBROUTINE com_covar(ndim,var1,var2,cov)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: ndim
    REAL(8),INTENT(IN) :: var1(ndim)
    REAL(8),INTENT(IN) :: var2(ndim)
    REAL(8),INTENT(OUT) :: cov

    REAL(8) :: amean1,amean2
    REAL(8) :: dev1(ndim),dev2(ndim)

    CALL com_mean(ndim,var1,amean1)
    CALL com_mean(ndim,var2,amean2)

    dev1(:) = var1(:) - amean1
    dev2(:) = var2(:) - amean2

    cov = SUM( dev1*dev2 ) / REAL(ndim-1,8)

    RETURN
END SUBROUTINE com_covar

subroutine rkbesl ( x, alpha, nb, ize, bk, ncalc )

    !*****************************************************************************80
    !
    !! RKBESL calculates K Bessel function with non-integer orders.
    !
    !  Discussion:
    !
    !    This routine calculates modified Bessel functions of the second
    !    kind, K SUB(N+ALPHA) (X), for non-negative argument X, and
    !    non-negative order N+ALPHA, with or without exponential scaling.
    !
    !    This program is based on a program written by J. B. Campbell
    !    that computes values of the Bessel functions K of real
    !    argument and real order.  Modifications include the addition
    !    of non-scaled functions, parameterization of machine
    !    dependencies, and the use of more accurate approximations
    !    for SINH and SIN.
    !
    !    In case of an error, NCALC .NE. NB, and not all K's are
    !    calculated to the desired accuracy.
    !
    !    NCALC < -1:  An argument is out of range. For example,
    !    NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
    !    XMAX.  In this case, the B-vector is not calculated,
    !    and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
    !
    !    NCALC = -1:  Either  K(ALPHA,X) .GE. XINF  or
    !    K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) .GE. XINF.  In this case,
    !    the B-vector is not calculated.  Note that again
    !    NCALC .NE. NB.
    !
    !    0 < NCALC < NB: Not all requested function values could
    !    be calculated accurately.  BK(I) contains correct function
    !    values for I <= NCALC, and contains the ratios
    !    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    JB Campbell,
    !    On Temme's Algorithm for the Modified Bessel Functions of the
    !    Third Kind,
    !    ACM Transactions on Mathematical Software,
    !    Volume 6, Number 4, December 1980, pages 581-586.
    !
    !    JB Campbell,
    !    A FORTRAN IV Subroutine for the Modified Bessel Functions of
    !    the Third Kind of Real Order and Real Argument,
    !    Report NRC/ERB-925,
    !    National Research Council, Canada.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the non-negative argument for which
    !    K's or exponentially scaled K's (K*EXP(X))
    !    are to be calculated.  If K's are to be calculated,
    !    X must not be greater than XMAX.
    !
    !    Input, real ( kind = 8 ) ALPHA, the fractional part of order for which
    !    K's or exponentially scaled K's (K*EXP(X)) are to be calculated.
    !    0 <= ALPHA < 1.0.
    !
    !    Input, integer ( kind = 4 ) NB, the number of functions to be calculated,
    !    NB .GT. 0.  The first function calculated is of order ALPHA, and the
    !    last is of order (NB - 1 + ALPHA).
    !
    !    Input, integer ( kind = 4 ) IZE, scaling option.
    !    1, unscaled functions are to calculated,
    !    2, exponentially scaled functions are to be calculated.
    !
    !    Output, real ( kind = 8 ) BK(NB), the results.  If the routine
    !    terminates normally, with NCALC = NB, the vector BK contains the
    !    functions K(ALPHA,X), ... , K(NB-1+ALPHA,X), or the corresponding
    !    exponentially scaled functions.
    !    If (0 < NCALC < NB), BK(I) contains correct function
    !    values for I <= NCALC, and contains the ratios
    !    K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
    !
    !    Output, integer ( kind = 4 ) NCALC, error indicator.  If NCALC = NB, then
    !    all the requested values were calculated to the desired accuracy.
    !
    implicit none

    real ( kind = 8 ) a
    real ( kind = 8 ) alpha
    real ( kind = 8 ) blpha
    real ( kind = 8 ),intent(inout) :: bk(nb)
    real ( kind = 8 ) bk1
    real ( kind = 8 ) bk2
    real ( kind = 8 ) c
    real ( kind = 8 ) d
    real ( kind = 8 ) dm
    real ( kind = 8 ) d1
    real ( kind = 8 ) d2
    real ( kind = 8 ) d3
    real ( kind = 8 ) enu
    real ( kind = 8 ) eps
    real ( kind = 8 ) estf(7)
    real ( kind = 8 ) estm(6)
    real ( kind = 8 ) ex
    real ( kind = 8 ) four
    real ( kind = 8 ) f0
    real ( kind = 8 ) f1
    real ( kind = 8 ) f2
    real ( kind = 8 ) half
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iend
    integer ( kind = 4 ) itemp
    integer ( kind = 4 ) ize
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mplus1
    integer ( kind = 4 ) nb
    integer ( kind = 4 ),intent(inout) ::  ncalc
    real ( kind = 8 ) one
    real ( kind = 8 ) p(8)
    real ( kind = 8 ) p0
    real ( kind = 8 ) q(7)
    real ( kind = 8 ) q0
    real ( kind = 8 ) r(5)
    real ( kind = 8 ) ratio
    real ( kind = 8 ) s(4)
    real ( kind = 8 ) sqxmin
    real ( kind = 8 ) t(6)
    real ( kind = 8 ) tinyx
    real ( kind = 8 ) two
    real ( kind = 8 ) twonu
    real ( kind = 8 ) twox
    real ( kind = 8 ) t1
    real ( kind = 8 ) t2
    real ( kind = 8 ) wminf
    real ( kind = 8 ) x
    real ( kind = 8 ) xinf
    real ( kind = 8 ) xmax
    real ( kind = 8 ) xmin
    real ( kind = 8 ) x2by4
    real ( kind = 8 ) zero
    !
    !  Mathematical constants
    !    A = LOG(2.D0) - Euler's constant
    !    D = SQRT(2.D0/PI)
    !
    data half / 0.5d0 /
    data one / 1.0d0 /
    data two / 2.0d0 /
    data zero /0.0d0/
    data four /4.0d0 /
    data tinyx / 1.0d-10/
    data a / 0.11593151565841244881d0/
    data d /0.797884560802865364d0/
    !
    !  Machine dependent parameters
    !
    data eps /2.22d-16/
    data sqxmin /1.49d-154/
    data xinf /1.79d+308/
    data xmin /2.23d-308/
    data xmax /705.342d0/
    !
    !  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
    !                                         + Euler's constant
    !  Coefficients converted from hex to decimal and modified
    !  by W. J. Cody, 2/26/82
    !
    !  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
    !  T    - Approximation for SINH(Y)/Y
    !
    data p/ 0.805629875690432845d00,    0.204045500205365151d02, &
        0.157705605106676174d03,    0.536671116469207504d03, &
        0.900382759291288778d03,    0.730923886650660393d03, &
        0.229299301509425145d03,    0.822467033424113231d00/
    data q/ 0.294601986247850434d02,    0.277577868510221208d03, &
        0.120670325591027438d04,    0.276291444159791519d04, &
        0.344374050506564618d04,    0.221063190113378647d04, &
        0.572267338359892221d03/
    data r/-0.48672575865218401848d+0,  0.13079485869097804016d+2, &
        -0.10196490580880537526d+3,  0.34765409106507813131d+3, &
        0.34958981245219347820d-3/
    data s/-0.25579105509976461286d+2,  0.21257260432226544008d+3, &
        -0.61069018684944109624d+3,  0.42269668805777760407d+3/
    data t/ 0.16125990452916363814d-9, 0.25051878502858255354d-7, &
        0.27557319615147964774d-5, 0.19841269840928373686d-3, &
        0.83333333333334751799d-2, 0.16666666666666666446d+0/
    data estm/5.20583d1, 5.7607d0, 2.7782d0, 1.44303d1, 1.853004d2, &
        9.3715d0/
    data estf/4.18341d1, 7.1075d0, 6.4306d0, 4.25110d1, 1.35633d0, &
        8.45096d1, 2.0d1/

    ex = x
    enu = alpha
    ncalc = min ( nb, 0 ) - 2

    if ( 0 < nb .and. &
        ( zero <= enu .and. enu < one ) .and. &
        ( 1 <= ize .and. ize <= 2 ) .and. &
        ( ize /= 1 .or. ex <= xmax ) .and. &
        zero < ex )  then

        k = 0
        if ( enu < sqxmin ) then
            enu = zero
        end if

        if ( half < enu ) then
            k = 1
            enu = enu - one
        end if

        twonu = enu + enu
        iend = nb + k - 1
        c = enu * enu
        d3 = -c

        if ( ex <= one ) then
            !
            !  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA,
            !                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA.
            !
            d1 = zero
            d2 = p(1)
            t1 = one
            t2 = q(1)

            do i = 2, 7, 2
                d1 = c * d1 + p(i)
                d2 = c * d2 + p(i+1)
                t1 = c * t1 + q(i)
                t2 = c * t2 + q(i+1)
            end do

            d1 = enu * d1
            t1 = enu * t1
            f1 = log ( ex )
            f0 = a + enu * ( p(8) &
                - enu * ( d1 + d2 ) / ( t1 + t2 ) ) - f1
            q0 = exp ( -enu * ( a - enu * &
                ( p(8) + enu * ( d1 - d2 ) / ( t1 - t2 ) ) - f1 ) )
            f1 = enu * f0
            p0 = exp ( f1 )
            !
            !  Calculation of F0.
            !
            d1 = r(5)
            t1 = one
            do i = 1, 4
                d1 = c * d1 + r(i)
                t1 = c * t1 + s(i)
            end do

            if ( abs ( f1 ) <= half ) then
                f1 = f1 * f1
                d2 = zero
                do i = 1, 6
                    d2 = f1 * d2 + t(i)
                end do
                d2 = f0 + f0 * f1 * d2
            else
                d2 = sinh ( f1 ) / enu
            end if

            f0 = d2 - enu * d1 / ( t1 * p0 )
            !
            !  X <= 1.0E-10.
            !
            !  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X).
            !
            if ( ex <= tinyx ) then

                bk(1) = f0 + ex * f0

                if ( ize == 1 ) then
                    bk(1) = bk(1) - ex * bk(1)
                end if

                ratio = p0 / f0
                c = ex * xinf
                !
                !  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
                !  1/2 <= ALPHA.
                !
                if ( k /= 0 ) then

                    ncalc = -1

                    if ( c / ratio <= bk(1) ) then
                        return
                    end if

                    bk(1) = ratio * bk(1) / ex
                    twonu = twonu + two
                    ratio = twonu

                end if

                ncalc = 1

                if ( nb == 1 ) then
                    return
                end if
                !
                !  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1.
                !
                ncalc = -1
                do i = 2, nb
                    if ( c <= ratio ) then
                        return
                    end if
                    bk(i) = ratio / ex
                    twonu = twonu + two
                    ratio = twonu
                end do

                ncalc = 1
                j = ncalc + 1

                do i = j, nb
                    if ( xinf / bk(i) <= bk(ncalc) ) then
                        return
                    end if
                    bk(i) = bk(ncalc) * bk(i)
                    ncalc = i
                end do

                return
            !
            !  1.0E-10 < X <= 1.0.
            !
            else

                c = one
                x2by4 = ex * ex / four
                p0 = half * p0
                q0 = half * q0
                d1 = -one
                d2 = zero
                bk1 = zero
                bk2 = zero
                f1 = f0
                f2 = p0

100         continue

            d1 = d1 + two
            d2 = d2 + one
            d3 = d1 + d3
            c = x2by4 * c / d2
            f0 = ( d2 * f0 + p0 + q0 ) / d3
            p0 = p0 / ( d2 - enu )
            q0 = q0 / ( d2 + enu )
            t1 = c * f0
            t2 = c * ( p0 - d2 * f0 )
            bk1 = bk1 + t1
            bk2 = bk2 + t2

            if ( eps < abs ( t1 / ( f1 + bk1 ) ) .or. &
                eps < abs ( t2 / ( f2 + bk2 ) ) )  then
                go to 100
            end if

            bk1 = f1 + bk1
            bk2 = two * ( f2 + bk2 ) / ex

            if ( ize == 2 ) then
                d1 = exp ( ex )
                bk1 = bk1 * d1
                bk2 = bk2 * d1
            end if

            wminf = estf(1) * ex + estf(2)

        end if
    !
    !  1/EPS < X.
    !
    else if ( one < eps * ex ) then

        ncalc = nb
        bk1 = one / ( d * sqrt ( ex ) )
        do i = 1, nb
            bk(i) = bk1
        end do

        return

    else
        !
        !  1 < X.
        !
        twox = ex + ex
        blpha = zero
        ratio = zero

        if ( ex <= four ) then
            !
            !  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0.
            !
            d2 = aint ( estm(1) / ex + estm(2) )
            m = int ( d2 )
            d1 = d2 + d2
            d2 = d2 - half
            d2 = d2 * d2
            do i = 2, m
                d1 = d1 - two
                d2 = d2 - d1
                ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
            end do
            !
            !  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
            !  recurrence and K(ALPHA,X) from the Wronskian.
            !
            d2 = aint ( estm(3) * ex + estm(4) )
            m = int ( d2 )
            c = abs ( enu )
            d3 = c + c
            d1 = d3 - one
            f1 = xmin
            f0 = ( two * ( c + d2 ) / ex &
                + half * ex / ( c + d2 + one ) ) * xmin

            do i = 3, m
                d2 = d2 - one
                f2 = ( d3 + d2 + d2 ) * f0
                blpha = ( one + d1 / d2 ) * ( f2 + blpha )
                f2 = f2 / ex + f1
                f1 = f0
                f0 = f2
            end do

            f1 = ( d3 + two ) * f0 / ex + f1
            d1 = zero
            t1 = one
            do i = 1, 7
                d1 = c * d1 + p(i)
                t1 = c * t1 + q(i)
            end do

            p0 = exp ( c * ( a + c * ( p(8) &
                - c * d1 / t1 ) - log ( ex ) ) ) / ex
            f2 = ( c + half - ratio ) * f1 / ex
            bk1 = p0 + ( d3 * f0 - f2 + f0 + blpha ) &
                / ( f2 + f1 + f0 ) * p0

            if ( ize == 1 ) then
                bk1 = bk1 * exp ( - ex )
            end if

            wminf = estf(3) * ex + estf(4)

        else
            !
            !  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
            !  recurrence, for 4 < X.
            !
            dm = aint ( estm(5) / ex + estm(6) )
            m = int ( dm )
            d2 = dm - half
            d2 = d2 * d2
            d1 = dm + dm

            do i = 2, m
                dm = dm - one
                d1 = d1 - two
                d2 = d2 - d1
                ratio = ( d3 + d2 ) / ( twox + d1 - ratio )
                blpha = ( ratio + ratio * blpha ) / dm
            end do

            bk1 = one / ( ( d + d * blpha ) * sqrt ( ex ) )

            if ( ize == 1 ) then
                bk1 = bk1 * exp ( - ex )
            end if

            wminf = estf(5) * ( ex - abs ( ex - estf(7) ) ) + estf(6)

        end if
        !
        !  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
        !  K(ALPHA+1,X)/K(ALPHA,X).
        !
        bk2 = bk1 + bk1 * ( enu + half - ratio ) / ex

    end if
    !
    !  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
    !  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1.
    !
    ncalc = nb
    bk(1) = bk1

    if ( iend == 0 ) then
        return
    end if

    j = 2 - k

    if ( 0 < j ) then
        bk(j) = bk2
    end if

    if ( iend == 1 ) then
        return
    end if

    m = min ( int ( wminf - enu ), iend )

    do i = 2, m

        t1 = bk1
        bk1 = bk2
        twonu = twonu + two

        if ( ex < one ) then

            if ( ( xinf / twonu ) * ex <= bk1 ) then
                exit
            end if

        else

            if ( xinf / twonu <= bk1 / ex ) then
                exit
            end if

        end if

        bk2 = twonu / ex * bk1 + t1
        itemp = i
        j = j + 1

        if ( 0 < j ) then
            bk(j) = bk2
        end if

    end do

    m = itemp

    if ( m == iend ) then
        return
    end if

    ratio = bk2 / bk1
    mplus1 = m + 1
    ncalc = -1

    do i = mplus1, iend

        twonu = twonu + two
        ratio = twonu / ex + one / ratio
        j = j + 1

        if ( 1 < j ) then
            bk(j) = ratio
        else
            if ( xinf / ratio <= bk2 ) then
                return
            end if
            bk2 = ratio * bk2
        end if

    end do

    ncalc = max ( mplus1 - k, 1 )

    if ( ncalc == 1 ) then
        bk(1) = bk2
    end if

    if ( nb == 1 ) then
        return
    end if

    j = ncalc + 1

    do i = j, nb
        if ( xinf / bk(i) <= bk(ncalc) ) then
            return
        end if
        bk(i) = bk(ncalc) * bk(i)
        ncalc = i
    end do

end if

return
end
