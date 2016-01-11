!$Id: main_offline.F90 1369 2013-04-24 16:38:17Z lnerger $
!BOP
!
! !ROUTINE: main --- Main program for example of PDAF offline implementation
!
! !INTERFACE:
PROGRAM MAIN_OFFLINE

    ! !DESCRIPTION:
    ! This is the main program for an example implementation of
    ! PDAF with domain-decomposition and offline configuration.
    !
    ! In the offline mode, we assume that the ensemble
    ! integrations are performed in a separate program (model)
    ! and the forecasted ensemble can be read from files. After
    ! initializing the ensemble information by reading model
    ! outputs, a single analysis step is performed. Subsequently,
    ! the analysis ensemble can be written to files that can be
    ! used to initialize another ensemble forecast.
    !
    ! Using PDAF for domain-decomposition, the offline
    ! mode can be used to perform assimilation with domain-
    ! decomposed models. If the models write results for each
    ! sub-domain, these can be read here using the same
    ! parallelization. Then, the filter analysis can be
    ! performed utilizing this parallelization. If the files
    ! contain the full model state, PDAF in offline mode
    ! can be used either on a single processor, or the
    ! fields can be distributed in this program to utilize
    ! the parallelization of the filters.
    !
    ! Parameters can be set in the code, or - preferably -
    ! by command line arguments that are parsed by the
    ! module PARSER. The format for this is
    ! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
    ! The handles are defined in the code before the calls
    ! to the routine PARSE.
    !
    ! !REVISION HISTORY:
    ! 2008-07 - Lars Nerger - Initial code
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_parallel, &     ! Parallelization
        ONLY: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
        init_parallel, finalize_parallel
    USE timer, &            ! Timing
        ONLY: timeit, time_tot
    USE mod_memcount, &     ! Counting allocated memory
        ONLY: memcount_ini, memcount_get
    USE mod_assimilation, & ! Model variables
        ONLY: STATE_DIM, OBS_DIM, ENSEMBLE_NUMBER, ncid, varid, Normal_Score_Trans, &
        XF_NC, XF_NC_Copy, HXF_NC, OBS_NC, XF_COORD_NC, OBS_COORD_NC, R_NC, H_NC, &
        FILE_NAME, STATE_DIM_NAME, OBS_DIM_NAME, ENSEMBLE_NUMBER_NAME, Normal_Score_Trans_NAME, &
        XF_NAME, HXF_NAME, H_NAME, OBS_NAME, XF_COORD_NAME, OBS_COORD_NAME, R_NAME, XA_NAME, XM_NAME, &
        State_DIM_Single_Layer, Parameter_DIM, Par_Sens_Dim, Par_Uniform_STD, Alpha_Inflation, &
        State_DIM_Single_Layer_NAME, Parameter_DIM_NAME, Par_Sens_Dim_NAME, Par_Uniform_STD_NAME, Alpha_Inflation_NAME, &
        Parameter_Optimization_Flag, Parameter_Optimization_Flag_NAME, &
        State_DIM_Single_Column, State_DIM_Single_Column_NAME, &
        X_Left, X_Right, Y_Lower, Y_Upper, &
        Bias_Model_Dim, Bias_Obs_Dim, Bias_Model_Dim_NAME, Bias_Obs_Dim_NAME,&
        Bias_Model_Uniform_STD,Bias_Model_Uniform_STD_NAME, &
        Bias_Obs_Uniform_STD, Bias_Obs_Uniform_STD_NAME, &
        Bias_Forecast_Model_Option, Bias_Forecast_Model_Option_NAME, &
        Bias_Observation_Model_Option, Bias_Observation_Model_Option_NAME, &
		Model_Inflation_Uniform_STD, Model_Inflation_Uniform_STD_NAME, &
		Correlation_Range, Correlation_Range_NAME, &
		lmbda_DIM, minimize_lbfgsb_n, minimize_lbfgsb_m, minimize_lbfgsb_iprint, &
		minimize_lbfgsb_factr, minimize_lbfgsb_pgtol, &
		lmbda_state, lmbda_parameter, lmbda_bias, hlmbda, minimize_lbfgsb_epsilon_in, &
		minimize_lbfgsb_n_NAME, minimize_lbfgsb_m_NAME, minimize_lbfgsb_iprint_NAME, &
		minimize_lbfgsb_factr_NAME, minimize_lbfgsb_pgtol_NAME, &
		minimize_lbfgsb_epsilon_in_NAME, lmbda_DIM_NAME, &
		GridSize_Sys, GridSize_Obs, GridSize_Sys_NAME, GridSize_Obs_NAME

    use netcdf

    IMPLICIT NONE
    !EOP


    ! Local variables
    INTEGER :: i, j                 ! Counter

    ! **********************
    ! *** Initialize MPI ***
    ! **********************

    CALL init_parallel() ! initializes MPI


    ! **********************************************************
    ! ***               PROGRAM CONTROL                      ***
    ! **********************************************************

    print *,"*** Reading file ", FILE_NAME, "! "
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
    ! the file.
    call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, STATE_DIM_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = STATE_DIM) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, OBS_DIM_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = OBS_DIM) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, ENSEMBLE_NUMBER_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = ENSEMBLE_NUMBER) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, Parameter_DIM_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = Parameter_DIM) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, Par_Sens_Dim_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = Par_Sens_Dim) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, Bias_Model_Dim_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = Bias_Model_Dim) )

    ! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, Bias_Obs_Dim_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = Bias_Obs_Dim) )
	
	! Get the varid of the Dimension ID, based on its name.
    call check( nf90_inq_dimid(ncid, lmbda_DIM_NAME, varid) )
    ! Get the varid of the Dimension, based on its name.
    call check( nf90_inquire_dimension(ncid, varid, len = lmbda_DIM) )
	
		
    print*,"STATE_DIM,OBS_DIM,ENSEMBLE_NUMBER,Parameter_DIM,Par_Sens_Dim"
	print*,STATE_DIM,OBS_DIM,ENSEMBLE_NUMBER,Parameter_DIM,Par_Sens_Dim
    print*,"Bias_Model_Dim,Bias_Obs_Dim",Bias_Model_Dim,Bias_Obs_Dim

    ALLOCATE(XF_NC(STATE_DIM, ENSEMBLE_NUMBER))
    ALLOCATE(XF_NC_Copy(STATE_DIM, ENSEMBLE_NUMBER))
    ALLOCATE(HXF_NC(STATE_DIM, ENSEMBLE_NUMBER))
    ALLOCATE(XF_COORD_NC(STATE_DIM, 2))
    ALLOCATE(H_NC(OBS_DIM, STATE_DIM))
    ALLOCATE(OBS_NC(STATE_DIM))
    ALLOCATE(R_NC(OBS_DIM,OBS_DIM))
    ALLOCATE(OBS_COORD_NC(OBS_DIM, 2))
    ALLOCATE(Par_Uniform_STD(Parameter_DIM))
    ALLOCATE(Bias_Model_Uniform_STD(Bias_Model_Dim))
    ALLOCATE(Bias_Obs_Uniform_STD(Bias_Obs_Dim))
	ALLOCATE(Model_Inflation_Uniform_STD(STATE_DIM))
	ALLOCATE(minimize_lbfgsb_epsilon_in(lmbda_DIM))
	
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, State_DIM_Single_Layer_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, State_DIM_Single_Layer) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, State_DIM_Single_Column_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, State_DIM_Single_Column) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Correlation_Range_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Correlation_Range) )
	
	! Get the varid of the data variable, based on its name.
	call check( nf90_inq_varid(ncid, GridSize_Sys_NAME, varid) )
	! Read the data.
	call check( nf90_get_var(ncid, varid, GridSize_Sys) )
	
	! Get the varid of the data variable, based on its name.
	call check( nf90_inq_varid(ncid, GridSize_Obs_NAME, varid) )
	! Read the data.
	call check( nf90_get_var(ncid, varid, GridSize_Obs) )
	
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Normal_Score_Trans_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Normal_Score_Trans) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Alpha_Inflation_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Alpha_Inflation) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Parameter_Optimization_Flag_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Parameter_Optimization_Flag) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Bias_Forecast_Model_Option_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Bias_Forecast_Model_Option) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Bias_Observation_Model_Option_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Bias_Observation_Model_Option) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_n_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_n) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_m_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_m) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_iprint_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_iprint) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_factr_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_factr) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_pgtol_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_pgtol) )
	
    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, XF_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, XF_NC, start = (/ 1, 1 /), count = (/ STATE_DIM, ENSEMBLE_NUMBER /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, HXF_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, HXF_NC, start = (/ 1, 1 /), count = (/ STATE_DIM, ENSEMBLE_NUMBER /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, XF_COORD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, XF_COORD_NC, start = (/ 1, 1 /), count = (/ STATE_DIM, 2 /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, H_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, H_NC, start = (/ 1, 1 /), count = (/ OBS_DIM, STATE_DIM  /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, OBS_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, OBS_NC, start = (/ 1 /), count = (/ STATE_DIM /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, R_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, R_NC, start = (/ 1, 1 /), count = (/ OBS_DIM, OBS_DIM /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, OBS_COORD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, OBS_COORD_NC, start = (/ 1, 1 /), count = (/ OBS_DIM, 2 /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Par_Uniform_STD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Par_Uniform_STD, start = (/ 1 /), count = (/ Parameter_DIM /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Bias_Model_Uniform_STD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Bias_Model_Uniform_STD, start = (/ 1 /), count = (/ Bias_Model_Dim /)) )

    ! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Bias_Obs_Uniform_STD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Bias_Obs_Uniform_STD, start = (/ 1 /), count = (/ Bias_Obs_Dim /)) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, Model_Inflation_Uniform_STD_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, Model_Inflation_Uniform_STD, start = (/ 1 /), count = (/ STATE_DIM /)) )
	
	! Get the varid of the data variable, based on its name.
    call check( nf90_inq_varid(ncid, minimize_lbfgsb_epsilon_in_NAME, varid) )
    ! Read the data.
    call check( nf90_get_var(ncid, varid, minimize_lbfgsb_epsilon_in, start = (/ 1 /), count = (/ lmbda_DIM /)) )
	
    !print*,"varid",varid

    print*,"Normal_Score_Trans",Normal_Score_Trans
    print*,"State_DIM_Single_Layer",State_DIM_Single_Layer
    print*,"Parameter_DIM",Parameter_DIM
    print*,"State_DIM_Single_Column",State_DIM_Single_Column
    print*,"Correlation_Range",Correlation_Range
    print*,"Alpha_Inflation",Alpha_Inflation
    print*,"Parameter_Optimization_Flag",Parameter_Optimization_Flag
    print*,"XF",XF_NC(:,1)
    print*,"HXF",HXF_NC(:,1)
    print*,"XF_COORD",XF_COORD_NC(:,1)
    print*,"H",H_NC(:,1)
    print*,"OBS",OBS_NC(1)
    print*,"R",R_NC(1,1)
    print*,"OBS_COORD",OBS_COORD_NC(:,1)
    print*,"Par_Uniform_STD",Par_Uniform_STD
    print*,"Bias_Model_Uniform_STD",Bias_Model_Uniform_STD
    print*,"Bias_Obs_Uniform_STD",Bias_Obs_Uniform_STD
	print*,"minimize_lbfgsb_n",minimize_lbfgsb_n
	print*,"minimize_lbfgsb_m",minimize_lbfgsb_m
	print*,"minimize_lbfgsb_iprint",minimize_lbfgsb_iprint
	print*,"minimize_lbfgsb_factr",minimize_lbfgsb_factr
	print*,"minimize_lbfgsb_pgtol",minimize_lbfgsb_pgtol
	
    !    print*,"XF",XF_NC(2346,:)
    !    print*,"HXF",HXF_NC(2346,:)
    !    print*,"XF_COORD",XF_COORD_NC(2346,:)
    !    print*,"H",H_NC(2346,:)
    !    print*,"OBS",OBS_NC(2346)
    !    print*,"R",R_NC(2346,2346)
    !    print*,"OBS_COORD",OBS_COORD_NC(2346,:)
    !    print*,""
    !    print*,"XF",XF_NC(2347,:)
    !    print*,"HXF",HXF_NC(2347,:)
    !    print*,"XF_COORD",XF_COORD_NC(2347,:)
    !    print*,"H",H_NC(2347,:)
    !    print*,"OBS",OBS_NC(2347)
    !    print*,"R",R_NC(2347,2347)
    !    print*,"OBS_COORD",OBS_COORD_NC(2347,:)
    print*,"minimize_lbfgsb_epsilon_in",minimize_lbfgsb_epsilon_in(:)
    ! Close the file. This frees up any internal netCDF resources
    ! associated with the file.
    call check( nf90_close(ncid) )

    print *,"*** SUCCESS reading file ", FILE_NAME, "! "
    print *,""

    XF_NC_Copy = XF_NC

    X_Left = minval(XF_COORD_NC(:,1))
    X_Right = maxval(XF_COORD_NC(:,1))
    Y_Lower = minval(XF_COORD_NC(:,2))
    Y_Upper = maxval(XF_COORD_NC(:,2))
    print *,"X_Left, X_Right, Y_Lower, Y_Upper",X_Left, X_Right, Y_Lower, Y_Upper
    print*,""

    ! ********************************
    ! ***      INITIALIZATION      ***
    ! ********************************

    ! *** Initial Screen output ***
    initscreen: IF (mype_world == 0) THEN

        WRITE (*, '(/8x, a/)') '+++++ PDAF tutorial - offline mode +++++'
        WRITE (*, '(16x, a)') 'Data assimilation with PDAF'

        IF (npes_world > 1) THEN
            WRITE (*, '(/21x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
        ELSE
            WRITE (*, '(/21x, a/)') 'Running on 1 PE'
        END IF
        WRITE (*, '(/)')
     
    END IF initscreen

    ! *** set number of timers ***
    CALL timeit(6, 'ini')

    ! *** set first timer ***
    CALL timeit(1, 'new')

    ! *** set number of memory counters ***
    CALL memcount_ini(3)

  
    ! *** Initialize MPI communicators for PDAF (model and filter) ***
    ! *** NOTE: It is always n_modeltasks=1 for offline mode       ***
    CALL init_parallel_pdaf(0, 1)

    ! *** Initialize model information ***
    ! *** This should only be information on the model dimension
    ! *** Generally, this could be joined with init_pdaf.
    CALL timeit(2, 'new')
    CALL initialize()
    CALL timeit(2, 'old')


    ! *******************************
    ! ***      ASSIMILATION       ***
    ! *******************************

    print*," *** Initialize PDAF ***"
    CALL timeit(4, 'new')
    CALL init_pdaf()
    CALL timeit(4, 'old')

    print*," *** Perform analysis ***"
    CALL timeit(3, 'new')
    IF (mype_world == 0) &
        WRITE (*, '(/2x, a)') 'PDAF test suite - offline mode: START ASSIMILATION'
    CALL assimilation_pdaf_offline()

    ! Syncronize at barrier for exit
    CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)
    WRITE (*,*) 'model PE exited: mype ', mype_world

    CALL timeit(3, 'old')


    ! ********************
    ! *** Finishing up ***
    ! ********************

    CALL timeit(1, 'old')

    ! *** Final screen output ***
    screen3: IF (mype_world == 0) THEN
        WRITE (*, '(/1x, a)') 'PDAF test suite: EXITED ASSIMILATION'

        ! *** Show allocated memory for the model ***
        WRITE (*, '(/18x, a)') 'Model - Memory overview'
        WRITE (*, '(10x, 45a)') ('-', i=1, 45)
        WRITE (*, '(21x, a, f10.3, a)') 'Allocated memory  (MB)'
        WRITE (*, '(14x, a, f10.5, a)') &
            'Model field:', memcount_get(1, 'M'), ' MB (persistent)'
        WRITE (*, '(12x, a, f10.5, a)') &
            'ensemble init:', memcount_get(2, 'M'), ' MB (temporary)'
        WRITE (*, '(13x, a, f10.5, a)') &
            'Pre-Poststep:', memcount_get(3, 'M'), ' MB (temporary)'

        ! Show allocated memory for PDAF
        CALL PDAF_print_info(2)

        ! *** Print timings onto screen ***

        ! Show timings for PDAF
        CALL PDAF_print_info(1)

        WRITE (*, '(/17x, a)') 'Offline - Timing information'
        WRITE (*, '(10x, 45a)') ('-', i=1, 45)
        ! Timing summary for assimilation
        WRITE (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
        WRITE (*, '(18x, a, F11.3, 1x, a)') 'initialize filter:', time_tot(4), 's'
        WRITE (*, '(23x, a, F11.3, 1x, a)') 'assimilation:', time_tot(3), 's'
        WRITE (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'

        WRITE (*, '(/1x, a)') 'PDAF test suite: END'
    END IF screen3

    ! *** deallocate timers ***
    CALL timeit(6, 'fin')

    ! *** Terminate MPI
    CALL finalize_parallel()

    DEALLOCATE(HXF_NC)
    DEALLOCATE(XF_COORD_NC)
    DEALLOCATE(Par_Uniform_STD)

contains
    subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

END PROGRAM MAIN_OFFLINE
