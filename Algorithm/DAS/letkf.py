import numpy, time, multiprocessing, gc, os

def letkf(letkf_common,nx,ny,nbv,num_local_obs,eps,Mask,Obs_Grid,h,R,Model_State,E0_SysModel,E0_ObsModel,
          Corelation_Par,GridSize_Sys,GridSize_Obs,Bias_Forecast_Model_Option, Bias_Observation_Model_Option,
          msw_infl,parm_infl,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,
          Par_Uniform_STD,Par_Sens_Dim,nx_copy,Def_Localization,Normal_Score_Trans,State_Layer_Num,
          Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,
          minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol):
    
#    numpy.savez("letkf.npz",nx,ny,nbv,num_local_obs,eps,Mask,Obs_Grid,h,R,Model_State,E0_SysModel,E0_ObsModel,
#          Corelation_Par,GridSize_Sys,GridSize_Obs,Bias_Forecast_Model_Option,
#          msw_infl,parm_infl,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,
#          Par_Uniform_STD,Par_Sens_Dim,nx_copy,Def_Localization)
#    os.abort()
    # 4D-LETKF
    
    bias_a = []
    xa = []
    innovation = []
    increments = []
    localization_map = []
    
    #xf = numpy.zeros((nx,nbv),dtype=numpy.float32)
    #bias_f = numpy.zeros(nx,dtype=numpy.float32)
    #bias_a = numpy.zeros(nx,dtype=numpy.float32)
    #xf_obsmodel = numpy.zeros((nx,nbv),dtype=numpy.float32)
    #xa = numpy.zeros((nx,nbv),dtype=numpy.float32)
    #localization_map = numpy.zeros(nx,dtype=numpy.float32)
    #dxf = numpy.zeros((nx,nbv),dtype=numpy.float32)
    #xm = numpy.zeros(nx,dtype=numpy.float32)
    #xm_ens = numpy.zeros(nx,dtype=numpy.float32)
    #y = numpy.zeros(ny,dtype=numpy.float32)
    #d = numpy.zeros(ny,dtype=numpy.float32)
    #h4d = numpy.zeros((ny,nx),dtype=numpy.float32)
    #hxf = numpy.zeros((ny,nbv),dtype=numpy.float32)
    #hxfm = numpy.zeros(ny,dtype=numpy.float32)
    #hdxf = numpy.zeros((ny,nbv),dtype=numpy.float32)
#    rdiag_loc = numpy.zeros(ny,dtype=numpy.float32)
#    d_loc = numpy.zeros(ny,dtype=numpy.float32)
    #hdxf_loc = numpy.zeros((ny,nbv),dtype=numpy.float32)
    
    # Record the model grid coordinates and the observation coordinates    
    #X_Coordiates_Mask = numpy.zeros(nx,dtype=numpy.float32)
    #Y_Coordiates_Mask = numpy.zeros(nx,dtype=numpy.float32)
    #X_Coordiates_Obs = numpy.zeros(ny,dtype=numpy.float32)
    #Y_Coordiates_Obs = numpy.zeros(ny,dtype=numpy.float32)
    
#    X_Coordiates_Mask = Mask[:,0] - GridSize_Sys/2.0      # Model x Coordinates
#    Y_Coordiates_Mask = Mask[:,1] - GridSize_Sys/2.0     # Model y Coordinates
#        
#    X_Coordiates_Obs = Obs_Grid[:,0] - GridSize_Obs/2.0     # Obs x Coordinates
#    Y_Coordiates_Obs = Obs_Grid[:,1] - GridSize_Obs/2.0      # Obs y Coordinates                   
    
    X_Coordiates_Mask = Mask[:,0]      # Model x Coordinates
    Y_Coordiates_Mask = Mask[:,1]     # Model y Coordinates
        
    X_Coordiates_Obs = Obs_Grid[:,0]     # Obs x Coordinates
    Y_Coordiates_Obs = Obs_Grid[:,1]      # Obs y Coordinates 
    
    xf = E0_SysModel
    hxf_whole = E0_ObsModel
    #bias_f = E0_Bias_Forecast[:]
    
    # read obs
    y = Obs_Grid[:,2]    # Observation Vector
    #print y
    
#    # ensemble mean -> xm
#    #print numpy.shape(xf)
#    for i in range(nx):
#        xm[i] = numpy.mean(xf[i,:])
#        #print xm[i,j]
#    
#    # ensemble ptb -> dxf
#    for i in range(nbv):
#        dxf[:,i] = xf[:,i] - xm
#        #print dxf[:,i]
            
    # analysis step
    
    # hxf = H xf    
#    for j in range(nbv):
#        h4d = h
#        hxf[:,j] = numpy.dot(h4d[:,0], E0_ObsModel[0,j])
#        for i in range(1,nx):
#            hxf[:,j] = hxf[:,j] + numpy.dot(h4d[:,i], E0_ObsModel[i,j])
##            if i == 1820:
##                print E0_ObsModel
#    
#    # hxfm = mean(H xf)
#    for i in range(ny):  
#        hxfm[i] = numpy.mean(hxf[i,:])
            
#    # d = y - hxfm
#    d = y - hxfm
#    #numpy.savetxt("d.txt", d)
#    
#    # hdxf
#    for i in range(nbv):     
#        hdxf[:,i] = hxf[:,i] - hxfm
            
    # LETKF
    
    pos_sys = numpy.zeros((nx, 2),dtype=numpy.float32)
    pos_sys[:,0] = X_Coordiates_Mask
    pos_sys[:,1] = Y_Coordiates_Mask
    pos_obs = numpy.zeros((ny, 2),dtype=numpy.float32)
    pos_obs[:,0] = X_Coordiates_Obs
    pos_obs[:,1] = Y_Coordiates_Obs
    
    #Judge the correlation type
    if Corelation_Par[0,0] == 3:
        vario_type = 'Spherical'
    elif Corelation_Par[0,0] == 2:
        vario_type = 'Exponential'
    elif Corelation_Par[0,0] == 4:
        vario_type = 'Gaussian'
    elif Corelation_Par[0,0] == 7:
        vario_type = 'SteMat'
    elif Corelation_Par[0,0] == 6:
        vario_type = 'Matern'
    elif Corelation_Par[0,0] == 11:
        vario_type = 'Gauss'
    elif Corelation_Par[0,0] == 12:
        vario_type = 'Gaspari_Cohn'
    elif Corelation_Par[0,0] == 13:
        vario_type = 'Cosine'
    elif Corelation_Par[0,0] == 14:
        vario_type = 'Cosine_Squared'
    elif Corelation_Par[0,0] == 15:
        vario_type = 'Exp3'
    elif Corelation_Par[0,0] == 16:
        vario_type = 'Lewitt'
    elif Corelation_Par[0,0] == 17:
        vario_type = 'Cubic' 
    elif Corelation_Par[0,0] == 18:
        vario_type = 'Quadro'
    elif Corelation_Par[0,0] == 19:
        vario_type = 'Step'
    
    #print "*****************************************LETKF Input Array Dimension Check*******************************************************"
#    print numpy.shape(pos_sys)
#    print numpy.shape(pos_obs),numpy.shape(R)
#    print nx,ny
#    print numpy.shape(xf),numpy.shape(xm),numpy.shape(d),numpy.shape(dxf),numpy.shape(hdxf),numpy.shape(hdxf_loc)
    
    if Def_Print:
        print "num_local_obs,msw_infl,eps,Corelation_Par,vario_type,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,Par_Uniform_STD,nx_copy,Def_Localization,GridSize_Sys"
        print num_local_obs,msw_infl,eps,Corelation_Par,vario_type,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,Par_Uniform_STD,nx_copy,Def_Localization,GridSize_Sys
        
        if Def_Print >= 2:
            print "numpy.min(R),numpy.min(pos_sys),numpy.min(pos_obs),numpy.min(xf),numpy.min(hxf_whole),numpy.min(h),numpy.min(y),numpy.min(parm_infl)"
            print numpy.min(R),numpy.min(pos_sys),numpy.min(pos_obs),numpy.min(xf),numpy.min(hxf_whole),numpy.min(h),numpy.min(y),numpy.min(parm_infl)
            print "numpy.max(R),numpy.max(pos_sys),numpy.max(pos_obs),numpy.max(xf),numpy.max(hxf_whole),numpy.max(h),numpy.max(y),numpy.max(parm_infl)"
            print numpy.max(R),numpy.max(pos_sys),numpy.max(pos_obs),numpy.max(xf),numpy.max(hxf_whole),numpy.max(h),numpy.max(y),numpy.max(parm_infl)
    print ""
    
    xa,innovation,increments,localization_map = letkf_common.call_letkf_2d(num_local_obs,msw_infl,eps,R,Corelation_Par,vario_type,pos_sys,pos_obs,xf,hxf_whole,h,y,\
                                                                           parm_infl,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,\
                                                                           Parameter_Regularization,Par_Uniform_STD,nx_copy,Def_Localization,GridSize_Sys,Normal_Score_Trans, State_Layer_Num,\
                                                                           Bias_Forecast_Model_Option,Bias_Observation_Model_Option,Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,
                                                                           minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol)
#     if Bias_Forecast_Model_Option == 0:
#         #print num_local_obs,msw_infl,eps,R,Corelation_Par,vario_type,pos_sys,pos_obs,xf,hxf_whole,h,y,parm_infl,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,Par_Uniform_STD,nx_copy,Def_Localization,GridSize_Sys
#         #print letkf_common.call_letkf_2d.__doc__
#         xa,innovation,increments,localization_map = letkf_common.call_letkf_2d(num_local_obs,msw_infl,eps,R,Corelation_Par,vario_type,pos_sys,pos_obs,xf,hxf_whole,h,y,parm_infl,Alpha_Inflation,nthreads,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,Par_Uniform_STD,nx_copy,Def_Localization,GridSize_Sys,Normal_Score_Trans, State_Layer_Num)
#     else:
#         xa,increments,localization_map,bias_a = letkf_common.call_letkf_2d_dynamic_bias(num_local_obs,msw_infl,eps,R,Corelation_Par,vario_type,pos_sys,pos_obs,xf,hxf_whole,h,y,parm_infl,Alpha_Inflation,bias_f,Gamma_Bias,nthreads,Def_Print,[nx,ny,nbv])
    
    #============== Memory Collection
    #del msw_infl,eps,R,Corelation_Par,vario_type,pos_sys,pos_obs,xf,hxf_whole
    del X_Coordiates_Mask, Y_Coordiates_Mask, X_Coordiates_Obs, Y_Coordiates_Obs
    del pos_sys, pos_obs
    
    gc.collect()
    del gc.garbage[:]
     
    return xa,innovation,increments,localization_map,bias_a
