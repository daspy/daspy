import os, sys, time, datetime, random, math, gc, subprocess, glob, string, shutil, warnings, commands, fileinput
import platform, numpy

def ReBEL(gssm_das_octave, letkf, letkf_common, octave, ftype,gssm_name,Gssm_model_tag,nx,ny,nwindow,nbv,num_local_obs,eps,Mask,Obs_Grid,h,B,R,Model_State,E0_SysModel,E0_ObsModel,
          Corelation_Par,GridSize_Sys,GridSize_Obs,bf,alpha_bias, Bias_Forecast_Model_Option,Bias_Observation_Model_Option,
          msw_infl,parm_infl,alpha_state,nthreads, U1_In, U2_In,Def_Print,Parameter_Optimization_Flag,Parameter_Regularization,Par_Uniform_STD,
          Par_Sens_Dim,nx_copy,Def_Localization,Normal_Score_Trans, State_Layer_Num,Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,
          minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol):
    
    xa = []
    bias_a = []
    localization_map = []
    innovation = []
    increments = []
    
    if ftype == 'letkf':
        #------------------- Local Ensemble Transform Kalman Filtering (LETKF) ---------------
        xa,innovation,increments,localization_map,bias_a = letkf.letkf(letkf_common,nx,ny,nbv,num_local_obs,eps,Mask,Obs_Grid,h,R,Model_State,E0_SysModel[:,0:nbv],E0_ObsModel,Corelation_Par,GridSize_Sys,GridSize_Obs,
                                                                 Bias_Forecast_Model_Option,Bias_Observation_Model_Option,msw_infl,parm_infl,alpha_state,nthreads,Def_Print,Parameter_Optimization_Flag,
                                                                 Parameter_Regularization,Par_Uniform_STD,Par_Sens_Dim,nx_copy,Def_Localization,Normal_Score_Trans,State_Layer_Num,Bias_Model_Uniform_STD,
                                                                 Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol)
        
        #print 'Analysis is:', numpy.mean(xa[0,:]), 'Model Ensemble Mean is:', numpy.mean(E0_SysModel[0,:]), 'Observation is:', Obs_Grid[0,2]
        return xa,innovation,increments,localization_map,bias_a
    
    
    
