# -*- coding: utf-8 -*- 
'''
Copyright of DasPy:
Author - Xujun Han (Forschungszentrum Jülich, Germany)
x.han@fz-juelich.de, xujunhan@gmail.com

DasPy was funded by:
1. Forschungszentrum Jülich, Agrosphere (IBG 3), Jülich, Germany
2. Cold and Arid Regions Environmental and Engineering Research Institute, Chinese Academy of Sciences, Lanzhou, PR China
3. Centre for High-Performance Scientific Computing in Terrestrial Systems: HPSC TerrSys, Geoverbund ABC/J, Jülich, Germany

Please include the following references related to DasPy:
1. Han, X., Li, X., He, G., Kumbhar, P., Montzka, C., Kollet, S., Miyoshi, T., Rosolem, R., Zhang, Y., Vereecken, H., and Franssen, H. J. H.: 
DasPy 1.0 : the Open Source Multivariate Land Data Assimilation Framework in combination with the Community Land Model 4.5, Geosci. Model Dev. Discuss., 8, 7395-7444, 2015.
2. Han, X., Franssen, H. J. H., Rosolem, R., Jin, R., Li, X., and Vereecken, H.: 
Correction of systematic model forcing bias of CLM using assimilation of cosmic-ray Neutrons and land surface temperature: a study in the Heihe Catchment, China, Hydrology and Earth System Sciences, 19, 615-629, 2015a.
3. Han, X., Franssen, H. J. H., Montzka, C., and Vereecken, H.: 
Soil moisture and soil properties estimation in the Community Land Model with synthetic brightness temperature observations, Water Resour Res, 50, 6081-6105, 2014a.
4. Han, X., Franssen, H. J. H., Li, X., Zhang, Y. L., Montzka, C., and Vereecken, H.: 
Joint Assimilation of Surface Temperature and L-Band Microwave Brightness Temperature in Land Data Assimilation, Vadose Zone J, 12, 0, 2013.
'''
import os, sys, time, datetime, random, math, gc, subprocess, glob, signal, string, shutil, warnings, multiprocessing, socket, getpass, ctypes, platform, functools, copy
import numpy, scipy, scipy.stats, scipy.signal, netCDF4, scipy.ndimage
import pp,imp

sys.path.append('SysModel/CLM')
sys.path.append('Utilities')
sys.path.append('Utilities/Soil')
sys.path.append('Algorithm')
sys.path.append('Algorithm/GSIF')
sys.path.append('Algorithm/ReBEL')
sys.path.append('Algorithm/Noise')
sys.path.append('Algorithm/MultiScale')
sys.path.append('Algorithm/Geostatistics/CorrelationModel')
sys.path.append('Algorithm/Geostatistics/Scripts')
sys.path.append('ForcingData')

from Call_CLM_CESM import *
from ParFor import *
from Read_Soil_Texture import *

from DAS_Assim_Common import *
from DAS_Driver_Common import *

    
def CLM_Assim_Common(Block_Index, Model_Driver, Def_PP, Def_First_Run, Def_Print, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs, eps, msw_infl, parm_infl, Post_Inflation_Alpha, Def_ParFor, Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict,  Call_Gstat_Flag, Assim_Algorithm_Name, Model_State, E0_SysModel, E0_ObsModel, \
                     Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Proj_String, Z_Resolution, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Variable_List,
                     Grid_Resolution_CEA, Prop_Grid_Array_Sys, Prop_Grid_Array_H_Trans, Model_Variance, Write_DA_File_Flag, Mask, Mask_Index, Land_Mask_Data, Observation_Variance, SensorQuantity, SensorQuantity_Index,
                     Observation_NLats, Observation_NLons, Observation_Longitude, Observation_Latitude, Observation_Matrix, DAS_Depends_Path,  DasPy_Path, CLM_NA, NAvalue, Soil_Layer_Index_DA, Soil_Layer_Num, ParFlow_Layer_Num, omp_get_num_procs_ParFor, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, NSLOTS, DAS_Output_Path, Region_Name,
                     Variable_Assimilation_Flag, Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point, SensorType, SensorVariable, SensorResolution, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial,
                     Observation_Corelation_Par,  Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, Low_Ratio_Par, High_Ratio_Par, 
                     Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                     CLM_Ground_Temperature_Ensemble_Mat,CLM_Vegetation_Temperature_Ensemble_Mat,CLM_Soil_Moisture_Ensemble_Mat,CLM_Soil_Temperature_Ensemble_Mat, PF_PRESSURE_Ensemble_Mat, PF_SATURATION_Ensemble_Mat,
                     Prop_Grid_Array_Sys_parm_infl, CLM_Latent_Heat_parm_infl, CLM_Surface_Temperature_parm_infl, CLM_Ground_Temperature_parm_infl,CLM_Vegetation_Temperature_parm_infl, CLM_Soil_Moisture_parm_infl,CLM_Soil_Temperature_parm_infl, PF_SATURATION_parm_infl,
                     CLM_Ground_Temperature_Ensemble_Mat_Bias,CLM_Vegetation_Temperature_Ensemble_Mat_Bias,CLM_Soil_Moisture_Ensemble_Mat_Bias,CLM_Soil_Temperature_Ensemble_Mat_Bias, 
                     CLM_Surface_Temperature_parm_infl_Bias, CLM_Ground_Temperature_parm_infl_Bias,CLM_Vegetation_Temperature_parm_infl_Bias, CLM_Soil_Moisture_parm_infl_Bias,CLM_Soil_Temperature_parm_infl_Bias,
                     Prop_Grid_Array_Bias, Observation_Bias, Prop_Grid_Array_Sys_parm_infl_Bias, Observation_parm_infl_Bias, Def_CDF_Matching, Plot_Analysis, Parameter_Optimization_Flag,
                     Start_Month, maxpft, Feedback_Assim, Dim_Soil_Par, Soil_Par_Sens, Dim_Veg_Par, Veg_Par_Sens, Dim_PFT_Par, PFT_Par_Sens, Dim_Hard_Par, Hard_Par_Sens, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, \
                     Parameter_Soil_Space_Ensemble, Parameter_Soil_Space_parm_infl, Parameter_Veg_Space_Ensemble, Parameter_Veg_Space_parm_infl, Parameter_PFT_Space_Ensemble, Parameter_PFT_Space_parm_infl, Parameter_Hard_Space_Ensemble, Parameter_Hard_Space_parm_infl, Parameter_Min_Max, \
                     Soil_Layer_Thickness_Ratio, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD, Par_Veg_Uniform_STD, Par_PFT_Uniform_STD, Par_Hard_Uniform_STD, 
                     Saturation_SSat, Saturation_SRes, Saturation_N, Saturation_Alpha, DateString_Plot, *vartuple):
    
    pyper = imp.load_source("pyper",DasPy_Path+"Utilities/pyper.py")
    Call_ReBEL_Octave = imp.load_source("Call_ReBEL_Octave",DasPy_Path+"Algorithm/ReBEL/Call_ReBEL_Octave.py")
    gssm_das_octave = []
    letkf = imp.load_source("letk",DasPy_Path+"Algorithm/DAS/letkf.py")
    letkf_common = imp.load_dynamic("letkf_common",DasPy_Path+"Algorithm/DAS/letkf_common.so")
    
    octave = vartuple[0]
    
    Analysis_Grid = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Localization_Map_Mask = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)    
    Analysis_Grid_Array = numpy.zeros((Ensemble_Number, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Innovation_State = numpy.zeros_like(Analysis_Grid_Array,dtype=numpy.float32)
    Increments_State = numpy.zeros_like(Analysis_Grid_Array,dtype=numpy.float32)
        
    #--------------Find Available Observations
    Observation_Matrix_Temp = numpy.copy(Observation_Matrix)
    Obs_Index = numpy.where(Observation_Matrix_Temp.flatten() != NAvalue)
    Obs_Index_Dim = numpy.size(Obs_Index)
    if Def_Print:
        print "Obs_Index_Dim", Obs_Index_Dim
    del Observation_Matrix_Temp
    
    #--------------- Save the Observation Grid for DA -----------------  
    Obs_Grid = numpy.zeros((Obs_Index_Dim, 3),dtype=numpy.float32)
    Obs_Grid[:, 0] = Observation_Longitude.flatten()[Obs_Index]
    Obs_Grid[:, 1] = Observation_Latitude.flatten()[Obs_Index]
    Obs_Grid[:, 2] = Observation_Matrix.flatten()[Obs_Index]
    
    if Write_DA_File_Flag:
        numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/X_Diff.txt', Mask[:, 0] - Observation_Longitude.flatten())
        numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Y_Diff.txt', Mask[:, 1] - Observation_Latitude.flatten())
        numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Obs_Grid.txt', Obs_Grid)
    #------------------------- Prepare DA --------------------
    
    State_DIM_Single_Layer = numpy.size(numpy.where(Mask_Index == False))
    
    # Expand the dimension to include the deep layers
    Mask_Copy = numpy.copy(Mask)
    Mask_Index_Single_Layer = numpy.copy(Mask_Index)
    
    Def_ReBEL_Temp = Def_ReBEL
    Def_Localization_Temp = Def_Localization
    
    if Def_ReBEL_Temp == 1:
        
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            # For State Estimation
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                if SensorType == "InSitu":
                    for i in range(Soil_Layer_Num - 5 -1):
                        Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                        Mask = numpy.vstack((Mask,Mask_Copy))
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        for i in range(Soil_Layer_Num + 2):
                            Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                            Mask = numpy.vstack((Mask,Mask_Copy))
                else:
                    for i in range(Soil_Layer_Num - 5):
                        Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                        Mask = numpy.vstack((Mask,Mask_Copy))
                            
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        for i in range(Soil_Layer_Num + 2):
                            Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                            Mask = numpy.vstack((Mask,Mask_Copy))
                            
            elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
                for i in range(Soil_Layer_Num + 2):
                    Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                    Mask = numpy.vstack((Mask,Mask_Copy))
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    for i in range(Soil_Layer_Num - 5):
                        Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                        Mask = numpy.vstack((Mask,Mask_Copy))
                
        else:
            # For Parameter Estimation
            if Soil_Par_Sens_Dim >= 1:
                for Par_Index in range(Dim_Soil_Par):
                    if Soil_Par_Sens[Par_Index]:
                        Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                        Mask = numpy.vstack((Mask,Mask_Copy))
            if PFT_Par_Sens_Dim >= 1:
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        Mask_Index = numpy.append(Mask_Index,Mask_Index_Single_Layer)
                        Mask = numpy.vstack((Mask,Mask_Copy))
            
            
        if Def_Print:
            print "numpy.shape(Mask_Index),numpy.shape(Mask)",numpy.shape(Mask_Index),numpy.shape(Mask)
        
        nx = numpy.size(numpy.where(Mask_Index == False))   # Number of Model Grids
        ny = Obs_Index_Dim   # Number of the Observations
        
        if Def_Print:
            if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
                # For State Estimatio
                if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                    if SensorType == "InSitu":
                        if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                            print "The Number of Model Grid is:", nx / (Soil_Layer_Num-5+17)
                        else:
                            print "The Number of Model Grid is:", nx / (Soil_Layer_Num-5)
                    else:
                        if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                            print "The Number of Model Grid is:", nx / (Soil_Layer_Num-5+1+17)
                        else:
                            print "The Number of Model Grid is:", nx / (Soil_Layer_Num-5+1)
                                
                elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        print "The Number of Model Grid is:", nx / (Soil_Layer_Num+3+10)
                    else:
                        print "The Number of Model Grid is:", nx / (Soil_Layer_Num+3)
                else:
                    print "The Number of Model Grid is:", nx / (1)
            
            else:
                # For Parameter Estimation
                if Soil_Par_Sens_Dim >= 1:
                    print "The Number of Model Grid is:", nx / (Soil_Par_Sens_Dim+1)
                if PFT_Par_Sens_Dim >= 1:
                    print "The Number of Model Grid is:", nx / (PFT_Par_Sens_Dim+1)
                
            print "The Number of Observation Grid is:", ny
            
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            # For State Estimation
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    if SensorType == "InSitu":
                        if (nx / (Soil_Layer_Num-5+17)) < ny:
                            print "******************nx / (Soil_Layer_Num-5+17) < ny************************",nx / (Soil_Layer_Num-5+17),ny
                            os.abort()
                    else:
                        if (nx / (Soil_Layer_Num-5+1+17)) < ny:
                            print "******************nx / (Soil_Layer_Num-5+1+17) < ny************************",nx / (Soil_Layer_Num-5+1+17),ny
                            os.abort()
                else:
                    if SensorType == "InSitu":
                        if (nx / (Soil_Layer_Num-5)) < ny:
                            print "******************nx / (Soil_Layer_Num-5) < ny************************",nx / (Soil_Layer_Num-5),ny
                            os.abort()
                    else:
                        if (nx / (Soil_Layer_Num-5+1)) < ny:
                            print "******************nx / (Soil_Layer_Num-5+1) < ny************************",nx / (Soil_Layer_Num-5+1),ny
                            os.abort()
                            
            elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    if (nx / (Soil_Layer_Num+3+10)) < ny:
                        print "******************nx / (Soil_Layer_Num+3+10) < ny************************",nx / (Soil_Layer_Num+3+10),ny
                        os.abort()
                else:
                    if (nx / (Soil_Layer_Num+3)) < ny:
                        print "******************nx / (Soil_Layer_Num+3) < ny************************",nx / (Soil_Layer_Num+3),ny
                        os.abort()
            else:
                if (nx / (1)) < ny:
                    print "******************nx / (1) < ny************************",nx / (1),ny
                    os.abort()
                
        else:
            # For Parameter Estimation
            if Soil_Par_Sens_Dim >= 1:
                if (nx / (Soil_Par_Sens_Dim + 1)) < ny:
                    print "******************nx / (numpy.size(numpy.where(Soil_Par_Sens == True)) + 1)************************",nx / (Soil_Par_Sens_Dim + 1),ny
                    os.abort()
            if PFT_Par_Sens_Dim >= 1:
                if (nx / (PFT_Par_Sens_Dim + 1)) < ny:
                    print "******************nx / (numpy.size(numpy.where(PFT_Par_Sens == True)) + 1)************************",nx / (PFT_Par_Sens_Dim + 1),ny
                    os.abort()
        
        #H Operator where the Value of H is 1.0 at the Location of the Observed Grid
        
        Mask_False = numpy.where(Mask_Index == False)
        Mask_False_Single_Layer = numpy.where(Mask_Index_Single_Layer == False)
        if Def_Print:
            print "Mask_False",Mask_False[0]
            print "Mask_False_Single_Layer",Mask_False_Single_Layer 
        
        h = numpy.zeros((ny, nx),dtype=numpy.integer)
        #h[~Mask_Index,~Mask_Index] = 1.0
        
        if Def_Print >= 3:
            print "Obs_Index",Obs_Index[0]
        
        if Def_ParFor:
            #print "State_DIM_Single_Layer",State_DIM_Single_Layer,numpy.shape(h[:,0:State_DIM_Single_Layer])
            #----------------------------------******* Using ParFor to Accelerate H_Operator"
            ParFor_H_Operator(h[:,0:State_DIM_Single_Layer],ny,Mask_False_Single_Layer[0],Obs_Index[0],Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")],\
                              SensorVariable,SensorType,Soil_Layer_Index_DA,State_DIM_Single_Layer,\
                              Parameter_Optimization_Flag,Bias_Estimation_Option_Model[Variable_List.index(SensorVariable)], Bias_Estimation_Option_Obs[Variable_List.index(SensorVariable)],DAS_Depends_Path,omp_get_num_procs_ParFor)
        else:
            for ny_index in range(ny):
                if numpy.size(numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])) > 0:
                    #print numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])
                    if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
                        # For State Estimation
                        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                            if SensorType == "InSitu":
                                #print Mask_False_Single_Layer[0],Obs_Index[0][ny_index],numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])
                                H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0] + Soil_Layer_Index_DA*State_DIM_Single_Layer
                            else:
                                #print Mask_False_Single_Layer[0],Obs_Index[0][ny_index],numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])
                                H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                        else:
                            H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                    elif (not Parameter_Optimization_Flag) and ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1)):
                        # For Bias Estimation
                        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                            H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                        else:
                            H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                    else:
                        # For Bias Estimation
                        H_Col_Index = numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                    #print "numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]",numpy.where(Mask_False_Single_Layer[0] == Obs_Index[0][ny_index])[0][0]
                    h[ny_index,H_Col_Index] = 1
    #    
        #-------------------------------------=========== Do Data Assimilation ========================================================'
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            # For State Estimation
            if Write_DA_File_Flag:
                #numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/h.txt", h)    # h is too large to save
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_SysModel_State.txt", E0_SysModel)
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_ObsModel_State.txt", E0_ObsModel)
        else:
            if Write_DA_File_Flag:
                #numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/h.txt", h)    # h is too large to save
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_SysModel_Parameter.txt", E0_SysModel)
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_ObsModel_Parameter.txt", E0_ObsModel)
        
        # 4D-LETKF
        nwindow = 1  # time window for 4D-LETKF
        
        
        R = numpy.diagflat(Observation_Variance.flatten()[Obs_Index])
        
        GridSize_Sys = abs((MODEL_X_Right - MODEL_X_Left) / Col_Numbers)
        GridSize_Obs = abs((MODEL_X_Right - MODEL_X_Left) / Col_Numbers)
        #------------------------- Call Assimilation Algorithm --------------------"
        
        ftype = Assim_Algorithm_Name
                
        bf = []
                            
        alpha_bias = 0.5     # LETKF Bias Correction Turning Parameter
        B = []
        
        gssm_name = 'CLM'
        Gssm_model_tag = 'GSSM_CLM'
        U1 = [0]
        U2 = [0]
        
        #print numpy.shape(Mask[~Mask_Index,:])
        #print numpy.shape(E0_SysModel[:,:]),numpy.shape(E0_ObsModel[:,:])
        #print nx,ny,numpy.shape(Mask[~Mask_Index,:])
        
        print 'There are ', nx, ' Grids need to be processed! ','nthreads',omp_get_num_procs_ParFor,' Num_Local_Obs ',Num_Local_Obs,"for Block",Block_Index
        if Def_Print:
            print "numpy.shape(parm_infl)",numpy.shape(parm_infl)
        
        #numpy.savetxt("E0_SysModel.txt",E0_SysModel[:,0:Ensemble_Number])
        #numpy.savetxt("E0_ObsModel.txt",E0_ObsModel)
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            # For State Estimation
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                State_Layer_Num_Single_Column = 11 + ParFlow_Layer_Num
                if Feedback_Assim:
                    State_Layer_Num_Single_Column = 11 + 15
            elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
                State_Layer_Num_Single_Column = 18
                if Feedback_Assim:
                    State_Layer_Num_Single_Column = 18 + 10
            else:
                State_Layer_Num_Single_Column = 1
                
        elif Parameter_Optimization_Flag:
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                State_Layer_Num_Single_Column = 1
                if Feedback_Assim:
                    State_Layer_Num_Single_Column = 1
            elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
                State_Layer_Num_Single_Column = 1
                if Feedback_Assim:
                    State_Layer_Num_Single_Column = 1
            else:
                State_Layer_Num_Single_Column = 1
            
        if PFT_Par_Sens_Dim >= 1:
            Par_Uniform_STD = numpy.asarray(Par_PFT_Uniform_STD,numpy.float32)
            Par_Sens_Dim = PFT_Par_Sens_Dim
        elif Soil_Par_Sens_Dim >= 1:
            Par_Uniform_STD = numpy.asarray(Par_Soil_Uniform_STD,numpy.float32)
            Par_Sens_Dim = Soil_Par_Sens_Dim
        else:
            Par_Uniform_STD = numpy.zeros(1,numpy.float32)
            Par_Sens_Dim = 1
        
        if Def_Print:
            print "Par_Uniform_STD",Par_Uniform_STD,"Par_Sens_Dim",Par_Sens_Dim
        
        Bias_Forecast_Model_Option = Bias_Estimation_Option_Model[Variable_List.index(SensorVariable)]
        Bias_Observation_Model_Option = Bias_Estimation_Option_Obs[Variable_List.index(SensorVariable)]
        
        Bias_Model_Dim = State_DIM_Single_Layer        
        Bias_Obs_Dim = State_DIM_Single_Layer
        Bias_Model_Uniform_STD = numpy.zeros(Bias_Model_Dim)
        Bias_Obs_Uniform_STD = numpy.zeros(Bias_Obs_Dim)
        
        Model_Inflation_Uniform_STD = numpy.zeros(nx)
        Model_Inflation_Uniform_STD[:] = Model_State_Inflation_Range_STD[Variable_List.index(SensorVariable)]
        
        if PFT_Par_Sens_Dim >= 1 or Soil_Par_Sens_Dim >= 1:
            if not numpy.size(Par_Uniform_STD) >= 1:
                print "not numpy.size(Par_Uniform_STD) >= 1 !!!!!!!!!!!!!!!!!!!!!!"
                os.abort()
        
        ########################NS++++++++++++++++++++++
        if (numpy.abs(numpy.min(Obs_Grid[:, 2]) - numpy.mean(E0_ObsModel[:, 0:Ensemble_Number])) > 3*numpy.std(E0_ObsModel[:, 0:Ensemble_Number])) or \
            (numpy.abs(numpy.max(Obs_Grid[:, 2]) - numpy.mean(E0_ObsModel[:, 0:Ensemble_Number])) > 3*numpy.std(E0_ObsModel[:, 0:Ensemble_Number])):
            Normal_Score_Trans_Temp = 0
        else:
            Normal_Score_Trans_Temp = Normal_Score_Trans
            
        ############################ BoxCox
        minimize_lbfgsb_m = 10
        minimize_lbfgsb_iprint = -1
        minimize_lbfgsb_factr = 1e1
        minimize_lbfgsb_pgtol = 1.0e-5
        minimize_lbfgsb_epsilon_in = numpy.asarray([1e-08,1e-08])
        
        if PDAF_Assim_Framework:
            
            print "********************************************** Using PDAF to Accelerate Assimilation"
            if PDAF_Filter_Type == 2 or PDAF_Filter_Type == 4 or PDAF_Filter_Type == 6:
                type_forget = 0
            else:
                type_forget = 0
            
            # Assing the processors for MPI
            if ny < NSLOTS:
                NSLOTS_PDAF = ny
            else:
                NSLOTS_PDAF = NSLOTS
            
            if PDAF_Assim_Framework == 2:
                PDAF_Path = "mpiexec -n "+str(NSLOTS_PDAF)+" "+DasPy_Path+"Algorithm/PDAF/bin/offline_2D_parallel/PDAF_offline -filtertype "+str(PDAF_Filter_Type)+" -type_forget "+str(type_forget)+" -locweight 2 -local_range "+str(Observation_Corelation_Par[3, 0])
                #PDAF_Path = "mpiexec -n 1 "+DasPy_Path+"Algorithm/PDAF/bin/offline_2D_parallel/PDAF_offline -filtertype "+str(PDAF_Filter_Type)+" -type_forget "+str(type_forget)+" -locweight 2 -local_range "+str(Observation_Corelation_Par[3, 0])
            
            else:
                PDAF_Path = DasPy_Path+"Algorithm/PDAF/bin/offline_2D_serial/PDAF_offline -filtertype "+str(PDAF_Filter_Type)+" -type_forget "+str(type_forget)+" -locweight 2 -local_range "+str(Observation_Corelation_Par[3, 0])
                #PDAF_Path = DasPy_Path+"Algorithm/PDAF/bin/offline_2D_serial/PDAF_offline -filtertype "+str(PDAF_Filter_Type)+" -type_forget "+str(type_forget)+" -locweight 2 -local_range 500"
            print "-----PDAF_Path-----",PDAF_Path
            
            os.chdir(DasPy_Path)
            
            PDAF_Work_Path = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_"+str(Block_Index+1)+"/"
            NC_FileName_PDAF = PDAF_Work_Path+"NC_to_PDAF.nc"
            
            NC_File_Out_PDAF = netCDF4.Dataset(NC_FileName_PDAF, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_PDAF.createDimension('STATE_DIM', nx)
            NC_File_Out_PDAF.createDimension('Parameter_DIM', Par_Sens_Dim)
            NC_File_Out_PDAF.createDimension('Par_Sens_Dim', Par_Sens_Dim)
            NC_File_Out_PDAF.createDimension('OBS_DIM', ny)
            NC_File_Out_PDAF.createDimension('SOIL_LAYER_NUM', Soil_Layer_Num)
            NC_File_Out_PDAF.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
            NC_File_Out_PDAF.createDimension('ENSEMBLE_NUMBER', Ensemble_Number)
            NC_File_Out_PDAF.createDimension('COORD_DIM', 2)
            NC_File_Out_PDAF.createDimension('Scalar', 1)
            NC_File_Out_PDAF.createDimension('Bias_Model_Dim',Bias_Model_Dim)
            NC_File_Out_PDAF.createDimension('Bias_Obs_Dim',Bias_Obs_Dim) 
            NC_File_Out_PDAF.createDimension('lmbda_DIM',2)
            
            NC_File_Out_PDAF.createVariable('Normal_Score_Trans','i4',('Scalar',),zlib=True)
            if Normal_Score_Trans_Temp:
                NC_File_Out_PDAF.variables['Normal_Score_Trans'][:] = 1
            else:
                NC_File_Out_PDAF.variables['Normal_Score_Trans'][:] = 0
            NC_File_Out_PDAF.createVariable('Alpha_Inflation','f4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['Alpha_Inflation'][:] = Post_Inflation_Alpha
            NC_File_Out_PDAF.createVariable('Parameter_Optimization_Flag','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['Parameter_Optimization_Flag'][:] = Parameter_Optimization_Flag
            NC_File_Out_PDAF.createVariable('Bias_Forecast_Model_Option','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['Bias_Forecast_Model_Option'][:] = Bias_Forecast_Model_Option
            NC_File_Out_PDAF.createVariable('Bias_Observation_Model_Option','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['Bias_Observation_Model_Option'][:] = Bias_Observation_Model_Option
            NC_File_Out_PDAF.createVariable('State_DIM_Single_Layer','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['State_DIM_Single_Layer'][:] = State_DIM_Single_Layer
            NC_File_Out_PDAF.createVariable('State_DIM_Single_Column','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['State_DIM_Single_Column'][:] = State_Layer_Num_Single_Column
            NC_File_Out_PDAF.createVariable('Correlation_Range','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['Correlation_Range'][:] = Observation_Corelation_Par[3, 0] / Grid_Resolution_CEA
            NC_File_Out_PDAF.createVariable('GridSize_Sys','f4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['GridSize_Sys'][:] = GridSize_Sys
            NC_File_Out_PDAF.createVariable('GridSize_Obs','f4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['GridSize_Obs'][:] = GridSize_Obs
            
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_n','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_n'][:] = 2
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_m','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_m'][:] = minimize_lbfgsb_m
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_iprint','i4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_iprint'][:] = minimize_lbfgsb_iprint
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_factr','f4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_factr'][:] = minimize_lbfgsb_factr
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_pgtol','f4',('Scalar',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_pgtol'][:] = minimize_lbfgsb_pgtol
            NC_File_Out_PDAF.createVariable('minimize_lbfgsb_epsilon_in','f4',('lmbda_DIM',),zlib=True)
            NC_File_Out_PDAF.variables['minimize_lbfgsb_epsilon_in'][:] = minimize_lbfgsb_epsilon_in
                                    
            NC_File_Out_PDAF.createVariable('XF_NC','f4',('ENSEMBLE_NUMBER','STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('HXF_NC','f4',('ENSEMBLE_NUMBER','STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('H_NC','i4',('STATE_DIM','OBS_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('OBS_NC','f4',('STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('XF_COORD_NC','f4',('COORD_DIM','STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('OBS_COORD_NC','f4',('COORD_DIM','OBS_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('R_NC','f4',('OBS_DIM','OBS_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('XA_NC','f4',('ENSEMBLE_NUMBER','STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('XM_NC','f4',('STATE_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('Par_Uniform_STD','f4',('Parameter_DIM',),zlib=True)
            NC_File_Out_PDAF.createVariable('Bias_Model_Uniform_STD','f4',('Bias_Model_Dim',),zlib=True)
            NC_File_Out_PDAF.createVariable('Bias_Obs_Uniform_STD','f4',('Bias_Obs_Dim',),zlib=True)
            NC_File_Out_PDAF.createVariable('Model_Inflation_Uniform_STD','f4',('STATE_DIM',),zlib=True)
            
            NC_File_Out_PDAF.close()
            
            NC_File_Out_PDAF = netCDF4.Dataset(NC_FileName_PDAF, 'r+', format='NETCDF4')
            NC_File_Out_PDAF.variables['XF_NC'][:,:] = numpy.transpose(E0_SysModel[:, 0:Ensemble_Number])
            #NC_File_Out_PDAF.variables['XF_NC_Init'][:,:] = numpy.transpose(E0_SysModel_Copy_Dual[:, 0:Ensemble_Number])
            NC_File_Out_PDAF.variables['HXF_NC'][:,:] = numpy.transpose(E0_ObsModel[:, 0:Ensemble_Number])
            NC_File_Out_PDAF.variables['H_NC'][:,:] = 0.0
            NC_File_Out_PDAF.variables['H_NC'][:,:] = numpy.transpose(h[:,:])
            NC_File_Out_PDAF.variables['OBS_NC'][:] = -9999.0
            NC_File_Out_PDAF.variables['OBS_NC'][Obs_Index] = Obs_Grid[:, 2]
            NC_File_Out_PDAF.variables['XF_COORD_NC'][0,:] = Mask[~Mask_Index, 0]
            NC_File_Out_PDAF.variables['XF_COORD_NC'][1,:] = Mask[~Mask_Index, 1]                         
            NC_File_Out_PDAF.variables['OBS_COORD_NC'][0,:] = Mask[~Mask_Index, 0][Obs_Index]
            NC_File_Out_PDAF.variables['OBS_COORD_NC'][1,:] = Mask[~Mask_Index, 1][Obs_Index]
            #NC_File_Out_PDAF.variables['OBS_COORD_NC'][0,:] = Observation_Longitude.flatten()[~Mask_Index_Single_Layer]
            #NC_File_Out_PDAF.variables['OBS_COORD_NC'][1,:] = Observation_Latitude.flatten()[~Mask_Index_Single_Layer]
            
            NC_File_Out_PDAF.variables['R_NC'][:,:] = numpy.transpose(R)
            NC_File_Out_PDAF.variables['Par_Uniform_STD'][:] = Par_Uniform_STD
            NC_File_Out_PDAF.variables['Bias_Model_Uniform_STD'][:] = Bias_Model_Uniform_STD[:]
            NC_File_Out_PDAF.variables['Bias_Obs_Uniform_STD'][:] = Bias_Obs_Uniform_STD[:]
            NC_File_Out_PDAF.variables['Model_Inflation_Uniform_STD'][:] = Model_Inflation_Uniform_STD[:]
            
            NC_File_Out_PDAF.sync()
            NC_File_Out_PDAF.close()
            
            os.chdir(PDAF_Work_Path)
            print "************Call PDAF"
            PDAF_Output = open("PDAF_Output.txt","w")
            subprocess.call(shlex.split(PDAF_Path), stdout=PDAF_Output, stderr=PDAF_Output, shell=False)
            #subprocess.call(shlex.split(PDAF_Path, shell=False)
            PDAF_Output.close()
            #subprocess.call(shlex.split("killall -9 -q -w PDAF_offline &> /dev/null"),shell=False)
            #subprocess.call(shlex.split("killall -9 -q -w psilogger PDAF_offline &> /dev/null"),shell=False)
            #os.abort()
            os.chdir(DasPy_Path)
            
            NC_File_Out_PDAF = netCDF4.Dataset(NC_FileName_PDAF, 'r')
            xa = numpy.transpose(NC_File_Out_PDAF.variables['XA_NC'][:,:])
            NC_File_Out_PDAF.close()
            
            innovation = numpy.zeros((nx,Ensemble_Number),dtype=numpy.float32)
            increments = numpy.zeros((nx,Ensemble_Number),dtype=numpy.float32)
            localization_map = numpy.zeros(nx,dtype=numpy.float32)
            bias_a = numpy.zeros(nx,dtype=numpy.float32)
            
            #os.abort()
                                
        else:
            try:
                xa,innovation,increments,localization_map,bias_a = Call_ReBEL_Octave.ReBEL(gssm_das_octave, letkf, letkf_common, octave,ftype,gssm_name,Gssm_model_tag,nx,ny,nwindow,Ensemble_Number,Num_Local_Obs,eps,Mask[~Mask_Index,:],Obs_Grid,h,B,R,Model_State,E0_SysModel[:,0:Ensemble_Number],E0_ObsModel,
                                              Observation_Corelation_Par,Grid_Resolution_CEA,Grid_Resolution_CEA,bf,alpha_bias, Bias_Forecast_Model_Option, Bias_Observation_Model_Option, msw_infl,parm_infl,Post_Inflation_Alpha,omp_get_num_procs_ParFor,U1,U2,Def_Print,Parameter_Optimization_Flag,
                                              Parameter_Regularization,Par_Uniform_STD,Par_Sens_Dim,State_DIM_Single_Layer,Def_Localization_Temp,Normal_Score_Trans,State_Layer_Num_Single_Column,Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,
                                              minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol)
            except:
                print "**************** User Default Correlation Parameters to Call ReBEL Again!"
                Observation_Corelation_Par[0, 0] = 6 # matern Model
                Observation_Corelation_Par[1, 0] = 0.0
                Observation_Corelation_Par[2, 0] = 1.0
                Observation_Corelation_Par[3, 0] = 4.0*Grid_Resolution_CEA
                Observation_Corelation_Par[4, 0] = 1.0
                xa,innovation,increments,localization_map,bias_a = Call_ReBEL_Octave.ReBEL(gssm_das_octave, letkf, letkf_common, octave,ftype,gssm_name,Gssm_model_tag,nx,ny,nwindow,Ensemble_Number,Num_Local_Obs,eps,Mask[~Mask_Index,:],Obs_Grid,h,B,R,Model_State,E0_SysModel[:,0:Ensemble_Number],E0_ObsModel,
                                              Observation_Corelation_Par,Grid_Resolution_CEA,Grid_Resolution_CEA,bf,alpha_bias, Bias_Forecast_Model_Option, Bias_Observation_Model_Option, msw_infl,parm_infl,Post_Inflation_Alpha,omp_get_num_procs_ParFor,U1,U2,Def_Print,Parameter_Optimization_Flag,
                                              Parameter_Regularization,Par_Uniform_STD,Par_Sens_Dim,State_DIM_Single_Layer,Def_Localization_Temp,Normal_Score_Trans,State_Layer_Num_Single_Column,Bias_Model_Uniform_STD,Bias_Obs_Uniform_STD,Model_Inflation_Uniform_STD,
                                              minimize_lbfgsb_m,minimize_lbfgsb_iprint,minimize_lbfgsb_epsilon_in,minimize_lbfgsb_factr,minimize_lbfgsb_pgtol)
        
        if Def_Print:
            print "********************************** Mean Innovation ************************************************"
            print "Mean Innovation Value is:",numpy.mean(innovation),"Max Innovation Value is:",numpy.max(innovation),"Min Innovation Value is:",numpy.min(innovation)
            print "********************************** Mean Innovation ************************************************"
            print "********************************** Mean Increments ************************************************"
            print "Mean Increments Value is:",numpy.mean(increments),"Max Increments Value is:",numpy.max(increments),"Min Increments Value is:",numpy.min(increments)
            print "********************************** Mean Increments ************************************************"
        #print xa[0,:]                
        #xa,parm_infl = ReBEL(ftype,gssm_name,msw_infl,nx,ny,nwindow,Ensemble_Number,eps,Mask,Obs_Grid,h,R,E0_SysModel,E0_ObsModel,Observation_Corelation_Par,GridSize_Sys,GridSize_Obs,1)
        
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            # State Estimation
            # ensemble mean
            if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature") or \
            (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Sensible_Heat"):
                if Def_Print:
                    print "******************************************************** Update CLM_Soil_Temperature_Ensemble_Mat"
                xm = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                localization_map_col = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                for i in range(State_DIM_Single_Layer):
                    xm[i] = numpy.mean(xa[i, :])
                    localization_map_col[i] = localization_map[i]
                
                for Ens_Index in range(Ensemble_Number):
                    Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+0*State_DIM_Single_Layer):(State_DIM_Single_Layer+(0+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                    CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+1*State_DIM_Single_Layer):(State_DIM_Single_Layer+(1+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                    CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                
                    for Soil_Layer_Index in range(Soil_Layer_Num):
                        Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
    #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Index+2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Index+2+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
    #                    #print Analysis_Grid_Col[~Mask_Index]
                        
                        CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                    Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(0*State_DIM_Single_Layer):((0+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                    #print Analysis_Grid_Col[~Mask_Index]
                    #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                    Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                
                
                if Def_Print:
                    print "******************************************************** Update CLM_Soil_Temperature_parm_infl"
                Analysis_Grid_Col = Analysis_Grid.flatten()
                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(0*State_DIM_Single_Layer):((0+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                CLM_Surface_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+0*State_DIM_Single_Layer):(State_DIM_Single_Layer+(0+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                CLM_Vegetation_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+1*State_DIM_Single_Layer):(State_DIM_Single_Layer+(1+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                CLM_Ground_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    Analysis_Grid_Col = Analysis_Grid.flatten()
#                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Index+2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Index+2+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
#                    #print Analysis_Grid_Col[~Mask_Index]
                    CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_Ensemble_Mat"
                    for Ens_Index in range(Ensemble_Number):
                        for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                            Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Index+Soil_Layer_Num+2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Index+Soil_Layer_Num+2+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
        #                    #print Analysis_Grid_Col[~Mask_Index]
                            
                            #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                            CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_parm_infl"
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        Analysis_Grid_Col = Analysis_Grid[::].flatten()
    #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Index+Soil_Layer_Num+2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Index+Soil_Layer_Num+2+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
    #                    #print Analysis_Grid_Col[~Mask_Index]
                        #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(parm_infl[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                        #print numpy.shape(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:]),numpy.shape(numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1)))
                        CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                                    
            elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
                
                if SensorType == "InSitu":
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_Ensemble_Mat"
                
                    xm = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                    localization_map_col = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                    for i in range(State_DIM_Single_Layer):
                        xm[i] = numpy.mean(xa[Soil_Layer_Index_DA*State_DIM_Single_Layer+i, :])
                        localization_map_col[i] = localization_map[Soil_Layer_Index_DA*State_DIM_Single_Layer+i]
                                
                    for Ens_Index in range(Ensemble_Number):
                        for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                            Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
        #                    #print Analysis_Grid_Col[~Mask_Index]
                            
                            #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                            CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
            #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(Soil_Layer_Index_DA*State_DIM_Single_Layer):((Soil_Layer_Index_DA+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                        #print Analysis_Grid_Col[~Mask_Index]
                        #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                        Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_parm_infl"
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        Analysis_Grid_Col = Analysis_Grid[::].flatten()
    #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
    #                    #print Analysis_Grid_Col[~Mask_Index]
                        #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(parm_infl[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                        #print numpy.shape(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:]),numpy.shape(numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1)))
                        CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                    
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        if Def_Print:
                            print "******************************************************** Update CLM_Soil_Temperature_Ensemble_Mat"
                        for Ens_Index in range(Ensemble_Number):
                            Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 0)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 0 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                            CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                            CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                            for Soil_Layer_Index in range(Soil_Layer_Num):
                                Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
            #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 2 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
            #                    #print Analysis_Grid_Col[~Mask_Index]
                                
                                CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        if Def_Print:
                            print "******************************************************** Update CLM_Soil_Temperature_parm_infl"
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 0)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 0 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                        CLM_Vegetation_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                        CLM_Ground_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            Analysis_Grid_Col = Analysis_Grid.flatten()
        #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 5 + Soil_Layer_Index + 2 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
        #                    #print Analysis_Grid_Col[~Mask_Index]
                            CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                                        
                else:
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_Ensemble_Mat"
                
                    xm = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                    localization_map_col = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                    for i in range(State_DIM_Single_Layer):
                        xm[i] = numpy.mean(xa[i, :])
                        localization_map_col[i] = localization_map[i]
                                
                    for Ens_Index in range(Ensemble_Number):
                        for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                            Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[((Soil_Layer_Index+1)*State_DIM_Single_Layer):((Soil_Layer_Index+2)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
        #                    #print Analysis_Grid_Col[~Mask_Index]
                            
                            #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                            CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
            #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[0:State_DIM_Single_Layer,Ens_Index],dtype=numpy.float32)
                        #print Analysis_Grid_Col[~Mask_Index]
                        #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                        Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                               
                        
                    if Def_Print:
                        print "******************************************************** Update CLM_Soil_Moisture_parm_infl"
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        Analysis_Grid_Col = Analysis_Grid[::].flatten()
    #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[((Soil_Layer_Index+1)*State_DIM_Single_Layer):((Soil_Layer_Index+2)*State_DIM_Single_Layer)],dtype=numpy.float32)
    #                    #print Analysis_Grid_Col[~Mask_Index]
                        #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(parm_infl[(Soil_Layer_Index*State_DIM_Single_Layer):((Soil_Layer_Index+1)*State_DIM_Single_Layer),:],dtype=numpy.float32),axis=1) * Random_Factor_Normal[Ens_Index]
                        #print numpy.shape(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:]),numpy.shape(numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1)))
                        CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                    
                    
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        if Def_Print:
                            print "******************************************************** Update CLM_Soil_Temperature_Ensemble_Mat"
                        for Ens_Index in range(Ensemble_Number):
                            Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                            CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 2 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                            CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                            for Soil_Layer_Index in range(Soil_Layer_Num):
                                Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
            #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 3)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 3 + 1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
            #                    #print Analysis_Grid_Col[~Mask_Index]
                                
                                CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        if Def_Print:
                            print "******************************************************** Update CLM_Soil_Temperature_parm_infl"
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 1 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                        CLM_Vegetation_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 2)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + 2 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                        CLM_Ground_Temperature_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                        
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            Analysis_Grid_Col = Analysis_Grid.flatten()
        #                    #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(parm_infl[:])
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 3)*State_DIM_Single_Layer):(State_DIM_Single_Layer+(Soil_Layer_Num - 6 + Soil_Layer_Index + 3 + 1)*State_DIM_Single_Layer)],dtype=numpy.float32)
        #                    #print Analysis_Grid_Col[~Mask_Index]
                            CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
            
            else:
                if Def_Print:
                    print "******************************************************** Update Prop_Grid_Array_Sys"
            
                xm = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                localization_map_col = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
                for i in range(State_DIM_Single_Layer):
                    xm[i] = numpy.mean(xa[i, :])
                    localization_map_col[i] = localization_map[i]
                            
                for Ens_Index in range(Ensemble_Number):
                    
                    Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[0:State_DIM_Single_Layer,Ens_Index],dtype=numpy.float32)
                    #print Analysis_Grid_Col[~Mask_Index]
                    #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                    Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                if Def_Print:
                    print "******************************************************** Update Prop_Grid_Array_Sys_parm_infl"
                Analysis_Grid_Col = Analysis_Grid.flatten()
                Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl[(0*State_DIM_Single_Layer):((0+1)*State_DIM_Single_Layer)],dtype=numpy.float32)
                Prop_Grid_Array_Sys_parm_infl[:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
        
        else:
            # Parameter Estimation
            xm = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
            localization_map_col = numpy.zeros(State_DIM_Single_Layer,dtype=numpy.float32)
            for i in range(State_DIM_Single_Layer):
                xm[i] = numpy.mean(xa[i, :])
                localization_map_col[i] = localization_map[i]
            
            if Soil_Par_Sens_Dim >= 1:
                
                # Assign the Optimized Parameters to Model Input
                Par_Index_Sub = 0
                for Par_Index in range(Dim_Soil_Par-1):
                    if Soil_Par_Sens[Par_Index]:
                        for Ens_Index in range(Ensemble_Number):
                            xa_temp = xa[((Par_Index_Sub+1)*State_DIM_Single_Layer):((Par_Index_Sub+2)*State_DIM_Single_Layer),Ens_Index]
                            Analysis_Grid_Col = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:].flatten()
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa_temp,dtype=numpy.float32)
                            Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            #Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = imadjust.imadjust(Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:],numpy.min(Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]),numpy.max(Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]),Parameter_Min_Max[Par_Index,0],Parameter_Min_Max[Par_Index,1])
                        parm_infl_temp = parm_infl[((Par_Index_Sub+1)*State_DIM_Single_Layer):((Par_Index_Sub+2)*State_DIM_Single_Layer)]
                        Analysis_Grid_Col = Parameter_Soil_Space_parm_infl[Par_Index,:,:].flatten()
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl_temp,dtype=numpy.float32)
                        Parameter_Soil_Space_parm_infl[Par_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            
                        Par_Index_Sub += 1
                               
                if Soil_Par_Sens[Dim_Soil_Par-1]:
                    xa_temp = xa[((Par_Index_Sub+1)*State_DIM_Single_Layer):((Par_Index_Sub+2)*State_DIM_Single_Layer),Ens_Index]
                    Analysis_Grid_Col = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:].flatten()
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa_temp,dtype=numpy.float32)
                    #print "data_matrix[Min_Index,:]",data_matrix[Min_Index,:]
                    #print "Low_Tail, High_Tail",Low_Tail, High_Tail, data_matrix[Min_Index,-2]+2
                    numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                    Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:] =  numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                # Soil Boundary Check
                Par_Index_Sub = 0
                for Par_Index in range(Dim_Soil_Par):
                    if Soil_Par_Sens[Par_Index]:
                        
                        # Remove the outliers based on the ensemble median
                        Parameter_Soil_Space_Ensemble_Median = numpy.median(Parameter_Soil_Space_Ensemble[:,Par_Index,:,:],axis=0)
                        Parameter_Soil_Space_Ensemble_Max = Parameter_Soil_Space_Ensemble_Median + Par_Soil_Uniform_STD[Par_Index_Sub]/(numpy.sqrt(1/12.0)*2.0)
                        Parameter_Soil_Space_Ensemble_Min = 2.0 * Parameter_Soil_Space_Ensemble_Median - Parameter_Soil_Space_Ensemble_Max
                        
                        for Ens_Index in range(Ensemble_Number):
                            numexpr_a = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]
                            numexpr_b = Parameter_Soil_Space_Ensemble_Min
                            numexpr_c = numpy.where(numexpr_a < numexpr_b)
                            Lower_Index = numexpr_c
                            numexpr_a = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]
                            numexpr_b = Parameter_Soil_Space_Ensemble_Max
                            numexpr_c = numpy.where(numexpr_a > numexpr_b)
                            Upper_Index = numexpr_c
                            
                            if numpy.size(Lower_Index) > 1:                            
                                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)+str(Ens_Index)))
                                Lower_Boundary_Ens = numpy.random.uniform(1.0,High_Ratio_Par,size=numpy.shape(Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]))
                                #print numpy.shape(Lower_Index),numpy.shape(Lower_Boundary_Ens)
                                Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:][Lower_Index] = numpy.multiply(Lower_Boundary_Ens[Lower_Index],Parameter_Soil_Space_Ensemble_Min[Lower_Index])
                                del Lower_Boundary_Ens
                            if numpy.size(Upper_Index) > 1:                            
                                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)+str(Ens_Index)))
                                Upper_Boundary_Ens = numpy.random.uniform(Low_Ratio_Par,1.0,size=numpy.shape(Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]))
                                Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:][Upper_Index] = numpy.multiply(Upper_Boundary_Ens[Upper_Index],Parameter_Soil_Space_Ensemble_Max[Upper_Index])
                                del Upper_Boundary_Ens
                                
                        del Parameter_Soil_Space_Ensemble_Median,Parameter_Soil_Space_Ensemble_Max,Parameter_Soil_Space_Ensemble_Min
                        
                        # Boundary Check
                        numexpr_a = Parameter_Soil_Space_Ensemble[:,Par_Index,:,:]
                        numexpr_b = Parameter_Range_Soil[0,Par_Index]
                        numexpr_c = numpy.where(numexpr_a < numexpr_b)
                        Lower_Index = numexpr_c
                        numexpr_a = Parameter_Soil_Space_Ensemble[:,Par_Index,:,:]
                        numexpr_b = Parameter_Range_Soil[1,Par_Index]
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        Upper_Index = numexpr_c
                        
                        if numpy.size(Lower_Index) > 1:                            
                            numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)))
                            #Upper = Parameter_Range_Soil[0,Par_Index] + Par_Soil_Uniform_STD[Par_Index_Sub]
                            Upper = Parameter_Range_Soil[0,Par_Index] + Par_Soil_Uniform_STD[Par_Index_Sub] / numpy.sqrt(1.0/12.0)
                            Lower_Boundary_Ens = numpy.random.uniform(Parameter_Range_Soil[0,Par_Index],Upper,size=numpy.shape(Parameter_Soil_Space_Ensemble[:,Par_Index,:,:]))
                            #print numpy.shape(Lower_Index),numpy.shape(Lower_Boundary_Ens)
                            Parameter_Soil_Space_Ensemble[:,Par_Index,:,:][Lower_Index] = Lower_Boundary_Ens[Lower_Index]
                        
                        if numpy.size(Upper_Index) > 1:                            
                            numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)))
                            #Lower = Parameter_Range_Soil[1,Par_Index] - Par_Soil_Uniform_STD[Par_Index_Sub]
                            Lower = Parameter_Range_Soil[1,Par_Index] - Par_Soil_Uniform_STD[Par_Index_Sub] / numpy.sqrt(1.0/12.0)
                            Upper_Boundary_Ens = numpy.random.uniform(Lower,Parameter_Range_Soil[1,Par_Index],size=numpy.shape(Parameter_Soil_Space_Ensemble[:,Par_Index,:,:]))
                            Parameter_Soil_Space_Ensemble[:,Par_Index,:,:][Upper_Index] = Upper_Boundary_Ens[Upper_Index]
                    
                        Par_Index_Sub = Par_Index_Sub + 1
                    
                if Soil_Par_Sens[Dim_Soil_Par-1]:
                    Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:] = numpy.asarray(numpy.round(Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:]),dtype=numpy.integer)
                    numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                    Lower_Boundary_Ens = numpy.random.randint(1,3,size=numpy.shape(Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:]))
                    numexpr_a = Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:]
                    numexpr_b = 1
                    numexpr_c = numpy.where(numexpr_a < numexpr_b)
                    Lower_Index = numexpr_c
                    Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:][Lower_Index] = Lower_Boundary_Ens[Lower_Index]
                    numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                    Upper_Boundary_Ens = numpy.random.randint(1,3,size=numpy.shape(Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:]))
                    numexpr_a = Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:]
                    numexpr_b = 20
                    numexpr_c = numpy.where(numexpr_a > numexpr_b)
                    Upper_Index = numexpr_c
                    Parameter_Soil_Space_Ensemble[:,Dim_Soil_Par-1,:,:][Upper_Index] = Upper_Boundary_Ens[Upper_Index]
                
                    del Lower_Boundary_Ens,Upper_Boundary_Ens
                #for Soil_Layer_Index_Sub in range(Dim_Soil_Par):
                #    print "*****************************2"
                #    print numpy.min(Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,:,:]),numpy.max(Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,:,:])
                
                if Def_ParFor:
                    Parameter_Soil_Space_Ensemble = ParFor_Texture_Check(Dim_Soil_Par, Ensemble_Number, Row_Numbers, Col_Numbers, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Soil_Space_Ensemble, DAS_Depends_Path, omp_get_num_procs_ParFor)
                else:
                    for Ens_Index in range(Ensemble_Number):
                        # Soil Texture Boundary Check
                        for Row_Index in range(Row_Numbers):
                            for Col_Index in range(Col_Numbers):
                                # if Sand + Clay is greater than their sum
                                for Soil_Layer_Index_Sub in range(Soil_Texture_Layer_Opt_Num):
                                    Texture_Sum = (Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index] + Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index])
                                    if Texture_Sum > 98.0:
                                        Ratio = Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index] / (Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index]+Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index])
                                        Diff = Texture_Sum - 98.0
                                        Diff_Part1 = Ratio*Diff
                                        Diff_Part2 = Diff - Diff_Part1
                                        Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index] -= Diff_Part1
                                        Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index] -= Diff_Part2
                    
                                        del Texture_Sum,Ratio,Diff,Diff_Part1,Diff_Part2
                #for Soil_Layer_Index_Sub in range(Dim_Soil_Par):
                #    print "*****************************3"
                #    print numpy.min(Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,:,:]),numpy.max(Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,:,:])
                numexpr_a = Parameter_Soil_Space_Ensemble[:,2*Soil_Texture_Layer_Opt_Num:3*Soil_Texture_Layer_Opt_Num,:,:]
                numexpr_b = 130.0
                numexpr_c = numpy.where(numexpr_a > numexpr_b)
                Parameter_Soil_Space_Ensemble[:,2*Soil_Texture_Layer_Opt_Num:3*Soil_Texture_Layer_Opt_Num,:,:][numexpr_c] = 130.0
                
                
                 
                for Ens_Index in range(Ensemble_Number):
                    Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(0*State_DIM_Single_Layer):((0+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                    #print Analysis_Grid_Col[~Mask_Index]
                    #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                    Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
            
            
            if PFT_Par_Sens_Dim >= 1:
                
                print "Assign the Optimized PFT Parameters to Model Input\n"
                #print numpy.shape(Parameter_PFT_Space_Ensemble)
                Par_Index_Sub = 0
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        for Ens_Index in range(Ensemble_Number):
                            xa_temp = xa[((Par_Index_Sub+1)*State_DIM_Single_Layer):((Par_Index_Sub+2)*State_DIM_Single_Layer),Ens_Index]
                            Analysis_Grid_Col = Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:].flatten()
                            Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa_temp,dtype=numpy.float32)
                            Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            #Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:] = imadjust.imadjust(Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:],numpy.min(Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]),numpy.max(Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]),Parameter_Min_Max[Par_Index,0],Parameter_Min_Max[Par_Index,1])
                        parm_infl_temp = parm_infl[((Par_Index_Sub+1)*State_DIM_Single_Layer):((Par_Index_Sub+2)*State_DIM_Single_Layer)]
                        Analysis_Grid_Col = Parameter_PFT_Space_parm_infl[Par_Index,:,:].flatten()
                        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(parm_infl_temp,dtype=numpy.float32)
                        Parameter_PFT_Space_parm_infl[Par_Index,:,:] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                            
                        Par_Index_Sub += 1
                    
                # Hard Boundary Check
                Par_Index_Sub = 0
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        
                        # Remove the outliers based on the ensemble median
                        Parameter_PFT_Space_Ensemble_Median = numpy.median(Parameter_PFT_Space_Ensemble[:,Par_Index,:,:],axis=0)
                        Parameter_PFT_Space_Ensemble_Max = Parameter_PFT_Space_Ensemble_Median + Par_PFT_Uniform_STD[Par_Index_Sub]/(numpy.sqrt(1/12.0)*2.0)
                        Parameter_PFT_Space_Ensemble_Min = 2.0 * Parameter_PFT_Space_Ensemble_Median - Parameter_PFT_Space_Ensemble_Max
                        
                        for Ens_Index in range(Ensemble_Number):
                            numexpr_a = Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]
                            numexpr_b = Parameter_PFT_Space_Ensemble_Min
                            numexpr_c = numpy.where(numexpr_a < numexpr_b)
                            Lower_Index = numexpr_c
                            numexpr_a = Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]
                            numexpr_b = Parameter_PFT_Space_Ensemble_Max
                            numexpr_c = numpy.where(numexpr_a > numexpr_b)
                            Upper_Index = numexpr_c
                            
                            if numpy.size(Lower_Index) > 1:                            
                                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)+str(Ens_Index)))
                                Lower_Boundary_Ens = numpy.random.uniform(1.0,High_Ratio_Par,size=numpy.shape(Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]))
                                #print numpy.shape(Lower_Index),numpy.shape(Lower_Boundary_Ens)
                                Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:][Lower_Index] = numpy.multiply(Lower_Boundary_Ens[Lower_Index],Parameter_PFT_Space_Ensemble_Min[Lower_Index])
                                del Lower_Boundary_Ens
                            if numpy.size(Upper_Index) > 1:                            
                                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)+str(Ens_Index)))
                                Upper_Boundary_Ens = numpy.random.uniform(Low_Ratio_Par,1.0,size=numpy.shape(Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:]))
                                Parameter_PFT_Space_Ensemble[Ens_Index,Par_Index,:,:][Upper_Index] = numpy.multiply(Upper_Boundary_Ens[Upper_Index],Parameter_PFT_Space_Ensemble_Max[Upper_Index])
                                del Upper_Boundary_Ens
                                
                        del Parameter_PFT_Space_Ensemble_Median,Parameter_PFT_Space_Ensemble_Max,Parameter_PFT_Space_Ensemble_Min
                        
                        # Boundary Check
                        numexpr_a = Parameter_PFT_Space_Ensemble[:,Par_Index,:,:]
                        numexpr_b = Parameter_Range_PFT[0,Par_Index]
                        numexpr_c = numpy.where(numexpr_a < numexpr_b)
                        Lower_Index = numexpr_c
                        numexpr_a = Parameter_PFT_Space_Ensemble[:,Par_Index,:,:]
                        numexpr_b = Parameter_Range_PFT[1,Par_Index]
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        Upper_Index = numexpr_c
                        
                        if numpy.size(Lower_Index) > 1:                            
                            numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)))
                            #Upper = Parameter_Range_PFT[0,Par_Index] + Par_PFT_Uniform_STD[Par_Index_Sub]
                            Upper = Parameter_Range_PFT[0,Par_Index] + Par_PFT_Uniform_STD[Par_Index_Sub] / numpy.sqrt(1.0/12.0)
                            Lower_Boundary_Ens = numpy.random.uniform(Parameter_Range_PFT[0,Par_Index],Upper,size=numpy.shape(Parameter_PFT_Space_Ensemble[:,Par_Index,:,:]))
                            #print numpy.shape(Lower_Index),numpy.shape(Lower_Boundary_Ens)
                            Parameter_PFT_Space_Ensemble[:,Par_Index,:,:][Lower_Index] = Lower_Boundary_Ens[Lower_Index]
                                                    
                        if numpy.size(Upper_Index) > 1:                            
                            numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)))
                            #Lower = Parameter_Range_PFT[1,Par_Index] - Par_PFT_Uniform_STD[Par_Index_Sub]
                            Lower = Parameter_Range_PFT[1,Par_Index] - Par_PFT_Uniform_STD[Par_Index_Sub] / numpy.sqrt(1.0/12.0)
                            Upper_Boundary_Ens = numpy.random.uniform(Lower,Parameter_Range_PFT[1,Par_Index],size=numpy.shape(Parameter_PFT_Space_Ensemble[:,Par_Index,:,:]))
                            Parameter_PFT_Space_Ensemble[:,Par_Index,:,:][Upper_Index] = Upper_Boundary_Ens[Upper_Index]
                    
                        Par_Index_Sub = Par_Index_Sub + 1
                    
                for Ens_Index in range(Ensemble_Number):
                    Analysis_Grid_Col = Analysis_Grid_Array[Ens_Index,::].flatten()
        #            #print numpy.shape(Analysis_Grid_Col[~Mask_Index_Single_Layer]),numpy.shape(xa[:,Ens_Index])
                    Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xa[(0*State_DIM_Single_Layer):((0+1)*State_DIM_Single_Layer),Ens_Index],dtype=numpy.float32)
                    #print Analysis_Grid_Col[~Mask_Index]
                    #Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.mean(numpy.asarray(xa[:,:]),axis=1) * Random_Factor_Normal[Ens_Index]
                    Analysis_Grid_Array[Ens_Index,::][::] = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
                    
                    
            if Write_DA_File_Flag:
                numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Analysis_Parameter.txt', xm)
                numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Analysis_Ens_Parameter.txt', xa)
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Innovation_Parameter.txt",numpy.mean(innovation,axis=1))
                numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Increments_Parameter.txt",numpy.mean(increments,axis=1))
                
            
        if Write_DA_File_Flag:
            numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Analysis_State.txt', xm)
            numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Localization_map_col.txt', localization_map_col)
            numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Analysis_Ens_State.txt', xa)
            numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Bias_Analysis_State.txt', bias_a)
            numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Innovation_State.txt",numpy.mean(innovation,axis=1))
            numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Increments_State.txt",numpy.mean(increments,axis=1))
        
        for Ens_Index in range(Ensemble_Number):
            Innovation_State_Col = Innovation_State[Ens_Index,:,:].flatten()
            Innovation_State_Col[Obs_Index] = numpy.asarray(innovation[0:ny,Ens_Index],dtype=numpy.float32)
            Innovation_State[Ens_Index,:,:] = numpy.reshape(Innovation_State_Col, (Row_Numbers, -1))
            
            Increments_State_Col = Increments_State[Ens_Index,:,:].flatten()
            Increments_State_Col[~Mask_Index_Single_Layer] = numpy.asarray(increments[0:State_DIM_Single_Layer,Ens_Index],dtype=numpy.float32)
            Increments_State[Ens_Index,:,:] = numpy.reshape(Increments_State_Col, (Row_Numbers, -1))
            
            del Innovation_State_Col,Increments_State_Col
        #os.abort()
        #Analysis_Grid[:,:] = (CLM_Soil_Moisture[:,:,0]+Observation_Matrix[:,:])/2.0
        #Analysis_Grid[:, :] = Prop_Grid_Array_Sys[:,:]
        Analysis_Grid_Col = Analysis_Grid.flatten()
        Analysis_Grid_Col[~Mask_Index_Single_Layer] = numpy.asarray(xm,dtype=numpy.float32)
        Analysis_Grid = numpy.reshape(Analysis_Grid_Col, (Row_Numbers, -1))
        del Analysis_Grid_Col
        
        Localization_Map_Mask = numpy.zeros_like(Analysis_Grid)
        Localization_Map_Mask_Col = Localization_Map_Mask.flatten()
        Localization_Map_Mask_Col[~Mask_Index_Single_Layer] = localization_map_col
        Localization_Map_Mask = numpy.reshape(Localization_Map_Mask_Col, (Row_Numbers, -1))
        del Localization_Map_Mask_Col
        
        xm = []
        xa_temp = []
        localization_map_col = []
        Lower_Index = []
        Upper_Index = []
        parm_infl_temp = []
        Lower_Boundary_Ens = []
        Upper_Boundary_Ens = []
        
    numexpr_a = Land_Mask_Data
    numexpr_b = NAvalue
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    NA_Index_Analysis_Grid = numexpr_c
    
        
    if Def_Print >= 2:
        print "NA_Index_Analysis_Grid",NA_Index_Analysis_Grid
    
#    print numpy.min(Analysis_Grid),numpy.max(Analysis_Grid)
#    print numpy.min(Teta_Residual[Soil_Layer_Index_DA,::]),numpy.max(Teta_Saturated[Soil_Layer_Index_DA,::])
    
    if Write_DA_File_Flag:
        numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Analysis_Grid.txt', Analysis_Grid)
        numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Localization_Map_Mask.txt', Localization_Map_Mask)
    
    if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
        # State Estimation
        if Def_Print:
            print "Check the Ourliers"
        NA_Flag = False
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, Analysis_Grid, NA_Flag, SensorVariable,  Variable_Assimilation_Flag, Variable_List,
                            Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue)    
    #print numpy.min(Analysis_Grid)
        if Def_Print:
            print "Check the Ourliers"
        NA_Flag = False
        for Ens_Index in range(Ensemble_Number):
            Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, Analysis_Grid_Array[Ens_Index,::], NA_Flag, SensorVariable, Variable_Assimilation_Flag, Variable_List,
                                Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue)  
        
        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
            if Def_Print:
                print "Check the Ourliers"
            NA_Flag = False
            for Ens_Index in range(Ensemble_Number):
                for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                    Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index], NA_Flag, SensorVariable, Variable_Assimilation_Flag, Variable_List,
                                Teta_Residual[Soil_Layer_Index,:,:], Teta_Saturated[Soil_Layer_Index,:,:], Teta_Field_Capacity[Soil_Layer_Index,:,:], Teta_Wilting_Point[Soil_Layer_Index,:,:], NAvalue) 
            
                
            if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                for Ens_Index in range(Ensemble_Number):
                    #CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = Analysis_Grid_Array[Ens_Index,::]
                    Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index], NA_Flag, 'Soil_Temperature', Variable_Assimilation_Flag, Variable_List,
                                    Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue) 
                    Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index], NA_Flag, 'Vegetation_Temperature', Variable_Assimilation_Flag, Variable_List,
                                    Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue) 
                    for Soil_Layer_Index in range(Soil_Layer_Num):
                        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index], NA_Flag, 'Soil_Temperature', Variable_Assimilation_Flag, Variable_List,
                                    Teta_Residual[Soil_Layer_Index,:,:], Teta_Saturated[Soil_Layer_Index,:,:], Teta_Field_Capacity[Soil_Layer_Index,:,:], Teta_Wilting_Point[Soil_Layer_Index,:,:], NAvalue) 
        
           
        if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature"):
            if Def_Print:
                print "Check the Ourliers"
            NA_Flag = False
            for Ens_Index in range(Ensemble_Number):
                #CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = Analysis_Grid_Array[Ens_Index,::]
                Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index], NA_Flag, 'Soil_Temperature', Variable_Assimilation_Flag, Variable_List,
                                Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue) 
                Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index], NA_Flag, 'Vegetation_Temperature', Variable_Assimilation_Flag, Variable_List,
                                Teta_Residual[Soil_Layer_Index_DA,:,:], Teta_Saturated[Soil_Layer_Index_DA,:,:], Teta_Field_Capacity[Soil_Layer_Index_DA,:,:], Teta_Wilting_Point[Soil_Layer_Index_DA,:,:], NAvalue) 
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index], NA_Flag, 'Soil_Temperature', Variable_Assimilation_Flag, Variable_List,
                                Teta_Residual[Soil_Layer_Index,:,:], Teta_Saturated[Soil_Layer_Index,:,:], Teta_Field_Capacity[Soil_Layer_Index,:,:], Teta_Wilting_Point[Soil_Layer_Index,:,:], NAvalue) 
                
            if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                for Ens_Index in range(Ensemble_Number):
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,:,:,Ens_Index], NA_Flag, "Soil_Moisture", Variable_Assimilation_Flag, Variable_List,
                                    Teta_Residual[Soil_Layer_Index,:,:], Teta_Saturated[Soil_Layer_Index,:,:], Teta_Field_Capacity[Soil_Layer_Index,:,:], Teta_Wilting_Point[Soil_Layer_Index,:,:], NAvalue) 
    
                    
    if Def_ReBEL_Temp == 1:
        #============== Memory Collection
        del Mask, Mask_Index, Mask_Copy, Mask_Index_Single_Layer, Mask_False, Obs_Grid, h, R
        del Obs_Index, xa, innovation, increments, localization_map, bias_a
        del xm, localization_map_col, xa_temp, parm_infl_temp
    
    Analysis_Grid_Col = []
    Innovation_State_Col = []
    Increments_State_Col = []
    Localization_Map_Mask_Col = []
    del numexpr_a,numexpr_b,numexpr_c,NA_Index_Analysis_Grid
    del pyper, Call_ReBEL_Octave, gssm_das_octave, letkf, letkf_common
    del Bias_Model_Uniform_STD, Bias_Obs_Uniform_STD, Model_Inflation_Uniform_STD
                
    gc.collect()
    del gc.garbage[:]
    
    return Analysis_Grid, Analysis_Grid_Array, Localization_Map_Mask, \
        CLM_Ground_Temperature_Ensemble_Mat, CLM_Vegetation_Temperature_Ensemble_Mat, CLM_Soil_Moisture_Ensemble_Mat, CLM_Soil_Temperature_Ensemble_Mat, PF_PRESSURE_Ensemble_Mat, PF_SATURATION_Ensemble_Mat, \
        Prop_Grid_Array_Sys_parm_infl, CLM_Latent_Heat_parm_infl, CLM_Surface_Temperature_parm_infl, CLM_Ground_Temperature_parm_infl, CLM_Vegetation_Temperature_parm_infl, CLM_Soil_Moisture_parm_infl, CLM_Soil_Temperature_parm_infl, PF_SATURATION_parm_infl,\
        CLM_Ground_Temperature_Ensemble_Mat_Bias, CLM_Vegetation_Temperature_Ensemble_Mat_Bias, CLM_Soil_Moisture_Ensemble_Mat_Bias, CLM_Soil_Temperature_Ensemble_Mat_Bias, \
        CLM_Surface_Temperature_parm_infl_Bias, CLM_Ground_Temperature_parm_infl_Bias, CLM_Vegetation_Temperature_parm_infl_Bias, CLM_Soil_Moisture_parm_infl_Bias, CLM_Soil_Temperature_parm_infl_Bias,\
        Prop_Grid_Array_Bias, Observation_Bias, Prop_Grid_Array_Sys_parm_infl_Bias, Observation_parm_infl_Bias, \
        Parameter_Soil_Space_Ensemble, Parameter_Soil_Space_parm_infl, Parameter_Veg_Space_Ensemble, Parameter_Veg_Space_parm_infl, Parameter_PFT_Space_Ensemble, Parameter_PFT_Space_parm_infl, \
        Parameter_Hard_Space_Ensemble, Parameter_Hard_Space_parm_infl, Innovation_State, Increments_State
        
def Block_Assim(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, 
                Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs, eps, msw_infl, Post_Inflation_Alpha, Def_ParFor, Ensemble_Number_Predict,
                Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                Observation_Corelation_Par, Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, Low_Ratio_Par, High_Ratio_Par,
                Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, 
                Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD, Par_Veg_Uniform_STD, Par_PFT_Uniform_STD, Par_Hard_Uniform_STD, DateString_Plot,
                DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path, *vartuple):
    
    NC_FileName_Block_Assim_Common = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_Common.nc"
    
    if Def_Print:
        print 'Read NetCDF File:',NC_FileName_Block_Assim_Common
    
    NC_File_Block_Assim_Common = netCDF4.Dataset(NC_FileName_Block_Assim_Common, 'r')
    Mask_Sub = NC_File_Block_Assim_Common.variables['Mask_Sub'][:,:,:]
    Mask_Index_Sub_NC = NC_File_Block_Assim_Common.variables['Mask_Index_Sub'][:,:]
    Model_Variance = NC_File_Block_Assim_Common.variables['Model_Variance'][:,:,:]
    Observation_Variance = NC_File_Block_Assim_Common.variables['Observation_Variance'][:,:,:]    
    Observation_Latitude = NC_File_Block_Assim_Common.variables['Observation_Latitude'][:,:,:]
    Observation_Longitude = NC_File_Block_Assim_Common.variables['Observation_Longitude'][:,:,:]    
    Observation_Matrix = NC_File_Block_Assim_Common.variables['Observation_Matrix'][:,:,:]    
    NC_File_Block_Assim_Common.close()
    
    Mask_Index_Sub = numpy.zeros_like(Mask_Index_Sub_NC,dtype=numpy.bool)
    
    numexpr_a = Mask_Index_Sub_NC
    numexpr_b = 0.0
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    Mask_Index_Sub[numexpr_c] = False
    
    numexpr_a = Mask_Index_Sub_NC
    numexpr_b = 1.0
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    Mask_Index_Sub[numexpr_c] = True
    
    #print Mask_Index_Sub
    
    #print "Teta_Residual",Teta_Residual
    #print "Teta_Saturated",Teta_Saturated
    #os.abort()
    
    Sub_Block_Index_Row = Sub_Block_Index_Row_Mat_Vector[Block_Index]
    Sub_Block_Index_Col = Sub_Block_Index_Col_Mat_Vector[Block_Index]
    
    # Run Data Assimilation
    Row_Numbers_SubBlock = Row_Numbers_SubBlock_Array[Block_Index]
    Col_Numbers_SubBlock = Col_Numbers_SubBlock_Array[Block_Index]
    
    # We define two boundary box for sub block assimilation: use big box to do assimilation, use small box to select the assimilation results to keep smooth
    # This is the box to select the assimilation results (small)
    Sub_Block_Row_Start = Sub_Block_Row_Start_Array[Block_Index]
    Sub_Block_Row_End = Sub_Block_Row_End_Array[Block_Index]
    Sub_Block_Col_Start = Sub_Block_Col_Start_Array[Block_Index]
    Sub_Block_Col_End = Sub_Block_Col_End_Array[Block_Index]
    
    if Def_Print:
        print "Sub_Block_Row_Start, Sub_Block_Row_End, Sub_Block_Col_Start, Sub_Block_Col_End",Sub_Block_Row_Start, Sub_Block_Row_End, Sub_Block_Col_Start, Sub_Block_Col_End
    # This is the box to select the states for assimilation (big)
    
    Observation_Box_Row_Index_Start = max([Sub_Block_Row_Start - Observation_Box, 0])
    Observation_Box_Row_Index_End = min([Sub_Block_Row_End + Observation_Box, Row_Numbers])
    Observation_Box_Col_Index_Start = max([Sub_Block_Col_Start - Observation_Box, 0])
    Observation_Box_Col_Index_End = min([Sub_Block_Col_End + Observation_Box, Col_Numbers])
    Observation_NLats_SubBlock = int(Observation_Box_Row_Index_End - Observation_Box_Row_Index_Start)
    Observation_NLons_SubBlock = int(Observation_Box_Col_Index_End - Observation_Box_Col_Index_Start)
    if Def_Print:
        print "Observation_Box_Row_Index_Start,Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start,Observation_Box_Col_Index_End",Observation_Box_Row_Index_Start,Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start,Observation_Box_Col_Index_End
    
    X_Left_SubBlock = MODEL_X_Left + Sub_Block_Index_Col * Observation_NLons_SubBlock * Grid_Resolution_CEA
    X_Right_SubBlock = MODEL_X_Right - (Sub_Block_Ratio_Col - Sub_Block_Index_Col - 1) * Observation_NLons_SubBlock * Grid_Resolution_CEA
    Y_Lower_SubBlock = MODEL_Y_Lower + (Sub_Block_Ratio_Row - Sub_Block_Index_Row - 1) * Observation_NLats_SubBlock * Grid_Resolution_CEA
    Y_Upper_SubBlock = MODEL_Y_Upper - Sub_Block_Index_Row * Observation_NLats_SubBlock * Grid_Resolution_CEA
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
    #print numpy.shape(Prop_Grid_Array_Sys[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End])
    Prop_Grid_Array_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
    Prop_Grid_Array_H_Trans_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, Prop_Grid_Array_Sys_Index, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
    
    Model_State_SubBlock = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0)[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
    Model_Variance_SubBlock = Model_Variance[Prop_Grid_Array_Sys_Index, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
    Mask_SubBlock = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, 3),dtype=numpy.float32)
    Mask_SubBlock[:, 0] = Mask_Sub[0, ::][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End].flatten()
    Mask_SubBlock[:, 1] = Mask_Sub[1, ::][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End].flatten()
    Mask_SubBlock[:, 2] = Mask_Sub[2, ::][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End].flatten()
    Mask_Index_SubBlock = Mask_Index_Sub[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End].flatten()
    Land_Mask_Data_SubBlock = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
    if Def_Print >= 2:
        print "numpy.size(numpy.where(Mask_Index_SubBlock == False))",numpy.size(numpy.where(Mask_Index_SubBlock == False))
    if SensorVariable_Sub == "Soil_Moisture":
        Soil_Layer_Thickness_Ratio = numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Moisture'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End])
    else:
        Soil_Layer_Thickness_Ratio = numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End])
    
    NC_File_Out_Assimilation_2_Diagnostic.close()
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Constant.close()
    
    Analysis_Grid_SubBlock = numpy.zeros((Observation_NLats_SubBlock, Observation_NLons_SubBlock),dtype=numpy.float32)
    Analysis_Grid_Array_SubBlock = numpy.zeros((Ensemble_Number, Observation_NLats_SubBlock, Observation_NLons_SubBlock),dtype=numpy.float32)
    Localization_Map_Mask_SubBlock = numpy.zeros((Observation_NLats_SubBlock, Observation_NLons_SubBlock),dtype=numpy.float32)
    
    Innovation_State_SubBlock = numpy.zeros((Ensemble_Number, Observation_NLats_SubBlock, Observation_NLons_SubBlock),dtype=numpy.float32)
    Increments_State_SubBlock = numpy.zeros((Ensemble_Number, Observation_NLats_SubBlock, Observation_NLons_SubBlock),dtype=numpy.float32)
    
    Mask_SubBlock_Row = numpy.shape(Mask_SubBlock)[0]
    
    if Mask_SubBlock_Row == 0:
        # Compose the Full Analysis Grid
        Analysis_Grid_SubBlock = numpy.mean(Prop_Grid_Array_SubBlock, axis=0)
        
        for Ens_Index in range(Ensemble_Number):
            Analysis_Grid_Array_SubBlock[Ens_Index,:,:] = Prop_Grid_Array_SubBlock[Ens_Index, :, :]
    else:
        Observation_Variance_SubBlock = Observation_Variance[Observation_Matrix_Index,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Observation_Longitude_SubBlock = Observation_Longitude[Observation_Matrix_Index,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Observation_Latitude_SubBlock = Observation_Latitude[Observation_Matrix_Index,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Observation_Matrix_SubBlock = Observation_Matrix[Observation_Matrix_Index,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        numexpr_a = Mask_Index_Sub[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        numexpr_b = True
        numexpr_c = numpy.where(numexpr_a == numexpr_b)
        Observation_Matrix_SubBlock[numexpr_c] = NAvalue
#                            Observation_Variance_SubBlock = Observation_Variance[Observation_Matrix_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End]
#                            Observation_Longitude_SubBlock = Observation_Longitude[Observation_Matrix_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End]
#                            Observation_Latitude_SubBlock = Observation_Latitude[Observation_Matrix_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End]
#                            Observation_Matrix_SubBlock = Observation_Matrix[Observation_Matrix_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End]
        if Plot_Analysis >= 2:
            import matplotlib
            # Force matplotlib to not use any Xwindows backend.
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import matplotlib.cm as cm
            import matplotlib.colors as colors
            from mpl_toolkits.axes_grid.inset_locator import inset_axes
            
            Observation_Matrix_SubBlock_Masked = numpy.ma.masked_values(Observation_Matrix_SubBlock, NAvalue)
            
            w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)            
            
            Variable_Min = numpy.min(Observation_Matrix_SubBlock_Masked)
            Variable_Max = numpy.max(Observation_Matrix_SubBlock_Masked)
            ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
            color_boun_list = []
            color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
            for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                color_bound[0] += color_bound[2]
                color_boun_list.append(color_bound[0])
            
            
            fig1 = plt.figure(figsize=(w*2, h*2), dpi=80)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax = fig1.add_subplot(1, 1, 1)
            im1 = ax.imshow(Observation_Matrix_SubBlock_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('Observation_Matrix_SubBlock')
            plt.grid(True)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Observation_Matrix_SubBlock.png")
            plt.close('all')
        
        NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
        Teta_Residual_SubBlock = NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Teta_Saturated_SubBlock = NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Teta_Field_Capacity_SubBlock = NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        Teta_Wilting_Point_SubBlock = NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        
        CLM_Ground_Temperature_Ensemble_Mat_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        CLM_Soil_Moisture_Ensemble_Mat_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        CLM_Soil_Temperature_Ensemble_Mat_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        
        NC_File_Out_Assimilation_2_Diagnostic.close()
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Constant.close()
        
        Prop_Grid_Array_Bias_SubBlock = []
        CLM_Ground_Temperature_Ensemble_Mat_Bias_SubBlock = []
        CLM_Vegetation_Temperature_Ensemble_Mat_Bias_SubBlock = []
        CLM_Soil_Moisture_Ensemble_Mat_Bias_SubBlock = []
        CLM_Soil_Temperature_Ensemble_Mat_Bias_SubBlock = []
        Observation_Bias_SubBlock = []
        
        NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
        if Soil_Par_Sens_Dim >= 1:
            Parameter_Soil_Space_Ensemble_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:, :, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
            Parameter_Soil_Space_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        else:
            Parameter_Soil_Space_Ensemble_SubBlock = []
            Parameter_Soil_Space_parm_infl_SubBlock = []
        Parameter_Veg_Space_Ensemble_SubBlock = []
        Parameter_Veg_Space_parm_infl_SubBlock = []
        if PFT_Par_Sens_Dim >= 1:
            Parameter_PFT_Space_Ensemble_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
            Parameter_PFT_Space_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        else:
            Parameter_PFT_Space_Ensemble_SubBlock = []
            Parameter_PFT_Space_parm_infl_SubBlock = []
        Parameter_Hard_Space_Ensemble_SubBlock = []
        Parameter_Hard_Space_parm_infl_SubBlock = []
        
        Saturation_SSat_SubBlock = []
        Saturation_SRes_SubBlock = []
        Saturation_N_SubBlock = []
        Saturation_Alpha_SubBlock = []
        
        NC_File_Out_Assimilation_2_Parameter.close()
        
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        Prop_Grid_Array_Sys_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys_parm_infl'][Prop_Grid_Array_Sys_Index,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        CLM_Surface_Temperature_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Surface_Temperature_parm_infl'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        CLM_Ground_Temperature_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_parm_infl'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        CLM_Vegetation_Temperature_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_parm_infl'][Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        CLM_Soil_Moisture_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_parm_infl'][:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        CLM_Soil_Temperature_parm_infl_SubBlock = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_parm_infl'][:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
           
        NC_File_Out_Assimilation_2_Initial.close()
        
        Prop_Grid_Array_Sys_parm_infl_Bias_SubBlock = []
        CLM_Surface_Temperature_parm_infl_Bias_SubBlock = []
        CLM_Ground_Temperature_parm_infl_Bias_SubBlock = []
        CLM_Vegetation_Temperature_parm_infl_Bias_SubBlock = []
        CLM_Soil_Moisture_parm_infl_Bias_SubBlock = []
        CLM_Soil_Temperature_parm_infl_Bias_SubBlock = []
        Observation_parm_infl_Bias_SubBlock = []
        
        if Def_Print:
            print "---------------------- Prepare the Model Ensemble Grid for Bayesian Filtering DA -------------------"
        E0_SysModel_SubBlock = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32) 
        E0_ObsModel_SubBlock = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
        for Ens_Index in range(Ensemble_Number):
            E0_SysModel_SubBlock[:, Ens_Index] = Prop_Grid_Array_SubBlock[Ens_Index, :, :].flatten()
            #print Prop_Grid_Array_H_Trans[:, Prop_Grid_Array_Sys_Index, Row_Index,Col_Index]
            E0_ObsModel_SubBlock[:, Ens_Index] = Prop_Grid_Array_H_Trans_SubBlock[Ens_Index, :, :].flatten()
        
        ########################################################################################################################################
        Mask_Index_Vector = ~Mask_Index_SubBlock
        E0_SysModel_SubBlock = E0_SysModel_SubBlock[Mask_Index_Vector, 0:Ensemble_Number]
        E0_ObsModel_SubBlock = E0_ObsModel_SubBlock[Mask_Index_Vector, 0:Ensemble_Number]
        #print E0_ObsMode
        
        Parameter_Min_Max = numpy.zeros((Dim_Soil_Par,2))
        
        parm_infl = numpy.zeros(Observation_NLats_SubBlock * Observation_NLons_SubBlock,dtype=numpy.float32)
        
        Prop_Grid_Array_Sys_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Soil_Moisture_parm_infl = numpy.zeros((Soil_Layer_Num, Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Soil_Temperature_parm_infl = numpy.zeros((Soil_Layer_Num, Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Surface_Temperature_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Ground_Temperature_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Vegetation_Temperature_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Latent_Heat_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        CLM_Sensible_Heat_parm_infl = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock),dtype=numpy.float32)
        
        CLM_Soil_Moisture_Ensemble = numpy.zeros((Soil_Layer_Num, Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
        CLM_Soil_Temperature_Ensemble = numpy.zeros((Soil_Layer_Num, Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
        CLM_Surface_Temperature_Ensemble = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
        CLM_Ground_Temperature_Ensemble = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
        CLM_Vegetation_Temperature_Ensemble = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)            
        
        Prop_Grid_Array_Sys_parm_infl[:] = Prop_Grid_Array_Sys_parm_infl_SubBlock[:, :].flatten()
        for Soil_Layer_Index in range(Soil_Layer_Num):
            CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,:] = CLM_Soil_Moisture_parm_infl_SubBlock[Soil_Layer_Index,:, :].flatten()
            CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,:] = CLM_Soil_Temperature_parm_infl_SubBlock[Soil_Layer_Index,:, :].flatten()
        CLM_Surface_Temperature_parm_infl[:] = CLM_Surface_Temperature_parm_infl_SubBlock[:, :].flatten()
        CLM_Ground_Temperature_parm_infl[:] = CLM_Ground_Temperature_parm_infl_SubBlock[:, :].flatten()
        CLM_Vegetation_Temperature_parm_infl[:] = CLM_Vegetation_Temperature_parm_infl_SubBlock[:, :].flatten() 
        
        parm_infl = []
        CLM_Soil_Moisture_parm_infl_Bias = []
        CLM_Soil_Temperature_parm_infl_Bias = []
        CLM_Surface_Temperature_parm_infl_Bias = []
        CLM_Ground_Temperature_parm_infl_Bias = []
        CLM_Vegetation_Temperature_parm_infl_Bias = []
        
        Prop_Grid_Array_Sys_parm_infl_Bias = []
        Observation_parm_infl_Bias = []
        
        Prop_Grid_Array_Bias = []
        Observation_Bias = []
        
        CLM_Soil_Moisture_Ensemble_Bias = []
        CLM_Soil_Temperature_Ensemble_Bias = []
        CLM_Surface_Temperature_Ensemble_Bias = []
        CLM_Ground_Temperature_Ensemble_Bias = []
        CLM_Vegetation_Temperature_Ensemble_Bias = []
            
        if (not Parameter_Optimization_Flag) and (not ((numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1))):
            
            print "------------------ State Definition"
            for Ens_Index in range(Ensemble_Number):
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,:,Ens_Index] = CLM_Soil_Moisture_Ensemble_Mat_SubBlock[Soil_Layer_Index,:, :,Ens_Index].flatten()
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    CLM_Soil_Temperature_Ensemble[Soil_Layer_Index,:,Ens_Index] = CLM_Soil_Temperature_Ensemble_Mat_SubBlock[Soil_Layer_Index,:, :,Ens_Index].flatten()
                CLM_Ground_Temperature_Ensemble[:,Ens_Index] = CLM_Ground_Temperature_Ensemble_Mat_SubBlock[:, :,Ens_Index].flatten()
                CLM_Vegetation_Temperature_Ensemble[:,Ens_Index] = CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock[:, :,Ens_Index].flatten()            
            
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Soil_Moisture":
                
                if SensorType_Sub == "InSitu":
                    if Soil_Layer_Index_DA == 1:
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[0,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        for Soil_Layer_Index in range(2,Soil_Layer_Num - 5,1):
                            CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    elif Soil_Layer_Index_DA == 2:
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[1,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[0,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        for Soil_Layer_Index in range(3,Soil_Layer_Num - 5,1):
                            CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    elif Soil_Layer_Index_DA == 3:
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[2,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[1,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[0,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        for Soil_Layer_Index in range(4,Soil_Layer_Num - 5,1):
                            CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    elif Soil_Layer_Index_DA == 4:
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[3,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[2,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[1,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[0,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_SysModel_SubBlock))
                        E0_ObsModel_SubBlock = numpy.vstack((CLM_Soil_Moisture_Col,E0_ObsModel_SubBlock))
                        for Soil_Layer_Index in range(5,Soil_Layer_Num - 5,1):
                            CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    
                    parm_infl = CLM_Soil_Moisture_parm_infl[0,Mask_Index_Vector]
                    for Soil_Layer_Index in range(1,Soil_Layer_Num - 5,1):
                        parm_infl = numpy.hstack((parm_infl,CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
                                        
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        # Couple the Soil Temperature to State Vector
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            CLM_Soil_Temperature_Col = CLM_Soil_Temperature_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Temperature_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Temperature_Col))
                        
                        parm_infl = numpy.hstack((parm_infl,CLM_Vegetation_Temperature_parm_infl[Mask_Index_Vector]))
                        parm_infl = numpy.hstack((parm_infl,CLM_Ground_Temperature_parm_infl[Mask_Index_Vector]))
                        
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            parm_infl = numpy.hstack((parm_infl,CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
                        
                else:
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    
                        
                    parm_infl = CLM_Soil_Moisture_parm_infl[0,Mask_Index_Vector]
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        parm_infl = numpy.hstack((parm_infl,CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
                    
                           
                    if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                        # Couple the Soil Temperature to State Vector
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            CLM_Soil_Temperature_Col = CLM_Soil_Temperature_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                            E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Temperature_Col))
                            E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Temperature_Col))
                        
                        parm_infl = numpy.hstack((parm_infl,CLM_Vegetation_Temperature_parm_infl[Mask_Index_Vector]))
                        parm_infl = numpy.hstack((parm_infl,CLM_Ground_Temperature_parm_infl[Mask_Index_Vector]))
                        
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            parm_infl = numpy.hstack((parm_infl,CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
            
            elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature"):
                E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Vegetation_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock, CLM_Ground_Temperature_Ensemble[Mask_Index_Vector, 0:Ensemble_Number]))
                
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    CLM_Soil_Temperature_Col = CLM_Soil_Temperature_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                    E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Temperature_Col))
                    E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Temperature_Col))
                
                if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature":
                    parm_infl = CLM_Surface_Temperature_parm_infl[Mask_Index_Vector]
                
                parm_infl = numpy.hstack((parm_infl,CLM_Vegetation_Temperature_parm_infl[Mask_Index_Vector]))
                parm_infl = numpy.hstack((parm_infl,CLM_Ground_Temperature_parm_infl[Mask_Index_Vector]))
                
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    parm_infl = numpy.hstack((parm_infl,CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
                
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    # Couple Soil Moisture to State Vector
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        CLM_Soil_Moisture_Col = CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,Mask_Index_Vector,0:Ensemble_Number]
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,CLM_Soil_Moisture_Col))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,CLM_Soil_Moisture_Col))
                    
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        parm_infl = numpy.hstack((parm_infl,CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Mask_Index_Vector]))
                        
        else:
            
            NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
            print "------------------------- Parameter Definition"
            if Soil_Par_Sens_Dim >= 1:
                if Def_Print:
                    print "**********************************************************************Optimize Soil Parameter"
                Parameter_Space_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Soil_Par_Sens, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                Par_Index_Sub = 0
                for Par_Index in range(Dim_Soil_Par):
                    if Soil_Par_Sens[Par_Index]:
                        Parameter_Col = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
                        for Ens_Index in range(Ensemble_Number):
                            Parameter_Col[:,Ens_Index] = Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:].flatten()
                            #Parameter_Min_Max[Par_Index,0] = numpy.min(Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:])
                            #Parameter_Min_Max[Par_Index,1] = numpy.max(Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:])
                            #print numpy.min(Parameter_Col[:,Ens_Index]),numpy.max(Parameter_Col[:,Ens_Index]),numpy.min(E0_ObsModel_SubBlock[:,Ens_Index]),numpy.max(E0_ObsModel_SubBlock[:,Ens_Index])
                            #Parameter_Col[:,Ens_Index] = imadjust.imadjust(Parameter_Col[:,Ens_Index],numpy.min(Parameter_Col[:,Ens_Index]),numpy.max(Parameter_Col[:,Ens_Index]),numpy.min(E0_ObsModel_SubBlock[:,Ens_Index]),numpy.max(E0_ObsModel_SubBlock[:,Ens_Index]))
                        Parameter_Col_Temp = Parameter_Col[Mask_Index_Vector,:]
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,Parameter_Col_Temp))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,Parameter_Col_Temp))
                        del Parameter_Col, Parameter_Col_Temp
                        Par_Index_Sub += 1
                
                parm_infl_Space_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][Soil_Par_Sens, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                
                if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Soil_Moisture":
                    parm_infl = CLM_Soil_Moisture_parm_infl[Soil_Layer_Index_DA,Mask_Index_Vector]
                elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature":
                    parm_infl = CLM_Surface_Temperature_parm_infl[Mask_Index_Vector]
                else:
                    parm_infl = Prop_Grid_Array_Sys_parm_infl[Mask_Index_Vector]
                    
                Par_Index_Sub = 0
                for Par_Index in range(Dim_Soil_Par):
                    if Soil_Par_Sens[Par_Index]:
                        parm_infl_Col = parm_infl_Space_SubBlock[Par_Index_Sub,:,:].flatten()
                        parm_infl = numpy.hstack((parm_infl,parm_infl_Col[Mask_Index_Vector]))
                        del parm_infl_Col
                        Par_Index_Sub += 1
                
                del Parameter_Space_SubBlock, parm_infl_Space_SubBlock
                
                
            if PFT_Par_Sens_Dim >= 1:
                if Def_Print:
                    print "**********************************************************************Optimize PFT Parameter"
                #print PFT_Par_Sens,numpy.shape(Parameter_PFT_Space_Ensemble)
                Parameter_Space_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,PFT_Par_Sens, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                Par_Index_Sub = 0
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        Parameter_Col = numpy.zeros((Observation_NLats_SubBlock * Observation_NLons_SubBlock, Ensemble_Number),dtype=numpy.float32)
                        for Ens_Index in range(Ensemble_Number):
                            Parameter_Col[:,Ens_Index] = Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:].flatten()
                            #Parameter_Min_Max[Par_Index,0] = numpy.min(Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:])
                            #Parameter_Min_Max[Par_Index,1] = numpy.max(Parameter_Space_SubBlock[Ens_Index,Par_Index_Sub,:,:])
                            #print numpy.min(Parameter_Col[:,Ens_Index]),numpy.max(Parameter_Col[:,Ens_Index]),numpy.min(E0_ObsModel_SubBlock[:,Ens_Index]),numpy.max(E0_ObsModel_SubBlock[:,Ens_Index])
                            #Parameter_Col[:,Ens_Index] = imadjust.imadjust(Parameter_Col[:,Ens_Index],numpy.min(Parameter_Col[:,Ens_Index]),numpy.max(Parameter_Col[:,Ens_Index]),numpy.min(E0_ObsModel_SubBlock[:,Ens_Index]),numpy.max(E0_ObsModel_SubBlock[:,Ens_Index]))
                        Parameter_Col_Temp = Parameter_Col[Mask_Index_Vector,:]
                        E0_SysModel_SubBlock = numpy.vstack((E0_SysModel_SubBlock,Parameter_Col_Temp))
                        E0_ObsModel_SubBlock = numpy.vstack((E0_ObsModel_SubBlock,Parameter_Col_Temp))
                        del Parameter_Col, Parameter_Col_Temp
                        Par_Index_Sub += 1
                
                #print "PFT_Par_Sens",PFT_Par_Sens
                parm_infl_Space_SubBlock = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][PFT_Par_Sens, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                #print numpy.shape(parm_infl_Space_SubBlock)
                
                if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "PFT_Moisture":
                    parm_infl = CLM_Soil_Moisture_parm_infl[Soil_Layer_Index_DA,Mask_Index_Vector]
                elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature":
                    parm_infl = CLM_Surface_Temperature_parm_infl[Mask_Index_Vector]
                else:
                    parm_infl = Prop_Grid_Array_Sys_parm_infl[Mask_Index_Vector]
                    
                Par_Index_Sub = 0
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        parm_infl_Col = parm_infl_Space_SubBlock[Par_Index_Sub,:,:].flatten()
                        parm_infl = numpy.hstack((parm_infl,parm_infl_Col[Mask_Index_Vector]))
                        del parm_infl_Col
                        Par_Index_Sub += 1
                
                del Parameter_Space_SubBlock, parm_infl_Space_SubBlock
            
            
            NC_File_Out_Assimilation_2_Parameter.close()
        #os.abort()
        #E0_ObsModel_All = E0_SysModel_All[:,0:Ensemble_Number]
        #########################################################################################################################################        
        
        PF_PRESSURE_Ensemble_Mat_SubBlock = []
        PF_SATURATION_Ensemble_Mat_SubBlock = []
        PF_SATURATION_parm_infl_SubBlock = []
        
        if Def_Print:
            start = time.time()
        
        print "numpy.size(numpy.where(Observation_Matrix_SubBlock != NAvalue)[0])",numpy.size(numpy.where(Observation_Matrix_SubBlock != NAvalue)[0])
        print "numpy.shape(E0_SysModel_SubBlock)[0]",numpy.shape(E0_SysModel_SubBlock)[0]
        if numpy.size(numpy.where(Observation_Matrix_SubBlock != NAvalue)[0]) > 0 and numpy.shape(E0_SysModel_SubBlock)[0] > 0:
            if Def_Print:
                print "numpy.shape(E0_SysModel_SubBlock),numpy.shape(E0_ObsModel_SubBlock),numpy.shape(parm_infl)",numpy.shape(E0_SysModel_SubBlock),numpy.shape(E0_ObsModel_SubBlock),numpy.shape(parm_infl)
                print "numpy.min(E0_SysModel_SubBlock),numpy.min(E0_ObsModel_SubBlock),numpy.min(parm_infl)",numpy.min(E0_SysModel_SubBlock),numpy.min(E0_ObsModel_SubBlock),numpy.min(parm_infl)
                print "numpy.max(E0_SysModel_SubBlock),numpy.max(E0_ObsModel_SubBlock),numpy.max(parm_infl)",numpy.max(E0_SysModel_SubBlock),numpy.max(E0_ObsModel_SubBlock),numpy.max(parm_infl)
                print "numpy.size(numpy.where(Observation_Matrix[Observation_Matrix_Index,::]!= NAvalue))",numpy.size(numpy.where(Observation_Matrix[Observation_Matrix_Index,::]!= NAvalue)[0])
            if Def_Print >= 2:
                print "E0_SysModel_SubBlock,E0_ObsModel_SubBlock,parm_infl",E0_SysModel_SubBlock,E0_ObsModel_SubBlock,parm_infl
            
            #-----------------------###################### Call CLM_Assim_Common"
            Analysis_Grid_SubBlock, Analysis_Grid_Array_SubBlock, Localization_Map_Mask_SubBlock, \
            CLM_Ground_Temperature_Ensemble_Mat_SubBlock, CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock, CLM_Soil_Moisture_Ensemble_Mat_SubBlock, CLM_Soil_Temperature_Ensemble_Mat_SubBlock, PF_PRESSURE_Ensemble_Mat_SubBlock, PF_SATURATION_Ensemble_Mat_SubBlock, \
            Prop_Grid_Array_Sys_parm_infl_SubBlock, CLM_Latent_Heat_parm_infl_SubBlock, CLM_Surface_Temperature_parm_infl_SubBlock, CLM_Ground_Temperature_parm_infl_SubBlock,CLM_Vegetation_Temperature_parm_infl_SubBlock, CLM_Soil_Moisture_parm_infl_SubBlock,CLM_Soil_Temperature_parm_infl_SubBlock, PF_SATURATION_parm_infl_SubBlock,\
            CLM_Ground_Temperature_Ensemble_Mat_Bias_SubBlock, CLM_Vegetation_Temperature_Ensemble_Mat_Bias_SubBlock, CLM_Soil_Moisture_Ensemble_Mat_Bias_SubBlock, CLM_Soil_Temperature_Ensemble_Mat_Bias_SubBlock, \
            CLM_Surface_Temperature_parm_infl_Bias_SubBlock, CLM_Ground_Temperature_parm_infl_Bias_SubBlock,CLM_Vegetation_Temperature_parm_infl_Bias_SubBlock, CLM_Soil_Moisture_parm_infl_Bias_SubBlock,CLM_Soil_Temperature_parm_infl_Bias_SubBlock, \
            Prop_Grid_Array_Bias_SubBlock, Observation_Bias_SubBlock, Prop_Grid_Array_Sys_parm_infl_Bias_SubBlock, Observation_parm_infl_Bias_SubBlock, \
            Parameter_Soil_Space_Ensemble_SubBlock, Parameter_Soil_Space_parm_infl_SubBlock, Parameter_Veg_Space_Ensemble_SubBlock, Parameter_Veg_Space_parm_infl_SubBlock, Parameter_PFT_Space_Ensemble_SubBlock, Parameter_PFT_Space_parm_infl_SubBlock, \
            Parameter_Hard_Space_Ensemble_SubBlock, Parameter_Hard_Space_parm_infl_SubBlock, Innovation_State_SubBlock, Increments_State_SubBlock = \
            CLM_Assim_Common(Block_Index, Model_Driver, Def_PP, Def_First_Run, Def_Print, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs, eps, msw_infl, parm_infl, Post_Inflation_Alpha, Def_ParFor, Observation_NLats_SubBlock, Observation_NLons_SubBlock, Ensemble_Number, Ensemble_Number_Predict,  
                            Call_Gstat_Flag, Assim_Algorithm_Name, Model_State_SubBlock, E0_SysModel_SubBlock, E0_ObsModel_SubBlock, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, X_Left_SubBlock, X_Right_SubBlock, Y_Lower_SubBlock, Y_Upper_SubBlock, Proj_String, Z_Resolution, 
                            Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Variable_List,
                             Grid_Resolution_CEA, Prop_Grid_Array_SubBlock, Prop_Grid_Array_H_Trans_SubBlock, Model_Variance_SubBlock, Write_DA_File_Flag, Mask_SubBlock, Mask_Index_SubBlock, Land_Mask_Data_SubBlock, Observation_Variance_SubBlock, SensorQuantity_Sub, SensorQuantity_Index,
                             Observation_NLats_SubBlock, Observation_NLons_SubBlock, Observation_Longitude_SubBlock, Observation_Latitude_SubBlock, Observation_Matrix_SubBlock, DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, Soil_Layer_Index_DA, Soil_Layer_Num, ParFlow_Layer_Num, omp_get_num_procs_ParFor, 
                             Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, NSLOTS, DAS_Output_Path, Region_Name,
                             Variable_Assimilation_Flag, Teta_Residual_SubBlock, Teta_Saturated_SubBlock, Teta_Field_Capacity_SubBlock, Teta_Wilting_Point_SubBlock, SensorType_Sub, SensorVariable_Sub, SensorResolution_Sub, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, 
                             Observation_Corelation_Par[Observation_Matrix_Index,::],  Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, Low_Ratio_Par, High_Ratio_Par, 
                             Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                             CLM_Ground_Temperature_Ensemble_Mat_SubBlock,CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock, CLM_Soil_Moisture_Ensemble_Mat_SubBlock,CLM_Soil_Temperature_Ensemble_Mat_SubBlock, PF_PRESSURE_Ensemble_Mat_SubBlock, PF_SATURATION_Ensemble_Mat_SubBlock,
                             Prop_Grid_Array_Sys_parm_infl_SubBlock, [], CLM_Surface_Temperature_parm_infl_SubBlock, CLM_Ground_Temperature_parm_infl_SubBlock,CLM_Vegetation_Temperature_parm_infl_SubBlock, CLM_Soil_Moisture_parm_infl_SubBlock,CLM_Soil_Temperature_parm_infl_SubBlock, PF_SATURATION_parm_infl_SubBlock,
                             CLM_Ground_Temperature_Ensemble_Mat_Bias_SubBlock,CLM_Vegetation_Temperature_Ensemble_Mat_Bias_SubBlock, CLM_Soil_Moisture_Ensemble_Mat_Bias_SubBlock,CLM_Soil_Temperature_Ensemble_Mat_Bias_SubBlock, 
                             CLM_Surface_Temperature_parm_infl_Bias_SubBlock, CLM_Ground_Temperature_parm_infl_Bias_SubBlock,CLM_Vegetation_Temperature_parm_infl_Bias_SubBlock, CLM_Soil_Moisture_parm_infl_Bias_SubBlock,CLM_Soil_Temperature_parm_infl_Bias_SubBlock,
                             Prop_Grid_Array_Bias_SubBlock, Observation_Bias_SubBlock, Prop_Grid_Array_Sys_parm_infl_Bias_SubBlock, Observation_parm_infl_Bias_SubBlock, Def_CDF_Matching, Plot_Analysis, Parameter_Optimization_Flag,
                             Start_Month, maxpft, Feedback_Assim, Dim_Soil_Par, Soil_Par_Sens, Dim_Veg_Par, Veg_Par_Sens, Dim_PFT_Par, PFT_Par_Sens, Dim_Hard_Par, Hard_Par_Sens, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Parameter_Soil_Space_Ensemble_SubBlock, Parameter_Soil_Space_parm_infl_SubBlock, 
                             Parameter_Veg_Space_Ensemble_SubBlock, Parameter_Veg_Space_parm_infl_SubBlock, Parameter_PFT_Space_Ensemble_SubBlock, Parameter_PFT_Space_parm_infl_SubBlock, Parameter_Hard_Space_Ensemble_SubBlock, Parameter_Hard_Space_parm_infl_SubBlock, Parameter_Min_Max, 
                             Soil_Layer_Thickness_Ratio, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, 
                             Par_Soil_Uniform_STD, Par_Veg_Uniform_STD, Par_PFT_Uniform_STD, Par_Hard_Uniform_STD, 
                             Saturation_SSat_SubBlock, Saturation_SRes_SubBlock, Saturation_N_SubBlock, Saturation_Alpha_SubBlock, DateString_Plot, *vartuple)
        else:
            Analysis_Grid_SubBlock = numpy.mean(Prop_Grid_Array_SubBlock[:, :, :],axis=0)
            for Ens_Index in range(Ensemble_Number):
                Analysis_Grid_Array_SubBlock[Ens_Index,:,:] = Prop_Grid_Array_SubBlock[Ens_Index, :, :]
#                if SensorType != "AMSR_E":
#                    Observation_Matrix[(Sub_Block_Index_Row*Row_Numbers_SubBlock):((Sub_Block_Index_Row+1)*Row_Numbers_SubBlock),(Sub_Block_Index_Col*Col_Numbers_SubBlock):((Sub_Block_Index_Col+1)*Col_Numbers_SubBlock)] = Observation_Matrix_SubBlock
#                else:
#                    Observation_Matrix = Observation_Matrix_SubBlock
        #Observation_Matrix = Observation_Matrix_SubBlock
        # Compose the Full Analysis Grid
        
        if Def_Print:
            end = time.time()
            print 'Time of Block',Block_Index,'is: ', (end - start), 'Seconds'
    
        if Def_Print:
            print "Before",Observation_Box_Row_Index_Start,Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start,Observation_Box_Col_Index_End
        if Observation_Box_Row_Index_Start == 0 and Sub_Block_Row_Start > 0:
            Observation_Box_Row_Index_Start = Sub_Block_Row_Start
        elif Observation_Box_Row_Index_Start > 0:
            Observation_Box_Row_Index_Start = Observation_Box
            
        if Observation_Box_Row_Index_End < Row_Numbers:
            Observation_Box_Row_Index_End = Observation_NLats_SubBlock - Observation_Box
        elif Observation_Box_Row_Index_End == Row_Numbers and Observation_NLats_SubBlock <= Row_Numbers:
            Observation_Box_Row_Index_End = Observation_Box_Row_Index_Start + Row_Numbers_SubBlock
            
        if Observation_Box_Col_Index_Start == 0 and Sub_Block_Col_Start > 0:
            Observation_Box_Col_Index_Start = Sub_Block_Col_Start
        elif Observation_Box_Col_Index_Start > 0:
            Observation_Box_Col_Index_Start = Observation_Box
            
        if Observation_Box_Col_Index_End < Col_Numbers:
            Observation_Box_Col_Index_End = Observation_NLons_SubBlock - Observation_Box
        elif Observation_Box_Col_Index_End == Col_Numbers and Observation_NLons_SubBlock <= Col_Numbers:
            Observation_Box_Col_Index_End = Observation_Box_Col_Index_Start + Col_Numbers_SubBlock
        
        if Def_Print:
            print "After",Observation_Box_Row_Index_Start,Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start,Observation_Box_Col_Index_End
        
        Sub_Block_Row_Start = Sub_Block_Row_Start_Array[Block_Index]
        Sub_Block_Row_End = Sub_Block_Row_End_Array[Block_Index]
        Sub_Block_Col_Start = Sub_Block_Col_Start_Array[Block_Index]
        Sub_Block_Col_End = Sub_Block_Col_End_Array[Block_Index]
        
        Row_Numbers_SubBlock_Patch = Sub_Block_Row_End - Sub_Block_Row_Start
        Col_Numbers_SubBlock_Patch = Sub_Block_Col_End - Sub_Block_Col_Start
        if Def_Print:
            print "Row_Numbers_SubBlock_Patch,Col_Numbers_SubBlock_Patch",Row_Numbers_SubBlock_Patch,Col_Numbers_SubBlock_Patch
            
        # Record the State Analysis
        NC_FileName_Out_Block_Assim = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_"+str(Block_Index+1)+".nc"
        if os.path.exists(NC_FileName_Out_Block_Assim):
            os.remove(NC_FileName_Out_Block_Assim)
        
        if Def_Print:
            print 'Write NetCDF File:',NC_FileName_Out_Block_Assim
        
        NC_File_Out_Block_Assim = netCDF4.Dataset(NC_FileName_Out_Block_Assim, 'w', diskless=True, persist=True, format='NETCDF4')
        # Dim the dimensions of NetCDF
        NC_File_Out_Block_Assim.createDimension('lon', Col_Numbers_SubBlock_Patch)
        NC_File_Out_Block_Assim.createDimension('lat', Row_Numbers_SubBlock_Patch)
        NC_File_Out_Block_Assim.createDimension('Soil_Layer_Num', Soil_Layer_Num)
        NC_File_Out_Block_Assim.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
        NC_File_Out_Block_Assim.createDimension('Ensemble_Number', Ensemble_Number)
        NC_File_Out_Block_Assim.createDimension('Dim_CLM_State', Dim_CLM_State)
        NC_File_Out_Block_Assim.createDimension('Dim_Soil_Par', Dim_Soil_Par)
        NC_File_Out_Block_Assim.createDimension('Dim_Veg_Par', Dim_Veg_Par)
        NC_File_Out_Block_Assim.createDimension('Dim_PFT_Par', Dim_PFT_Par)
        NC_File_Out_Block_Assim.createDimension('Dim_Hard_Par', Dim_Hard_Par)
        NC_File_Out_Block_Assim.createDimension('maxpft', maxpft)
        NC_File_Out_Block_Assim.createDimension('Dim_Observation_Quantity', Dim_Observation_Quantity)
        
        NC_File_Out_Block_Assim.createVariable('Prop_Grid_Array_Sys','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Prop_Grid_Array_H_Trans','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True)
                
        if Parameter_Optimization_Flag:
            if Soil_Par_Sens_Dim >= 1:
                NC_File_Out_Block_Assim.createVariable('Parameter_Soil_Space_Ensemble','f4',('Ensemble_Number','Dim_Soil_Par','lat','lon',),zlib=True)
                NC_File_Out_Block_Assim.createVariable('Parameter_Soil_Space_parm_infl','f4',('Dim_Soil_Par','lat','lon',),zlib=True)
            if PFT_Par_Sens_Dim >= 1:
                NC_File_Out_Block_Assim.createVariable('Parameter_PFT_Space_Ensemble','f4',('Ensemble_Number','Dim_PFT_Par','lat','lon',),zlib=True)
                NC_File_Out_Block_Assim.createVariable('Parameter_PFT_Space_parm_infl','f4',('Dim_PFT_Par','lat','lon',),zlib=True)
            
        NC_File_Out_Block_Assim.createVariable('Analysis_Grid','f4',('lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Localization_Map_Mask','f4',('lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Analysis_Grid_Array','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Innovation_State','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Increments_State','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True)
        
        NC_File_Out_Block_Assim.createVariable('CLM_Soil_Moisture_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Soil_Temperature_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_Soil_Ice_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Vegetation_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Ground_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_2m_Air_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_Snow_Depth_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_Snow_Water_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_ROOTFR_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Prop_Grid_Array_Sys_parm_infl','f4',('Dim_CLM_State','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Soil_Moisture_parm_infl','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Soil_Temperature_parm_infl','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
        #NC_File_Out_Block_Assim.createVariable('CLM_Soil_Ice_parm_infl','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Vegetation_Temperature_parm_infl','f4',('lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Ground_Temperature_parm_infl','f4',('lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('CLM_Surface_Temperature_parm_infl','f4',('lat','lon',),zlib=True)
        NC_File_Out_Block_Assim.createVariable('Observation','f4',('Dim_CLM_State','lat','lon',),zlib=True)
                   
        if Def_Print:
            print "Row_Numbers_SubBlock_Patch,Col_Numbers_SubBlock_Patch",Row_Numbers_SubBlock_Patch,Col_Numbers_SubBlock_Patch
        if Def_Print:
            print 'Write NetCDF File:',NC_FileName_Out_Block_Assim
        
        if Parameter_Optimization_Flag:
            if Soil_Par_Sens_Dim >= 1:
                NC_File_Out_Block_Assim.variables['Parameter_Soil_Space_Ensemble'][:, :, :, :] = Parameter_Soil_Space_Ensemble_SubBlock[:, :, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                NC_File_Out_Block_Assim.variables['Parameter_Soil_Space_parm_infl'][:, :, :] = Parameter_Soil_Space_parm_infl_SubBlock[:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
            if PFT_Par_Sens_Dim >= 1:
                NC_File_Out_Block_Assim.variables['Parameter_PFT_Space_Ensemble'][:, :, :, :] = Parameter_PFT_Space_Ensemble_SubBlock[:, :, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                NC_File_Out_Block_Assim.variables['Parameter_PFT_Space_parm_infl'][:, :, :] = Parameter_PFT_Space_parm_infl_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
            
        NC_File_Out_Block_Assim.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :] = Prop_Grid_Array_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Prop_Grid_Array_H_Trans'][:, Prop_Grid_Array_Sys_Index, :, :] = Prop_Grid_Array_H_Trans_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Innovation_State'][:,Prop_Grid_Array_Sys_Index,:,:] = Innovation_State_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Increments_State'][:,Prop_Grid_Array_Sys_Index,:,:] = Increments_State_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Analysis_Grid'][:, :] = Analysis_Grid_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Localization_Map_Mask'][:, :] = Localization_Map_Mask_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Analysis_Grid_Array'][:, Prop_Grid_Array_Sys_Index, :, :] = Analysis_Grid_Array_SubBlock[:,Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Soil_Moisture_Ensemble_Mat'][:, :, :, :] = CLM_Soil_Moisture_Ensemble_Mat_SubBlock[:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        NC_File_Out_Block_Assim.variables['CLM_Soil_Temperature_Ensemble_Mat'][:, :, :, :] = CLM_Soil_Temperature_Ensemble_Mat_SubBlock[:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        NC_File_Out_Block_Assim.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:, :, :] = CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        NC_File_Out_Block_Assim.variables['CLM_Ground_Temperature_Ensemble_Mat'][:, :, :] = CLM_Ground_Temperature_Ensemble_Mat_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End, :]
        NC_File_Out_Block_Assim.variables['Prop_Grid_Array_Sys_parm_infl'][Prop_Grid_Array_Sys_Index,:, :] = Prop_Grid_Array_Sys_parm_infl_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Soil_Moisture_parm_infl'][:, :, :] = CLM_Soil_Moisture_parm_infl_SubBlock[:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Soil_Temperature_parm_infl'][:, :, :] = CLM_Soil_Temperature_parm_infl_SubBlock[:, Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Vegetation_Temperature_parm_infl'][:, :] = CLM_Vegetation_Temperature_parm_infl_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Ground_Temperature_parm_infl'][:, :] = CLM_Ground_Temperature_parm_infl_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['CLM_Surface_Temperature_parm_infl'][:, :] = CLM_Surface_Temperature_parm_infl_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
        NC_File_Out_Block_Assim.variables['Observation'][Prop_Grid_Array_Sys_Index,:,:] = Observation_Matrix_SubBlock[Observation_Box_Row_Index_Start:Observation_Box_Row_Index_End, Observation_Box_Col_Index_Start:Observation_Box_Col_Index_End]
                   
        NC_File_Out_Block_Assim.sync()
        NC_File_Out_Block_Assim.close()
    
    del Observation_Variance_SubBlock,Observation_Longitude_SubBlock,Observation_Latitude_SubBlock,Observation_Matrix_SubBlock
    del E0_SysModel_SubBlock,E0_ObsModel_SubBlock, parm_infl
    del Parameter_Soil_Space_Ensemble_SubBlock, Parameter_Soil_Space_parm_infl_SubBlock, Parameter_Hard_Space_Ensemble_SubBlock, Parameter_Hard_Space_parm_infl_SubBlock, Parameter_Veg_Space_Ensemble_SubBlock
    del Parameter_Veg_Space_parm_infl_SubBlock, Parameter_PFT_Space_Ensemble_SubBlock, Parameter_PFT_Space_parm_infl_SubBlock
    del Prop_Grid_Array_SubBlock, Prop_Grid_Array_H_Trans_SubBlock, Model_State_SubBlock, Model_Variance_SubBlock, Mask_SubBlock, Mask_Index_SubBlock, Land_Mask_Data_SubBlock
    del Teta_Residual_SubBlock, Teta_Saturated_SubBlock, Teta_Field_Capacity_SubBlock, Teta_Wilting_Point_SubBlock, Analysis_Grid_SubBlock, Localization_Map_Mask_SubBlock, Analysis_Grid_Array_SubBlock
    del Innovation_State_SubBlock, Increments_State_SubBlock, Soil_Layer_Thickness_Ratio
    del CLM_Soil_Moisture_Ensemble_Mat_SubBlock, CLM_Soil_Temperature_Ensemble_Mat_SubBlock, CLM_Vegetation_Temperature_Ensemble_Mat_SubBlock
    del CLM_Ground_Temperature_Ensemble_Mat_SubBlock, CLM_Soil_Moisture_parm_infl_SubBlock, CLM_Soil_Temperature_parm_infl_SubBlock
    del Prop_Grid_Array_Sys_parm_infl_SubBlock, CLM_Vegetation_Temperature_parm_infl_SubBlock, CLM_Ground_Temperature_parm_infl_SubBlock, CLM_Surface_Temperature_parm_infl_SubBlock, CLM_Latent_Heat_parm_infl_SubBlock
    del Prop_Grid_Array_Sys_parm_infl, CLM_Soil_Moisture_parm_infl, CLM_Soil_Temperature_parm_infl, CLM_Surface_Temperature_parm_infl, CLM_Ground_Temperature_parm_infl, CLM_Vegetation_Temperature_parm_infl,CLM_Latent_Heat_parm_infl, CLM_Sensible_Heat_parm_infl
    del CLM_Soil_Moisture_Ensemble, CLM_Soil_Temperature_Ensemble, CLM_Surface_Temperature_Ensemble, CLM_Ground_Temperature_Ensemble, CLM_Vegetation_Temperature_Ensemble
    del Mask_Sub, Mask_Index_Sub, Mask_Index_Sub_NC, Model_Variance, Observation_Variance, Observation_Longitude, Observation_Latitude, Observation_Matrix,
                           
    del numexpr_a,numexpr_b,numexpr_c
    CLM_Soil_Moisture_Col = []
    
    del CLM_Soil_Moisture_Ensemble_Mat_Bias_SubBlock, CLM_Soil_Temperature_Ensemble_Mat_Bias_SubBlock, CLM_Vegetation_Temperature_Ensemble_Mat_Bias_SubBlock
    del CLM_Ground_Temperature_Ensemble_Mat_Bias_SubBlock, CLM_Soil_Moisture_parm_infl_Bias_SubBlock, CLM_Soil_Temperature_parm_infl_Bias_SubBlock
    del CLM_Vegetation_Temperature_parm_infl_Bias_SubBlock, CLM_Ground_Temperature_parm_infl_Bias_SubBlock, CLM_Surface_Temperature_parm_infl_Bias_SubBlock
    del CLM_Soil_Moisture_parm_infl_Bias, CLM_Soil_Temperature_parm_infl_Bias, CLM_Surface_Temperature_parm_infl_Bias, CLM_Ground_Temperature_parm_infl_Bias, CLM_Vegetation_Temperature_parm_infl_Bias,
    del CLM_Soil_Moisture_Ensemble_Bias, CLM_Soil_Temperature_Ensemble_Bias, CLM_Surface_Temperature_Ensemble_Bias, CLM_Ground_Temperature_Ensemble_Bias, CLM_Vegetation_Temperature_Ensemble_Bias
    del Prop_Grid_Array_Bias, Observation_Bias, Prop_Grid_Array_Bias_SubBlock, Observation_Bias_SubBlock, Prop_Grid_Array_Sys_parm_infl_Bias_SubBlock, Observation_parm_infl_Bias_SubBlock
    del Prop_Grid_Array_Sys_parm_infl_Bias, Observation_parm_infl_Bias
    
    gc.collect()
    del gc.garbage[:]
    
    
    return

#"******************************************************************Assimilation*******************************************************************************************"
def Assimilation_Update(mpi4py_comm, mpi4py_rank, mpi4py_name, Model_Driver, NSLOTS, finidat_initial_Array, Def_ParFor, Def_Region, Def_Initial, Irrig_Scheduling, Irrigation_Hours,  Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Model_Path, CLM_Flag, Def_PP, job_server_node_array, active_nodes_server,
                  Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, Weather_Forecast_Days, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water, N0, nlyr,
                  DAS_Data_Path, DAS_Depends_Path, DasPy_Path, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs, eps, msw_infl, Plot_Analysis, Def_Figure_Output, DateString_Plot, 
                  Def_Write_Initial,DA_Flag, Write_DA_File_Flag, Mask, Mask_Index, COSMOS_Circle_Array, COSMOS_Circle_Index_Array, COSMOS_Circle_Num_Array, Call_Gstat_Flag, mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Station_XY_Index, Station_XY, Observation_Box,
                  Variable_Assimilation_Flag, Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, Model_Variance,  Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, PP_Servers_Per_Node, Def_CESM_Multi_Instance, PP_Port,
                  Z_Resolution, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Proj_String, MODEL_CEA_X, MODEL_CEA_Y, Hydraulic_File_Name, Assim_Algorithm_Name, Low_Ratio_Par, High_Ratio_Par, Post_Inflation_Alpha_State, irrig_nsteps_per_day, PFT_Num, 
                  Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset, fpftcon_name, Crop_Sum, 
                  Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Dim_Observation_Quantity, 
                  Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array, finidat_name,
                  Ensemble_Number, Ensemble_Number_Predict,  Soil_Layer_Num, Snow_Layer_Num, maxpft, Forcing_File_Path, dtime, Observation_Path, Dim_CLM_State, Dim_Obs_Type, CLM_NA, NAvalue, Variable_List, ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor,
                  Grid_Resolution_CEA, Grid_Resolution_GEO, SensorQuantity_Sub, SensorType_Sub, SensorVariable_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub,Analysis_Variable_Name, Soil_Layer_Index_DA,
                  Observation_Matrix, Observation_Variance, Observation_Latitude, Observation_Longitude, Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Corelation_Par,
                  octave,r,Def_CDF_Matching, numrad, cols1d_ixy, cols1d_jxy, pfts1d_ixy, pfts1d_jxy, cols1d_ityplun, pfts1d_ityplun, column_len, pft_len, pfts1d_itypveg, pfts1d_ci,
                  diskless_flag, persist_flag, Forcing_File_Path_Home, Forcing_File_Path_Array, history_file_name, Constant_File_Name, Run_Dir_Array, Feedback_Assim,
                  Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, \
                  Analysis_Grid, Localization_Map_Mask, ObsModel_Mat, ObsModel_Variance_Mat, Prop_Grid_Array_Sys_Index, Observation_Matrix_Index, Mask_Sub, Mask_Index_Sub, 
                  SensorQuantity_Index, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, ParFlow_Layer_Num,
                  NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, \
                  NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Bias_Copy, NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy, NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy, 
                  NC_FileName_Parameter_Space_Single, DAS_Output_Path, COSMIC_Py, window, memory_profiler, COSMIC, finidat_name_string, Observation_Time_File_Path):
            
    ##################################################################################################
    
    if True:
        
        Parameter_Optimization_Flag = 0
        Bias_Estimation_Option_Model_Assim = numpy.asarray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])   # Model Bias
        Bias_Estimation_Option_Obs_Assim = numpy.asarray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])   # Observation Bias
        Soil_Par_Sens = []
        Veg_Par_Sens = []
        PFT_Par_Sens = []
        Hard_Par_Sens = []
        
        if Def_PP and (not PDAF_Assim_Framework == 2) and (Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) > 1 and len(active_nodes_server) > 1:
            print "********************************************** Using PP to Accelerate Block_Assim"
            
            Job_Num_Per_Node = int(numpy.ceil(float(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) / len(active_nodes_server)))
            
            if Job_Num_Per_Node == 0:
                Job_Num_Per_Node = 1
            job_server_node_results = []
            job_server_node_results_wise = [[] for i in range(len(active_nodes_server))]
            
            print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"
            Block_Index = 0
            
            Node_Status = numpy.zeros(len(active_nodes_server),dtype=numpy.bool)
            Node_Status[:] = True
            
            while Block_Index < Sub_Block_Ratio_Row*Sub_Block_Ratio_Col:
                if numpy.size(numpy.where(Node_Status==True)) > 0:
                    Node_Index = numpy.where(Node_Status==True)[0][0]
                    print "***********************Node_Index",Node_Index,"Block_Index",Block_Index,"is submitted!"
                    job_server_node = job_server_node_array[numpy.min([Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]                      
                                        
                    job_server_node_results.append(job_server_node.submit(Block_Assim, args=(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                            Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                            Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                            Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                            Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                            Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                            Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_State, Def_ParFor, Ensemble_Number_Predict,
                            Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                            Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                            Observation_Corelation_Par, Bias_Estimation_Option_Model_Assim, Bias_Estimation_Option_Obs_Assim, Low_Ratio_Par, High_Ratio_Par,
                            Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                            DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                            NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path),
                          depfuncs=(CLM_Assim_Common, Check_Outliers, ParFor_Fusion, ParFor_H_Operator, ParFor_Texture_Check, ParFor_Check_Outliers, ParFor_Check_Outliers_NA,),
                          modules=("numpy", "netCDF4", "sys", "os", "re", "gc", "imp", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "signal", "gc", "imp", "getpass", "calendar", "glob","scipy.stats", 'scipy.weave'), group='Block_Assim'))
                    
                    job_server_node_results_wise[Node_Index] = job_server_node_results[Block_Index]
                    
                    Node_Status[Node_Index] = False
                    
                    Block_Index = Block_Index + 1
                
                if Block_Index >= len(active_nodes_server):
                    for job in job_server_node_results_wise:
                        if job != [] and job.finished:
                            Node_Index = job_server_node_results_wise.index(job)
                            print "*********************************************************************Node_Index",Node_Index,"is finished!"
                            Node_Status[Node_Index] = True
                            job_server_node_results_wise[Node_Index] = []
            
            for job_server_node in job_server_node_array:
                job_server_node.wait()
                if Def_Print >= 2:
                    job_server_node.print_stats()
            
            if Def_Print:
                if len(job_server_node_results) > 0:
                    for job in job_server_node_results:
                        job_index = job_server_node_results.index(job)
                        if job_index > (Ensemble_Number - 1):
                            break
                        print "Results of ",job_index,"is", job()
                    
        else:
            print "********* Run Block_Assim Sequentially"
            if PDAF_Assim_Framework == 2:
                DAS_Driver_Common.Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                
            for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
                
                Block_Assim(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                            Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                            Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                            Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                            Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                            Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                            Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_State, Def_ParFor, Ensemble_Number_Predict,
                            Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                            Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                            Observation_Corelation_Par, Bias_Estimation_Option_Model_Assim, Bias_Estimation_Option_Obs_Assim, Low_Ratio_Par, High_Ratio_Par,
                            Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, 
                            Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                            DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, 
                            NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                            NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path, octave, r)
        
        print "Write NC_File_Out_Assimilation_2_Initial.nc"
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
        
        Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, :, :, :]
        Prop_Grid_Array_Sys_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys_parm_infl'][:, :, :]
        CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:, :, :, :]
        CLM_Soil_Moisture_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_parm_infl'][:, :, :]
        CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:, :, :, :]
        CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:, :, :]
        CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:, :, :]
        CLM_Soil_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_parm_infl'][:, :, :]
        CLM_Vegetation_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_parm_infl'][:, :]
        CLM_Ground_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_parm_infl'][:, :]
        CLM_Surface_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Surface_Temperature_parm_infl'][:, :]
        Innovation_State = NC_File_Out_Assimilation_2_Diagnostic.variables['Innovation_State'][:,:, :, :]
        Increments_State = NC_File_Out_Assimilation_2_Diagnostic.variables['Increments_State'][:,:, :, :]
        Observation = NC_File_Out_Assimilation_2_Diagnostic.variables['Observation'][:,:, :]
        
        for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
            print "Block_Index",Block_Index
            Sub_Block_Row_Start = Sub_Block_Row_Start_Array[Block_Index]
            Sub_Block_Row_End = Sub_Block_Row_End_Array[Block_Index]
            Sub_Block_Col_Start = Sub_Block_Col_Start_Array[Block_Index]
            Sub_Block_Col_End = Sub_Block_Col_End_Array[Block_Index]
                                    
            NC_FileName_Out_Block = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_"+str(Block_Index+1)+".nc"
            
            NC_File_Out_Block = netCDF4.Dataset(NC_FileName_Out_Block, 'r')
            
            Analysis_Grid[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid'][:,:]
            Localization_Map_Mask[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Localization_Map_Mask'][:,:]
            
            Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid_Array'][:,Prop_Grid_Array_Sys_Index,:,:]
            Prop_Grid_Array_Sys_parm_infl[Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Prop_Grid_Array_Sys_parm_infl'][Prop_Grid_Array_Sys_Index,:,:]
            CLM_Soil_Moisture_Ensemble_Mat[:, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End, :] = NC_File_Out_Block.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:]
            CLM_Soil_Moisture_parm_infl[:, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['CLM_Soil_Moisture_parm_infl'][:,:,:]
            CLM_Soil_Temperature_Ensemble_Mat[:, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End, :] = NC_File_Out_Block.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:]
            CLM_Vegetation_Temperature_Ensemble_Mat[Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End, :] = NC_File_Out_Block.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:]
            CLM_Ground_Temperature_Ensemble_Mat[Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End, :] = NC_File_Out_Block.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:]
            CLM_Soil_Temperature_parm_infl[:, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['CLM_Soil_Temperature_parm_infl'][:,:,:]
            CLM_Vegetation_Temperature_parm_infl[Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['CLM_Vegetation_Temperature_parm_infl'][:,:]
            CLM_Ground_Temperature_parm_infl[Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['CLM_Ground_Temperature_parm_infl'][:,:]
            CLM_Surface_Temperature_parm_infl[Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['CLM_Surface_Temperature_parm_infl'][:,:]
            Innovation_State[:,Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Innovation_State'][:,Prop_Grid_Array_Sys_Index,:,:]
            Increments_State[:,Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Increments_State'][:,Prop_Grid_Array_Sys_Index,:,:]
            Observation[:,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Observation'][:,:,:]
            
            NC_File_Out_Block.close()
        
        
        NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, Prop_Grid_Array_Sys_Index, :, :] = Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, :, :]
        NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys_parm_infl'][Prop_Grid_Array_Sys_Index, :, :] = Prop_Grid_Array_Sys_parm_infl[Prop_Grid_Array_Sys_Index, :, :]
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:, :, :, :] = CLM_Soil_Moisture_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_parm_infl'][:, :, :] = CLM_Soil_Moisture_parm_infl
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:, :, :, :] = CLM_Soil_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:, :, :] = CLM_Vegetation_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:, :, :] = CLM_Ground_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_parm_infl'][:, :, :] = CLM_Soil_Temperature_parm_infl
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_parm_infl'][:, :] = CLM_Vegetation_Temperature_parm_infl
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_parm_infl'][:, :] = CLM_Ground_Temperature_parm_infl
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Surface_Temperature_parm_infl'][:, :] = CLM_Surface_Temperature_parm_infl
        NC_File_Out_Assimilation_2_Diagnostic.variables['Innovation_State'][:,Prop_Grid_Array_Sys_Index, :, :] = Innovation_State[:,Prop_Grid_Array_Sys_Index, :, :]
        NC_File_Out_Assimilation_2_Diagnostic.variables['Increments_State'][:,Prop_Grid_Array_Sys_Index, :, :] = Increments_State[:,Prop_Grid_Array_Sys_Index, :, :]
        NC_File_Out_Assimilation_2_Diagnostic.variables['Observation'][:,:, :] = Observation
            
        del Analysis_Grid_Array,Prop_Grid_Array_Sys_parm_infl,CLM_Soil_Moisture_Ensemble_Mat,CLM_Soil_Moisture_parm_infl
        del CLM_Soil_Temperature_Ensemble_Mat,CLM_Vegetation_Temperature_Ensemble_Mat,CLM_Ground_Temperature_Ensemble_Mat
        del CLM_Soil_Temperature_parm_infl,CLM_Vegetation_Temperature_parm_infl,CLM_Ground_Temperature_parm_infl
        del CLM_Surface_Temperature_parm_infl,Innovation_State,Increments_State,Observation
            
        NC_File_Out_Assimilation_2_Initial.sync()
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Diagnostic.sync()
        NC_File_Out_Assimilation_2_Diagnostic.close()
                    
    #--------------------------------Finish Assimilation
    
    NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
    Observation_Matrix_Copy = Observation_Matrix[Observation_Matrix_Index,::]
    Observation_Matrix_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, Observation_Matrix_Copy)
    Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid[Prop_Grid_Array_Sys_Index,::])
    Analysis_Grid_Temp = numpy.ma.masked_where(Analysis_Grid_Temp == NAvalue, Analysis_Grid_Temp)
    Model_State = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], numpy.mean(NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0))
    ObsModel_Mat_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, ObsModel_Mat)
    NC_File_Out_Assimilation_2_Initial_Copy.close()
    
    if Def_Print:
        print "numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)",numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)
        print "Min Observation Value is:", Observation_Matrix_Copy.min(), "Maximum Observation Value is:", Observation_Matrix_Copy.max()
        print "Min Model_State Value is:", Model_State.min(), "Maximum Model_State Value is:", Model_State.max()
        print "Min Analysis_Grid Value is:", Analysis_Grid_Temp.min(), "Maximum Analysis_Grid Value is:", Analysis_Grid_Temp.max()
    print "Analysis Mean is:", numpy.mean(Analysis_Grid_Temp), "Model Ensemble Mean is:", numpy.mean(Model_State), "(Analysis - Model_State) Mean is:", numpy.mean(Analysis_Grid_Temp.flatten() - Model_State.flatten())
    print "ObsModel_Mat Mean is:", numpy.mean(ObsModel_Mat_Copy),"Observation Mean is:", numpy.mean(Observation_Matrix_Copy), "(ObsModel_Mat - Observation) Mean is:", numpy.mean(ObsModel_Mat.flatten() - Observation_Matrix_Copy.flatten())
        
    if Def_Print != 0:
        
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:,:,:,:]
        CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:]
        CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:]
        #CLM_Soil_Ice_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_Ensemble_Mat'][:,:,:,:]
        CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:]
        CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:]
        #CLM_2m_Air_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:]
        #CLM_Snow_Depth_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:]
        #CLM_Snow_Water_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:]
        #CLM_ROOTFR_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,:]
        
        CLM_Soil_Moisture_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_parm_infl'][:,:,:]
        CLM_Soil_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_parm_infl'][:,:,:]
        #CLM_Soil_Ice_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_parm_infl'][:,:,:]
        CLM_Vegetation_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_parm_infl'][:,:]
        CLM_Ground_Temperature_parm_infl = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_parm_infl'][:,:]
        
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        
        NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
        CLM_Soil_Moisture_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:]
        CLM_Soil_Temperature_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:]
        #CLM_Soil_Ice_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Ice_Ensemble_Mat'][:,:,:,:]
        CLM_Vegetation_Temperature_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:]
        CLM_Ground_Temperature_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:]
        #CLM_2m_Air_Temperature_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:]
        #CLM_Snow_Depth_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:]
        #CLM_Snow_Water_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:]
        #CLM_ROOTFR_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,:]
        
        CLM_Soil_Moisture_parm_infl_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Moisture_parm_infl'][:,:,:]
        CLM_Soil_Temperature_parm_infl_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Temperature_parm_infl'][:,:,:]
        #CLM_Soil_Ice_parm_infl_Copy = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_parm_infl'][:,:,:]
        CLM_Vegetation_Temperature_parm_infl_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Vegetation_Temperature_parm_infl'][:,:]
        CLM_Ground_Temperature_parm_infl_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Ground_Temperature_parm_infl'][:,:]
        
        NC_File_Out_Assimilation_2_Initial_Copy.close()
        
        print "******************************************************** Station Statistics"
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
        for Station_Index in range(numpy.size(Station_XY)/2):
            print "Station_"+str(Station_Index+1),"Analysis:",Analysis_Grid_Temp[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Model Value:",Model_State[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Observation_Value:",Observation_Matrix_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
            print "ObsModel_Variance:",ObsModel_Variance_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",ObsModel_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
            
            Prop_Grid_Array_Sys_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, :, :, :]
            Prop_Grid_Array_H_Trans = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, :, :, :]
            
            for Ens_Index in range(Ensemble_Number):
                SysModel_Mat_Ens = Prop_Grid_Array_Sys_Copy[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                ObsModel_Mat_Ens = Prop_Grid_Array_H_Trans[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                #ObsModel_Mat = numpy.ma.masked_where(ObsModel_Mat == 0, ObsModel_Mat)
                print "Ens_Index",Ens_Index,"SysModel:",SysModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",ObsModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                        "Analysis:",Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                                
            del Prop_Grid_Array_Sys_Copy,Prop_Grid_Array_H_Trans
            
                
            if Def_Print >= 2:
                if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature":
                    for Soil_Layer_Index in range(Soil_Layer_Num):
                        for Ens_Index in range(Ensemble_Number):
                            print "Soil_Layer_Index",Soil_Layer_Index,"Ens_Index",Ens_Index,"SysModel:",CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],Ens_Index],\
                                "Analysis:",CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],Ens_Index]
                            if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                                print "Soil_Layer_Index",Soil_Layer_Index,"Ens_Index",Ens_Index,"SysModel:",CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],Ens_Index],\
                                "Analysis:",CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],Ens_Index]
                        
            if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Soil_Moisture"):
                for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                    print "Soil_Layer_Index",Soil_Layer_Index,"Model",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis-Model",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
                    if msw_infl < 0.0:
                        print "Soil_Layer_Index",Soil_Layer_Index,"Model_parm_infl",numpy.mean(CLM_Soil_Moisture_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis_parm_infl",numpy.mean(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis-Model",numpy.mean(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Soil_Moisture_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
                
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    print "Vegetation_Temperature","Model",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis-Model",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
                    
                    print "Ground_Temperature","Model",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis-Model",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Ground_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
        
                    for Soil_Layer_Index in range(Soil_Layer_Num):
                        print "Soil_Layer_Index",Soil_Layer_Index,"Model",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis-Model",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
                    
                    if msw_infl < 0.0:
                        print "Vegetation_Temperature","Model_parm_infl",numpy.mean(CLM_Vegetation_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis_parm_infl",numpy.mean(CLM_Vegetation_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis-Model",numpy.mean(CLM_Vegetation_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Vegetation_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
        
                        print "Ground_Temperature","Model_parm_infl",numpy.mean(CLM_Ground_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis_parm_infl",numpy.mean(CLM_Ground_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis-Model",numpy.mean(CLM_Ground_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Ground_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
            
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            print "Soil_Layer_Index",Soil_Layer_Index,"Model_parm_infl",numpy.mean(CLM_Soil_Temperature_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis_parm_infl",numpy.mean(CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis-Model",numpy.mean(CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Soil_Temperature_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])

            
            elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature":
                print "Vegetation_Temperature","Model",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis-Model",numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Vegetation_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
    
                print "Ground_Temperature","Model",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis-Model",numpy.mean(CLM_Ground_Temperature_Ensemble_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Ground_Temperature_Ensemble_Mat_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
    
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    print "Soil_Layer_Index",Soil_Layer_Index,"Model",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                          "Analysis-Model",numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
                
                if msw_infl < 0.0:
                    print "Vegetation_Temperature","Model_parm_infl",numpy.mean(CLM_Vegetation_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                          "Analysis_parm_infl",numpy.mean(CLM_Vegetation_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                          "Analysis-Model",numpy.mean(CLM_Vegetation_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Vegetation_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
    
                    print "Ground_Temperature","Model_parm_infl",numpy.mean(CLM_Ground_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis_parm_infl",numpy.mean(CLM_Ground_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis-Model",numpy.mean(CLM_Ground_Temperature_parm_infl[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Ground_Temperature_parm_infl_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
        
                    for Soil_Layer_Index in range(Soil_Layer_Num):
                        print "Soil_Layer_Index",Soil_Layer_Index,"Model_parm_infl",numpy.mean(CLM_Soil_Temperature_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis_parm_infl",numpy.mean(CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                              "Analysis-Model",numpy.mean(CLM_Soil_Temperature_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Soil_Temperature_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])

                
                if Feedback_Assim: # and (string.atoi(Stop_Month) >= 4) and (string.atoi(Stop_Month) <= 10):
                    Prop_Grid_Array_Sys_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, :, :, :]
                    Prop_Grid_Array_H_Trans = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, :, :, :]
            
                    for Ens_Index in range(Ensemble_Number):
                        SysModel_Mat_Ens = Prop_Grid_Array_Sys_Copy[Ens_Index, Variable_List.index("Soil_Moisture"), :, :]
                        ObsModel_Mat_Ens = Prop_Grid_Array_H_Trans[Ens_Index, Variable_List.index("Surface_Temperature"), :, :]
                        #ObsModel_Mat = numpy.ma.masked_where(ObsModel_Mat == 0, ObsModel_Mat)
                        print "Ens_Index",Ens_Index,"SysModel:",SysModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",ObsModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                                "Analysis:",Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                    
                    del Prop_Grid_Array_Sys_Copy,Prop_Grid_Array_H_Trans
                    
                    for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                        print "Soil_Layer_Index",Soil_Layer_Index,"Model",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:]),\
                              "Analysis-Model",numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])-numpy.mean(CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0],:])
                        if msw_infl < 0.0:
                            print "Soil_Layer_Index",Soil_Layer_Index,"Model_parm_infl",numpy.mean(CLM_Soil_Moisture_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis_parm_infl",numpy.mean(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]),\
                                  "Analysis-Model",numpy.mean(CLM_Soil_Moisture_parm_infl[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])-numpy.mean(CLM_Soil_Moisture_parm_infl_Copy[Soil_Layer_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]])
        
        print "******************************************************** Station Statistics"
        Prop_Grid_Array_Sys_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, :, :, :]
        
        for Ens_Index in range(Ensemble_Number):
            Model_State = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Prop_Grid_Array_Sys_Copy[Ens_Index, Prop_Grid_Array_Sys_Index, :, :])
            Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,::])
            Analysis_Grid_Temp = numpy.ma.masked_where(Analysis_Grid_Temp == NAvalue, Analysis_Grid_Temp)
            print "Min Model:", Model_State.min(), "Max Model:", Model_State.max(), "Min Analysis:", Analysis_Grid_Temp.min(), "Max Analysis:", Analysis_Grid_Temp.max()
            
            #Analysis_Grid[numpy.where(Analysis_Grid[Prop_Grid_Array_Sys_Index,::] == NAvalue)] = CLM_NA
            numexpr_a = Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,::]
            numexpr_b = NAvalue
            numexpr_c = numpy.where(numexpr_a == numexpr_b)
            NA_Index_Analysis_Grid = numexpr_c
            Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,::][NA_Index_Analysis_Grid] = numpy.mean(Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,::][NA_Index_Analysis_Grid])
        
        print "Finish the Analysis of", Analysis_Variable_Name[Prop_Grid_Array_Sys_Index]
        print "Dim NA_Value of Analysis", numpy.size(numpy.where(Analysis_Grid[Prop_Grid_Array_Sys_Index,::] == NAvalue))
         
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Initial_Copy.close()
        del Prop_Grid_Array_Sys_Copy
        
    #os.abort()
    numexpr_a = Analysis_Grid[Prop_Grid_Array_Sys_Index,::]
    numexpr_b = NAvalue
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    NA_Index_Analysis_Grid = numexpr_c      
    
    OutputDate=Stop_Year+ Stop_Month+Stop_Day
    
    if Write_DA_File_Flag:
        if not os.path.exists(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Localization_Map_Mask"):
            os.makedirs(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Localization_Map_Mask")
        Localization_Map_Mask_File_Name = DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Localization_Map_Mask/Localization_Map_Mask_"+SensorVariable_Sub+"_"+OutputDate+".txt"
        numpy.savetxt(Localization_Map_Mask_File_Name,Localization_Map_Mask[Prop_Grid_Array_Sys_Index,::])
    
    if Def_Write_Initial:
        if Def_PP and Ensemble_Number > 1:
            print "********************************************** Using PP to Accelerate Write_Initial_File"
            if PDAF_Assim_Framework == 2:    # Restart PP sever after PDAF MPI
                job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = DAS_Driver_Common.Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
                while len(job_server_node_array) < 1:
                    job_server_node_array = DAS_Driver_Common.Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                    job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = DAS_Driver_Common.Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
            
            Job_Num_Per_Node = int(numpy.ceil(float(Ensemble_Number) / len(active_nodes_server)))
            print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"                    
            if Job_Num_Per_Node == 0:
                Job_Num_Per_Node = 1
            job_server_node_results = []
            Ens_Index = 0
            
            Job_Num_Per_Node_Index = 0
            while Job_Num_Per_Node_Index < Job_Num_Per_Node and Ens_Index < Ensemble_Number:
                for Node_Index in range(len(active_nodes_server)):
                    job_server_node = job_server_node_array[numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]
                    #job_server_node = job_server_node_array[Node_Index]
                    if Ens_Index > Ensemble_Number - 1:
                        break
                    
                    job_server_node_results.append(job_server_node.submit(Write_Initial_File, args=(Ens_Index, Model_Driver, Def_PP, DasPy_Path, Run_Dir_Array, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Row_Numbers, Col_Numbers, finidat_name, SensorVariable_Sub, Variable_ID_Sub, CLM_NA, Feedback_Assim, Stop_Month, Stop_Hour, UTC_Zone, \
                                                   pft_len,maxpft, Def_Region, DAS_Data_Path, Region_Name, Crop_Sum, SensorType_Sub, Row_Numbers_String, Col_Numbers_String, \
                                                   Snow_Layer_Num, column_len, Mask_Index_Sub, Def_Print, Def_ParFor, dtime, irrig_nsteps_per_day, PFT_Num, NAvalue, Density_of_liquid_water, Freezing_temperature_of_fresh_water, Density_of_ice, \
                                                   DAS_Depends_Path, omp_get_num_procs_ParFor, Soil_Layer_Index_DA, Variable_Assimilation_Flag, Variable_List, \
                                                   NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy, fpftcon_name, finidat_name_string),
                                      depfuncs=(Check_Outliers,),
                                      modules=("numpy", "netCDF4", "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "gc", "imp", "getpass", "calendar","scipy.stats", 'scipy.weave'), group='Write_Initial_File'))
                    
                    Ens_Index = Ens_Index + 1         
                Job_Num_Per_Node_Index = Job_Num_Per_Node_Index + 1                   
            
            for job_server_node in job_server_node_array:
                job_server_node.wait()
                if Def_Print >= 2:
                    job_server_node.print_stats() 
            
            if Def_Print:
                if len(job_server_node_results) > 0:
                    for job in job_server_node_results:
                        job_index = job_server_node_results.index(job)
                        if job_index > (Ensemble_Number - 1):
                            break
                        print "Results of ",job_index,"is", job()
                        
        else:   # ********* Run Read_History_File Sequentially
            for Ens_Index in range(Ensemble_Number):
                Write_Initial_File(Ens_Index, Model_Driver, Def_PP, DasPy_Path, Run_Dir_Array, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Row_Numbers, Col_Numbers, finidat_name, SensorVariable_Sub, Variable_ID_Sub, CLM_NA, Feedback_Assim, Stop_Month, Stop_Hour, UTC_Zone, \
                   pft_len, maxpft, Def_Region, DAS_Data_Path, Region_Name, Crop_Sum, SensorType_Sub, Row_Numbers_String, Col_Numbers_String, \
                   Snow_Layer_Num, column_len, Mask_Index_Sub, Def_Print, Def_ParFor, dtime, irrig_nsteps_per_day, PFT_Num, NAvalue, Density_of_liquid_water, Freezing_temperature_of_fresh_water, Density_of_ice, \
                   DAS_Depends_Path, omp_get_num_procs_ParFor, Soil_Layer_Index_DA, Variable_Assimilation_Flag, Variable_List, 
                   NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy, fpftcon_name, finidat_name_string)
                    
    gc.collect()
    del gc.garbage[:]
    
        
    return Analysis_Grid, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, job_server_node_array, active_nodes_server

def Parameter_Update(mpi4py_comm, mpi4py_rank, mpi4py_name, gelmna_threshold, Optimized_Parameter_Index, Model_Driver, NSLOTS,Def_PP, Def_First_Run, Def_Print, Feedback_Assim, Def_Par_Optimized, Parameter_Optimization, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, 
                     Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, SensorQuantity_Sub, SensorType_Sub, SensorVariable_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub, Variable_List, maxpft, \
                              Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict, Dim_Obs_Type, Observation_Matrix, Observation_Longitude, Observation_Latitude, job_server_node_array, active_nodes_server,  ntasks_CLM,
                              Mask, Mask_Index, NAvalue, COSMOS_Circle_Array, COSMOS_Circle_Index_Array, COSMOS_Circle_Num_Array, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Par_Index_Increment_Soil_Par, DasPy_Path, \
                              Variable_Assimilation_Flag, DAS_Depends_Path, Def_ParFor, omp_get_num_procs_ParFor, Def_CDF_Matching, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, PP_Servers_Per_Node, Def_CESM_Multi_Instance, PP_Port, \
                             Plot_Analysis, Soil_Layer_Index_DA, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, Post_Inflation_Alpha_Par, \
                              Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Datetime_Start, Datetime_Initial, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, Write_DA_File_Flag, 
                              r, Observation_Box, Def_Region, Dim_CLM_State, Num_Local_Obs, Model_Variance, DateString_Plot,
                            Def_Multiresolution, Def_ReBEL, Def_Localization, Assim_Algorithm_Name, eps, msw_infl, Region_Name, Call_Gstat_Flag, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Proj_String, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution,
                            dtime, Irrigation_Hours, column_len, Weather_Forecast_Days, Datetime_End, Hydraulic_File_Name, fpftcon_name, Run_Dir_Array, 
                            Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Dim_Observation_Quantity, 
                            Snow_Layer_Num, Def_Write_Initial, cols1d_ixy, cols1d_jxy, cols1d_ityplun, pfts1d_ityplun, Freezing_temperature_of_fresh_water, Density_of_ice, N0, nlyr,
                            Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,
                            Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                            diskless_flag, persist_flag, Irrig_Scheduling, Run_Dir_Home, Start_Month, Stop_Year, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, finidat_name, Density_of_liquid_water, Irrigation_Grid_Flag_Array,
                            mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Datetime_Stop, Datetime_Stop_Init, CLM_NA,
                            Observation_Variance, Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Corelation_Par, octave, Station_XY, Station_XY_Index, Soil_Layer_Num, Analysis_Variable_Name,
                            Analysis_Grid, Localization_Map_Mask, ObsModel_Mat, ObsModel_Variance_Mat, Mask_Sub, Mask_Index_Sub, Mask_Index_Vector, Observation_Matrix_Index, Prop_Grid_Array_Sys_Index, Model_State,
                            SensorQuantity_Index,  E0_ObsModel_Mask, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Soil_Par_Accum_Dim, Veg_Par_Accum_Dim, PFT_Par_Accum_Dim, Hard_Par_Accum_Dim, ParFlow_Layer_Num,
                            Forcing_File_Path, Observation_Path, DAS_Data_Path, Grid_Resolution_CEA, Grid_Resolution_GEO, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                            NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy, NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy,
                            NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy, NC_FileName_Assimilation_2_Parameter_Obs_Dim, 
                            NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy, NC_FileName_Parameter_Space_Single, DAS_Output_Path, \
                            COSMIC_Py, window, memory_profiler, COSMIC, Observation_Time_File_Path):
    
    
    NC_File_Parameter_Space_Single = netCDF4.Dataset(NC_FileName_Parameter_Space_Single,'r')
    
    Bias_Estimation_Option_Model_Par = numpy.asarray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])   # Model Bias
    Bias_Estimation_Option_Obs_Par = numpy.asarray([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])   # Observation Bias
    
    print "Optimized_Parameter_Index",Optimized_Parameter_Index
    
    if Def_Print >= 4:
        NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
        print NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][:,:,:], NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:]
        NC_File_Out_Assimilation_2_Parameter.close()
    
    # Split the Block into SubBlocks to do data assimilation   
    
    Soil_Par_Sens_Dim = numpy.size(numpy.where(Soil_Par_Sens == True))
    Veg_Par_Sens_Dim = 0
    PFT_Par_Sens_Dim = 0
    Hard_Par_Sens_Dim = 0
    
    if Soil_Par_Sens_Dim >= 1:
        Soil_Par_Accum_Dim = Soil_Par_Accum_Dim + 1
        Optimized_Parameter_Index[0] = Optimized_Parameter_Index[0] + 1
        print "**********************************************************************Optimize Soil Parameter"
        if Parameter_Optimization == 2:
            print "############################## Parameter Estimation using Augmentation"
            
            Parameter_Optimization_Flag = 1
                            
            if Def_PP and (not PDAF_Assim_Framework == 2) and (Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) > 1 and len(active_nodes_server) > 1:
                print "********************************************** Using PP to Accelerate Block_Assim"
                Job_Num_Per_Node = int(numpy.ceil(float(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) / len(active_nodes_server)))
                print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"
                if Job_Num_Per_Node == 0:
                    Job_Num_Per_Node = 1
                job_server_node_results = []
                job_server_node_results_wise = [[] for i in range(len(active_nodes_server))]
                
                # The following submits 1 job to 1 node and then retrieves the results
                print "+++++++++++++++++ The following submits",Job_Num_Per_Node,"jobs to 1 node and then retrieves the results"
                Block_Index = 0
                
                Node_Status = numpy.zeros(len(active_nodes_server),dtype=numpy.bool)
                Node_Status[:] = True
                
                while Block_Index < Sub_Block_Ratio_Row*Sub_Block_Ratio_Col:
                    if numpy.size(numpy.where(Node_Status==True)) > 0:
                        Node_Index = numpy.where(Node_Status==True)[0][0]
                        print "***********************Node_Index",Node_Index,"Block_Index",Block_Index,"is submitted!"                         
                        job_server_node = job_server_node_array[numpy.min([Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]                      
                        
                        job_server_node_results.append(job_server_node.submit(Block_Assim, args=(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                                Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                                Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                                Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                                Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                                Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                                Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_Par, Def_ParFor, Ensemble_Number_Predict,
                                Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                                Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                                Observation_Corelation_Par, Bias_Estimation_Option_Model_Par, Bias_Estimation_Option_Obs_Par, Low_Ratio_Par, High_Ratio_Par,
                                Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                                DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                                NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path),
                              depfuncs=(CLM_Assim_Common, Check_Outliers, ParFor_PFT, ParFor_PFT_Block_Assim, ParFor_Fusion, ParFor_H_Operator, ParFor_Texture_Check, ParFor_Check_Outliers, ParFor_Check_Outliers_NA,),
                              modules=("numpy", "netCDF4", "sys", "os", "re", "gc", "imp", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "getpass", "calendar", "glob","scipy.stats",'scipy.weave'), group='Block_Assim'))
                        
                        job_server_node_results_wise[Node_Index] = job_server_node_results[Block_Index]
                        
                        Node_Status[Node_Index] = False
                        
                        Block_Index = Block_Index + 1
                    
                    if Block_Index >= len(active_nodes_server):
                        for job in job_server_node_results_wise:
                            if job != [] and job.finished:
                                Node_Index = job_server_node_results_wise.index(job)
                                print "*********************************************************************Node_Index",Node_Index,"is finished!"
                                Node_Status[Node_Index] = True
                                job_server_node_results_wise[Node_Index] = []
                
                for job_server_node in job_server_node_array:
                    job_server_node.wait()
                    if Def_Print >= 2:
                        job_server_node.print_stats()
                
                if Def_Print:
                    if len(job_server_node_results) > 0:
                        for job in job_server_node_results:
                            job_index = job_server_node_results.index(job)
                            if job_index > (Ensemble_Number - 1):
                                break
                            print "Results of ",job_index,"is", job()
                        
            else:
                print "********* Run Block_Assim Sequentially"
                if PDAF_Assim_Framework == 2:
                    DAS_Driver_Common.Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                
                for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
                    
                    Block_Assim(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                                Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                                Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                                Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                                Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                                Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                                Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_Par, Def_ParFor, Ensemble_Number_Predict,
                                Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                                Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                                Observation_Corelation_Par, Bias_Estimation_Option_Model_Par, Bias_Estimation_Option_Obs_Par, Low_Ratio_Par, High_Ratio_Par,
                                Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                                DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                                NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path, octave, r)
            
            
            print "Write NC_File_Out_Assimilation_2_Initial.nc"
            NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')
            NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
            NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
            Parameter_Soil_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:, :, :, :]
            Parameter_Soil_Space_parm_infl = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][:, :, :]
            Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, :, :, :]
            
            
            for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
                print "Block_Index",Block_Index
                Sub_Block_Row_Start = Sub_Block_Row_Start_Array[Block_Index]
                Sub_Block_Row_End = Sub_Block_Row_End_Array[Block_Index]
                Sub_Block_Col_Start = Sub_Block_Col_Start_Array[Block_Index]
                Sub_Block_Col_End = Sub_Block_Col_End_Array[Block_Index]
                                    
                NC_FileName_Out_Block = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_"+str(Block_Index+1)+".nc"
                
                NC_File_Out_Block = netCDF4.Dataset(NC_FileName_Out_Block, 'r')
                
                Analysis_Grid[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid'][:,:]
                Localization_Map_Mask[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Localization_Map_Mask'][:,:]
                
                Parameter_Soil_Space_Ensemble[:, :, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:]
                Parameter_Soil_Space_parm_infl[:, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Parameter_Soil_Space_parm_infl'][:,:,:]
                Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid_Array'][:,Prop_Grid_Array_Sys_Index,:,:]
                
                NC_File_Out_Block.close()
            
           
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:, :, :, :] = Parameter_Soil_Space_Ensemble
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][:, :, :] = Parameter_Soil_Space_parm_infl
            NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, Prop_Grid_Array_Sys_Index, :, :] = Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, :, :]
            del Parameter_Soil_Space_Ensemble,Parameter_Soil_Space_parm_infl,Analysis_Grid_Array
            
            NC_File_Out_Assimilation_2_Initial.sync()
            NC_File_Out_Assimilation_2_Initial.close()
            
            NC_File_Out_Assimilation_2_Parameter.sync()
            NC_File_Out_Assimilation_2_Parameter.close()
            
            NC_File_Out_Assimilation_2_Diagnostic.sync()
            NC_File_Out_Assimilation_2_Diagnostic.close()
            
            #------------------------------------------Finish Assimilation
            NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
            Observation_Matrix_Copy = numpy.copy(Observation_Matrix[Observation_Matrix_Index,::])
            Observation_Matrix_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, Observation_Matrix_Copy)
            Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid[Prop_Grid_Array_Sys_Index,::])
            Analysis_Grid_Temp = numpy.ma.masked_where(Analysis_Grid_Temp == NAvalue, Analysis_Grid_Temp)
            Model_State = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], numpy.mean(NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0))
            ObsModel_Mat_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, ObsModel_Mat)
            NC_File_Out_Assimilation_2_Initial_Copy.close()
            
            if Def_Print:
                print "numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)",numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)
                print "Min Observation Value is:", Observation_Matrix_Copy.min(), "Maximum Observation Value is:", Observation_Matrix_Copy.max()
                print "Min Model_State Value is:", Model_State.min(), "Maximum Model_State Value is:", Model_State.max()
                print "Min Analysis_Grid Value is:", Analysis_Grid_Temp.min(), "Maximum Analysis_Grid Value is:", Analysis_Grid_Temp.max()
            print "Analysis Mean is:", numpy.mean(Analysis_Grid_Temp), "Model Ensemble Mean is:", numpy.mean(Model_State), "(Analysis - Model_State) Mean is:", numpy.mean(Analysis_Grid_Temp.flatten() - Model_State.flatten())
            print "ObsModel_Mat Mean is:", numpy.mean(ObsModel_Mat_Copy),"Observation Mean is:", numpy.mean(Observation_Matrix_Copy), "(ObsModel_Mat - Observation) Mean is:", numpy.mean(ObsModel_Mat.flatten() - Observation_Matrix_Copy.flatten())
        
            if Def_Print:
                print "******************************************************** Station Statistics"
                for Station_Index in range(numpy.size(Station_XY)/2):
                    print "Station_"+str(Station_Index+1),"Analysis:",Analysis_Grid_Temp[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Model Value:",Model_State[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Observation_Value:",Observation_Matrix_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                    print "ObsModel_Variance:",ObsModel_Variance_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",ObsModel_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
            
            del Observation_Matrix_Copy,Analysis_Grid_Temp,Model_State,ObsModel_Mat_Copy
        
        else:
            print "############################## Wrong Parameter_Optimization Value, Should be 1 or 2"
            os.abort()
        
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
        NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Obs_Dim, 'r+')
        NC_File_Out_Assimilation_2_Parameter_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Copy, 'r')
            
        if Def_Print:
            Parameter_Soil_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:, :, :,:]
            Parameter_Soil_Space_Ensemble_Copy = NC_File_Out_Assimilation_2_Parameter_Copy.variables['Parameter_Soil_Space_Ensemble'][:, :, :,:]
                
            print "******************************************************** Station Statistics"
            for Station_Index in range(numpy.size(Station_XY)/2):
                print "-------------------------------------------Results of Station_"+str(Station_Index+1)
                Par_Index_Sub = 0
                for Par_Index in range(Dim_Soil_Par):
                    if Soil_Par_Sens[Par_Index]:
                        print "Soil_Par_Sens[Par_Index]",Par_Index
                        print "Parameter_Soil_Space_Ensemble_Copy:",numpy.mean(Parameter_Soil_Space_Ensemble_Copy[:, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0), \
                        "Parameter_Soil_Space_Ensemble:",numpy.mean(Parameter_Soil_Space_Ensemble[:, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0)
                        if msw_infl < 0.0:
                            print "Parameter_Soil_Space_parm_infl_Copy:",numpy.mean(NC_File_Out_Assimilation_2_Parameter_Copy.variables['Parameter_Soil_Space_parm_infl'][Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0), \
                            "Parameter_Soil_Space_parm_infl:",numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0)
                        
                        Prop_Grid_Array_Sys_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, :, :, :]
                        Prop_Grid_Array_H_Trans = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, :, :, :]
                        Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:,:,:,:]
                        
                        for Ens_Index in range(Ensemble_Number):
                            SysModel_Mat_Ens = Prop_Grid_Array_Sys_Copy[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                            ObsModel_Mat_Ens = Prop_Grid_Array_H_Trans[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                            print "Ens_Index",Ens_Index,"Parameter_Soil_Space_Ensemble_Copy:",Parameter_Soil_Space_Ensemble_Copy[Ens_Index, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                            "Parameter_Soil_Space_Ensemble:",Parameter_Soil_Space_Ensemble[Ens_Index, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                                "SysModel:",SysModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                                "ObsModel:",ObsModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Analysis:",Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                            del SysModel_Mat_Ens,ObsModel_Mat_Ens
                            
                        del Prop_Grid_Array_Sys_Copy,Prop_Grid_Array_H_Trans,Analysis_Grid_Array
                        
                        Par_Index_Sub += 1
            
            del Parameter_Soil_Space_Ensemble,Parameter_Soil_Space_Ensemble_Copy
        for Dim_Soil_Par_Index in range(Dim_Soil_Par):
            #print numpy.shape(Parameter_Soil_Space_Ensemble_Obs_Dim[:,Dim_Soil_Par_Index,:,:]),numpy.shape(Parameter_Soil_Space_Ensemble_Temp_Copy[:,Dim_Soil_Par_Index,:,:])
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Obs_Dim'][:,Dim_Soil_Par_Index,:,:] += NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par_Index,:,:]
    
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Initial_Copy.close()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        NC_File_Out_Assimilation_2_Parameter.close()
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.close()
        NC_File_Out_Assimilation_2_Parameter_Copy.close()
        
    Soil_Par_Sens_Dim = 0
    Veg_Par_Sens_Dim = 0
    PFT_Par_Sens_Dim = numpy.size(numpy.where(PFT_Par_Sens == True))
    Hard_Par_Sens_Dim = 0
    
    if PFT_Par_Sens_Dim >= 1:
        PFT_Par_Accum_Dim = PFT_Par_Accum_Dim + 1
        Optimized_Parameter_Index[2] = Optimized_Parameter_Index[2] + 1
        print "**********************************************************************Optimize PFT Parameter"
        if Parameter_Optimization == 2:
            print "############################## Parameter Estimation using Augmentation"
    
            Parameter_Optimization_Flag = 1
            
            if Def_PP and (not PDAF_Assim_Framework == 2) and (Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) > 1 and len(active_nodes_server) > 1:
                print "********************************************** Using PP to Accelerate Block_Assim"
                Job_Num_Per_Node = int(numpy.ceil(float(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col) / len(active_nodes_server)))
                print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"
                if Job_Num_Per_Node == 0:
                    Job_Num_Per_Node = 1
                job_server_node_results = []
                job_server_node_results_wise = [[] for i in range(len(active_nodes_server))]
                
                # The following submits 1 job to 1 node and then retrieves the results
                print "+++++++++++++++++ The following submits",Job_Num_Per_Node,"jobs to 1 node and then retrieves the results"
                Block_Index = 0
                
                Node_Status = numpy.zeros(len(active_nodes_server),dtype=numpy.bool)
                Node_Status[:] = True
                
                while Block_Index < Sub_Block_Ratio_Row*Sub_Block_Ratio_Col:
                    if numpy.size(numpy.where(Node_Status==True)) > 0:
                        Node_Index = numpy.where(Node_Status==True)[0][0]
                        print "***********************Node_Index",Node_Index,"Block_Index",Block_Index,"is submitted!"
                        job_server_node = job_server_node_array[numpy.min([Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]                      
                        
                        job_server_node_results.append(job_server_node.submit(Block_Assim, args=(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                                Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                                Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                                Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                                Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                                Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                                Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_Par, Def_ParFor, Ensemble_Number_Predict,
                                Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                                Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                                Observation_Corelation_Par, Bias_Estimation_Option_Model_Par, Bias_Estimation_Option_Obs_Par, Low_Ratio_Par, High_Ratio_Par,
                                Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                                DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                                NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path),
                              depfuncs=(CLM_Assim_Common, Check_Outliers, ParFor_PFT, ParFor_PFT_Block_Assim, ParFor_Fusion, ParFor_H_Operator, ParFor_Texture_Check, ParFor_Check_Outliers, ParFor_Check_Outliers_NA),
                              modules=("numpy", "netCDF4", "sys", "os", "re", "gc", "imp", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "getpass", "calendar", "glob","scipy.stats",'scipy.weave'), group='Block_Assim'))
                        
                        job_server_node_results_wise[Node_Index] = job_server_node_results[Block_Index]
                        
                        Node_Status[Node_Index] = False
                        
                        Block_Index = Block_Index + 1
                    
                    if Block_Index >= len(active_nodes_server):
                        for job in job_server_node_results_wise:
                            if job != [] and job.finished:
                                Node_Index = job_server_node_results_wise.index(job)
                                print "*********************************************************************Node_Index",Node_Index,"is finished!"
                                Node_Status[Node_Index] = True
                                job_server_node_results_wise[Node_Index] = []
                
                for job_server_node in job_server_node_array:
                    job_server_node.wait()
                    if Def_Print >= 2:
                        job_server_node.print_stats()
                
                if Def_Print:
                    if len(job_server_node_results) > 0:
                        for job in job_server_node_results:
                            job_index = job_server_node_results.index(job)
                            if job_index > (Ensemble_Number - 1):
                                break
                            print "Results of ",job_index,"is", job()
                        
            else:
                print "********* Run Block_Assim Sequentially"
                if PDAF_Assim_Framework == 2:
                    DAS_Driver_Common.Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                
                for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
                    
                    Block_Assim(Block_Index, Model_Driver, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Numbers, Col_Numbers, Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Row_Offset, Col_Offset,
                                Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                                Start_Month, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, Ensemble_Number, Prop_Grid_Array_Sys_Index, 
                                Dim_Observation_Quantity, SensorQuantity_Index, Observation_Box, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD,
                                Variable_List, Observation_Matrix_Index, Soil_Layer_Num, ParFlow_Layer_Num, SensorVariable_Sub, SensorType_Sub, SensorQuantity_Sub, SensorResolution_Sub, 
                                Variable_Assimilation_Flag, Soil_Layer_Index_DA, Feedback_Assim, Parameter_Optimization_Flag, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Dim_CLM_State, maxpft, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type,
                                Def_First_Run, Def_Print, Def_PP, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs[Prop_Grid_Array_Sys_Index], eps, msw_infl, Post_Inflation_Alpha_Par, Def_ParFor, Ensemble_Number_Predict,
                                Call_Gstat_Flag, diskless_flag, persist_flag, Assim_Algorithm_Name, Proj_String, Z_Resolution, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper,
                                Grid_Resolution_CEA, Write_DA_File_Flag, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial, Region_Name, NSLOTS,
                                Observation_Corelation_Par, Bias_Estimation_Option_Model_Par, Bias_Estimation_Option_Obs_Par, Low_Ratio_Par, High_Ratio_Par,
                                Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, Par_Veg_Uniform_STD_Sub, Par_PFT_Uniform_STD_Sub, Par_Hard_Uniform_STD_Sub, DateString_Plot,
                                DAS_Depends_Path, DasPy_Path, CLM_NA, NAvalue, omp_get_num_procs_ParFor, Def_CDF_Matching, Plot_Analysis, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, 
                                NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, DAS_Output_Path, octave, r)
            
            print "Write NC_File_Out_Assimilation_2_Initial.nc"
            NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')
            NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
            NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
            
            Parameter_PFT_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:]
            Parameter_PFT_Space_parm_infl = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][:,:,:]
            Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, :, :, :]
            
            for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
                print "Block_Index",Block_Index
                Sub_Block_Row_Start = Sub_Block_Row_Start_Array[Block_Index]
                Sub_Block_Row_End = Sub_Block_Row_End_Array[Block_Index]
                Sub_Block_Col_Start = Sub_Block_Col_Start_Array[Block_Index]
                Sub_Block_Col_End = Sub_Block_Col_End_Array[Block_Index]
                                
                NC_FileName_Out_Block = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_"+str(Block_Index+1)+".nc"
                
                NC_File_Out_Block = netCDF4.Dataset(NC_FileName_Out_Block, 'r')
                
                Analysis_Grid[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid'][:,:]
                Localization_Map_Mask[Prop_Grid_Array_Sys_Index,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Localization_Map_Mask'][:,:]
                
                Parameter_PFT_Space_Ensemble[:,:,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:]
                Parameter_PFT_Space_parm_infl[:,Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Parameter_PFT_Space_parm_infl'][:,:,:]
                Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, Sub_Block_Row_Start:Sub_Block_Row_End,Sub_Block_Col_Start:Sub_Block_Col_End] = NC_File_Out_Block.variables['Analysis_Grid_Array'][:,Prop_Grid_Array_Sys_Index,:,:]
                
                NC_File_Out_Block.close()
            
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:] = Parameter_PFT_Space_Ensemble
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][:,:,:] = Parameter_PFT_Space_parm_infl
            NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:, Prop_Grid_Array_Sys_Index, :, :] = Analysis_Grid_Array[:, Prop_Grid_Array_Sys_Index, :, :]
            del Parameter_PFT_Space_Ensemble,Parameter_PFT_Space_parm_infl,Analysis_Grid_Array
            
            NC_File_Out_Assimilation_2_Initial.sync()
            NC_File_Out_Assimilation_2_Initial.close()
            NC_File_Out_Assimilation_2_Diagnostic.sync()
            NC_File_Out_Assimilation_2_Diagnostic.close()
            NC_File_Out_Assimilation_2_Parameter.sync()
            NC_File_Out_Assimilation_2_Parameter.close()
    
            #------------------------------------------Finish Assimilation
            NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
            Observation_Matrix_Copy = numpy.copy(Observation_Matrix[Observation_Matrix_Index,::])
            Observation_Matrix_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, Observation_Matrix_Copy)
            Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid[Prop_Grid_Array_Sys_Index,::])
            Analysis_Grid_Temp = numpy.ma.masked_where(Analysis_Grid_Temp == NAvalue, Analysis_Grid_Temp)
            Model_State = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], numpy.mean(NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0))
            ObsModel_Mat_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, ObsModel_Mat)
            NC_File_Out_Assimilation_2_Initial_Copy.close()
            
            if Def_Print:
                print "numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)",numpy.shape(Analysis_Grid_Temp),numpy.shape(Model_State)
                print "Min Observation Value is:", Observation_Matrix_Copy.min(), "Maximum Observation Value is:", Observation_Matrix_Copy.max()
                print "Min Model_State Value is:", Model_State.min(), "Maximum Model_State Value is:", Model_State.max()
                print "Min Analysis_Grid Value is:", Analysis_Grid_Temp.min(), "Maximum Analysis_Grid Value is:", Analysis_Grid_Temp.max()
            print "Analysis Mean is:", numpy.mean(Analysis_Grid_Temp), "Model Ensemble Mean is:", numpy.mean(Model_State), "(Analysis - Model_State) Mean is:", numpy.mean(Analysis_Grid_Temp.flatten() - Model_State.flatten())
            print "ObsModel_Mat Mean is:", numpy.mean(ObsModel_Mat_Copy),"Observation Mean is:", numpy.mean(Observation_Matrix_Copy), "(ObsModel_Mat - Observation) Mean is:", numpy.mean(ObsModel_Mat.flatten() - Observation_Matrix_Copy.flatten())
        
            if Def_Print:
                print "******************************************************** Station Statistics"
                for Station_Index in range(numpy.size(Station_XY)/2):
                    print "Station_"+str(Station_Index+1),"Analysis:",Analysis_Grid_Temp[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Model Value:",Model_State[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Observation_Value:",Observation_Matrix_Copy[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                    print "ObsModel_Variance:",ObsModel_Variance_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",ObsModel_Mat[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
            
            del Observation_Matrix_Copy,Analysis_Grid_Temp,Model_State,ObsModel_Mat_Copy
            
        else:
            print "############################## Wrong Parameter_Optimization Value, Should be 1 or 2"
            os.abort()
        
        
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
        NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Obs_Dim, 'r+')
        NC_File_Out_Assimilation_2_Parameter_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Copy, 'r')
        
        if Def_Print:
            
            print "******************************************************** Station Statistics"
            for Station_Index in range(numpy.size(Station_XY)/2):
                print "-------------------------------------------Results of Station_"+str(Station_Index+1)
                Par_Index_Sub = 0
                for Par_Index in range(Dim_PFT_Par):
                    if PFT_Par_Sens[Par_Index]:
                        print "PFT_Par_Sens[Par_Index]",Par_Index
                        print "Parameter_PFT_Space_Ensemble_Copy:",numpy.mean(NC_File_Out_Assimilation_2_Parameter_Copy.variables['Parameter_PFT_Space_Ensemble'][:, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0), \
                        "Parameter_PFT_Space_Ensemble:",numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0)
                        if msw_infl < 0.0:
                            print "Parameter_PFT_Space_parm_infl_Copy:",numpy.mean(NC_File_Out_Assimilation_2_Parameter_Copy.variables['Parameter_PFT_Space_parm_infl'][Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0), \
                            "Parameter_PFT_Space_parm_infl:",numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]], axis=0)
                        
                        Prop_Grid_Array_Sys_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:, :, :, :]
                        Prop_Grid_Array_H_Trans = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, :, :, :]
                        Parameter_PFT_Space_Ensemble_Copy = NC_File_Out_Assimilation_2_Parameter_Copy.variables['Parameter_PFT_Space_Ensemble'][:, :, :, :]
                        Parameter_PFT_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:, :, :,:]
                        Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:,:,:,:]
                    
                        for Ens_Index in range(Ensemble_Number):
                            SysModel_Mat_Ens = Prop_Grid_Array_Sys_Copy[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                            ObsModel_Mat_Ens = Prop_Grid_Array_H_Trans[Ens_Index, Prop_Grid_Array_Sys_Index, :, :]
                            print "Ens_Index",Ens_Index,"Parameter_PFT_Space_Ensemble_Copy:",Parameter_PFT_Space_Ensemble_Copy[Ens_Index, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                            "Parameter_PFT_Space_Ensemble:",Parameter_PFT_Space_Ensemble[Ens_Index, Par_Index, Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                                "SysModel:",SysModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],\
                                "ObsModel:",ObsModel_Mat_Ens[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"Analysis:",Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                            del SysModel_Mat_Ens,ObsModel_Mat_Ens
                        
                        del Prop_Grid_Array_Sys_Copy,Prop_Grid_Array_H_Trans,Parameter_PFT_Space_Ensemble_Copy,Parameter_PFT_Space_Ensemble,Analysis_Grid_Array
                        Par_Index_Sub += 1
                #os.abort()
        
        for Dim_PFT_Par_Index in range(Dim_PFT_Par):
            #print numpy.shape(Parameter_PFT_Space_Ensemble_Obs_Dim[:,Dim_PFT_Par_Index,:,:]),numpy.shape(Parameter_PFT_Space_Ensemble_Temp_Copy[:,Dim_PFT_Par_Index,:,:])
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Obs_Dim'][:,Dim_PFT_Par_Index,:,:] += NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,Dim_PFT_Par_Index,:,:]
        
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Initial_Copy.close()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        NC_File_Out_Assimilation_2_Parameter.close()
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.close()
        NC_File_Out_Assimilation_2_Parameter_Copy.close()
    
    print "Soil_Par_Accum_Dim",Soil_Par_Accum_Dim,"Veg_Par_Accum_Dim",Veg_Par_Accum_Dim,"PFT_Par_Accum_Dim",PFT_Par_Accum_Dim
    print ""
    print "Optimized_Parameter_Index",Optimized_Parameter_Index
    
    numexpr_a = []
    numexpr_b = []
    numexpr_c = []
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    NC_File_Out_Assimilation_2_Parameter_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Copy, 'r')
    NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Obs_Dim, 'r+')
            
    if Soil_Par_Accum_Dim > 0:
        Parameter_Soil_Space_Ensemble_Obs_Dim = NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Obs_Dim'][:,:,:,:]
        
        for Dim_Soil_Par_Index in range(Dim_Soil_Par):
            if Soil_Par_Sens_Array[0][Dim_Soil_Par_Index] or Soil_Par_Sens_Array[1][Dim_Soil_Par_Index]:
                #print numpy.shape(Parameter_Soil_Space_Ensemble_Obs_Dim[:,Dim_Soil_Par_Index,:,:]),numpy.shape(Parameter_Soil_Space_Ensemble_Temp_Copy[:,Dim_Soil_Par_Index,:,:])
                Parameter_Soil_Space_Ensemble_Obs_Dim[:,Dim_Soil_Par_Index,:,:] = Parameter_Soil_Space_Ensemble_Obs_Dim[:,Dim_Soil_Par_Index,:,:] / float(Soil_Par_Accum_Dim)
                if Parameter_Optimization == 1:
                    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Predict_Obs_Dim'][:,Dim_Soil_Par_Index,:,:] = NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Predict_Obs_Dim'][:,Dim_Soil_Par_Index,:,:] / float(Soil_Par_Accum_Dim)                    
            else:
                for Ens_Index in range(Ensemble_Number):
                    Parameter_Soil_Space_Ensemble_Obs_Dim[Ens_Index,Dim_Soil_Par_Index,:,:] = NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][Dim_Soil_Par_Index,:,:]
            if Parameter_Optimization == 1:
                for Ens_Index in range(Ensemble_Number_Predict):
                    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Predict_Obs_Dim'][Ens_Index,Dim_Soil_Par_Index,:,:] = NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][Dim_Soil_Par_Index,:,:]
         
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Obs_Dim'][:,:,:,:] = Parameter_Soil_Space_Ensemble_Obs_Dim
        del Parameter_Soil_Space_Ensemble_Obs_Dim
    
    if PFT_Par_Accum_Dim > 0:
        Parameter_PFT_Space_Ensemble_Obs_Dim = NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Obs_Dim'][:,:,:,:]
        
        for Dim_PFT_Par_Index in range(Dim_PFT_Par):
            if numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array)[:,Dim_PFT_Par_Index] == True)) >= 1:
                Parameter_PFT_Space_Ensemble_Obs_Dim[:,Dim_PFT_Par_Index,:,:] = Parameter_PFT_Space_Ensemble_Obs_Dim[:,Dim_PFT_Par_Index,:,:] / float(PFT_Par_Accum_Dim)
                if Parameter_Optimization == 1:
                    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Predict_Obs_Dim'][Ens_Index,Dim_PFT_Par_Index,:,:] = NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Predict_Obs_Dim'][Ens_Index,Dim_PFT_Par_Index,:,:] / float(PFT_Par_Accum_Dim)
            else:
                for Ens_Index in range(Ensemble_Number):
                    Parameter_PFT_Space_Ensemble_Obs_Dim[Ens_Index,Dim_PFT_Par_Index,:,:] = NC_File_Parameter_Space_Single.variables['Parameter_PFT_Space_Single'][Dim_PFT_Par_Index,:,:]
            if Parameter_Optimization == 1:
                for Ens_Index in range(Ensemble_Number_Predict):
                    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Predict_Obs_Dim'][Ens_Index,Dim_PFT_Par_Index,:,:] = NC_File_Parameter_Space_Single.variables['Parameter_PFT_Space_Single'][Dim_PFT_Par_Index,:,:]
        
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Obs_Dim'][:,:,:,:] = Parameter_PFT_Space_Ensemble_Obs_Dim
        del Parameter_PFT_Space_Ensemble_Obs_Dim
               
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Veg_Space_Ensemble_Obs_Dim'][:,:,:] = 0.0
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Hard_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
    
    if (numpy.size(numpy.where(numpy.asarray(Veg_Par_Sens_Array) == True)) >= 1):
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Veg_Space_Ensemble_Matrix_Obs_Dim'][:,:,:,:] = 0.0
    if (numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) >= 1):
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
    
    if Parameter_Optimization == 1:
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Predict_Obs_Dim'][:,:,:,:] = 0.0
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Veg_Space_Ensemble_Predict_Obs_Dim'][:,:,:] = 0.0
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Predict_Obs_Dim'][:,:,:,:] = 0.0
        NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Hard_Space_Ensemble_Predict_Obs_Dim'][:,:,:,:] = 0.0
        
    NC_File_Out_Assimilation_2_Parameter.sync()
    NC_File_Out_Assimilation_2_Parameter.close()
    NC_File_Out_Assimilation_2_Parameter_Copy.close()
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.sync()
    NC_File_Out_Assimilation_2_Parameter_Obs_Dim.close()
    NC_File_Out_Assimilation_2_Initial.close()
    #Parameter_Soil_Space_Single = numpy.mean(Parameter_Soil_Space_Ensemble_Obs_Dim,axis=0)
    #Parameter_Veg_Space_Single = numpy.mean(Parameter_Veg_Space_Ensemble_Obs_Dim,axis=0)
      
    NC_File_Parameter_Space_Single.close()
    
    del numexpr_a,numexpr_b,numexpr_c
        
    gc.collect()
    del gc.garbage[:]
        
    if PDAF_Assim_Framework == 2:    # Restart PP sever after PDAF MPI
        job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = DAS_Driver_Common.Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
        while len(job_server_node_array) < 1:
            job_server_node_array = DAS_Driver_Common.Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
            job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = DAS_Driver_Common.Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
                    
                    
    return Def_Par_Optimized, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, job_server_node_array, active_nodes_server, Optimized_Parameter_Index


    