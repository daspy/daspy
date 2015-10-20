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
import numpy, string, smtplib, sys, imp, math, multiprocessing, shutil, fnmatch, shlex
from mpi4py import MPI
from concurrent import futures
sys.path.append('Utilities')

from DAS_Observation_Operator import *
from DAS_Assim import *
from DAS_Assim_Common import *
from DAS_Utilities import *

def DAS_Driver_Common(mpi4py_comm, mpi4py_null, mpi4py_rank, mpi4py_name, mpi4py_comm_split, mpipy_comm_decomposition, Model_Driver, PDAF_Assim_Framework, PDAF_Filter_Type, Def_PP, Def_CESM_Multi_Instance, Def_Par_Optimized, Def_Par_Sensitivity, Def_Par_Correlation, Do_DA_Flag, CLM_NA, NAvalue, finidat_initial_CLM, finidat_initial_CLM_Copy, Def_ParFor, Def_Region, 
                    Def_Initial, Irrig_Scheduling, Irrigation_Hours, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                    Start_Year, Start_Month, Start_Day, Start_Hour, Stop_Year, Stop_Month, Stop_Day, Stop_Hour, Datetime_Initial, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End,
                    DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime, N_Steps, Ensemble_Number, Row_Numbers, Col_Numbers, Ensemble_Number_Predict, 
                    Row_Numbers_String, Col_Numbers_String, Mask_Index, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Two_Step_Bias_Estimation_Flag, Two_Step_Bias_Estimation_Active, Mean_Dir,
                    PP_Port, NSLOTS, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, Dim_Observation_Quantity,
                    Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, ParFlow_Layer_Num, Density_of_liquid_water, Freezing_temperature_of_fresh_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, diskless_flag, persist_flag,
                    CLM_File_Name_List, Parameter_Range_Soil, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Observation_Bias_Initialization_Flag,
                    stop_tod_string, history_file_name, finidat_name, Observation_Matrix, Observation_Longitude, Observation_Latitude, Observation_Variance, \
                    COSMOS_Circle_Index_Array, COSMOS_Circle_Num_Array, COSMOS_Circle_Array, Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, \
                    Observation_Y_Lower, Observation_Y_Upper, Observation_Misc, Observation_View_Zenith_Angle, Observation_View_Time, Observation_Corelation_Par, \
                    Initial_Perturbation_ST_Flag, Def_First_Run, Def_SpinUp,  active_nodes_server, job_server_node_array, PP_Servers_Per_Node,
                    Forcing_File_Path_Home, SensorType, SensorVariable, SensorQuantity, SensorResolution, Variable_ID, QC_ID, Observation_File_Name, Dim_Obs_Type,
                    Write_DA_File_Flag, Use_Mask_Flag, Mask_File, Def_ReBEL, Def_Localization, Call_Gstat_Flag, Plot_Analysis, Num_Local_Obs_State,
                    Observation_Path,  Grid_Resolution_CEA, Grid_Resolution_GEO, mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen,
                    LATIXY_Mat, LONGXY_Mat, MODEL_CEA_X, MODEL_CEA_Y, MODEL_X_Left,MODEL_X_Right, MODEL_Y_Lower,MODEL_Y_Upper, LAI_Year_String, Month_String, 
                    DAS_Output_Path, Dim_CLM_State, Parameter_Optimization, Def_Debug, PCT_PFT_High, PCT_PFT_Low, PCT_PFT_WATER, Soil_Layer_Index_DA,
                    ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat, column_len, pft_len, cols1d_jxy, cols1d_ixy, pfts1d_jxy, pfts1d_ixy, 
                    Constant_File_Name_Header, finidat_name_string, Feedback_Assim, numrad, 
                    Density_of_ice, omp_get_num_procs_ParFor, Model_Variance, Additive_Noise_SM_Par, Additive_Noise_SM, Additive_Noise_ST,
                    Variable_List, Variable_Assimilation_Flag, Analysis_Variable_Name, Mask, Mask_X, Mask_Y, N0, nlyr, Station_XY, Station_XY_Index, Def_First_Run_RTM, Soil_Density, CMEM_Work_Path_Array,
                    DAS_File_Name_List, COUP_OAS_PFL, CESM_Init_Flag,
                    MODIS_LAI_Data_ID, Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, PFT_Par_Sens_Array, UTC_Zone, plt, cm, colors, DateString_Plot, octave, r, COSMIC_Py, window, memory_profiler, COSMIC, Observation_Time_File_Path):
    
    NC_FileName_Assimilation_2_Constant = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Constant.nc")[0]
    NC_FileName_Assimilation_2_Diagnostic = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Diagnostic.nc")[0]
    NC_FileName_Assimilation_2_Initial = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Initial.nc")[0]
    NC_FileName_Assimilation_2_Initial_Copy = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Initial_Copy.nc")[0]
    NC_FileName_Assimilation_2_Bias = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Bias.nc")[0]
    NC_FileName_Assimilation_2_Bias_Copy = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Bias_Copy.nc")[0]
    NC_FileName_Assimilation_2_Bias_Monthly = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Bias_Monthly.nc")[0]
    NC_FileName_Assimilation_2_Bias_Monthly_Copy = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Bias_Monthly_Copy.nc")[0]
    NC_FileName_Estimated_Bias = fnmatch.filter(DAS_File_Name_List,"*Estimated_Bias.nc")[0]
    NC_FileName_Assimilation_2_Parameter = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Parameter.nc")[0]
    NC_FileName_Assimilation_2_Parameter_Copy = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Parameter_Copy.nc")[0]
    NC_FileName_Assimilation_2_Parameter_Obs_Dim = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Parameter_Obs_Dim.nc")[0]
    NC_FileName_Assimilation_2_Parameter_Monthly = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Parameter_Monthly.nc")[0]
    NC_FileName_Assimilation_2_Parameter_Monthly_Copy = fnmatch.filter(DAS_File_Name_List,"*Assimilation_2_Parameter_Monthly_Copy.nc")[0]
    NC_FileName_Optimized_Parameter = fnmatch.filter(DAS_File_Name_List,"*Optimized_Parameter.nc")[0]
    NC_FileName_Soil_Moisture_Difference = fnmatch.filter(DAS_File_Name_List,"*Soil_Moisture_Difference.nc")[0]
    NC_FileName_Parameter_Space_Single = fnmatch.filter(DAS_File_Name_List,"*Parameter_Space_Single.nc")[0]
    
    job_tid = 0
    
    if Datetime_Start != Datetime_Stop and (Two_Step_Bias_Estimation_Flag == 0): # If There are severl Observation at the Same Time, CLM Only Need to Be Run Once.
                
        Datetime_Start_Mean = Datetime_Start
        stop_tod_string_final = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
        
        if mpi4py_rank == 0:
                        
            if Initial_Perturbation and ((numpy.sum(Initial_Perturbation_SM_Flag) or numpy.sum(Initial_Perturbation_ST_Flag)) or (Def_First_Run and Def_SpinUp)):
                NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
                NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')
            
                print "**************************** Initial Perturbation"
                
                r.assign("Ensemble_Number",Ensemble_Number)                    
                r.assign("seed_value",string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                
                if numpy.sum(Initial_Perturbation_SM_Flag) or (Def_First_Run and Def_SpinUp):
                    r('set.seed(seed_value)')
                    r('grf_data <- grf(n = Row_Numbers*Col_Numbers, grid = "reg", nx = Col_Numbers, ny = Row_Numbers, xlims = c(1, Col_Numbers), ylims = c(1, Row_Numbers), nsim = Ensemble_Number, cov.model = "exponential", \
                        cov.pars = c(0.0016, 10.0), \
                        kappa = 0.5, nugget = 0.0, mean = 0, RF=TRUE)')
                    #print r('dim(grf_data$data)')
                    #print numpy.shape(numpy.array(r('grf_data$data')))
                    GaussRF_Array = numpy.asarray(r['grf_data$data'])
                    
                    Initial_SM_Noise = numpy.asarray(NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_SM_Noise'][:,:,:])
                    
                    for Ens_Index in range(Ensemble_Number):
                        Initial_SM_Noise[Ens_Index,:,:] = numpy.reshape(GaussRF_Array[:,Ens_Index],(Row_Numbers,Col_Numbers))
                        #print numpy.min(Initial_SM_Noise[Ens_Index,:,:]),numpy.max(Initial_SM_Noise[Ens_Index,:,:])
                        numexpr_a = Initial_SM_Noise[Ens_Index,:,:]
                        numexpr_b = 0.02
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_SM_Noise[Ens_Index,:,:] = numexpr_a
                        numexpr_a = Initial_SM_Noise[Ens_Index,:,:]
                        numexpr_b = -0.02
                        numexpr_c = numpy.where(numexpr_a < numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_SM_Noise[Ens_Index,:,:] = numexpr_a
                        numexpr_a = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][::]
                        numexpr_b = NAvalue
                        numexpr_c = numpy.where(numexpr_a == numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_SM_Noise[Ens_Index,:,:] = numexpr_a
                    
                    NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_SM_Noise'][:,:,:] = Initial_SM_Noise
                    del Initial_SM_Noise,GaussRF_Array
                    
                if numpy.sum(Initial_Perturbation_ST_Flag) or (Def_First_Run and Def_SpinUp):
                    r('set.seed(seed_value)')
                    r('grf_data <- grf(n = Row_Numbers*Col_Numbers, grid = "reg", nx = Col_Numbers, ny = Row_Numbers, xlims = c(1, Col_Numbers), ylims = c(1, Row_Numbers), nsim = Ensemble_Number, cov.model = "exponential", \
                        cov.pars = c(1.0, 5.0), \
                        kappa = 0.5, nugget = 0.0, mean = 0, RF=TRUE)')
                    #print r('dim(grf_data$data)')
                    #print numpy.shape(numpy.asarray(r('grf_data$data')))
                    GaussRF_Array = numpy.asarray(r['grf_data$data'])
                    
                    Initial_ST_Noise = numpy.asarray(NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_ST_Noise'][:,:,:])
                    
                    for Ens_Index in range(Ensemble_Number):
                        Initial_ST_Noise[Ens_Index,:,:] = numpy.reshape(GaussRF_Array[:,Ens_Index],(Row_Numbers,Col_Numbers))
                        #print numpy.min(Initial_ST_Noise[Ens_Index,:,:]),numpy.max(Initial_ST_Noise[Ens_Index,:,:])
                        numexpr_a = Initial_ST_Noise[Ens_Index,:,:]
                        numexpr_b = 1.0
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_ST_Noise[Ens_Index,:,:] = numexpr_a
                        numexpr_a = Initial_ST_Noise[Ens_Index,:,:]
                        numexpr_b = -1.0
                        numexpr_c = numpy.where(numexpr_a < numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_ST_Noise[Ens_Index,:,:] = numexpr_a
                        numexpr_a = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][::]
                        numexpr_b = NAvalue
                        numexpr_c = numpy.where(numexpr_a == numexpr_b)
                        numexpr_a[numexpr_c] = numexpr_b
                        Initial_ST_Noise[Ens_Index,:,:] = numexpr_a
                    
                    NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_ST_Noise'][:,:,:] = Initial_ST_Noise
                    del Initial_ST_Noise,GaussRF_Array
                print ""
                
                NC_File_Out_Assimilation_2_Constant.close()
                NC_File_Out_Assimilation_2_Diagnostic.sync()
                NC_File_Out_Assimilation_2_Diagnostic.close()
        
        if Def_PP == 2:
            mpi4py_comm.barrier()
            mpi4py_comm.Barrier()
            
        if mpi4py_rank == 0:
            print "*************************************************Start Model Simulation***************************************************************"
        
            if Def_PP and Ensemble_Number > 1:
                
                print "********************************************** Using PP to Accelerate Prepare_Model_Operator"
                Job_Num_Per_Node = int(numpy.ceil(float(Ensemble_Number) / len(active_nodes_server)))
                print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"
                if Job_Num_Per_Node == 0:
                    Job_Num_Per_Node = 1
                job_server_node_results = []
                Ens_Index = 0
                
                Job_Num_Per_Node_Index = 0
                while Job_Num_Per_Node_Index < Job_Num_Per_Node and Ens_Index < Ensemble_Number:
                    for Node_Index in range(len(active_nodes_server)):
                        #print "Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server)",Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),"len(job_server_node_array)-1",len(job_server_node_array)-1,numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])
                        job_server_node = job_server_node_array[numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]
                        #print "job_server_node.get_active_nodes(),Job_Num_Per_Node_Index,Node_Index",job_server_node.get_active_nodes(),Job_Num_Per_Node_Index,Node_Index
                        if Ens_Index > Ensemble_Number - 1:
                            break
                        job_server_node_results.append(job_server_node.submit(Prepare_Model_Operator, args=(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                    Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                    CLM_File_Name_List, Parameter_Range_Soil,
                                                  Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, 
                                                  DAS_Data_Path, DasPy_Path, DAS_Output_Path, Forcing_File_Path_Array, dtime, Variable_Assimilation_Flag, Variable_List,\
                                                  Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, 0,\
                                                  omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                  Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                                                  NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single),
                                      depfuncs=(Run_CLM, Call_CLM_3D, Write_datm_atm_in, Write_datm_streams_txt, Write_presaero_stream_txt, Write_lnd_in, Write_rof_in, Write_Config_Files, Write_drv_in, Write_seq_maps, copyLargeFile,),
                                      modules=("numpy", "netCDF4", "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "signal", "gc", "imp", "getpass", "calendar", "glob","scipy.stats","scipy.signal",'scipy.weave'), group='Prepare_Model_Operator'))
                        
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
                
                
            else:
                print "*************************************************** Run Prepare_Model_Operator Sequentially"
                for Ens_Index in range(Ensemble_Number):
                    if Def_Print:
                        print "Ens_Index",Ens_Index
                    Prepare_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                    Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                    CLM_File_Name_List, Parameter_Range_Soil,
                                                  Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, 
                                                  DAS_Data_Path, DasPy_Path, DAS_Output_Path, Forcing_File_Path_Array, dtime, Variable_Assimilation_Flag, Variable_List,\
                                                  Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, 0,\
                                                  omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                  Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                                                  NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial,  NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single)                        
        
        
        if Def_PP == 2:
            mpi4py_comm.barrier()
            mpi4py_comm.Barrier()
            
        ##############################################################################
        if True:
            if Def_PP == 1:
                
                print "********************************************** Using PP to Accelerate Call_Model_Operator"
                Job_Num_Per_Node = int(numpy.ceil(float(Ensemble_Number) / len(active_nodes_server)))
                print "The following submits",Job_Num_Per_Node,"jobs on each node and then retrieves the results"
                if Job_Num_Per_Node == 0:
                    Job_Num_Per_Node = 1
                job_server_node_results = []
                Ens_Index = 0
                
                Job_Num_Per_Node_Index = 0
                while Job_Num_Per_Node_Index < Job_Num_Per_Node and Ens_Index < Ensemble_Number:
                    for Node_Index in range(len(active_nodes_server)):
                        #print "Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server)",Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),"len(job_server_node_array)-1",len(job_server_node_array)-1,numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])
                        job_server_node = job_server_node_array[numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]
                        #job_server_node = job_server_node_array[Node_Index]
                        if Ens_Index > Ensemble_Number - 1:
                            break
                        job_server_node_results.append(job_server_node.submit(Call_Model_Operator, args=(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized,  Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                                                                             Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                                                                             CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                                                                                  Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, 0,\
                                                                                  Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                                                  Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                                                                  NC_FileName_Assimilation_2_Initial,  NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null),
                                      depfuncs=(Run_CLM, Call_CLM_3D, Write_datm_atm_in, Write_datm_streams_txt, Write_presaero_stream_txt, Write_lnd_in, Write_rof_in, Write_Config_Files, Write_drv_in, Write_seq_maps),
                                      modules=("numpy", "netCDF4", "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "signal", "gc", "imp", "getpass", "calendar", "glob","scipy.stats","scipy.signal",'scipy.weave'), group='DAS'))
                        
                        Ens_Index = Ens_Index + 1
                    Job_Num_Per_Node_Index = Job_Num_Per_Node_Index + 1
                
            elif Def_PP == 2:
                if mpi4py_rank == 0:
                    print "********************************************** Using Mpi4Py to Accelerate Call_Model_Operator"
                Ens_Index = mpi4py_rank/mpipy_comm_decomposition
                if Def_Print:
                    print "mpi4py_rank",mpi4py_rank,"Ens_Index",Ens_Index
                if Ens_Index < Ensemble_Number:
                    Call_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                            CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                                          Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, 0,\
                                          Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                          Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                          NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null)    
                
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
        
            else:
                print "*************************************************** Run DAS Sequentially"
                for Ens_Index in range(Ensemble_Number):
                    Call_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                            CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                                          Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, 0,\
                                          Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                          Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                          NC_FileName_Assimilation_2_Initial,  NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null)                          
                  
        if Def_PP == 1 and (not Def_CESM_Multi_Instance):
            print "***************************************** Collect CLM and Observation Information"
            
            job_tid = job_server_node_results[0].tid
            # We append the job of read observation to the end of Run CLM
            for job_server_node in job_server_node_array:
                job_server_node.wait()
                if Def_Print >= 2:
                    job_server_node.print_stats()
    
    if mpi4py_rank == 0:
        if Two_Step_Bias_Estimation_Flag == 0:
            
            print "*************************************************** Run Observation_Blocks Sequentially"
            for Observation_Matrix_Index in range(Dim_Obs_Type):
                print "Read",str(Observation_Matrix_Index+1)+"th","Observation Matrix!"
                SensorType_Sub = SensorType[Observation_Matrix_Index]
                SensorVariable_Sub = SensorVariable[Observation_Matrix_Index]
                SensorQuantity_Sub = SensorQuantity[Observation_Matrix_Index]
                SensorResolution_Sub = SensorResolution[Observation_Matrix_Index]
                Variable_ID_Sub = Variable_ID[Observation_Matrix_Index]
                QC_ID_Sub = QC_ID[Observation_Matrix_Index]
                Observation_File_Name_Sub = Observation_File_Name[Observation_Matrix_Index]
                if Observation_File_Name_Sub != "None":
                    Observation_NLons[Observation_Matrix_Index], Observation_NLats[Observation_Matrix_Index], \
                    Observation_X_Left[Observation_Matrix_Index], Observation_X_Right[Observation_Matrix_Index], \
                    Observation_Y_Lower[Observation_Matrix_Index], Observation_Y_Upper[Observation_Matrix_Index], \
                    Observation_Corelation_Par[Observation_Matrix_Index,::], Observation_Bias_Range, \
                    Observation_Bias_Range_STD, Observation_Bias_Initialization_Flag = \
                    Observation_Blocks(Observation_Matrix_Index, Def_PP, Def_CESM_Multi_Instance, Def_Region, Def_ReBEL, Def_Localization, DasPy_Path, DAS_Depends_Path, DAS_Output_Path, Region_Name, Num_Local_Obs_State,\
                                       Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict,  Call_Gstat_Flag, Plot_Analysis, \
                                       Write_DA_File_Flag, Use_Mask_Flag, Mask_File, Forcing_File_Path_Home, dtime, Observation_Path, Observation_File_Name_Sub, \
                                       DAS_Data_Path, Grid_Resolution_CEA, Grid_Resolution_GEO, NAvalue, Variable_List, Variable_Assimilation_Flag, \
                                       SensorType_Sub, SensorVariable_Sub, SensorQuantity_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub, PDAF_Assim_Framework, PDAF_Filter_Type, \
                                       mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, \
                                       Stop_Year, Stop_Month, Stop_Day, Def_Print, Observation_Bias_Range, Observation_Bias_Range_STD, Observation_Bias_Initialization_Flag, plt, cm, colors, octave, r)
                    
                    NC_FileName_Observation = DAS_Output_Path+"Analysis/"+Region_Name+"/Observation_"+str(Observation_Matrix_Index+1)+".nc"
                    NC_File_Observation = netCDF4.Dataset(NC_FileName_Observation, 'r')
                    Observation_Matrix[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_Matrix'][::]
                    Observation_Variance[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_Variance'][::]
                    Observation_Latitude[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_Latitude'][::]
                    Observation_Longitude[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_Longitude'][::]
                    Observation_View_Zenith_Angle[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_View_Zenith_Angle'][::]
                    Observation_View_Time[Observation_Matrix_Index,::] = NC_File_Observation.variables['Observation_View_Time'][::]
                    Observation_Misc[Observation_Matrix_Index,:,::] = NC_File_Observation.variables['Observation_Misc'][:,:,:]
                    NC_File_Observation.close()
            
            # Re-Initialize
            Def_Par_Optimized = 0
            
            
            print "******************* job_server_node.tid",job_tid
            print ""
            if Def_PP and job_tid >= 500:
                print "Restart PP becuase of too many jobs or too many open files, will make PP Hang on!!"
                job_server_node_array = Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
                while len(job_server_node_array) < 1:
                    job_server_node_array = Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                    job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
                
                
                
            print "************************************************Read Model Results****************************************************************************"
            
        
            
            #os.abort()
            Variable_Min = numpy.zeros(Dim_Obs_Type)
            Variable_Max = numpy.zeros(Dim_Obs_Type)
            # There is a bug in the mismatch between model xy and obs xy because of the project funcion in radal (On JUROPA)
            for Observation_Matrix_Index in range(Dim_Obs_Type):
                SensorType_Sub = SensorType[Observation_Matrix_Index]
                SensorVariable_Sub = SensorVariable[Observation_Matrix_Index]
                SensorQuantity_Sub = SensorQuantity[Observation_Matrix_Index]
                SensorResolution_Sub = SensorResolution[Observation_Matrix_Index]
                Variable_ID_Sub = Variable_ID[Observation_Matrix_Index]
                QC_ID_Sub = QC_ID[Observation_Matrix_Index]
                Observation_File_Name_Sub = Observation_File_Name[Observation_Matrix_Index]
                if Observation_File_Name_Sub != "None":
        
                    Observation_Matrix_Masked = numpy.ma.masked_values(Observation_Matrix[Observation_Matrix_Index,:,:], NAvalue)
                    print "numpy.size(numpy.where(numpy.ma.getmask(Observation_Matrix_Masked)==True))",numpy.size(numpy.where(numpy.ma.getmask(Observation_Matrix_Masked)==True)),numpy.shape(Observation_Matrix_Masked)
                    if (numpy.size(numpy.where(numpy.ma.getmask(Observation_Matrix_Masked)==True)) / 2) == Row_Numbers*Col_Numbers:
                        continue
                    
                    Variable_Min[Observation_Matrix_Index] = numpy.min(Observation_Matrix_Masked)
                    Variable_Max[Observation_Matrix_Index] = numpy.max(Observation_Matrix_Masked)
                    print Observation_Matrix_Index,SensorType_Sub, SensorVariable_Sub, "Variable_Min_Obs",Variable_Min[Observation_Matrix_Index],"Variable_Max_Obs",Variable_Max[Observation_Matrix_Index]
                    
                    if Plot_Analysis and (SensorType_Sub != "COSMOS") and (SensorType_Sub != "InSitu"):
                                
                        w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)       
                        
                        ticks = numpy.arange(Variable_Min[Observation_Matrix_Index], Variable_Max[Observation_Matrix_Index], (Variable_Max[Observation_Matrix_Index] - Variable_Min[Observation_Matrix_Index]) / 5.0)
                        color_boun_list = []
                        color_bound = [Variable_Min[Observation_Matrix_Index], Variable_Max[Observation_Matrix_Index], (Variable_Max[Observation_Matrix_Index] - Variable_Min[Observation_Matrix_Index]) / 100.0]
                        for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                            color_bound[0] += color_bound[2]
                            color_boun_list.append(color_bound[0])
                        
                        
                        fig1 = plt.figure(figsize=(w*2, h*2), dpi=80)
                        fig1.suptitle(DateString_Plot, fontsize=16)
                        ax = fig1.add_subplot(1, 1, 1)
                        im1 = ax.imshow(Observation_Matrix_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                        ax.set_title('Observation_Matrix')
                        plt.grid(True)
                        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Observation_Matrix_"+str(Observation_Matrix_Index)+".png")
                        plt.close('all')
                    #print "******************** Re-assign the Observation Coordinates using the Model Coordinates, Becuase of  the version bug of PROJ.4"
                    #Observation_Longitude[Observation_Matrix_Index,:,:] = MODEL_CEA_X
                    #Observation_Latitude[Observation_Matrix_Index,:,:] = MODEL_CEA_Y
                    Obs_X_Left = numpy.min(Observation_Longitude[Observation_Matrix_Index,:,:])
                    Obs_X_Right = numpy.max(Observation_Longitude[Observation_Matrix_Index,:,:])
                    Obs_Y_Lower = numpy.min(Observation_Latitude[Observation_Matrix_Index,:,:])
                    Obs_Y_Upper = numpy.max(Observation_Latitude[Observation_Matrix_Index,:,:])
                    
                    if Def_Print:
                        print "MODEL_X_Left,MODEL_X_Right",MODEL_X_Left,MODEL_X_Right
                        print "MODEL_Y_Lower,MODEL_Y_Upper",MODEL_Y_Lower,MODEL_Y_Upper
                        print "Obs_X_Left,Obs_X_Right",Obs_X_Left,Obs_X_Right
                        print "Obs_Y_Lower,Obs_Y_Upper",Obs_Y_Lower,Obs_Y_Upper
                    if Grid_Resolution_CEA == SensorResolution[Observation_Matrix_Index]:
                        if numpy.abs(MODEL_X_Left - Obs_X_Left) > 2*Grid_Resolution_CEA:
                            print "The mismatch of model and obs X coordinates is very huge, and the mismatch is", numpy.abs(MODEL_X_Left - Obs_X_Left)
                            os.abort()
                        elif numpy.abs(MODEL_Y_Lower - Obs_Y_Lower) > 2*Grid_Resolution_CEA:
                            print "The mismatch of model and obs Y coordinates is very huge, and the mismatch is", numpy.abs(MODEL_Y_Lower - Obs_Y_Lower)
                            os.abort()        
            
            
            if Initial_Perturbation:
                Additive_Noise_SM_Cov = numpy.square(Additive_Noise_SM_Par[:,0]) * Additive_Noise_SM_Par[:,1:11]
                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                Additive_Noise_SM = numpy.random.multivariate_normal(numpy.zeros(Soil_Layer_Num-5),Additive_Noise_SM_Cov, size=((Ensemble_Number)))
                numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)))
                Additive_Noise_ST = numpy.random.multivariate_normal(numpy.zeros(2),[[0.04,0.64],[0.64,0.04]], size=((Ensemble_Number)))
            
        #                    if Def_Print >= 3:
        #                        print Additive_Noise_SM,Additive_Noise_ST
            
            if Parameter_Optimization or Def_Debug:
                NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
                Sand_Top_Region = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,0,::],axis=0)
                Clay_Top_Region = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,0+Soil_Texture_Layer_Opt_Num,::],axis=0)
                NC_File_Out_Assimilation_2_Parameter.close()
                
            if Def_First_Run == 1:
                print "Open Constant File:", Run_Dir_Home+"_Ens1/"+Constant_File_Name_Header
                NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
                NC_File_Out_Assimilation_2_Initial.sync()
                NC_File_Out_Assimilation_2_Initial.close()
                
            
            if Def_PP==1 and Def_CESM_Multi_Instance and (socket.gethostname()[0:4] != 'node'):
                job_server_node_array = Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
                while len(job_server_node_array) < 1:
                    job_server_node_array = Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                    job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
            
            
        
               
        if Def_PP and Ensemble_Number > 1:
            print "********************************************** Using PP to Accelerate Read_History_File"
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
                    job_server_node_results.append(job_server_node.submit(Read_History_File, args=(Ens_Index, Model_Driver, Def_First_Run, Def_Region, Def_PP, Dim_CLM_State, Row_Numbers, Col_Numbers, column_len, pft_len, \
                                                                                                   Datetime_Start, Datetime_Stop, Start_Year, Stop_Month, Region_Name, Run_Dir_Home, history_file_name, Ensemble_Number, \
                                                                                                   Constant_File_Name_Header, Variable_Assimilation_Flag, Variable_List, Additive_Noise_SM, Additive_Noise_ST, N0, nlyr,\
                                                                 Mean_Dir, Two_Step_Bias_Estimation_Active, Feedback_Assim, Soil_Texture_Layer_Opt_Num, Grid_Resolution_CEA, Grid_Resolution_GEO, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,\
                                                                 SensorType, SensorVariable, Variable_ID, Analysis_Variable_Name, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Def_ParFor, DAS_Depends_Path, DAS_Data_Path, Def_Print,\
                                                                 Irrig_Scheduling, Density_of_liquid_water, Density_of_ice, NAvalue, CLM_NA, omp_get_num_procs_ParFor, DAS_Output_Path, \
                                                                 NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial,  \
                                                                 NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, finidat_name_string, Plot_Analysis, DasPy_Path),
                                                                depfuncs=(Check_Outliers,),
                                  modules=("numpy", "netCDF4", "datetime", "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "gc", "imp",\
                                            "getpass", "calendar","scipy.stats","scipy.signal",'scipy.weave', 'scipy.ndimage'), group='Read_History_File'))
                    
                    Ens_Index = Ens_Index + 1
                Job_Num_Per_Node_Index = Job_Num_Per_Node_Index + 1
                
            for job_server_node in job_server_node_array:
                job_server_node.wait()
                if Def_Print >= 2:
                    job_server_node.print_stats()
            
            if len(job_server_node_results) > 0:
                for job in job_server_node_results:
                    job_index = job_server_node_results.index(job)
                    if job_index > (Ensemble_Number - 1):
                        break
                    if Def_Print:
                        print "Results of ",job_index,"is", job()
                    Soil_Layer_Index_DA = job()[0]
                    Analysis_Variable_Name = job()[1]
                    Constant_File_Name = job()[2]
                              
        else:   
            print "********* Run Read_History_File Sequentially"
            for Ens_Index in range(Ensemble_Number):
                if Def_Print:
                    print "Ens_Index",Ens_Index
                Soil_Layer_Index_DA, Analysis_Variable_Name, Constant_File_Name = \
                    Read_History_File(Ens_Index, Model_Driver, Def_First_Run, Def_Region, Def_PP, Dim_CLM_State, Row_Numbers, Col_Numbers, column_len, pft_len, \
                                      Datetime_Start, Datetime_Stop, Start_Year, Stop_Month, Region_Name, Run_Dir_Home, history_file_name, Ensemble_Number, \
                                      Constant_File_Name_Header, Variable_Assimilation_Flag, Variable_List, Additive_Noise_SM, Additive_Noise_ST, N0, nlyr,
                                     Mean_Dir, Two_Step_Bias_Estimation_Active, Feedback_Assim, Soil_Texture_Layer_Opt_Num, Grid_Resolution_CEA, Grid_Resolution_GEO, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,
                                     SensorType, SensorVariable, Variable_ID, Analysis_Variable_Name, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Def_ParFor, DAS_Depends_Path, DAS_Data_Path, Def_Print,
                                     Irrig_Scheduling, Density_of_liquid_water, Density_of_ice, NAvalue, CLM_NA, omp_get_num_procs_ParFor, DAS_Output_Path, \
                                     NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial,  \
                                     NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, finidat_name_string, Plot_Analysis, DasPy_Path, r)
        
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')        
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')        
        
        Prop_Grid_Array_Sys = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:,:,:,:]
        CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:]
        CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:]
        CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:]
        CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:]
        CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat'][:,:,:]
        CLM_Snow_Depth_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:]
        CLM_Snow_Water_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:]
        CLM_INT_SNOW_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_INT_SNOW_Ensemble_Mat'][:,:,:]
        CLM_FH2OSFC_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_FH2OSFC_Ensemble_Mat'][:,:,:]
        
        CLM_2m_Air_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:]
        CLM_Air_Pressure_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Air_Pressure_Ensemble_Mat'][:,:,:]
        
           
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:] =  0.0
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:,:] =  0.0
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:] = 0.0
        CLM_Soil_Temperature_Ratio_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:]
        CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:,:]
        CLM_Soil_Moisture_Ratio_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:]
        
        print "Collect the Read_History Results."
        for Ens_Index in range(Ensemble_Number):
            if Def_Print:
                print "Ens_Index",Ens_Index
            
            NC_FileName_Out_Ens = DAS_Output_Path+"Analysis/"+Region_Name+"/CLM_History_Ens_"+str(Ens_Index+1)+".nc"
            
            NC_File_Out_Ens = netCDF4.Dataset(NC_FileName_Out_Ens, 'r')
            
            Prop_Grid_Array_Sys[Ens_Index,:,:,:] = NC_File_Out_Ens.variables['Prop_Grid_Array_Sys'][:,:,:]
            CLM_Soil_Moisture_Ensemble_Mat[:,:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:]
            CLM_Soil_Temperature_Ensemble_Mat[:,:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:]
            CLM_Vegetation_Temperature_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:]
            CLM_Ground_Temperature_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:]
            CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat'][:,:]
            CLM_Snow_Depth_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:]
            CLM_Snow_Water_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Snow_Water_Ensemble_Mat'][:,:]
            CLM_INT_SNOW_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_INT_SNOW_Ensemble_Mat'][:,:]
            CLM_FH2OSFC_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_FH2OSFC_Ensemble_Mat'][:,:]
            
            CLM_2m_Air_Temperature_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:]
            CLM_Air_Pressure_Ensemble_Mat[:,:,Ens_Index] = NC_File_Out_Ens.variables['CLM_Air_Pressure_Ensemble_Mat'][:,:]
            
            CLM_Soil_Temperature_Ratio_Ensemble_Mat +=  NC_File_Out_Ens.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:] / numpy.float(Ensemble_Number)
            CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale[:,:,:,Ens_Index] =  NC_File_Out_Ens.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:]
            CLM_Soil_Moisture_Ratio_Ensemble_Mat += NC_File_Out_Ens.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:] / numpy.float(Ensemble_Number)
               
            NC_File_Out_Ens.close()
        
        NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:,:,:,:] = Prop_Grid_Array_Sys
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:] = CLM_Soil_Moisture_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:] = CLM_Soil_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:] = CLM_Vegetation_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:] = CLM_Ground_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat'][:,:,:] = CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:] = CLM_Snow_Depth_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:] = CLM_Snow_Water_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_INT_SNOW_Ensemble_Mat'][:,:,:] = CLM_INT_SNOW_Ensemble_Mat
        NC_File_Out_Assimilation_2_Initial.variables['CLM_FH2OSFC_Ensemble_Mat'][:,:,:] = CLM_FH2OSFC_Ensemble_Mat
        
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:] = CLM_Soil_Temperature_Ratio_Ensemble_Mat
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:,:] = CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:] = CLM_Soil_Moisture_Ratio_Ensemble_Mat
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:] = CLM_2m_Air_Temperature_Ensemble_Mat
        NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Air_Pressure_Ensemble_Mat'][:,:,:] = CLM_Air_Pressure_Ensemble_Mat
        del Prop_Grid_Array_Sys,CLM_Soil_Moisture_Ensemble_Mat,CLM_Soil_Temperature_Ensemble_Mat,CLM_Vegetation_Temperature_Ensemble_Mat
        del CLM_Ground_Temperature_Ensemble_Mat,CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat
        del CLM_Snow_Depth_Ensemble_Mat,CLM_Snow_Water_Ensemble_Mat,CLM_INT_SNOW_Ensemble_Mat,CLM_FH2OSFC_Ensemble_Mat
        del CLM_Soil_Temperature_Ratio_Ensemble_Mat,CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale,CLM_Soil_Moisture_Ratio_Ensemble_Mat
        del CLM_2m_Air_Temperature_Ensemble_Mat,CLM_Air_Pressure_Ensemble_Mat
        
        
        NC_File_Out_Assimilation_2_Initial.sync()
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Diagnostic.sync()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        
            
        print "***************************** Call H Operator"
            
        if Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")] == 1 or \
            Variable_Assimilation_Flag[Variable_List.index("Sensible_Heat")] == 1:        
            
            print "################### Using TSM as Observation Operator for Surface Temperature Assimilation and Sensible Heat Assimilation"
            if Def_PP and Ensemble_Number > 1:
                NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
                Prop_Grid_Array_Sys = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, :, :, :]
                
                print "********************************************** Using PP to Accelerate TSM"             
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
                        job_server_node_results.append(job_server_node.submit(Run_TSM, args=(Ens_Index, \
                                Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, NAvalue, Variable_Assimilation_Flag, Variable_List,
                                DasPy_Path, Def_Print, SensorVariable, Observation_View_Zenith_Angle, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias),
                                depfuncs=(Check_Outliers,),
                        modules=("numpy", "netCDF4", "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "gc", "imp", "getpass", "calendar", "scipy.signal",'scipy.weave'), group='TSM'))
                        
                        Ens_Index = Ens_Index + 1
                    Job_Num_Per_Node_Index = Job_Num_Per_Node_Index + 1
                        
                for job_server_node in job_server_node_array:
                    job_server_node.wait()
                    if Def_Print >= 2:
                        job_server_node.print_stats()
                
                if len(job_server_node_results) > 0:
                    
                    for job in job_server_node_results:
                        job_index = job_server_node_results.index(job)
                        if job_index > (Ensemble_Number - 1):
                            break
                        if Def_Print:
                            print "Results of ",job_index,"is", job()[Station_XY_Index[0][1],Station_XY_Index[0][0]]
                        Prop_Grid_Array_Sys[job_index, Variable_List.index("Surface_Temperature"), :, :] = job()
                        Prop_Grid_Array_Sys[job_index, Variable_List.index("Sensible_Heat"), :, :] = job()
                
                NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, :, :, :] = Prop_Grid_Array_Sys
                del Prop_Grid_Array_Sys
                NC_File_Out_Assimilation_2_Initial.sync()
                NC_File_Out_Assimilation_2_Initial.close()
                
            else:
                print "********************************************** Run TSM Sequentially"
                Prop_Grid_Array_Sys_Temp = numpy.zeros((Ensemble_Number,Row_Numbers,Col_Numbers),dtype=numpy.float32)
                for Ens_Index in range(Ensemble_Number):
                    Prop_Grid_Array_Sys_Temp[Ens_Index,:,:] =  Run_TSM(Ens_Index, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, NAvalue, Variable_Assimilation_Flag, Variable_List,
                                                                        DasPy_Path, Def_Print, SensorVariable, Observation_View_Zenith_Angle, NC_FileName_Assimilation_2_Constant, 
                                                                        NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias)
                
                NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
                NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Variable_List.index("Surface_Temperature"), :, :] = Prop_Grid_Array_Sys_Temp
                #NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Variable_List.index("Sensible_Heat"), :, :] = Prop_Grid_Array_Sys_Temp
                NC_File_Out_Assimilation_2_Initial.sync()
                NC_File_Out_Assimilation_2_Initial.close()
                del Prop_Grid_Array_Sys_Temp
    
        print '=========================== Generate Model Ensembles ========================================================'  
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')                          
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        CLM_2m_Air_Temperature_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:],axis=2)
        
        Mean_Index_Prop_Grid_Array_Sys = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'],axis=0)
        
        NC_File_Out_Assimilation_2_Initial.close()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+')
        for Prop_Grid_Array_Sys_Index in range(Dim_CLM_State):
            if Variable_Assimilation_Flag[Prop_Grid_Array_Sys_Index]:
                if Def_Print >= 2:
                    print numpy.where(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :] == NAvalue)
                    #print Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :]
                    #Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :][numpy.where(numpy.isnan(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :]))] = numpy.min(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :][not numpy.where(numpy.isnan(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :]))])
                    #Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :][numpy.where(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :] == NAvalue)] = numpy.mean(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :][numpy.where(Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index, :, :] != NAvalue)])
                    print numpy.min(numpy.ma.masked_where(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :] == NAvalue, \
                                                          NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :])), \
                                                          numpy.max(numpy.ma.masked_where(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :] \
                                                                                          == NAvalue, NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :]))
                
                NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, Prop_Grid_Array_Sys_Index, :, :] = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :]
                #print "SysModel:",Prop_Grid_Array_Sys[:,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],"ObsModel:",\
                #NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:,Prop_Grid_Array_Sys_Index,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                #os.abort()
        NC_File_Out_Assimilation_2_Initial.sync()
        NC_File_Out_Assimilation_2_Initial.close()
                
        print '=========================== Finish Generating Model Ensembles ========================================================'     
        
        
        if Def_Debug and Def_First_Run:
            NC_FileName_Assimilation_2_Parameter_First = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_First.nc"
            copyLargeFile(NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_First)
        
        NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
        
        for Prop_Grid_Array_Sys_Index in range(Dim_CLM_State):
            if Variable_Assimilation_Flag[Prop_Grid_Array_Sys_Index]:
                if Use_Mask_Flag:
                    # find the model grids which are need to be assiilated.
                    # Read the Mask Grid
                    
                    Mask_Value = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][:,:]                    
                    
                    Mask[Prop_Grid_Array_Sys_Index, 0, ::] = Mask_X
                    Mask[Prop_Grid_Array_Sys_Index, 1, ::] = Mask_Y
                    Mask[Prop_Grid_Array_Sys_Index, 2, ::] = Mask_Value
                    
                    Mask_Index[Prop_Grid_Array_Sys_Index, ::] = False
                    numexpr_a = Mask[Prop_Grid_Array_Sys_Index, 2, ::]
                    numexpr_b = NAvalue
                    numexpr_c = numpy.where(numexpr_a == numexpr_b)
                    NA_Index = numexpr_c
                    if Def_Print:
                        print numpy.shape(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys']), numpy.shape(Mask_Index), numpy.shape(NA_Index)
                        print "numpy.where(Mask[Prop_Grid_Array_Sys_Index, 2, ::] == NAvalue)",numpy.size(NA_Index),numpy.size(Mask_Index[Prop_Grid_Array_Sys_Index,::])
                    if numpy.shape(NA_Index)[1] > 0:
                        Mask_Index[Prop_Grid_Array_Sys_Index,::][NA_Index] = True
                    numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][0, Prop_Grid_Array_Sys_Index, :, :])
                    numexpr_b = NAvalue
                    numexpr_c = numpy.where(numexpr_a == numexpr_b)
                    NA_Index = numexpr_c
                    if Def_Print:
                        print numpy.shape(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys']), numpy.shape(Mask_Index), numpy.shape(NA_Index)
                        print "numpy.where(Prop_Grid_Array_Sys[0, Prop_Grid_Array_Sys_Index, :, :] == NAvalue)",numpy.size(NA_Index)
                    if numpy.shape(NA_Index)[1] > 0:
                        Mask_Index[Prop_Grid_Array_Sys_Index,::][NA_Index] = True
                        
                    if Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")] == 1 or Feedback_Assim == 1:
                        # Exclude the Fronzen Grid Cells
                        CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat'][:,:,:],axis=2)
                        
                        numexpr_a = CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat
                        numexpr_b = 0.0
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        NA_Index = numexpr_c
                        print "numpy.where(CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat > 0)",numpy.size(NA_Index)
                        if numpy.shape(NA_Index)[1] > 0:
                            Mask_Index[Prop_Grid_Array_Sys_Index,::][NA_Index] = True
                     
                    Valid_Grid_Num = numpy.size(numpy.where(Mask_Index[Prop_Grid_Array_Sys_Index,::] == False))
                    print "numpy.size(numpy.where(Mask_Index[Prop_Grid_Array_Sys_Index,::] == False))",Valid_Grid_Num
                    if Valid_Grid_Num == 0:
                        break
                        continue
                
                else:
                    
                    Mask[Prop_Grid_Array_Sys_Index, 0, ::] = Mask_X
                    Mask[Prop_Grid_Array_Sys_Index, 1, ::] = Mask_Y
                    Mask[Prop_Grid_Array_Sys_Index, 2, ::] = 1
                    
                    Mask_Index[Prop_Grid_Array_Sys_Index, ::] = False
                    
                    numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][0, Prop_Grid_Array_Sys_Index, :, :])
                    numexpr_b = NAvalue
                    numexpr_c = numpy.where(numexpr_a == numexpr_b)
                    
                    NA_Index = numexpr_c
                    if Def_Print:
                        print numpy.shape(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][0, Prop_Grid_Array_Sys_Index, :, :]), numpy.shape(Mask_Index[Prop_Grid_Array_Sys_Index, :, :]), numpy.shape(NA_Index)
                    Mask_Index[Prop_Grid_Array_Sys_Index,::][NA_Index] = True
                                  
                    #print Mask[1,::]
                if Write_DA_File_Flag:
                    Mask_Mat = numpy.zeros((numpy.size(Mask[Prop_Grid_Array_Sys_Index, 2, ::]), 3),dtype=numpy.float32)
                    Mask_Mat[:, 0] = Mask[Prop_Grid_Array_Sys_Index, 0, ::].flatten()
                    Mask_Mat[:, 1] = Mask[Prop_Grid_Array_Sys_Index, 1, ::].flatten()
                    Mask_Mat[:, 2] = Mask[Prop_Grid_Array_Sys_Index, 2, ::].flatten()
                    numpy.savetxt(DasPy_Path+'Analysis/DAS_Temp/Mask_'+str(Prop_Grid_Array_Sys_Index+1)+'.txt', Mask_Mat)
        
    #                w,h = plt.figaspect(float(Row_Numbers)/Col_Numbers)
    #                fig1 = plt.figure(figsize=(w,h))
    #                ax1 = fig1.add_subplot(1,1,1)
    #                im1 = ax1.imshow(Mask[1, 2, ::], cmap=cm.jet, interpolation='bilinear')
    #                plt.colorbar(im1)
    #                plt.show()
        NC_File_Out_Assimilation_2_Constant.close()
        NC_File_Out_Assimilation_2_Initial.close()
        
        print "------------------ Call DAS_Observation_Operator"
        DAS_Observation_Operator(Ensemble_Number,Def_PP,Def_Print, active_nodes_server, job_server_node_array, DAS_Depends_Path, Def_Region, Soil_Layer_Num, Grid_Resolution_CEA, DasPy_Path, \
            Start_Year, Start_Month, Start_Day, Start_Hour, LATIXY_Mat, LONGXY_Mat, Mask_Index, Station_XY, Station_XY_Index, N0, nlyr,
            Def_First_Run_RTM, DAS_Output_Path, CMEM_Work_Path_Array, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, Datetime_Initial,
            maxpft, PCT_PFT_High, PCT_PFT_Low, PCT_PFT_WATER, Density_of_liquid_water, Run_Dir_Array, Row_Numbers_String, Col_Numbers_String,
            Forcing_File_Path_Array, history_file_name, Constant_File_Name, NAvalue, CLM_NA,
            Row_Numbers, Col_Numbers, Variable_Assimilation_Flag,Variable_List,Variable_ID,SensorVariable,SensorType,
            Soil_Density, DAS_Data_Path, Region_Name, mksrf_edgen, mksrf_edges, mksrf_edgew, mksrf_edgee, Grid_Resolution_GEO, Soil_Layer_Index_DA, Soil_Texture_Layer_Opt_Num,
            NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Parameter,
            NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, omp_get_num_procs_ParFor, octave, r)
        
        os.chdir(DasPy_Path)                                   
                
        
        NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')
        if Def_Print:
            print numpy.shape(NC_File_Out_Assimilation_2_Diagnostic.variables['Mask_Index'][:,:,:]),numpy.shape(Mask_Index)
        NC_File_Out_Assimilation_2_Diagnostic.variables['Mask_Index'][:,:,:] = Mask_Index
        NC_File_Out_Assimilation_2_Diagnostic.sync()
        NC_File_Out_Assimilation_2_Diagnostic.close()
        
        #os.abort()
        
        try:
            if Parameter_Optimization or Def_Debug:
                del Sand_Top_Region, Clay_Top_Region
        except:
            pass
    
        Backup_Initial_Files(NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy,
                             NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy,
                             NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy,
                             NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy,
                             NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy)
    
    else:
        Observation_Matrix = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_Variance = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_Latitude = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_Longitude = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_View_Zenith_Angle = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_View_Time = numpy.empty((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_NLons = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_NLats = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_X_Left = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_X_Right = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_Y_Lower = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_Y_Upper = numpy.empty(len(SensorType),dtype=numpy.float32)
        Observation_Misc = numpy.empty((len(SensorType), 10, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Observation_Corelation_Par = numpy.empty((len(SensorType), 5, 2),dtype=numpy.float32)   
        Analysis_Variable_Name = None
        Constant_File_Name = None
        Def_Par_Optimized = None
        Soil_Layer_Index_DA = None
        Mask = numpy.empty((Dim_CLM_State, 3, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask_Index = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.bool)
        Model_Variance = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Def_First_Run_RTM = None
        ECOCVL_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOCVH_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVL_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVH_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOWAT_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        Mean_Index_Prop_Grid_Array_Sys = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Model_State_Inflation_Range = numpy.empty((Dim_CLM_State,2), dtype=numpy.float32)
        Model_State_Inflation_Range_STD = numpy.empty(Dim_CLM_State, dtype=numpy.float32)
        Model_Bias_Range = numpy.empty((Dim_CLM_State, 2), dtype=numpy.float32)
        Observation_Bias_Range = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, 2), dtype=numpy.float32)
        Model_Bias_Range_STD = numpy.empty((Dim_CLM_State, 2), dtype=numpy.float32)
        Observation_Bias_Range_STD = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, 2), dtype=numpy.float32)
        Model_Bias_STD = numpy.empty(Dim_CLM_State, dtype=numpy.float32)
        Observation_Bias_STD = numpy.empty((Dim_CLM_State,Dim_Observation_Quantity), dtype=numpy.float32)
        Observation_Bias_Initialization_Flag = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, Ensemble_Number), dtype=numpy.float32)
        job_server_node_array = None
        active_nodes_server = None
        
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
    
    gc.collect()
    del gc.garbage[:]
    
    return Observation_Matrix, Observation_Longitude, Observation_Latitude, Observation_Variance, \
        Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, \
        Observation_Y_Lower, Observation_Y_Upper, Observation_Misc, Observation_View_Zenith_Angle, Observation_View_Time, Observation_Corelation_Par, \
        Analysis_Variable_Name, Constant_File_Name, Def_Par_Optimized, Soil_Layer_Index_DA,\
        Mask, Mask_Index, Model_Variance, Def_First_Run_RTM, ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat, \
        Mean_Index_Prop_Grid_Array_Sys, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, \
        Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, \
        Observation_Bias_Initialization_Flag, job_server_node_array, active_nodes_server
    




def DAS_Config(mpi4py_rank, Model_Driver, Start_Year, Start_Month, Start_Day, Start_Hour, Start_Minute, Datetime_Start, \
                Def_CESM_Multi_Instance, Ensemble_Number, Region_Name, Def_Initial, Run_Dir_Home,\
                Datetime_Start_Init, Def_ParFor, DAS_Data_Path, Def_Region,\
                Def_PP, Row_Numbers, Col_Numbers, Def_Print, DAS_Depends_Path, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type):
    #--------------------------------- UTC Time ----------------------------------
    Current_Day = datetime.date.today().day
    if Current_Day < 10:
        Current_Day_String = "0"+ str(Current_Day)
    else:
        Current_Day_String = str(Current_Day)
    
    # Model Initial Time
    Initial_Year = string.atoi(Start_Year)
    Initial_Month = 01
    Initial_Day = 01
    Initial_Hour = 00
    Initial_Minute = 00
    Datetime_Initial = datetime.datetime(Initial_Year, Initial_Month, Initial_Day, Initial_Hour, Initial_Minute)
    
    # Get the Number of Days in Each Day
    Num_of_Days_Monthly = numpy.zeros(12)
    Num_of_Days_Monthly = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    Num_of_Days_Monthly[1] = (datetime.datetime(Initial_Year,3,1,0)-datetime.datetime(Initial_Year,2,1,0)).days
    if mpi4py_rank == 0:
        print Num_of_Days_Monthly
    
    Rapid_delta = datetime.timedelta(days=5)
    Datetime_Start_Rapid = datetime.datetime(string.atoi(Start_Year), string.atoi(Start_Month), string.atoi(Start_Day), string.atoi(Start_Hour), string.atoi(Start_Minute))
    
    Datetime_Start_Initial_File = Datetime_Start - datetime.timedelta(hours=1)
    Datetime_Start_Initial_File_Init = datetime.datetime(Datetime_Start_Initial_File.year, Datetime_Start_Initial_File.month, Datetime_Start_Initial_File.day, 00, 00)
    
    Month_String = str(Datetime_Start_Initial_File.month).zfill(2)
    Day_String = str(Datetime_Start_Initial_File.day).zfill(2)
    Hour_Tod_string = str((Datetime_Start_Initial_File - Datetime_Start_Initial_File_Init).seconds).zfill(5)

    finidat_initial_CLM = Region_Name + ".clm2.r."+str(Datetime_Start_Initial_File.year)+"-"+Month_String+"-"+Day_String+"-"+Hour_Tod_string+".nc"
    finidat_initial_PFCLM = Region_Name + ".clm2.r."+str(Datetime_Start_Initial_File.year)+"-"+Month_String+"-"+Day_String+"-"+Hour_Tod_string+".nc"
    
    if not Def_Initial:
        finidat_initial_CLM = ""
        finidat_initial_PFCLM = ""
    
    if mpi4py_rank == 0:
        print "============================================================================================================================"
        print "finidat_initial_CLM",finidat_initial_CLM
        print "finidat_initial_PFCLM",finidat_initial_PFCLM
        print "============================================================================================================================"   

    Constant_File_Name_Header = Region_Name + ".clm2.h0."+Start_Year+"-"+Start_Month+"-"+Start_Day+"-"+str((Datetime_Start - Datetime_Start_Init).seconds).zfill(5)+".nc"
    
    if mpi4py_rank == 0:
        print "============================================================================================================================"
        print "Constant_File_Name_Header",Constant_File_Name_Header
        print "============================================================================================================================"
    
    Soil_Layer_Num = 15
    Snow_Layer_Num = 5
    ParFlow_Layer_Num = 30
    
    
    CLM_Flag = "3D" # Define whether to use the 3D version or 1D version
    maxpft = 17
    
    if Model_Driver == "CLM_CROP":
        maxpft = 21
    
    numrad = 2 #(1 = visible, 2 = NIR)
    Density_of_liquid_water = 1000.0 # Units (km/m^3)
    Density_of_ice = 917.0 # Units (km/m^3)
    Freezing_temperature_of_fresh_water = 273.15  # Units (K)
    
    num_processors = multiprocessing.cpu_count()
    if mpi4py_rank == 0:
        print 'cpu_count = %d\n' % num_processors
    
    if Def_ParFor:
        if num_processors <= 4:
            omp_get_num_procs_ParFor = 2
        else:
            omp_get_num_procs_ParFor = num_processors
    else:
        omp_get_num_procs_ParFor = 2
    
    if socket.gethostname()[0:4] == 'node':
        num_processors = 16
        omp_get_num_procs_ParFor = 16
    elif socket.gethostname()[0:2] == 'jr' or socket.gethostname()[0:2] == 'j3':
        num_processors = 48
        omp_get_num_procs_ParFor = 48
    elif socket.gethostname()[0] == 'j':
        num_processors = 16
        omp_get_num_procs_ParFor = 16
    elif socket.gethostname()[0] == 'n':
        num_processors = 16
        omp_get_num_procs_ParFor = 16
    
    omp_get_num_procs_ParFor = max([1,num_processors/Ensemble_Number])
    NSLOTS = int(sys.argv[1])   # How many processors assigned
    Node_Num = NSLOTS / num_processors
    
    ntasks_CLM = numpy.ones(8,dtype=numpy.integer)
    rootpe_CLM = numpy.zeros(8,dtype=numpy.integer)
    ntasks_CLM[:] = 1
    #nthreads_CLM = num_processors
    nthreads_CLM = 1    # OpenMP does not work for CLM4.5
    
    Model_Path = DAS_Depends_Path + "bin/cesm_sp_serial.exe"
    
    if mpi4py_rank == 0:
        print "**********Model_Path",Model_Path
    #os.abort()
    #---------------------------Split the Block into SubBlocks to do data assimilation
    # Default Values
    Sub_Block_Ratio_Row = 1
    Sub_Block_Ratio_Col = 1
       
    if PDAF_Assim_Framework == 2:
        Sub_Block_Ratio_Row = 1
        Sub_Block_Ratio_Col = 1
    else:
        if Def_Region == 3:
            Sub_Block_Ratio_Row = 4
            Sub_Block_Ratio_Col = 1
    
    Row_Offset = numpy.zeros((Sub_Block_Ratio_Row*Sub_Block_Ratio_Col,2),dtype=numpy.integer)
    Col_Offset = numpy.zeros((Sub_Block_Ratio_Row*Sub_Block_Ratio_Col,2),dtype=numpy.integer)
    
    Row_Numbers_SubBlock_Array_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Col_Numbers_SubBlock_Array_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Row_Numbers_SubBlock_Array_Mat[:,:] = numpy.floor(Row_Numbers / Sub_Block_Ratio_Row)
    Col_Numbers_SubBlock_Array_Mat[:,:] = numpy.floor(Col_Numbers / Sub_Block_Ratio_Col)
    
    for Sub_Block_Index_Col in range(Sub_Block_Ratio_Col):
        Row_Numbers_SubBlock_Array_Mat[-1,Sub_Block_Index_Col] += Row_Numbers % Sub_Block_Ratio_Row
    for Sub_Block_Index_Row in range(Sub_Block_Ratio_Row):
        Col_Numbers_SubBlock_Array_Mat[Sub_Block_Index_Row,-1] += Col_Numbers % Sub_Block_Ratio_Col
    Row_Numbers_SubBlock_Array = Row_Numbers_SubBlock_Array_Mat.flatten()
    Col_Numbers_SubBlock_Array = Col_Numbers_SubBlock_Array_Mat.flatten()
    if mpi4py_rank == 0:
        if Def_Print:
            print "Row_Numbers_SubBlock_Array",Row_Numbers_SubBlock_Array
            print "Col_Numbers_SubBlock_Array",Col_Numbers_SubBlock_Array
    
    Row_Numbers_SubBlock_Array_Mat_Cumsum = numpy.cumsum(Row_Numbers_SubBlock_Array_Mat,axis=0)
    Col_Numbers_SubBlock_Array_Mat_Cumsum = numpy.cumsum(Col_Numbers_SubBlock_Array_Mat,axis=1)
    if mpi4py_rank == 0:
        if Def_Print:
            print "Row_Numbers_SubBlock_Array_Mat_Cumsum",Row_Numbers_SubBlock_Array_Mat_Cumsum
            print "Col_Numbers_SubBlock_Array_Mat_Cumsum",Col_Numbers_SubBlock_Array_Mat_Cumsum
        
    Sub_Block_Row_Start_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Sub_Block_Row_End_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Sub_Block_Col_Start_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Sub_Block_Col_End_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Sub_Block_Row_Start_Mat[0,:] = 0
    Sub_Block_Row_End_Mat[0,:] = Row_Numbers_SubBlock_Array_Mat_Cumsum[0,:]
    Sub_Block_Col_Start_Mat[:,0] = 0
    Sub_Block_Col_End_Mat[:,0] = Col_Numbers_SubBlock_Array_Mat_Cumsum[:,0]
    
    for Sub_Block_Index_Row in range(1,Sub_Block_Ratio_Row,1):
        Sub_Block_Row_Start_Mat[Sub_Block_Index_Row,:] = Row_Numbers_SubBlock_Array_Mat_Cumsum[Sub_Block_Index_Row-1,:]
        Sub_Block_Row_End_Mat[Sub_Block_Index_Row,:] = Row_Numbers_SubBlock_Array_Mat_Cumsum[Sub_Block_Index_Row,:]
            
    for Sub_Block_Index_Col in range(1,Sub_Block_Ratio_Col,1):
        Sub_Block_Col_Start_Mat[:,Sub_Block_Index_Col] = Col_Numbers_SubBlock_Array_Mat_Cumsum[:,Sub_Block_Index_Col-1]
        Sub_Block_Col_End_Mat[:,Sub_Block_Index_Col] = Col_Numbers_SubBlock_Array_Mat_Cumsum[:,Sub_Block_Index_Col]
    
    Sub_Block_Index_Row_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    Sub_Block_Index_Col_Mat = numpy.zeros((Sub_Block_Ratio_Row,Sub_Block_Ratio_Col),dtype=numpy.integer)
    
    for Sub_Block_Index_Row in range(Sub_Block_Ratio_Row):
        Sub_Block_Index_Col_Mat[Sub_Block_Index_Row,:] = range(Sub_Block_Ratio_Col)
    for Sub_Block_Index_Col in range(Sub_Block_Ratio_Col):
        Sub_Block_Index_Row_Mat[:,Sub_Block_Index_Col] = range(Sub_Block_Ratio_Row)
    Sub_Block_Index_Row_Mat_Vector = Sub_Block_Index_Row_Mat.flatten()
    Sub_Block_Index_Col_Mat_Vector = Sub_Block_Index_Col_Mat.flatten()
    
    for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
        Sub_Block_Index_Row = Sub_Block_Index_Row_Mat_Vector[Block_Index]
        Sub_Block_Index_Col = Sub_Block_Index_Col_Mat_Vector[Block_Index]
        Sub_Block_Row_Start_Mat[Sub_Block_Index_Row,Sub_Block_Index_Col] += Row_Offset[Block_Index,0]
        Sub_Block_Row_End_Mat[Sub_Block_Index_Row,Sub_Block_Index_Col] += Row_Offset[Block_Index,1]
        Sub_Block_Col_Start_Mat[Sub_Block_Index_Row,Sub_Block_Index_Col] += Col_Offset[Block_Index,0]
        Sub_Block_Col_End_Mat[Sub_Block_Index_Row,Sub_Block_Index_Col] += Col_Offset[Block_Index,1]
            
    Sub_Block_Row_Start_Array = Sub_Block_Row_Start_Mat.flatten()
    Sub_Block_Row_End_Array = Sub_Block_Row_End_Mat.flatten()
    Sub_Block_Col_Start_Array = Sub_Block_Col_Start_Mat.flatten()
    Sub_Block_Col_End_Array = Sub_Block_Col_End_Mat.flatten()
    Row_Numbers_SubBlock_Array = Sub_Block_Row_End_Array - Sub_Block_Row_Start_Array
    Col_Numbers_SubBlock_Array = Sub_Block_Col_End_Array - Sub_Block_Col_Start_Array
    if mpi4py_rank == 0:
        if Def_Print:
            print "Sub_Block_Row_Start_Array", Sub_Block_Row_Start_Array, "Sub_Block_Row_End_Array",Sub_Block_Row_End_Array
            print "Sub_Block_Col_Start_Array", Sub_Block_Col_Start_Array, "Sub_Block_Col_End_Array",Sub_Block_Col_End_Array
            print "Row_Numbers_SubBlock_Array", Row_Numbers_SubBlock_Array, "Col_Numbers_SubBlock_Array",Col_Numbers_SubBlock_Array
    
    return Num_of_Days_Monthly, Datetime_Initial, NSLOTS, Constant_File_Name_Header, finidat_initial_CLM, finidat_initial_PFCLM, \
            Soil_Layer_Num, Snow_Layer_Num, ParFlow_Layer_Num, maxpft, numrad, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water,\
            ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor, Model_Path, CLM_Flag, \
            Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,\
            Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,


def PP_Worker(Node, Port):
    active_nodes_server_index = Node
    ppserver_port = Port
    ppservers_node=(active_nodes_server_index+":"+str(ppserver_port),)
    print "ppservers_node",ppservers_node
    job_server_node = pp.Server(ncpus=0, ppservers=ppservers_node,secret="123456",socket_timeout=999999999)
    time.sleep(1)
    print "Starting pp with",job_server_node.get_active_nodes(), "workers on node",active_nodes_server_index,ppservers_node
    active_nodes = job_server_node.get_active_nodes()
    del active_nodes['local']
    print active_nodes
    while len(active_nodes.keys()) == 0:
        job_server_node = pp.Server(ncpus=0, ppservers=ppservers_node,secret="123456",socket_timeout=999999999)
        time.sleep(1)
        print "Starting pp with",job_server_node.get_ncpus(),job_server_node.get_active_nodes(), "workers on node",active_nodes_server_index,ppservers_node
        active_nodes = job_server_node.get_active_nodes()
        del active_nodes['local']
    
    return job_server_node


def Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port):
    
    if mpi4py_rank == 0:
        print "Starting ppserver on ",NSLOTS," processors......" 
        print "one ppserver will start per processor , log messages in the file "+DAS_Output_Path+"ppserver_log.txt"
    
    # Should be changed for specific application
    PROCS_PER_NODE = 16
    PP_Servers_Per_Node = 1
    
    if mpi4py_rank == 0:
        try:
            os.remove(DAS_Output_Path+"ppserver_log.txt")
        except:
            pass
    
    if mpi4py_rank == 0:
        start = time.time()
    
    job_server = []
    job_server_node_array = []
    nodelist = []
    
    if mpi4py_rank == 0:
        print "len(sys.argv)",len(sys.argv)
    
    if Def_PP and Ensemble_Number > 1:
        
        # Parallel Python
        import pp
        # tuple of all parallel python servers to connect with
        
        if True:
            if True:
                
                if mpi4py_rank == 0:
                    print "Open Nodelist File",DAS_Output_Path+"nodefile.txt"
                    subprocess.call("rm -rf "+DAS_Output_Path+"nodefile.txt",shell=True)
                
                if mpi4py_rank == 0:
                    CMD_String = "echo `cat $PBS_NODEFILE` > "+DAS_Output_Path+"nodefile.txt"
                    print "CMD_String",CMD_String
                    pipe = subprocess.Popen(CMD_String, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,stderr=subprocess.PIPE, close_fds=True)
                    
                    stdout_value, stderr_value = pipe.communicate()
                    print 'stdout_value:', repr(stdout_value)
                    print 'stderr_value:', repr(stderr_value)
                    
                    nodelist_file = open(DAS_Output_Path+"nodefile.txt",'r')
                    nodelist = nodelist_file.readline()
                    nodelist_file.close()
                    active_nodes_server = list(set(string.split(nodelist)))
                else:
                    active_nodes_server = None
                
                if Def_PP == 2:
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
                    active_nodes_server = mpi4py_comm.bcast(active_nodes_server)
                    if mpi4py_rank == 0:
                        print "type(active_nodes_server)",type(active_nodes_server)
                
                if mpi4py_rank == 0:
                    print "active_nodes_server", active_nodes_server
                
                if mpi4py_rank == 0:
                    print "------------Calculate the PP_Servers_Per_Node"
                PP_Servers_Per_Node = int(numpy.ceil(float(Ensemble_Number+1.0) / len(active_nodes_server)))
                #PP_Servers_Per_Node = PROCS_PER_NODE
                if mpi4py_rank == 0:
                    print "----------------------PP_Servers_Per_Node is",PP_Servers_Per_Node
                    print ""
                
                if Def_PP == 2:
                    
                    print "mpi4py_rank",mpi4py_rank,"mpi4py_name",mpi4py_name
                    
                    Command_String = "ppserver.py -p "+str(PP_Port+mpi4py_rank)+" -w 1 -t 120 -k 9999999 -s '123456' &"
                    
                    if mpi4py_rank == 0:
                        PP_Port_array = ['' for Index in range(NSLOTS)]
                        Node_array = ['' for Index in range(NSLOTS)]
                    else:
                        PP_Port_array= None
                        Node_array = None
                    
                    PP_Port_array = mpi4py_comm.scatter(PP_Port_array)
                    Node_array = mpi4py_comm.scatter(Node_array)
                    
                    if mpi4py_rank%PROCS_PER_NODE < PP_Servers_Per_Node:
                        pipe = subprocess.Popen(Command_String, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,stderr=subprocess.PIPE, close_fds=True)
                        active_nodes_server_index = (mpi4py_rank+1)/PROCS_PER_NODE
                        ppserver_port = PP_Port+mpi4py_rank
                        print "mpi4py_rank",mpi4py_rank,"ppserver_port",ppserver_port,"mpi4py_name",mpi4py_name
                        PP_Port_array = str(ppserver_port)
                        Node_array = mpi4py_name
                    
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
                    
                    PP_Port_array = mpi4py_comm.gather(PP_Port_array)
                    Node_array = mpi4py_comm.gather(Node_array)
                    
                    if mpi4py_rank == 0:
                        print "PP_Port_array",PP_Port_array
                        print "Node_array",Node_array
                        ppserver_file = open(DAS_Output_Path+"ppserver_log.txt",'w')
                        for Index in range(NSLOTS):
                            if PP_Port_array[Index] != '':
                                ppserver_file.write("Starting ppserver on "+Node_array[Index]+" "+PP_Port_array[Index]+"\n")
                        ppserver_file.close()
                    
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
                    
                else:
                    
                    # use mpiexec to call start_pp_server
                    ppserver_CMD_String = DAS_Depends_Path+"bin/mpiexec -n "+str(NSLOTS)+" "+DasPy_Path+"Utilities/Start_PP/start_pp_server "+str(PP_Port)+" "+str(PP_Servers_Per_Node)+" "+str(PROCS_PER_NODE)+" | tee "+DAS_Output_Path+"/ppserver_log.txt &"
                    #ppserver_CMD_String = DAS_Depends_Path+"bin/mpiexec --machinefile="+MPI_Machine_File_Name_PP+" -n "+str(PP_Servers_Per_Node*len(active_nodes_server_without_duplicates))+" "+DasPy_Path+"Utilities/Start_PP/start_pp_server "+str(PP_Port)+" | tee "+DAS_Output_Path+"/ppserver_log.txt &"
                    print "ppserver_CMD_String",ppserver_CMD_String
                    subprocess.call("killall -9 -q -w psilogger &> /dev/null",shell=True)
                    subprocess.call(ppserver_CMD_String,shell=True)
                    
                    print "active_nodes_server",active_nodes_server
                    
                    Node_ID = numpy.zeros(len(active_nodes_server),dtype=numpy.integer)
                    for Node_Index in range(len(active_nodes_server)):
                        # Get the node number ID
                        if socket.gethostname()[0] == 'n':
                            Node_ID[Node_Index] = int(active_nodes_server[Node_Index][1:])
                        elif socket.gethostname()[0] == 'j':
                            Node_ID[Node_Index] = int(active_nodes_server[Node_Index][5:])
                        elif socket.gethostname()[0:4] == 'node':
                            Node_ID[Node_Index] = int(active_nodes_server[Node_Index][4:])
                    Node_ID = numpy.sort(Node_ID)
                    
                    for Node_Index in range(len(active_nodes_server)):
                        if socket.gethostname()[0:4] == 'node':
                            active_nodes_server[Node_Index] = "node"+str(Node_ID[Node_Index])
                    print "active_nodes_server",active_nodes_server
                
                if mpi4py_rank == 0:
                    print "Wait untile "+DAS_Output_Path+"/ppserver_log.txt written....."
                    ppserver_lines_length = 0
                    while ppserver_lines_length != len(active_nodes_server) * PP_Servers_Per_Node:
                        print "ppserver_lines_length, len(active_nodes_server)*PP_Servers_Per_Node",ppserver_lines_length, len(active_nodes_server)*PP_Servers_Per_Node
                        time.sleep(1)
                        ppserver_file = open(DAS_Output_Path+"ppserver_log.txt",'r')
                        ppserver_lines = ppserver_file.readlines()
                        ppserver_file.close()
                        #print ppserver_lines
                        ppserver_lines_length = len(ppserver_lines)
                    
                        end = time.time()
                        
                        #print end-start
                        
                        if (end-start) > 30:
                            print "ppserver not started, please restart................"
                            job_server_node_array = Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
                            subprocess.call(ppserver_CMD_String,shell=True)
                            start = time.time()
                    
                
                    Node_List_Prior = []
                    Port_List_Prior = []
                    for ppserver_line_index in range(len(ppserver_lines)):
                        Node_List_Prior.append(string.split(ppserver_lines[ppserver_line_index])[3])
                        Port_List_Prior.append(string.split(ppserver_lines[ppserver_line_index])[4])
                    
                    # Get the Node Processor Number
                    Node_Count = numpy.zeros(len(Node_List_Prior))
                    nodelist = list(string.split(nodelist))
                    print "nodelist",nodelist
                    
                    #find all occurrences of an element in a list
                    def count(my_string, my_list):
                        indices = []
                        for idx, elem in enumerate(my_list):
                            if elem==my_string:
                                indices.append(idx)
                        return indices

                    for Node_List_Elem in Node_List_Prior:
                        print "Node_List_Elem",Node_List_Elem,nodelist.count(Node_List_Elem),count(Node_List_Elem,Node_List_Prior)
                        Node_Count[count(Node_List_Elem,Node_List_Prior)] = nodelist.count(Node_List_Elem)
                    
                    Node_Count_Index = numpy.argsort(Node_Count)[::-1]
                    
                    print "--------------------------------------Node_Count",Node_Count
                    print "--------------------------------------Node_List_Prior",Node_List_Prior
                    print "--------------------------------------Port_List_Prior",Port_List_Prior
                    Node_List = list(Node_List_Prior)
                    Port_List = list(Port_List_Prior)
                    for Node_Index in range(len(Node_Count_Index)):
                        Node_List[Node_Index] = Node_List_Prior[Node_Count_Index[Node_Index]]
                        Port_List[Node_Index] = Port_List_Prior[Node_Count_Index[Node_Index]]
                    print "--------------------------------------Node_List", Node_List
                    print "--------------------------------------Port_List", Port_List
                
                    start = time.time()
                    job_server_node_array = ['' for Node_List_Index in range(len(Node_List))]
                    
                    print "Check whether the ppserver is started on the node"
                    active_nodes_server_index = Node_List[0]
                    ppserver_port = Port_List[0]
                    ppservers_node=(active_nodes_server_index+":"+str(ppserver_port),)
                    print "ppservers_node",ppservers_node
                    job_server_node = pp.Server(ncpus=0, ppservers=ppservers_node,secret="123456",socket_timeout=999999999)
                    time.sleep(1)
                    print "Starting pp with",job_server_node.get_active_nodes(), "workers on node",active_nodes_server_index,ppservers_node
                    active_nodes = job_server_node.get_active_nodes()
                    del active_nodes['local']
                    while len(active_nodes.keys()) == 0:
                        job_server_node = pp.Server(ncpus=0, ppservers=ppservers_node,secret="123456",socket_timeout=999999999)
                        time.sleep(1)
                        print "Starting pp with",job_server_node.get_ncpus(),job_server_node.get_active_nodes(), "workers on node",active_nodes_server_index,ppservers_node
                        active_nodes = job_server_node.get_active_nodes()
                        del active_nodes['local']
                        end = time.time()
                        if (end-start) > 60:
                            print "Failed to start ppserver"
                            subprocess.call("killall -9 -q -w psilogger &> /dev/null",shell=True)
                            os.abort()
                    
                    print "Prepare the PP"                  
                    
                    #num_processors = multiprocessing.cpu_count()
                    print "Spawn the PP Workers using Multithread"
                    executor = futures.ThreadPoolExecutor(max_workers=len(Node_List))
                    for Node_List_Index in range(len(Node_List)):
                        #Thread Version
                        job_server_node_array[Node_List_Index] = executor.submit(PP_Worker,Node_List[Node_List_Index], Port_List[Node_List_Index]).result()
                        # Serial Version
                        #job_server_node_array[Node_List_Index] = PP_Worker(Node_List[Node_List_Index], Port_List[Node_List_Index])
                    
                    print "job_server_node_array",job_server_node_array
                
                
                    if len(job_server_node_array) == 0:
                        print "len(job_server_node_array) == 0................"
                        print "No enough node server"
                        subprocess.call("killall -9 -q -w psilogger &> /dev/null",shell=True)
                        os.abort()
                
                if Def_PP == 2:
                    # Wait for all ranks
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
    
    if mpi4py_rank == 0:
        print "ppservers are now started, now submit python jobs....."
        print "**************************len(job_server_node_array)",len(job_server_node_array)
        print "*******************************************************************"
    return job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node


def Stop_ppserver(mpi4py_rank, Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node):
    
    print "Stop ppserver to realease mpi for CLM........"
    for job_server_node in job_server_node_array:
        job_server_node.destroy()
    
    subprocess.call("killall -9 -q -w psilogger ppserver.py &> /dev/null",shell=True)

    time.sleep(1)

    job_server_node_array = []
    
    return job_server_node_array

def Backup_Initial_Files(NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy,
                         NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy,
                         NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy,
                         NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy,
                         NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy):
    
    print "************************************* Make a Copy of Assimilation_2_Initial.nc"
    copyLargeFile(NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy)
    copyLargeFile(NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy)
    
            