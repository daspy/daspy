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
from mpi4py import MPI
import os, sys, time, datetime, calendar, random, math, gc, imp, subprocess, glob, signal, string, shutil, fnmatch, warnings, multiprocessing, socket, getpass, ctypes, platform, functools, copy
import numpy, scipy, scipy.stats, scipy.signal, netCDF4, scipy.ndimage

sys.path.append('SysModel')
sys.path.append('SysModel/CLM')
sys.path.append('Utilities')
sys.path.append('Utilities/Soil')
sys.path.append('Algorithm')
sys.path.append('Algorithm/DAS')
sys.path.append('ForcingData')

from Call_CLM_CESM import *

from numpy import min, max
from Read_Soil_Texture import *

import ParFor

import pp, pyper
#from IPython.parallel import depend, require, dependent


#********************************************************************Call Model************************************************************************

def Prepare_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                            CLM_File_Name_List, Parameter_Range_Soil, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, 
                          DAS_Data_Path, DasPy_Path, DAS_Output_Path, Forcing_File_Path_Array, dtime, Variable_Assimilation_Flag, Variable_List,\
                          Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, 
                          DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Irrigation_Scheduling_Flag,\
                          omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                          Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                          NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single):
    
    pyper = imp.load_source("pyper",DasPy_Path+"Utilities/pyper.py")
    DAS_Utilities = imp.load_source("pyper",DasPy_Path+"DAS_Utilities.py")
    ParFor = imp.load_source("ParFor",DasPy_Path+"ParFor.py")
    
    fndepdat_name = fnmatch.filter(CLM_File_Name_List,"fndep*")[0]
    fatmgrid_name = fnmatch.filter(CLM_File_Name_List,"griddata_*")[0]
    fatmlndfrc_name = fnmatch.filter(CLM_File_Name_List,"domain*")[0]
    fsurdat_name = fnmatch.filter(CLM_File_Name_List,"surfdata_*")[0]
    fglcmask_name = fnmatch.filter(CLM_File_Name_List,"glcmaskdata_*")[0]
    flndtopo_name = fnmatch.filter(CLM_File_Name_List,"topodata_*")[0]
    fsnowoptics_name = fnmatch.filter(CLM_File_Name_List,"snicar_optics_5bnd_c090915*")[0]
    fsnowaging_name = fnmatch.filter(CLM_File_Name_List,"snicar_drdt_bst_fit_60_c070416*")[0]
    fpftcon_name = fnmatch.filter(CLM_File_Name_List,"pft*")[0]
    domain_name = fnmatch.filter(CLM_File_Name_List,"domain*")[0]
    rdirc_name = fnmatch.filter(CLM_File_Name_List,"rdirc_0*")[0]
    popd_streams_name = fnmatch.filter(CLM_File_Name_List,"*clmforc.Li_2012_hdm*")[0]
    light_streams_name = fnmatch.filter(CLM_File_Name_List,"*clmforc.Li_2012_climo*")[0]
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')    
    Land_Mask_Data = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][:,:]
    Teta_Saturated = NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][:,:,:]
    
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')    
    Initial_SM_Noise = NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_SM_Noise'][:,:,:]
    Initial_ST_Noise = NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_ST_Noise'][:,:,:]
    NC_File_Out_Assimilation_2_Diagnostic.close()
    
    NC_File_Parameter_Space_Single = netCDF4.Dataset(NC_FileName_Parameter_Space_Single,'r')
    Parameter_ParFlow_Space_Single = NC_File_Parameter_Space_Single.variables['Parameter_ParFlow_Space_Single'][:,:,:]
    Parameter_Soil_Space_Single = NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][:,:,:]
    Parameter_Veg_Space_Single = NC_File_Parameter_Space_Single.variables['Parameter_Veg_Space_Single'][:,:]
    GaussRF_Array = NC_File_Parameter_Space_Single.variables['GaussRF_Array'][:,:]
    Sand_Ratio = NC_File_Parameter_Space_Single.variables['Sand_Ratio'][:,:,:]
    Clay_Ratio = NC_File_Parameter_Space_Single.variables['Clay_Ratio'][:,:,:]
    Organic_Ratio = NC_File_Parameter_Space_Single.variables['Organic_Ratio'][:,:,:]
    MONTHLY_LAI_Empirical = NC_File_Parameter_Space_Single.variables['MONTHLY_LAI_Empirical'][:,:]
    NC_File_Parameter_Space_Single.close()
        
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
    
    Run_Dir = Run_Dir_Array[Ens_Index]
    Forcing_File_Path = Forcing_File_Path_Array[Ens_Index]
    
    if Def_Print:
        print "Processing " + str(Ens_Index + 1) + 'th Ensemble under',Run_Dir
        
        
    if Def_First_Run == 1:
        if Def_Print:
            print "*******************************Copy Surfdata File**********************"
         
        #print DAS_Data_Path + "SysModel/CLM/tools/mksurfdata/surfdata_"+Row_Numbers_String+"x"+Col_Numbers_String+".nc",Run_Dir+"/"+"surfdata_"+Row_Numbers_String+"x"+Col_Numbers_String+".nc"
        #DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/tools/mkgriddata/griddata_"+Row_Numbers_String+"x"+Col_Numbers_String+".nc",Run_Dir+"/"+"griddata_"+Row_Numbers_String+"x"+Col_Numbers_String+".nc")
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/tools/"+fatmlndfrc_name,Run_Dir+fatmlndfrc_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name,Run_Dir+fsurdat_name)
        #DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/glc/cism/griddata/"+fglcmask_name,Run_Dir+fglcmask_name)
        #DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/tools/"+flndtopo_name,Run_Dir+flndtopo_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name,Run_Dir+fpftcon_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"+rdirc_name,Run_Dir+rdirc_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/snicardata/"+fsnowoptics_name,Run_Dir+fsnowoptics_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/snicardata/"+fsnowaging_name,Run_Dir+fsnowaging_name)
        DAS_Utilities.copyLargeFile(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/ndepdata/"+fndepdat_name,Run_Dir+fndepdat_name)
        
    
    mksurfdata_NC_FileName_In = Run_Dir+"/"+fsurdat_name
         
    pft_physiology_file_name = Run_Dir+"/"+fpftcon_name            
        
    if Ensemble_Number > 1 and Def_Par_Optimized:
        if Def_Print:
            print "******************************* Parameter Perturbation and Modify the surfdata.nc *******************************"         
                
        mksurfdata_NC_File_Orig = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name, 'r')
        
        mksurfdata_NC_File_In = netCDF4.Dataset(mksurfdata_NC_FileName_In, 'r+')
        
        Sand_Mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Clay_Mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
         
        Parameter_Soil_Space_Ensemble = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:])
          
        for Soil_Layer_Index in range(10):
            Sand_Mat[Soil_Layer_Index,:,:] = numpy.flipud(Parameter_Soil_Space_Ensemble[Ens_Index,0,::]) * Sand_Ratio[Soil_Layer_Index,:,:]
            numexpr_a = Sand_Mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[0,0]
            numexpr_c = numpy.where(numexpr_a < numexpr_b)
            Sand_Mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[0,0]
            numexpr_a = Sand_Mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[1,0]
            numexpr_c = numpy.where(numexpr_a > numexpr_b)
            Sand_Mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[1,0]
            Clay_Mat[Soil_Layer_Index,:,:] = numpy.flipud(Parameter_Soil_Space_Ensemble[Ens_Index,0+Soil_Texture_Layer_Opt_Num,::]) * Clay_Ratio[Soil_Layer_Index,:,:]
            numexpr_a = Clay_Mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[0,1]
            numexpr_c = numpy.where(numexpr_a < numexpr_b)
            Clay_Mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[0,1]
            numexpr_a = Clay_Mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[1,1]
            numexpr_c = numpy.where(numexpr_a > numexpr_b)
            Clay_Mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[1,1]
              
        mksurfdata_NC_File_In.variables["PCT_SAND"][:,:,:] = Sand_Mat
        mksurfdata_NC_File_In.variables["PCT_CLAY"][:,:,:] = Clay_Mat
        del Sand_Mat,Clay_Mat
          
        Organic_mat = numpy.zeros((10,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        for Soil_Layer_Index in range(8):
            Organic_mat[Soil_Layer_Index,:,:] = numpy.flipud(Parameter_Soil_Space_Ensemble[Ens_Index,0+2*Soil_Texture_Layer_Opt_Num,::]) * Organic_Ratio[Soil_Layer_Index,:,:]
            numexpr_a = Organic_mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[0,2]
            numexpr_c = numpy.where(numexpr_a < numexpr_b)
            Organic_mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[0,2]
            numexpr_a = Organic_mat[Soil_Layer_Index,:,:]
            numexpr_b = Parameter_Range_Soil[1,2]
            numexpr_c = numpy.where(numexpr_a > numexpr_b)
            Organic_mat[Soil_Layer_Index,:,:][numexpr_c] = Parameter_Range_Soil[1,2]
        mksurfdata_NC_File_In.variables["ORGANIC"][:,:,:] = Organic_mat
        del Organic_mat
          
        del Parameter_Soil_Space_Ensemble
        del numexpr_a, numexpr_b, numexpr_c
        
        if (numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) >= 1):
            mksurfdata_NC_File_In.variables["MONTHLY_LAI"][:,:,:,:] = \
            numpy.multiply(mksurfdata_NC_File_Orig.variables["MONTHLY_LAI"][:,:,:,:], numpy.flipud(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][Ens_Index,0,::]))
            mksurfdata_NC_File_In.variables["MONTHLY_SAI"][:,:,:,:] = \
            numpy.multiply(mksurfdata_NC_File_Orig.variables["MONTHLY_SAI"][:,:,:,:], numpy.flipud(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][Ens_Index,1,::]))
            mksurfdata_NC_File_In.variables["MONTHLY_HEIGHT_TOP"][:,:,:,:] = \
            numpy.multiply(mksurfdata_NC_File_Orig.variables["MONTHLY_HEIGHT_TOP"][:,:,:,:], numpy.flipud(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][Ens_Index,2,::]))            
            
            del MONTHLY_LAI_Empirical
            
        mksurfdata_NC_File_In.sync()
        mksurfdata_NC_File_In.close()
        mksurfdata_NC_File_Orig.close()
                
        if Initial_Perturbation and ((numpy.sum(Initial_Perturbation_SM_Flag) or numpy.sum(Initial_Perturbation_ST_Flag)) or Def_First_Run) and Def_Initial:
            if Def_Print:
                print "**************************Initial Value Perturbation, Open Initial File:", Run_Dir+"/"+ finidat_initial_CLM
            CLM_Initial_File = netCDF4.Dataset(Run_Dir+"/"+ finidat_initial_CLM, "r+")
            column_len = len(CLM_Initial_File.dimensions['column'])
             
            for Soil_Layer_Index in range(Soil_Layer_Num-5):
                if Initial_Perturbation_SM_Flag[Soil_Layer_Index]:
                    #print numpy.shape(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index]),numpy.shape(Initial_SM_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten())
                    CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index] += numpy.repeat(Initial_SM_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten(),column_len/(Row_Numbers*Col_Numbers))  * (Soil_Thickness[Soil_Layer_Index] * Density_of_liquid_water)
                    #print numpy.max(Initial_SM_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten()  * (Soil_Thickness[Soil_Layer_Index] * Density_of_liquid_water)),numpy.min(Initial_SM_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten()  * (Soil_Thickness[Soil_Layer_Index] * Density_of_liquid_water))
                    CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index][numpy.where(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index] < 0.0)] = 0.2 * (Soil_Thickness[Soil_Layer_Index] * Density_of_liquid_water)
             
            CLM_Initial_File.variables["T_GRND"][:] += numpy.repeat(Initial_ST_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten(),column_len/(Row_Numbers*Col_Numbers))
            for Soil_Layer_Index in range(Soil_Layer_Num):
                if Initial_Perturbation_ST_Flag[Soil_Layer_Index]:
                    CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index][:] += Initial_ST_Noise[Ens_Index,:,:][numpy.where(Land_Mask_Data != NAvalue)].flatten()
            CLM_Initial_File.sync()
            CLM_Initial_File.close()
    
    NC_File_Out_Assimilation_2_Parameter.sync()
    NC_File_Out_Assimilation_2_Parameter.close()
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Constant.close()
    #NC_File_Assimilation_2_Parameter_Ens.sync()
    #NC_File_Assimilation_2_Parameter_Ens.close()
    
    del Initial_SM_Noise,Initial_ST_Noise,Land_Mask_Data, Teta_Saturated
    del Parameter_ParFlow_Space_Single,Parameter_Soil_Space_Single,Parameter_Veg_Space_Single
    del Sand_Ratio, Clay_Ratio, Organic_Ratio, GaussRF_Array
    
    gc.collect()
    del gc.garbage[:]
    
    return 0


def Call_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                            CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                          Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Irrigation_Scheduling_Flag,\
                          Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                          Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                          NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null):
    
    if Def_PP == 2:
        fcomm_rank = mpi4py_comm_split.Get_rank()
        fcomm = mpi4py_comm_split
        fcomm_null = mpi4py_null
    else:
        fcomm_rank = 0
        fcomm = 0
        fcomm_null = 0
    
    #if Def_Print:
    #    print "fcomm,fcomm_null,fcomm_rank",fcomm,fcomm_null,fcomm_rank
        
    fndepdat_name = fnmatch.filter(CLM_File_Name_List,"fndep*")[0]
    fatmgrid_name = fnmatch.filter(CLM_File_Name_List,"griddata_*")[0]
    fatmlndfrc_name = fnmatch.filter(CLM_File_Name_List,"domain*")[0]
    fsurdat_name = fnmatch.filter(CLM_File_Name_List,"surfdata_*")[0]
    fglcmask_name = fnmatch.filter(CLM_File_Name_List,"glcmaskdata_*")[0]
    flndtopo_name = fnmatch.filter(CLM_File_Name_List,"topodata_*")[0]
    fsnowoptics_name = fnmatch.filter(CLM_File_Name_List,"snicar_optics_5bnd_c090915*")[0]
    fsnowaging_name = fnmatch.filter(CLM_File_Name_List,"snicar_drdt_bst_fit_60_c070416*")[0]
    fpftcon_name = fnmatch.filter(CLM_File_Name_List,"pft*")[0]
    domain_name = fnmatch.filter(CLM_File_Name_List,"domain*")[0]
    rdirc_name = fnmatch.filter(CLM_File_Name_List,"rdirc_0*")[0]
    popd_streams_name = fnmatch.filter(CLM_File_Name_List,"*clmforc.Li_2012_hdm*")[0]
    light_streams_name = fnmatch.filter(CLM_File_Name_List,"*clmforc.Li_2012_climo*")[0]
    
    Run_Dir = Run_Dir_Array[Ens_Index]
    Forcing_File_Path = Forcing_File_Path_Array[Ens_Index]
    
    if Def_Print:
        print "Processing " + str(Ens_Index + 1) + 'th Ensemble under',Run_Dir
    
    if Def_Initial and os.path.isfile(Run_Dir + finidat_initial_CLM):
        finidat_name = finidat_initial_CLM
    else:
        finidat_name = ""
    if Def_Print:
        print "Initial File Information", Def_Initial, finidat_name
    
    if Def_Print:
        print "---------------------------------------------- Run CLM Model ---------------------------"
    
    Run_CLM(Model_Driver, Def_PP, Do_DA_Flag, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized,   Def_Region, Def_Initial, \
            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, 
            Run_Dir, Run_Dir_Array, Model_Path, CLM_Flag, domain_name, rdirc_name, \
            fatmgrid_name, fatmlndfrc_name, fndepdat_name, fsurdat_name, fsnowoptics_name, fsnowaging_name, fglcmask_name, finidat_name, flndtopo_name, \
            fpftcon_name, popd_streams_name, light_streams_name, N_Steps, Row_Numbers_String, Col_Numbers_String, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, \
            Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init,  DasPy_Path, DAS_Data_Path, Forcing_File_Path, Forcing_File_Path_Array, 
            dtime, ntasks_CLM, rootpe_CLM, nthreads_CLM, Ensemble_Number, num_processors,
            COUP_OAS_PFL, CESM_Init_Flag, fcomm, fcomm_null, fcomm_rank)
    
    stop_tod_string = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
    history_file_name = Region_Name + '.clm2.h0.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
    restart_file_name = Region_Name + '.clm2.r.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
        
    gc.collect()
    del gc.garbage[:]
    
    if Def_Print:
        print history_file_name, restart_file_name
    
    return 0



#*******************************************************************Obervation********************************************************************************************

def Observation_Blocks(Observation_Matrix_Index, Def_PP, Def_CESM_Multi_Instance, Def_Region, Def_ReBEL, Def_Localization, DasPy_Path, DAS_Depends_Path, DAS_Output_Path, Region_Name, Num_Local_Obs,\
                       Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict,  Call_Gstat_Flag, Plot_Analysis, \
                       Write_DA_File_Flag, Use_Mask_Flag, Mask_File, Forcing_File_Path, dtime, Observation_Path, Observation_File_Name, \
                       DAS_Data_Path, Grid_Resolution_CEA, Grid_Resolution_GEO, NAvalue, Variable_List, Variable_Assimilation_Flag, \
                       SensorType, SensorVariable, SensorQuantity, SensorResolution, Variable_ID, QC_ID, PDAF_Assim_Framework, PDAF_Filter_Type, \
                       mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Datetime_Start, Datetime_Stop, Datetime_Stop_Init, \
                       Stop_Year, Stop_Month, Stop_Day, Def_Print, Observation_Bias_Range, Observation_Bias_Range_STD, Observation_Bias_Initialization_Flag, plt, cm, colors, *vartuple):
    
    octave = vartuple[0]
    r = vartuple[1]
    
    # Matrix to store the observation supplementary information
    Observation_Misc = numpy.zeros((10, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    print "######################### Open Observation File: ", Observation_File_Name
    print "\n"    
    
    if numpy.sum(Variable_Assimilation_Flag) > 0:
        Observation_File = netCDF4.Dataset(Observation_File_Name, "r")
        #----------------- Get the Observation Matrix Dimension
        #print Observation_File.dimensions["lat"]
        Observation_NLats = len(Observation_File.dimensions["lat"])
        Observation_NLons = len(Observation_File.dimensions["lon"])
        Observation_Latitude = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_Longitude = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_View_Zenith_Angle = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_View_Time = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_Latitude = Observation_File.variables["CEA_Y"][::]
        Observation_Longitude = Observation_File.variables["CEA_X"][::]
        Observation_Latitude = numpy.flipud(Observation_Latitude)
        Observation_Matrix = Observation_File.variables[Variable_ID][::]
        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature" or\
            Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Sensible_Heat":
            if SensorType == "Terra" or SensorType == "Aqua":
                Observation_View_Zenith_Angle = numpy.abs(numpy.flipud(Observation_File.variables["View_angl"][::]) - 65.0)
                numexpr_a = Observation_View_Zenith_Angle
                numexpr_b = 255.0-65.0
                numexpr_c = numpy.where(numexpr_a == numexpr_b)
                Observation_View_Zenith_Angle[numexpr_c] = 0.0
            else:
                Observation_View_Zenith_Angle = numpy.abs(numpy.flipud(Observation_File.variables["View_angl"][::]))
            if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature":
                Observation_View_Time = numpy.abs(numpy.flipud(Observation_File.variables["View_time"][::]))    # MODIS is local time, but changed in MRT_Py to UTC
            
    else:
        Observation_NLats = Row_Numbers
        Observation_NLons = Col_Numbers
        Observation_Latitude = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_Longitude = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_Matrix = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_View_Zenith_Angle = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
        Observation_View_Time = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
    
    
    if Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Soil_Moisture":
        if SensorType == "AMSR_E" and SensorQuantity != "K":
            if Variable_ID == "A_Soil_Moisture" or Variable_ID == "D_Soil_Moisture":
                AMSR_E_QA = Observation_File.variables[QC_ID][::]
                AMSR_E_QA_Index = numpy.zeros((Observation_NLats, Observation_NLons), dtype=numpy.bool)
                for Lat_Index in range(Observation_NLats):
                    for Long_Index in range(Observation_NLons):                
                        AMSR_E_QA_String = numpy.binary_repr(int(AMSR_E_QA[Lat_Index, Long_Index]), width=16)
                        #print AMSR_E_QA[Lat_Index,Long_Index],AMSR_E_QA_String[6],AMSR_E_QA_String[-10]
                        if AMSR_E_QA_String[-10] == "0":
                            AMSR_E_QA_Index[Lat_Index, Long_Index] = True
                        else:
                            AMSR_E_QA_Index[Lat_Index, Long_Index] = False
                    
                Observation_Matrix = Observation_Matrix / 1000.0
                numexpr_a = Observation_Matrix
                numexpr_b = 0.0
                numexpr_c = numpy.where(numexpr_a < numexpr_b)
                Observation_Matrix[numexpr_c] = NAvalue
                numexpr_a = Observation_Matrix
                numexpr_b = 0.8
                numexpr_c = numpy.where(numexpr_a >= numexpr_b)
                Observation_Matrix[numexpr_c] = NAvalue
            Grid_Resolution_GEO_Global = (179.999999415 + 179.999995782) / 1383.0
            mksrf_edgew_temp = mksrf_edgew + Grid_Resolution_GEO[0] / 2.0
            mksrf_edgen_temp = mksrf_edgen - Grid_Resolution_GEO[1] / 2.0
            Corner_Row_Index = int(numpy.floor((86.716744081 - mksrf_edgen_temp) / Grid_Resolution_GEO_Global)) - 5
            Corner_Col_Index = int(numpy.floor((mksrf_edgew_temp + 179.999999415) / Grid_Resolution_GEO_Global)) - 5
            Observation_Latitude = Observation_Latitude[Corner_Row_Index:(Corner_Row_Index + numpy.ceil(Row_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO[1])) + 10), \
                                                                                  Corner_Col_Index:(Corner_Col_Index + numpy.ceil(Col_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO[0])) + 10)]
            Observation_Longitude = Observation_Longitude[Corner_Row_Index:(Corner_Row_Index + numpy.ceil(Row_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO)) + 10), \
                                                                                  Corner_Col_Index:(Corner_Col_Index + numpy.ceil(Col_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO[0])) + 10)]
            Observation_Matrix = numpy.flipud(numpy.flipud(Observation_Matrix)[Corner_Row_Index:(Corner_Row_Index + numpy.ceil(Row_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO[1])) + 10), \
                                                                                  Corner_Col_Index:(Corner_Col_Index + numpy.ceil(Col_Numbers / (Grid_Resolution_GEO_Global / Grid_Resolution_GEO[0])) + 10)])
            Observation_Matrix = numpy.flipud(Observation_Matrix)
            
        elif SensorType == "SMOS" and SensorQuantity == "K":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "SMOS" and SensorQuantity == "m3/m3":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "ASCAT" and SensorQuantity == "m3/m3":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "ASAR" and SensorQuantity == "DB":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "PALSAR" and SensorQuantity == "DB":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif (SensorType == "Terra" or SensorType == "Aqua") and SensorQuantity == "m3/m3":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "ECV_SM" and SensorQuantity == "m3/m3":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        elif SensorType == "COSMOS" and SensorQuantity == "Neutron_Count":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
            
            if Def_Region == 3 or Def_Region == 8:
                Observation_Misc[0,:,:] = numpy.abs(numpy.flipud(Observation_File.variables["bd"][::]))
                Observation_Misc[1,:,:] = numpy.abs(numpy.flipud(Observation_File.variables["lw"][::]))
                Observation_Misc[2,:,:] = numpy.abs(numpy.flipud(Observation_File.variables["Ncosmic"][::]))
            
        elif SensorType == "InSitu":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        else:
            Observation_Matrix = numpy.flipud(Observation_Matrix)
    
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable)] and SensorVariable == "Surface_Temperature":
        if SensorType == "Terra" or SensorType == "Aqua":
            LST_QA = Observation_File.variables[QC_ID][::]
            LST_QA_Index = numpy.zeros((Observation_NLats, Observation_NLons), dtype=numpy.bool)
            for Lat_Index in range(Observation_NLats):
                for Long_Index in range(Observation_NLons):
                    LST_QA_String = numpy.binary_repr(int(LST_QA[Lat_Index, Long_Index]), width=8)
                    if LST_QA_String[6:8] == "00":
                        LST_QA_Index[Lat_Index, Long_Index] = True
                    else:
                        LST_QA_Index[Lat_Index, Long_Index] = False
            
            numexpr_a = Observation_Matrix
            numexpr_b = 65535
            numexpr_c = numpy.where(numexpr_a > numexpr_b)
            Observation_Matrix[numexpr_c] = NAvalue / 0.02
            numexpr_a = Observation_Matrix
            numexpr_b = 7500
            numexpr_c = numpy.where(numexpr_a < numexpr_b)
            Observation_Matrix[numexpr_c] = NAvalue / 0.02
            Observation_Matrix[~LST_QA_Index] = NAvalue / 0.02
            # Because there are different observation time for one catchment, we need 1 most  common time clock.
            Observation_Matrix[numpy.where(Observation_View_Time == -9999)] = NAvalue / 0.02
            Observation_Matrix = numpy.flipud(Observation_Matrix * 0.02)
            
        elif SensorType == "InSitu":
            Observation_Matrix = numpy.flipud(Observation_Matrix)
        
    else:
        Observation_Matrix = numpy.flipud(Observation_Matrix)
    
    # Make sure the Observation Bias is perturbed only once
    Observation_Bias_Initialization_Flag[:,:,:] += 1
    
    # Observation Corner Location
    Observation_X_Left = float(Observation_Longitude.min())
    Observation_X_Right = float(Observation_Longitude.max())
    Observation_Y_Lower = float(Observation_Latitude.min())
    Observation_Y_Upper = float(Observation_Latitude.max())
    
    
    Observation_Longitude = Observation_Longitude - Observation_X_Left
    Observation_Latitude = Observation_Latitude - Observation_Y_Lower
    Observation_X_Right = Observation_X_Right - Observation_X_Left
    Observation_X_Left = Observation_X_Left - Observation_X_Left    
    Observation_Y_Upper = Observation_Y_Upper - Observation_Y_Lower
    Observation_Y_Lower = Observation_Y_Lower- Observation_Y_Lower
    
    print "--------------Observation_X_Left,Observation_X_Right,Observation_Y_Lower,Observation_Y_Upper-----------------"
    print Observation_X_Left,Observation_X_Right,Observation_Y_Lower,Observation_Y_Upper
    
    r.assign("Observation_X_Left", Observation_X_Left)
    r.assign("Observation_X_Right", Observation_X_Right)
    r.assign("Observation_Y_Lower", Observation_Y_Lower)
    r.assign("Observation_Y_Upper", Observation_Y_Upper)
    
    Observation_Variance = numpy.zeros((Observation_NLats, Observation_NLons),dtype=numpy.float32)
    print "****************************** Pre-Defined Observation Variance!"
    Observation_Matrix[numpy.isnan(Observation_Matrix)] = NAvalue
    Observation_Matrix_None_NA_Index = numpy.where(Observation_Matrix != NAvalue)
    print "numpy.size(Observation_Matrix_None_NA_Index) / 2",numpy.size(Observation_Matrix_None_NA_Index) / 2
    print "numpy.min(Observation_Matrix[Observation_Matrix_None_NA_Index]),numpy.max(Observation_Matrix[Observation_Matrix_None_NA_Index])"
    print numpy.min(Observation_Matrix[Observation_Matrix_None_NA_Index]),numpy.max(Observation_Matrix[Observation_Matrix_None_NA_Index])
    
    if SensorVariable == "Soil_Moisture":
        if SensorQuantity == "K":
            Observation_Variance[Observation_Matrix_None_NA_Index] = numpy.square(Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.01)
            Observation_Variance[Observation_Matrix_None_NA_Index] = 4.0
        elif SensorQuantity == "DB":
            Observation_Variance[Observation_Matrix_None_NA_Index] = numpy.square(Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.05)
        elif SensorQuantity == "m3/m3":
            Observation_Variance[Observation_Matrix_None_NA_Index] = numpy.square(Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.05)
            if SensorType == "ECV_SM" or SensorType == "ASCAT":
                # To avoid the missing values of variance
                Variance_Max = numpy.max(numpy.asarray(Observation_File.variables[Variable_ID+"_noise"][::]))
                print "-------------- Variance_Max",Variance_Max
                if Variance_Max > 0.0:
                    Observation_Variance = numpy.square(numpy.flipud(Observation_File.variables[Variable_ID+"_noise"][::]))
                else:
                    Observation_Variance[Observation_Matrix_None_NA_Index] = 0.0016
            else:
                if Def_Region == -1:
                    Observation_Variance[Observation_Matrix_None_NA_Index] = 0.0016
                elif Def_Region == 0:
                    Observation_Variance[Observation_Matrix_None_NA_Index] = 0.0009
                else:
                    Observation_Variance[Observation_Matrix_None_NA_Index] = 0.0016
        elif SensorQuantity == "Neutron_Count":
            print "-----------------------------Assign COSMOS Observation Error"
            Observation_Variance[Observation_Matrix_None_NA_Index] = Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.25
            if Def_Region == -2:
                Observation_Variance[Observation_Matrix_None_NA_Index] = Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.25
        
    elif SensorVariable == "Surface_Temperature" or SensorVariable == "Vegetation_Temperature":
        Observation_Variance[Observation_Matrix_None_NA_Index] = numpy.square(Observation_Matrix[Observation_Matrix_None_NA_Index] * 0.01)
        Observation_Variance[Observation_Matrix_None_NA_Index] = 1.0
        #if SensorType == "Terra" or SensorType == "Aqua":
        #    Observation_Variance[Observation_Matrix_None_NA_Index] = 0.25   # http://landval.gsfc.nasa.gov/ProductStatus.php?ProductID=MOD11
       
    if numpy.sum(Variable_Assimilation_Flag) > 0:
        Observation_File.close()
    
    Grid_Resolution_CEA_Local = (Observation_X_Right-Observation_X_Left)/Col_Numbers
    
    Observation_Corelation_Par = numpy.zeros((5, 2))
    Observation_Corelation_Par[0, 0] = 2 # exponential Model
    #Observation_Corelation_Par[0, 0] = 12 # Gaspari_Cohn Model
    Observation_Corelation_Par[1, 0] = 0.0
    Observation_Corelation_Par[2, 0] = 1.0
    Observation_Corelation_Par[3, 0] = 10.0*Grid_Resolution_CEA_Local
    Observation_Corelation_Par[4, 0] = 1.0    
                
    print "Observation Variance is:", numpy.mean(Observation_Variance[Observation_Matrix_None_NA_Index])
        
    if Write_DA_File_Flag:
        numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Observation_Corelation_Par.txt", Observation_Corelation_Par)
        numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/R.txt", Observation_Variance)
        numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/Observation_Matrix.txt",Observation_Matrix)
    
               

    print " -------------Define Observation Collection NC Files"
    NC_FileName_Observation = DAS_Output_Path+"Analysis/"+Region_Name+"/Observation_"+str(Observation_Matrix_Index+1)+".nc"
    if os.path.exists(NC_FileName_Observation):
        os.remove(NC_FileName_Observation)
        
    if Def_Print:
        print 'Write NetCDF File:',NC_FileName_Observation
        
    NC_File_Observation = netCDF4.Dataset(NC_FileName_Observation, 'w', diskless=True, persist=True, format='NETCDF4')
    # Dim the dimensions of NetCDF
    NC_File_Observation.createDimension('lon', Observation_NLons)
    NC_File_Observation.createDimension('lat', Observation_NLats)
    NC_File_Observation.createDimension('misc', 10)
        
    NC_File_Observation.createVariable('Observation_Variance','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_Variance'][:,:] = Observation_Variance
    
    NC_File_Observation.createVariable('Observation_Latitude','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_Latitude'][:,:] = Observation_Latitude
    
    NC_File_Observation.createVariable('Observation_Longitude','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_Longitude'][:,:] = Observation_Longitude
    
    NC_File_Observation.createVariable('Observation_Matrix','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_Matrix'][:,:] = Observation_Matrix
    
    NC_File_Observation.createVariable('Observation_View_Zenith_Angle','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_View_Zenith_Angle'][:,:] = Observation_View_Zenith_Angle
    
    NC_File_Observation.createVariable('Observation_View_Time','f4',('lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_View_Time'][:,:] = Observation_View_Time
    
    NC_File_Observation.createVariable('Observation_Misc','f4',('misc','lat','lon',),zlib=True)
    NC_File_Observation.variables['Observation_Misc'][:,:,:] = Observation_Misc
    
    
    NC_File_Observation.sync()
    NC_File_Observation.close()
    
    numexpr_a = []
    numexpr_b = []
    numexpr_c = []
    del Observation_Variance, Observation_Latitude, Observation_Longitude, Observation_Matrix, Observation_View_Zenith_Angle, Observation_View_Time
    
    r('gc(TRUE)')
    
    gc.collect()
    del gc.garbage[:]
    
    return Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, \
            Observation_Corelation_Par, Observation_Bias_Range, Observation_Bias_Range_STD, Observation_Bias_Initialization_Flag


def Read_History_File(Ens_Index, Model_Driver, Def_First_Run, Def_Region, Def_PP, Dim_CLM_State, Row_Numbers, Col_Numbers, column_len, pft_len, \
                      Datetime_Start, Datetime_Stop, Start_Year, Stop_Month, Region_Name, Run_Dir_Home, history_file_name, Ensemble_Number, Constant_File_Name_Header, 
                      Variable_Assimilation_Flag, Variable_List, Additive_Noise_SM, Additive_Noise_ST, N0, nlyr,
                     Mean_Dir, Two_Step_Bias_Estimation_Active, Feedback_Assim, Soil_Texture_Layer_Opt_Num,                 
                     Grid_Resolution_CEA, Grid_Resolution_GEO, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,
                     SensorType, SensorVariable, Variable_ID, Analysis_Variable_Name, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Def_ParFor, DAS_Depends_Path, DAS_Data_Path, Def_Print,
                     Irrig_Scheduling, Density_of_liquid_water, Density_of_ice, NAvalue, CLM_NA, omp_get_num_procs_ParFor, DAS_Output_Path, \
                     NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, \
                     NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, finidat_name_string, Plot_Analysis, DasPy_Path, *vartuple):
    
#    numpy.savez("Read_History_File.npz",Ens_Index, Dim_CLM_State, Row_Numbers, Col_Numbers, Run_Dir_Home, history_file_name, Ensemble_Number, Constant_File_Name_Header, Bias_Forecast, Additive_Noise_SM, \
#                     Soil_Moisture_DA_Flag, Surface_Temperature_DA_Flag, Vegetation_Temperature_DA_Flag, Albedo_DA_Flag, Snow_Depth_DA_Flag, Snow_Cover_Fraction_DA_Flag,Canopy_Water_DA_Flag,
#                     Snow_Water_Equivalent_DA_Flag, Crop_Planting_DA_Flag, Crop_Harvest_DA_Flag, LAI_DA_Flag, Latent_Heat_DA_Flag, Latent_Heat_Daily_DA_Flag, Sensible_Heat_DA_Flag, Water_Storage_DA_Flag, Water_Table_Depth_DA_Flag, Emissivity_DA_Flag, \
#                     SensorType, SensorVariable, Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point, Analysis_Variable_Name, Soil_Layer_Num, numrad, Def_ParFor, DAS_Depends_Path, \
#                     Irrigation_Rate_DA_Flag, Irrig_Scheduling, Irrigation_Scheduling_Flag, Model_Driver, Density_of_liquid_water, Density_of_ice, NAvalue, CLM_NA, omp_get_num_procs_ParFor)
    
    COSMIC_Py = imp.load_source("COSMIC_Py",DasPy_Path+"ObsModel/COSMOS/COSMIC_Py.py")
    COSMIC = imp.load_dynamic("COSMIC",DasPy_Path+"ObsModel/COSMOS/COSMIC.so")
    Clumping_Index = imp.load_source("Clumping_Index",DasPy_Path+"ObsModel/LST/Clumping_Index.py")
    ParFor = imp.load_source("ParFor",DasPy_Path+"ParFor.py")
    
    CLM_Initial_File = netCDF4.Dataset(finidat_name_string, 'r')
    cols1d_ixy = CLM_Initial_File.variables['cols1d_ixy'][:]
    cols1d_jxy = CLM_Initial_File.variables['cols1d_jxy'][:]
    cols1d_ityplun = CLM_Initial_File.variables['cols1d_ityplun'][:]
    
    #numpy.savetxt('cols1d_ixy',cols1d_ixy)
    #numpy.savetxt('cols1d_jxy',cols1d_jxy)
    pfts1d_ixy = CLM_Initial_File.variables['pfts1d_ixy'][:]
    pfts1d_jxy = CLM_Initial_File.variables['pfts1d_jxy'][:]
    pfts1d_itypveg = CLM_Initial_File.variables['pfts1d_itypveg'][:]
    pfts1d_ci = CLM_Initial_File.variables['pfts1d_ci'][:]
    pfts1d_ityplun = CLM_Initial_File.variables['pfts1d_ityplun'][:]
    
    CLM_Initial_File.close()
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    PCT_Veg = NC_File_Out_Assimilation_2_Constant.variables['PCT_Veg'][:,:]
    STD_ELEV = NC_File_Out_Assimilation_2_Constant.variables['STD_ELEV'][:,:]
    DEM_Data = NC_File_Out_Assimilation_2_Constant.variables['DEM_Data'][:,:]
    Land_Mask_Data = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][:,:]
    Bulk_Density_Top_Region = NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Top_Region'][:,:]
    Bulk_Density_Sub_Region = NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Sub_Region'][:,:]
    PFT_Dominant_Index = NC_File_Out_Assimilation_2_Constant.variables['PFT_Dominant_Index'][:,:]
    Bare_Grid_Index = numpy.where(PFT_Dominant_Index == 0)
    
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    Soil_Layer_Thickness_Ratio_Moisture = NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Moisture'][:, :, :]
    Soil_Layer_Thickness_Ratio_Moisture_Sum = numpy.sum(Soil_Layer_Thickness_Ratio_Moisture,axis=0)
    for Soil_Layer_Index in range(Soil_Layer_Num):
        Soil_Layer_Thickness_Ratio_Moisture[Soil_Layer_Index,:,:] = Soil_Layer_Thickness_Ratio_Moisture[Soil_Layer_Index,:,:]/Soil_Layer_Thickness_Ratio_Moisture_Sum
    
    PCT_LAKE = NC_File_Out_Assimilation_2_Constant.variables['PCT_LAKE'][:,:]
    PCT_LAKE_Index = numpy.where(PCT_LAKE > 50)
        
    NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
    Parameter_Soil_Space_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:],axis=0)
    #Parameter_Soil_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][Ens_Index,:,:,:]
    NC_File_Out_Assimilation_2_Parameter.close()
    
    Sand_Top_Region = Parameter_Soil_Space_Ensemble[0,::]
    Sand_Sub_Region = Parameter_Soil_Space_Ensemble[1,::]
    Clay_Top_Region = Parameter_Soil_Space_Ensemble[0+Soil_Texture_Layer_Opt_Num,::]
    Clay_Sub_Region = Parameter_Soil_Space_Ensemble[1+Soil_Texture_Layer_Opt_Num,::]
    
    
    Prop_Grid_Array_Sys_Sub = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Soil_Moisture_Ensemble_Mat_Sub = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Soil_Temperature_Ensemble_Mat_Sub = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_ROOTFR_Ensemble_Mat_Sub = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_ROOTR_Ensemble_Mat_Sub = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    CLM_Soil_Moisture = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32) # Units: m3/m3
    CLM_Soil_Temperature = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)  # soil-snow temperature Units: K
    CLM_ROOTFR = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32) # fraction of roots in each soil layer
    CLM_ROOTR = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32) # effective fraction of roots in each soil layer
    CLM_Ground_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_Onset_Freezing_Degree_Days_Counter = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: C degree-days
    CLM_Vegetation_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_Urban_Ground_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)  # Urban Ground Temperature Units: K
    CLM_Rural_Ground_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)  # Rural Ground Temperature Units: K
    CLM_Lake_Temperature = numpy.zeros((10,Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_LAI = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Albedo = numpy.zeros((numrad*2, Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units None
    CLM_Snow_Depth = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units Meters
    CLM_INT_SNOW = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units mm
    CLM_FH2OSFC = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # surface water mm
    CLM_Snow_Cover_Fraction = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units Meters
    CLM_Snow_Water = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_QIRRIG = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Soil_Moisture_Ratio = numpy.ones((Soil_Layer_Num, Row_Numbers,Col_Numbers),dtype=numpy.float32) # Ration Among the different columns
    CLM_Soil_Temperature_Ratio = numpy.ones((Soil_Layer_Num, Row_Numbers,Col_Numbers),dtype=numpy.float32) # Ration Among the different columns
    CLM_Soil_Moisture_Ratio_MultiScale = numpy.ones((Soil_Layer_Num, Row_Numbers,Col_Numbers),dtype=numpy.float32) # Ration Among the different columns
    CLM_HTOP = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_HBOT = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Water_Storage = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Water_Table_Depth = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Irrigation_Rate = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Latent_Heat_Daily = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: W/m^2
    CLM_Sensible_Heat = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: W/m^2
    CLM_Latent_Heat = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: W/m^2       
    CLM_Vcmax_SCOPE = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: umol co2/m**2/s
    CLM_m_SCOPE = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: 
    CLM_SLA_SCOPE = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: m2 g-1 C)
    
    CLM_2m_Air_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_Air_Pressure = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Shortwave_Radiation = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Longwave_Radiation = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Specific_Humdity = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Wind_Speed = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Z0MV = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_DISPLA = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    
    Analysis_Variable_Name = ['' for i in range(Dim_CLM_State)]
    Soil_Layer_Index_DA = 0
    
    # For soil moisture, soil temperature bias estimation, we need the mean model simulation, for latent and sensbile heat bias estimation (observation bias)
    # We use the ensemble simulation
    if Two_Step_Bias_Estimation_Active and (not Variable_Assimilation_Flag[Variable_List.index("Latent_Heat")] == 1) and \
        (not Variable_Assimilation_Flag[Variable_List.index("Sensible_Heat")] == 1):
        Run_Dir = Mean_Dir
    else:
        Run_Dir = Run_Dir_Home+"_Ens"+str(Ens_Index+1)
    
    if Def_Print:
        print "Open History File:", Run_Dir+"/"+history_file_name
    try:
        CLM_History_File = netCDF4.Dataset(Run_Dir+"/"+history_file_name, 'r')
    except:
        print Run_Dir+"/"+history_file_name,"does not exist!...\n"
        os.abort()
#                            if history_file_name[-20:-1] == finidat_name[-20:-1]:
#                                Constant_File_Name = finidat_name
#                            else:
    #Constant_File_Name = []
#                        if Def_First_Run == 1:
#                            Constant_File_Name = Run_Dir +"/"+ CLM_History_File.Time_constant_3Dvars_filename[2:]
    #Constant_File_Name = Run_Dir +"/"+ Region_Name+'.clm2.h0.'+Start_Year+'-'+Start_Month+'-'+Start_Day+'-'+str((Datetime_Start - Datetime_Start_Init).seconds)+'.nc'
    if Ensemble_Number == 1:
        Constant_File_Name = Run_Dir_Home +"/"+ Constant_File_Name_Header
    else:
        if Def_Region == 4 or Def_Region == 8 or Def_Region == 9:
            if Def_First_Run == 1:
                Constant_File_Name = Run_Dir+ Constant_File_Name_Header
            else:
                Mean_Dir = Run_Dir_Home+"_Ens_Mean"
                Constant_File_Name = Mean_Dir +"/"+ Constant_File_Name_Header
        else:
            Constant_File_Name = Run_Dir +"/"+ Constant_File_Name_Header
    
    CLM_Onset_Freezing_Degree_Days_Counter = numpy.flipud(CLM_History_File.variables['ONSET_FDD'][0, :, :])  # onset freezing degree days counter
    
    if Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")] == 1 or Feedback_Assim == 1:
        
        Analysis_Variable_Name[0] = "Soil Moisture"
        if Def_Print:
            print Analysis_Variable_Name[0]
        #------------------------------------------- Read the Soil Moisture Data
        
        if True:
            
            CLM_Ground_Temperature = numpy.flipud(CLM_History_File.variables['TG'][0, :, :])  # Ground Temperature Units: K
            CLM_Ground_Temperature[numpy.isnan(CLM_Ground_Temperature)] = 283.15
            
            CLM_Vegetation_Temperature = numpy.flipud(CLM_History_File.variables['TV'][0, :, :]) # Units: K
            CLM_Vegetation_Temperature[numpy.isnan(CLM_Vegetation_Temperature)] = 283.15
            
            CLM_Snow_Depth = numpy.flipud(CLM_History_File.variables['SNOWDP'][0, :, :])
            #CLM_Snow_Depth[numpy.ma.getmask(CLM_Snow_Depth)] = NAvalue
            CLM_Snow_Water = numpy.flipud(CLM_History_File.variables['H2OSNO'][0, :, :]) # Snow Water: Units: (mm)
            #CLM_Snow_Water[numpy.ma.getmask(CLM_Snow_Water)] = NAvalue
            CLM_2m_Air_Temperature = numpy.flipud(CLM_History_File.variables['TSA'][0, :, :])
            #CLM_2m_Air_Temperature[numpy.ma.getmask(CLM_2m_Air_Temperature)] = 273.15
            CLM_Air_Pressure = numpy.flipud(CLM_History_File.variables['PSurf'][0, :, :])
            
            for Soil_Layer_Index in range(Soil_Layer_Num-5):
                #CLM_Soil_Moisture[Soil_Layer_Index, :, :] = numpy.flipud(CLM_History_File.variables['H2OSOI'][0, Soil_Layer_Index, :, :]) 
                CLM_Soil_Moisture[Soil_Layer_Index, :, :] = numpy.flipud(CLM_History_File.variables['H2OSOI'][0, Soil_Layer_Index, :, :])
                CLM_Soil_Temperature[Soil_Layer_Index, :, :] = numpy.flipud(CLM_History_File.variables['TSOI'][0, Soil_Layer_Index, :, :])
                
                CLM_Soil_Moisture_Ratio[Soil_Layer_Index,:,:] = CLM_Soil_Moisture[Soil_Layer_Index, :, :] / CLM_Soil_Moisture[0, :, :]
                
                
                if Def_Print >= 2:
                    print "Check the Ourliers"
                NA_Flag = False
                numexpr_a = CLM_Soil_Temperature[Soil_Layer_Index,:,:]
                numexpr_b = CLM_NA
                numexpr_c = numpy.where(numexpr_a == numexpr_b)
                CLM_Soil_Moisture[Soil_Layer_Index,:,:][numexpr_c] = 0.2
                CLM_Soil_Moisture[Soil_Layer_Index,:,:][numpy.isnan(CLM_Soil_Moisture[Soil_Layer_Index,:,:])] = 0.2
                Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Soil_Moisture[Soil_Layer_Index,:,:],NA_Flag, 'Soil_Moisture',  Variable_Assimilation_Flag, Variable_List,
                            numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][Soil_Layer_Index,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][Soil_Layer_Index,:,:]), 
                            numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][Soil_Layer_Index,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][Soil_Layer_Index,:,:]), NAvalue)
                
                
                CLM_Soil_Moisture_Ensemble_Mat_Sub[Soil_Layer_Index,:,:] = CLM_Soil_Moisture[Soil_Layer_Index, :, :]
                
                numexpr_a = CLM_Soil_Temperature[Soil_Layer_Index,:,:]
                numexpr_b = CLM_NA
                numexpr_c = numpy.where(numexpr_a == numexpr_b)
                CLM_Soil_Temperature[Soil_Layer_Index,:,:][numexpr_c] = 283.15
                CLM_Soil_Temperature[Soil_Layer_Index,:,:][numpy.isnan(CLM_Soil_Temperature[Soil_Layer_Index,:,:])] = 283.15
                Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Soil_Temperature[Soil_Layer_Index, :, :],NA_Flag,'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
                CLM_Soil_Temperature_Ensemble_Mat_Sub[Soil_Layer_Index,:,:] = CLM_Soil_Temperature[Soil_Layer_Index, :, :]
            
            try:
                SensorType_Name = SensorType[SensorVariable.index("Soil_Moisture")]
                
                if SensorType_Name == 'ECV_SM' or SensorType_Name == 'AMSR_E' or SensorType_Name == 'SMOS' or SensorType_Name == 'ASCAT' or SensorType_Name == 'ASAR' or SensorType_Name == 'PALSAR':
                    Soil_Layer_Index_DA = 1
                elif SensorType_Name == 'Terra' or SensorType_Name == 'Aqua':
                    Soil_Layer_Index_DA = 2
                elif SensorType_Name == "InSitu":
                    Soil_Layer_Index_DA = 1
                elif SensorType_Name == "COSMOS":
                    Soil_Layer_Index_DA = 0
            except:
                SensorType_Name = "InSitu"
                Soil_Layer_Index_DA = 1
            
            Prop_Grid_Array_Sys_Sub[0, :, :] = CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :] + Additive_Noise_SM[Ens_Index, Soil_Layer_Index_DA]
                    
            if SensorType_Name == "COSMOS":
                CLM_Soil_Moisture_Ratio_MultiScale[:, :, :] = 1.0
                CLM_Soil_Moisture_Ensemble_Mat_Sub[:,:,:] = CLM_Soil_Moisture[:, :, :]
                
                Air_Pressure = CLM_Air_Pressure
                Observation_Frequency = 1 # Hours
                #Sensor_Layer, Sensor_Depth, Sensor_Soil_Moisture, N_uncorr = COSMOS(N0, nlyr, Def_Print, Observation_Frequency, DEM_Data, CLM_Soil_Moisture, CLM_Soil_Layer_Thickness, Air_Pressure, Bulk_Density_Top_Region, Bulk_Density_Sub_Region, Row_Numbers, Col_Numbers, Soil_Layer_Num)
                Mask = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.bool)
                Mask[::] = True
                numexpr_a = CLM_Soil_Moisture[0,:,:]
                numexpr_b = NAvalue
                numexpr_c = numpy.where(numexpr_a == numexpr_b)
                Mask[numexpr_c] = False
                Mask_Col = Mask.flatten()
                Mask_Index = numpy.where(Mask_Col == True)
                Mask_Size = numpy.size(Mask_Index)
                
                
                CLM_Soil_Layer_Thickness_Cumsum = numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness_Cumsum'][:,:,:]) * 100
                
                soil_moisture = numpy.zeros((Mask_Size,Soil_Layer_Num-5),dtype=numpy.float32)
                layerz = numpy.zeros((Mask_Size,Soil_Layer_Num-5),dtype=numpy.float32)
                
                for Soil_Layer_Index in range(Soil_Layer_Num-5):
                    soil_moisture[:,Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,:,:].flatten()[Mask_Index]
                    layerz[:,Soil_Layer_Index] = CLM_Soil_Layer_Thickness_Cumsum[Soil_Layer_Index,:,:].flatten()[Mask_Index]
                
                bd = Bulk_Density_Top_Region.flatten()[Mask_Index]
                
                n = numpy.zeros(Mask_Size,dtype=numpy.float32)
                lattwat = numpy.zeros(Mask_Size,dtype=numpy.float32)
                
                n[:] = N0
                lattwat[:] = 0.03
                nthreads = omp_get_num_procs_ParFor
                
                Neutron_COSMOS, Sensor_Soil_Moisture_COSMOS, Sensor_Depth_COSMOS = COSMIC_Py.COSMIC_Py(COSMIC,n,nlyr,soil_moisture,layerz,bd,lattwat,nthreads)
                
                Sensor_Soil_Moisture = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Sensor_Soil_Moisture_Col = Sensor_Soil_Moisture.flatten()
                Sensor_Soil_Moisture_Col[Mask_Index] = Sensor_Soil_Moisture_COSMOS
                Sensor_Soil_Moisture = numpy.reshape(Sensor_Soil_Moisture_Col,(Row_Numbers, Col_Numbers))
                
                Prop_Grid_Array_Sys_Sub[0, :, :] = Sensor_Soil_Moisture + Additive_Noise_SM[Ens_Index, Soil_Layer_Index_DA]
                    #Prop_Grid_Array_Sys_Sub[0, :, :] = CLM_Soil_Moisture[1, :, :]
                
                if Plot_Analysis >= 1:
                    import matplotlib
                    # Force matplotlib to not use any Xwindows backend.
                    matplotlib.use('Agg')
                    import matplotlib.pyplot as plt
                    import matplotlib.cm as cm
                    import matplotlib.colors as colors
    
                    Mask_Valid = numpy.zeros((Row_Numbers,Col_Numbers),dtype=numpy.bool)
                    Mask_Valid[::] = False
                    numexpr_a = Land_Mask_Data
                    numexpr_b = NAvalue
                    numexpr_c = numpy.where(numexpr_a == numexpr_b)
                    Mask_Valid[numexpr_c] = True
                    Data1 = numpy.ma.masked_where(Mask_Valid, CLM_Soil_Moisture[1, :, :])
                    Data2 = numpy.ma.masked_where(Mask_Valid, Sensor_Soil_Moisture[:, :])
                    Data3 = numpy.ma.masked_where(Mask_Valid, Sensor_Soil_Moisture)
                    
                    Variable_Min = 0.05
                    Variable_Max = 0.55
                    ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
                    color_boun_list = []
                    color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
                    for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                        color_bound[0] += color_bound[2]
                        color_boun_list.append(color_bound[0])
                    
                    fig1 = plt.figure(figsize=(15, 10), dpi=80)
                    ax = fig1.add_subplot(1, 3, 1)
                    im1 = ax.imshow(Data1, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                    plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                    ax.set_title('CLM_Soil_Moisture')
                    ax = fig1.add_subplot(1, 3, 2)
                    im1 = ax.imshow(Data2, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                    plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                    ax.set_title('CLM_Soil_Moisture_Smoothed')
                    ax = fig1.add_subplot(1, 3, 3)
                    im1 = ax.imshow(Data3, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                    plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                    ax.set_title('Sensor_Soil_Moisture')
                    plt.grid(True)
                    plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/COSMOS_SM_Operator_XF.png")
                    plt.close('all')
                
                del Mask, numexpr_a, numexpr_b, numexpr_c, Mask_Col, Mask_Index, CLM_Soil_Layer_Thickness_Cumsum, soil_moisture, layerz, n, lattwat
                del Neutron_COSMOS, Sensor_Soil_Moisture_COSMOS, Sensor_Depth_COSMOS, Sensor_Soil_Moisture
            
            else:
                Prop_Grid_Array_Sys_Sub[0, :, :] = CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :] + Additive_Noise_SM[Ens_Index, Soil_Layer_Index_DA]
        
        if Def_Print >= 2:
            print "Check the Ourliers"
        NA_Flag = False
        numexpr_a = Prop_Grid_Array_Sys_Sub[0, :, :]
        numexpr_b = CLM_NA
        numexpr_c = numpy.where(numexpr_a == numexpr_b)
        Prop_Grid_Array_Sys_Sub[0, :, :][numexpr_c] = 0.2
        Prop_Grid_Array_Sys_Sub[0, :, :][numpy.isnan(Prop_Grid_Array_Sys_Sub[0, :, :])] = 0.2
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,Prop_Grid_Array_Sys_Sub[0, :, :],NA_Flag, 'Soil_Moisture',  Variable_Assimilation_Flag, Variable_List,
                    numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][Soil_Layer_Index_DA,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][Soil_Layer_Index_DA,:,:]), 
                    numpy.asarray(numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][Soil_Layer_Index_DA,:,:])), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][Soil_Layer_Index_DA,:,:]), NAvalue)
        
        CLM_LAI = numpy.flipud(CLM_History_File.variables['TLAI'][0, ::])   # Units: m2/m2
        
        
        if Def_Print >= 2:
            print "Check the Ourliers"
        NA_Flag = False
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_LAI,NA_Flag, "LAI",  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        Prop_Grid_Array_Sys_Sub[12, :, :] = CLM_LAI[:, :]
                
    if Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")] == 1 or Feedback_Assim == 1:
        Analysis_Variable_Name[1] = "Soil Temperature"
        if Def_Print:
            print Analysis_Variable_Name[1]
        #------------------------------------------- Read the Soil Temperature Data
        CLM_Ground_Temperature = numpy.flipud(CLM_History_File.variables['TG'][0, :, :])  # Ground Temperature Units: K
        CLM_Ground_Temperature[numpy.isnan(CLM_Ground_Temperature)] = 283.15
        
        CLM_Vegetation_Temperature[numpy.isnan(CLM_Vegetation_Temperature)] = 283.15
        #print numpy.max(CLM_Vegetation_Temperature),numpy.min(CLM_Vegetation_Temperature)
        CLM_Urban_Ground_Temperature = numpy.flipud(CLM_History_File.variables['TG_U'][0, :, :])  # Urban Ground Temperature Units: K
        CLM_Rural_Ground_Temperature = numpy.flipud(CLM_History_File.variables['TG_R'][0, :, :])  # Rural Ground Temperature Units: K
        
        #CLM_Snow_Temperature =  # T_SOISNO
        for Soil_Layer_Index in range(Soil_Layer_Num):
            if Def_Print >= 2:
                print "Check the Ourliers"
            NA_Flag = False
            CLM_Soil_Temperature[Soil_Layer_Index, :, :] = numpy.flipud(CLM_History_File.variables['TSOI'][0, Soil_Layer_Index, :, :])
            
            CLM_Soil_Temperature_Ratio[Soil_Layer_Index,:,:] = CLM_Soil_Temperature[Soil_Layer_Index, :, :] / CLM_Soil_Temperature[0, :, :]
            
            
            numexpr_a = CLM_Soil_Temperature[Soil_Layer_Index,:,:]
            numexpr_b = CLM_NA
            numexpr_c = numpy.where(numexpr_a == numexpr_b)
            CLM_Soil_Temperature[Soil_Layer_Index,:,:][numexpr_c] = 283.15
            CLM_Soil_Temperature[Soil_Layer_Index,:,:][numpy.isnan(CLM_Soil_Temperature[Soil_Layer_Index,:,:])] = 283.15
            Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Soil_Temperature[Soil_Layer_Index, :, :],NA_Flag,'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                    numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
            
            CLM_Soil_Temperature_Ensemble_Mat_Sub[Soil_Layer_Index,:,:] = CLM_Soil_Temperature[Soil_Layer_Index, :, :]
            
            if Soil_Layer_Index <= 9:
                CLM_Lake_Temperature[Soil_Layer_Index, :, :] = numpy.flipud(CLM_History_File.variables['TLAKE'][0, Soil_Layer_Index, :, :])  # Lake Ground Temperature Units: K
                
                # Assign the Lake temperature to lake grid cells
                Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Lake_Temperature[Soil_Layer_Index,::],NA_Flag, 'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                                numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
                CLM_Soil_Temperature[Soil_Layer_Index, :, :][PCT_LAKE_Index] = CLM_Lake_Temperature[Soil_Layer_Index, :, :][PCT_LAKE_Index]
                CLM_Soil_Temperature_Ensemble_Mat_Sub[Soil_Layer_Index, :, :][PCT_LAKE_Index] = CLM_Soil_Temperature[Soil_Layer_Index, :, :][PCT_LAKE_Index]
            
        #print CLM_Ground_Temperature
        if Def_Print >= 2:
            print "Check the Ourliers"
        NA_Flag = False
        #print CLM_Ground_Temperature
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Ground_Temperature,NA_Flag, 'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Vegetation_Temperature,NA_Flag, "Vegetation_Temperature",  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Urban_Ground_Temperature,NA_Flag, 'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_Rural_Ground_Temperature,NA_Flag, 'Soil_Temperature',  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        
        CLM_Ground_Temperature[PCT_LAKE_Index] = CLM_Lake_Temperature[0,::][PCT_LAKE_Index]
        CLM_Vegetation_Temperature[PCT_LAKE_Index] = CLM_Lake_Temperature[0,::][PCT_LAKE_Index]
        
        Prop_Grid_Array_Sys_Sub[1, :, :] = CLM_Ground_Temperature[::]
        
        #Prop_Grid_Array_Sys_Sub[1, :, :] = numpy.power(Fractional_Vegetation_Cover * numpy.power(CLM_Vegetation_Temperature[:, :], 4) + (1.0 - Fractional_Vegetation_Cover) * numpy.power(CLM_Ground_Temperature[:, :], 4), 0.25)
        
        if Def_Print >= 2:
            print numpy.max(CLM_Vegetation_Temperature),numpy.min(CLM_Vegetation_Temperature)
            print numpy.max(CLM_Ground_Temperature),numpy.min(CLM_Ground_Temperature)
            print numpy.max(Prop_Grid_Array_Sys_Sub[1, :, :]),numpy.min(Prop_Grid_Array_Sys_Sub[1, :, :])
            print numpy.where(numpy.isnan(Prop_Grid_Array_Sys_Sub[1, :, :]))
        
        #------------------------------------------- Read the LAI Data
        # Record the PFT Based LAI
        CLM_LAI = numpy.flipud(CLM_History_File.variables['TLAI'][0, ::])   # Units: m2/m2
        
        if Def_Print >= 2:
            print "Check the Ourliers"
        NA_Flag = False
        Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,CLM_LAI,NA_Flag, "LAI",  Variable_Assimilation_Flag, Variable_List,
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][0,:,:]), 
                        numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][0,:,:]), numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][0,:,:]), NAvalue)
        Prop_Grid_Array_Sys_Sub[Variable_List.index('LAI'), :, :] = CLM_LAI[:, :]
        
    
    CLM_History_File.close()
    
    NC_FileName_Out_CLM_History_Ens = DAS_Output_Path+"Analysis/"+Region_Name+"/CLM_History_Ens_"+str(Ens_Index+1)+".nc"
    if os.path.exists(NC_FileName_Out_CLM_History_Ens):
        os.remove(NC_FileName_Out_CLM_History_Ens)
        
    if Def_Print:
        print 'Write NetCDF File:',NC_FileName_Out_CLM_History_Ens

    NC_File_Out_CLM_History_Ens = netCDF4.Dataset(NC_FileName_Out_CLM_History_Ens, 'w', diskless=True, persist=True, format='NETCDF4')
    # Dim the dimensions of NetCDF
    NC_File_Out_CLM_History_Ens.createDimension('lon', Col_Numbers)
    NC_File_Out_CLM_History_Ens.createDimension('lat', Row_Numbers)
    NC_File_Out_CLM_History_Ens.createDimension('Soil_Layer_Num', Soil_Layer_Num)
    NC_File_Out_CLM_History_Ens.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
    NC_File_Out_CLM_History_Ens.createDimension('Dim_CLM_State', Dim_CLM_State)
    NC_File_Out_CLM_History_Ens.createVariable('Prop_Grid_Array_Sys','f4',('Dim_CLM_State','lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Soil_Moisture_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Soil_Temperature_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Vegetation_Temperature_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Ground_Temperature_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Snow_Depth_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Snow_Water_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_INT_SNOW_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_FH2OSFC_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    
    NC_File_Out_CLM_History_Ens.createVariable('CLM_2m_Air_Temperature_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Air_Pressure_Ensemble_Mat','f4',('lat','lon',),zlib=True)
    
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Soil_Temperature_Ratio_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
    
    NC_File_Out_CLM_History_Ens.createVariable('CLM_Soil_Moisture_Ratio_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True)
    
    NC_File_Out_CLM_History_Ens.variables['Prop_Grid_Array_Sys'][:,:,:] = Prop_Grid_Array_Sys_Sub
    NC_File_Out_CLM_History_Ens.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:] = CLM_Soil_Moisture_Ensemble_Mat_Sub
    NC_File_Out_CLM_History_Ens.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:] = CLM_Soil_Temperature_Ensemble_Mat_Sub
    NC_File_Out_CLM_History_Ens.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:] = CLM_Vegetation_Temperature
    NC_File_Out_CLM_History_Ens.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:] = CLM_Ground_Temperature
    NC_File_Out_CLM_History_Ens.variables['CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat'][:,:] = CLM_Onset_Freezing_Degree_Days_Counter
    NC_File_Out_CLM_History_Ens.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:] = CLM_Snow_Depth
    NC_File_Out_CLM_History_Ens.variables['CLM_Snow_Water_Ensemble_Mat'][:,:] = CLM_Snow_Water
    NC_File_Out_CLM_History_Ens.variables['CLM_INT_SNOW_Ensemble_Mat'][:,:] = CLM_INT_SNOW
    NC_File_Out_CLM_History_Ens.variables['CLM_FH2OSFC_Ensemble_Mat'][:,:] = CLM_FH2OSFC
    
    NC_File_Out_CLM_History_Ens.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:] = CLM_2m_Air_Temperature
    NC_File_Out_CLM_History_Ens.variables['CLM_Air_Pressure_Ensemble_Mat'][:,:] = CLM_Air_Pressure
    
            
    NC_File_Out_CLM_History_Ens.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:] = CLM_Soil_Temperature_Ratio
    NC_File_Out_CLM_History_Ens.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:] = CLM_Soil_Moisture_Ratio_MultiScale
    NC_File_Out_CLM_History_Ens.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:] = CLM_Soil_Moisture_Ratio    
        
    NC_File_Out_CLM_History_Ens.sync()
    NC_File_Out_CLM_History_Ens.close()
        
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Constant.close()
    
    del Parameter_Soil_Space_Ensemble,Sand_Top_Region,Sand_Sub_Region,Clay_Top_Region,Clay_Sub_Region
    del Prop_Grid_Array_Sys_Sub, CLM_Soil_Moisture_Ensemble_Mat_Sub, CLM_Soil_Temperature_Ensemble_Mat_Sub, CLM_ROOTFR_Ensemble_Mat_Sub, CLM_ROOTR_Ensemble_Mat_Sub
    del CLM_Soil_Moisture, CLM_Soil_Temperature, CLM_ROOTFR, CLM_ROOTR, CLM_Ground_Temperature, CLM_Urban_Ground_Temperature, CLM_Rural_Ground_Temperature, CLM_Lake_Temperature
    del CLM_2m_Air_Temperature, CLM_Air_Pressure, CLM_Shortwave_Radiation, CLM_Longwave_Radiation, CLM_Specific_Humdity, CLM_Wind_Speed, CLM_Z0MV, CLM_DISPLA
    del CLM_LAI, CLM_Vegetation_Temperature, CLM_Albedo, CLM_Snow_Depth, CLM_Snow_Cover_Fraction,CLM_Onset_Freezing_Degree_Days_Counter
    del CLM_INT_SNOW, CLM_FH2OSFC, CLM_Snow_Water, CLM_QIRRIG, CLM_Soil_Moisture_Ratio, CLM_Soil_Temperature_Ratio, CLM_Soil_Moisture_Ratio_MultiScale
    del CLM_HTOP, CLM_HBOT, CLM_Water_Storage, CLM_Water_Table_Depth, CLM_Irrigation_Rate, CLM_Vcmax_SCOPE, CLM_m_SCOPE, CLM_SLA_SCOPE
    del CLM_Latent_Heat_Daily, CLM_Latent_Heat, CLM_Sensible_Heat, Soil_Layer_Thickness_Ratio_Moisture, Soil_Layer_Thickness_Ratio_Moisture_Sum
    del cols1d_ixy,cols1d_jxy,cols1d_ityplun,pfts1d_ixy,pfts1d_jxy,pfts1d_itypveg,pfts1d_ci,pfts1d_ityplun
    del PCT_Veg,STD_ELEV,DEM_Data,Land_Mask_Data,Bulk_Density_Top_Region,Bulk_Density_Sub_Region, PCT_LAKE, PCT_LAKE_Index,PFT_Dominant_Index,Bare_Grid_Index
    
    del COSMIC_Py, COSMIC 
    
    numexpr_a = []
    numexpr_b = []
    numexpr_c = []

    gc.collect()
    del gc.garbage[:]
    
    
    return Soil_Layer_Index_DA, Analysis_Variable_Name, Constant_File_Name


def Write_Initial_File(Ens_Index, Model_Driver, Def_PP, DasPy_Path, Run_Dir_Array, Soil_Layer_Num, ParFlow_Layer_Num, numrad, Row_Numbers, Col_Numbers, finidat_name, SensorVariable_Sub, Variable_ID_Sub, CLM_NA, Feedback_Assim, Stop_Month, Stop_Hour, UTC_Zone, \
                       pft_len,maxpft, Def_Region, DAS_Data_Path, Region_Name, Crop_Sum, SensorType_Sub, Row_Numbers_String, Col_Numbers_String, \
                       Snow_Layer_Num, column_len, Mask_Index_Sub, Def_Print, Def_ParFor, dtime, irrig_nsteps_per_day, PFT_Num, NAvalue, Density_of_liquid_water, Freezing_temperature_of_fresh_water, Density_of_ice, \
                       DAS_Depends_Path, omp_get_num_procs_ParFor, Soil_Layer_Index_DA,  Variable_Assimilation_Flag, Variable_List, 
                       NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy, \
                       NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy, fpftcon_name, finidat_name_string):
        
#    numpy.savez("Write_Initial_File.npz",Ens_Index, Soil_Layer_Num, numrad, Row_Numbers, Col_Numbers, finidat_name_Array, SensorVariable_Sub, Variable_ID_Sub, CLM_NA, Teta_Saturated, \
#                       Model_Driver, pft_len,pfts1d_ixy,pfts1d_jxy,pfts1d_ci,pfts1d_itypveg,maxpft, Def_Region, DAS_Data_Path, Region_Name, PCT_PFT, Land_Mask_Data, SensorType_Sub, Row_Numbers_String, Col_Numbers_String, \
#                       Snow_Layer_Num, column_len, cols1d_jxy, cols1d_ixy, cols1d_ityplun, Mask_Index_Sub, Def_Print, Def_ParFor, dtime, irrig_nsteps_per_day, PFT_Num, NAvalue, Density_of_liquid_water, Freezing_temperature_of_fresh_water, Density_of_ice, \
#                       DAS_Depends_Path, omp_get_num_procs_ParFor, Soil_Layer_Index_DA, Analysis_Grid_Array, Variable_List, CLM_Soil_Moisture_Ensemble_Mat, CLM_Soil_Temperature_Ensemble_Mat, CLM_Soil_Ice_Ensemble_Mat, CLM_2m_Air_Temperature_Ensemble_Mat, \
#                       Soil_Moisture_DA_Flag, Surface_Temperature_DA_Flag, Vegetation_Temperature_DA_Flag, Albedo_DA_Flag, Snow_Depth_DA_Flag, Snow_Cover_Fraction_DA_Flag, Canopy_Water_DA_Flag,\
#                       Snow_Water_Equivalent_DA_Flag, Crop_Planting_DA_Flag, Crop_Harvest_DA_Flag, LAI_DA_Flag, Latent_Heat_DA_Flag, Latent_Heat_Daily_DA_Flag, Sensible_Heat_DA_Flag, Water_Storage_DA_Flag, Water_Table_Depth_DA_Flag, Emissivity_DA_Flag, Irrigation_Rate_DA_Flag, Irrigation_Scheduling_Flag, CLM_Soil_Layer_Thickness)
    
    ParFor = imp.load_source("ParFor",DasPy_Path+"ParFor.py")
    
    CLM_Initial_File = netCDF4.Dataset(finidat_name_string, 'r')
    cols1d_ixy = CLM_Initial_File.variables['cols1d_ixy'][:]
    cols1d_jxy = CLM_Initial_File.variables['cols1d_jxy'][:]
    cols1d_ityplun = CLM_Initial_File.variables['cols1d_ityplun'][:]
    
    #numpy.savetxt('cols1d_ixy',cols1d_ixy)
    #numpy.savetxt('cols1d_jxy',cols1d_jxy)
    pfts1d_ixy = CLM_Initial_File.variables['pfts1d_ixy'][:]
    pfts1d_jxy = CLM_Initial_File.variables['pfts1d_jxy'][:]
    pfts1d_itypveg = CLM_Initial_File.variables['pfts1d_itypveg'][:]
    pfts1d_ci = CLM_Initial_File.variables['pfts1d_ci'][:]
    pfts1d_ityplun = CLM_Initial_File.variables['pfts1d_ityplun'][:]
    
    CLM_Initial_File.close()
    
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][Ens_Index,:,:,:]
    CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Snow_Depth_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Snow_Water_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,Ens_Index]
    CLM_INT_SNOW_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_INT_SNOW_Ensemble_Mat'][:,:,Ens_Index]
    CLM_FH2OSFC_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_FH2OSFC_Ensemble_Mat'][:,:,Ens_Index]
    #CLM_ROOTFR_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Soil_Moisture_Ratio_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:]
    
    CLM_2m_Air_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    CLM_Soil_Layer_Thickness = numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][:,:,:])
    Teta_Residual = NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][:,:,:]    
    Teta_Saturated = NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][:,:,:]
    Teta_Field_Capacity = NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][:,:,:]    
    Teta_Wilting_Point = NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][:,:,:]
    STD_ELEV = NC_File_Out_Assimilation_2_Constant.variables['STD_ELEV'][:,:]
            
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Constant.close()
    NC_File_Out_Assimilation_2_Diagnostic.close()
    
    NC_File_Out_Assimilation_2_Initial_Copy = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial_Copy, 'r')
    CLM_Soil_Moisture_Ensemble_Mat_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,Ens_Index]
    NC_File_Out_Assimilation_2_Initial_Copy.close()
        
    CLM_Soil_Moisture = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32) # Units: m3/m3
    CLM_Soil_Temperature = numpy.zeros((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)  # soil-snow temperature Units: K
    CLM_Vegetation_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_Ground_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_Albedo = numpy.zeros((numrad*2, Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units None
    CLM_Snow_Depth = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units Meters
    CLM_Snow_Cover_Fraction = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units Meters
    CLM_Snow_Water = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_INT_SNOW = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_FH2OSFC = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Canopy_Water = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Emissivity = numpy.zeros((numrad, Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units None
    CLM_2m_Air_Temperature = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)   # Units: K
    CLM_LAI = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_HTOP = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_HBOT = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Water_Storage = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Water_Table_Depth = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Irrigation_Rate = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Planting_Date = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    CLM_Harvest_Date = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    CLM_Snow_Density = numpy.zeros((Snow_Layer_Num, column_len),dtype=numpy.float32)   # Snow Water: Units: (kg/m^3)
    CLM_Snow_Density_Whole = numpy.zeros(column_len,dtype=numpy.float32)   # Snow Water: Units: (kg/m^3)
    
    CLM_Soil_Temperature = CLM_Soil_Temperature_Ensemble_Mat[:,:,:]
    CLM_2m_Air_Temperature = CLM_2m_Air_Temperature_Ensemble_Mat[:,:]
    CLM_Soil_Moisture_Ratio = CLM_Soil_Moisture_Ratio_Ensemble_Mat[:,:,:]
    
    CLM_FH2OSFC = CLM_FH2OSFC_Ensemble_Mat[::]
    CLM_INT_SNOW = CLM_INT_SNOW_Ensemble_Mat[::]
    
    print "============================= Write the Model Initial Data ==========================================="
    finidat_name_string = Run_Dir_Array[Ens_Index]+finidat_name
    #------------------------------------------- Read the CLM Initial File
    try:
        CLM_Initial_File = netCDF4.Dataset(finidat_name_string, "r+")
    except:
        print "Open Initial File:", finidat_name_string, "Failed!!"
        os.abort()
        
    if Def_Print:
        print "SensorVariable_Sub",SensorVariable_Sub
        print "Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)",Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)]
#                            print "============================= Update the Deep Layers (Need to be studied in the future) ================================================="
    if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Soil_Moisture") or \
    (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature") or \
        (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Latent_Heat") or \
        (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Latent_Heat_Daily") or  \
        (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Sensible_Heat"):
        #Prop_Grid_Array_Sys[0, :, :] = CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :]
        #CLM_Soil_Moisture[:, :, :] = 0.2
        #CLM_Soil_Moisture[:, :, :] = CLM_Soil_Moisture[:, :, :] * Analysis_Grid[Prop_Grid_Array_Sys_Index,::] / CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :]
        #print CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :],Analysis_Grid_Array[Ens_Index,Variable_List.index("Soil_Moisture"),::]
        #CLM_Soil_Moisture[Soil_Layer_Index_DA, :, :] = Analysis_Grid_Array[Variable_List.index("Soil_Moisture"),::]
#                                for Soil_Layer_Index in range(1,Soil_Layer_Num):
#                                    CLM_Soil_Moisture[Soil_Layer_Index,:,:] = CLM_Soil_Moisture[Soil_Layer_Index,:,:] * (CLM_Soil_Moisture[0,:,:] / Prop_Grid_Array_Sys[:,:,0])
        CLM_Soil_Moisture = CLM_Soil_Moisture_Ensemble_Mat[:,:,:]
        #CLM_Soil_Moisture_Copy = CLM_Soil_Moisture_Ensemble_Mat_Copy[:,:,:]
        #CLM_Ground_Temperature = Analysis_Grid_Array[Variable_List.index("Surface_Temperature"),::]
        CLM_Vegetation_Temperature = CLM_Vegetation_Temperature_Ensemble_Mat[:,:]
        CLM_Ground_Temperature = CLM_Ground_Temperature_Ensemble_Mat[:,:]
        #print "Minimum Analysis Value is:", CLM_Ground_Temperature.min(), "Maximum Analysis Value is:", CLM_Ground_Temperature.max()
#                                for Soil_Layer_Index in range(Soil_Layer_Num):
#                                    CLM_Soil_Temperature[Soil_Layer_Index,:,:] = CLM_Soil_Temperature[Soil_Layer_Index,:,:] * (CLM_Ground_Temperature[:,:] / Prop_Grid_Array_Sys[:,:,0])
            
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Vegetation_Temperature":
        CLM_Vegetation_Temperature[:, :] = Analysis_Grid_Array[Variable_List.index("Vegetation_Temperature"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Albedo":
        if Variable_ID_Sub == "Albedo_BSA_Band_vis":
            CLM_Albedo[0, :, :] = Analysis_Grid_Array[Variable_List.index("Albedo_BSA_Band_vis"),::]
        if Variable_ID_Sub == "Albedo_BSA_Band_nir":
            CLM_Albedo[1, :, :] = Analysis_Grid_Array[Variable_List.index("Albedo_BSA_Band_nir"),::]
        if Variable_ID_Sub == "Albedo_WSA_Band_vis":
            CLM_Albedo[2, :, :] = Analysis_Grid_Array[Variable_List.index("Albedo_WSA_Band_vis"),::]
        if Variable_ID_Sub == "Albedo_WSA_Band_nir":
            CLM_Albedo[3, :, :] = Analysis_Grid_Array[Variable_List.index("Albedo_WSA_Band_nir"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Snow_Cover_Fraction":
        CLM_Snow_Cover_Fraction[:, :] = Analysis_Grid_Array[Variable_List.index("Snow_Cover_Fraction"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Snow_Depth":
        CLM_Snow_Depth[:, :] = Analysis_Grid_Array[Variable_List.index("Snow_Depth"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Snow_Water_Equivalent":
        CLM_Snow_Water[:, :] = Analysis_Grid_Array[Variable_List.index("Snow_Water_Equivalent"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Crop_Planting_Date":
        CLM_Planting_Date[:, :] = Analysis_Grid_Array[Variable_List.index("Crop_Planting_Date"),::]
        CLM_Planting_Date[numpy.where(CLM_Planting_Date < 0.0)] = 0.0
        CLM_Planting_Date[numpy.where(CLM_Planting_Date > 366)] = 0.0
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Crop_Harvest_Date":
        CLM_Harvest_Date[:, :] = Analysis_Grid_Array[Variable_List.index("Crop_Harvest_Date"),::]
        CLM_Harvest_Date[numpy.where(CLM_Harvest_Date < 0.0)] = 0.0
        CLM_Harvest_Date[numpy.where(CLM_Harvest_Date > 366)] = 0.0
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "LAI":
        CLM_LAI[:, :] = Analysis_Grid_Array[Variable_List.index("LAI"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Water_Storage":
        CLM_Water_Storage[:, :] = Analysis_Grid_Array[Variable_List.index("Water_Storage"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Water_Table":
        CLM_Water_Table_Depth[:, :] = Analysis_Grid_Array[Variable_List.index("Water_Table"),::]
    elif Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Irrigation_Rate":
        CLM_Irrigation_Rate[:, :] = Analysis_Grid_Array[Variable_List.index("Irrigation_Rate"),::]
    else:
        print "Wrong Assimilation Variable in Write Initial\n"
#        print numpy.ma.count(Observation_Matrix[Observation_Matrix_Index]),numpy.shape(Mask_Index_Sub)
#        print numpy.ma.count(Observation_Matrix[Observation_Matrix_Index][numpy.where(Mask_Index_Sub==False)])
#        print numpy.size(numpy.where(Mask_Index_Sub==False)[0]) * 0.01
    if Def_Print:
        print "Write Initial File:", finidat_name_string
    
    #------------------------------------------- Update the Soil Moisture Data
    if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Soil_Moisture"):
    
        print "------------------------------------------- Update the Soil Moisture Data"
        if Def_Print >= 2:
            print "Soil Liquid Water: ", numpy.min(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA]), numpy.max(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA])
        
        
        if Def_ParFor:
            H2OSOI_LIQ_Temp = ParFor.ParFor_Soil_Moisture(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,Soil_Layer_Index_DA,\
                                    Snow_Layer_Num,Soil_Layer_Num-Snow_Layer_Num,CLM_Soil_Moisture,CLM_Soil_Layer_Thickness,Density_of_liquid_water,Density_of_ice,\
                                    Teta_Saturated,cols1d_ityplun,CLM_Soil_Temperature,Freezing_temperature_of_fresh_water,\
                                    numpy.asarray(CLM_Initial_File.variables["H2OSOI_LIQ"][::]),Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
             
            # Only update the first 7 layers to avoid unreasonable subsurface flow
            CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num - 5)] = H2OSOI_LIQ_Temp[:,0:(Soil_Layer_Num - 5)]
            
            del H2OSOI_LIQ_Temp
            if Feedback_Assim:# and (string.atoi(Stop_Month) >= 5) and (string.atoi(Stop_Month) <= 9):
                CLM_Initial_File.variables["T_GRND"][:] = ParFor.ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_GRND"][:]),CLM_Ground_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                CLM_Initial_File.variables["T_GRND_U"][:] = CLM_Initial_File.variables["T_GRND"][:]
                CLM_Initial_File.variables["T_GRND_R"][:] = CLM_Initial_File.variables["T_GRND"][:]
                CLM_Initial_File.variables["TH2OSFC"][:] = CLM_Initial_File.variables["T_GRND"][:]
                CLM_Initial_File.variables["T_VEG"][:] = ParFor.ParFor_Common(CLM_NA,pft_len,Row_Numbers,Col_Numbers,pfts1d_ixy,pfts1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_VEG"][:]),CLM_Vegetation_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                #T_GRND_Post = CLM_Initial_File.variables["T_GRND"][:]
                T_SOISNO_Temp = numpy.asarray(CLM_Initial_File.variables["T_SOISNO"][:,:])
                for Soil_Layer_Index in range(Soil_Layer_Num):
                    #print numpy.shape(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Index)])
                    #print numpy.shape(ParFor_Ratio(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Index)],T_GRND_Post,T_GRND_Original,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue))
                    #CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index] = ParFor_Ratio(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index],T_GRND_Post,T_GRND_Original,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                    T_SOISNO_Temp[:,Snow_Layer_Num+Soil_Layer_Index] = ParFor.ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index][:],CLM_Soil_Temperature[Soil_Layer_Index,:,:],Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                
                CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Num] = T_SOISNO_Temp[:,Snow_Layer_Num+Soil_Layer_Num]
                del T_SOISNO_Temp
                
                Lake_Layer_Num = 10
                T_LAKE_Temp = numpy.asarray(CLM_Initial_File.variables["T_LAKE"][:,:])
                for Lake_Layer_Index in range(Lake_Layer_Num):
                    T_LAKE_Temp[:,Lake_Layer_Index] = ParFor.ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_LAKE"][:,Lake_Layer_Index]),CLM_Soil_Temperature[Lake_Layer_Index,:,:],Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                
                CLM_Initial_File.variables["T_LAKE"][:,:] = T_LAKE_Temp
                del T_LAKE_Temp
                    
        else:
            H2OSOI_LIQ_Temp = numpy.asarray(CLM_Initial_File.variables["H2OSOI_LIQ"][:,:])
            for Column_Index in range(column_len):
                Row_Index = Row_Numbers - cols1d_jxy[Column_Index]
                Col_Index = cols1d_ixy[Column_Index] - 1
                cols1d_ityplun_index = cols1d_ityplun[Column_Index]
                
                if cols1d_ityplun_index != 3:
                
                    if not Mask_Index_Sub[Row_Index,Col_Index]:
                        if cols1d_ityplun_index == 1 or cols1d_ityplun_index == 5 or cols1d_ityplun_index == 8:
                            
                            # Transform the CLM Soil Moisture Unit (m3/m3) to (kg/m2) after data assimilation
                            #for Soil_Layer_Index in range(Soil_Layer_Index_DA,Soil_Layer_Index_DA+1,1):
                            for Soil_Layer_Index in range(Soil_Layer_Num - 5):  # Only update the first 7 layers to avoid unreasonable subsurface flow
                                if CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] != NAvalue and CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] < CLM_NA:
                                    #H2OSOI_LIQ_Original = CLM_Initial_File.variables["H2OSOI_LIQ"][Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                    #H2OSOI_ICE_Original = CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                    if CLM_Soil_Temperature[Soil_Layer_Index,Row_Index, Col_Index] <= Freezing_temperature_of_fresh_water:
                                        
                                        H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water)
                                        #print H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                        # Update the Soil Ice
                                        #if CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] + CLM_Soil_Ice[Soil_Layer_Index,Row_Index,Col_Index] > Teta_Saturated[Soil_Layer_Index,Row_Index,Col_Index]:
                                        #    CLM_Soil_Ice[Snow_Layer_Num+Soil_Layer_Index,Row_Index,Col_Index] = max([Teta_Saturated[Soil_Layer_Index,Row_Index,Col_Index] - CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index],0.0])
                                        #    CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Ice[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_ice)
    
                                        #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / H2OSOI_LIQ_Temp[Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                        #CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                    else:
                                        H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water)
                                        #CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = 0.0
                                        #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / H2OSOI_LIQ_Temp[Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                        #for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                                        #    if CLM_Soil_Temperature[Soil_Layer_Index,Row_Index, Col_Index] > Freezing_temperature_of_fresh_water:
                                        #        CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = 0.0
                        else:
                            pass
    #                            if H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] > (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water) or H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] < 0.0:
    #                                print Soil_Layer_Index,H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index],CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index],CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index], Density_of_liquid_water 
                            #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = 0.2 * (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
                            #H2OSOI_LIQ_Temp[Column_Index,range(Soil_Layer_Num)] = CLM_Soil_Moisture[range(Soil_Layer_Num),Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
                            #print CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,:]
                            #print CLM_Soil_Moisture[range(Soil_Layer_Num),Row_Index,Col_Index], (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
            CLM_Initial_File.variables["H2OSOI_LIQ"][:,:] = H2OSOI_LIQ_Temp
            del H2OSOI_LIQ_Temp
            
        if Def_Print >= 2:
            print "Soil Liquid Water: ", numpy.min(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA]), numpy.max(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA])
        
             
    #------------------------------------------- Update the Soil Temperature Data
    elif (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Surface_Temperature") or (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] and SensorVariable_Sub == "Sensible_Heat"):
        print "------------------------------------------- Update the Soil Temperature Data"
        if Def_Print >= 2:
            print "T_GRND: ", numpy.min(CLM_Initial_File.variables["T_GRND"][:]), numpy.max(CLM_Initial_File.variables["T_GRND"][:])
            print "T_GRND_U: ", numpy.min(CLM_Initial_File.variables["T_GRND_U"][:]), numpy.max(CLM_Initial_File.variables["T_GRND_U"][:])
            print "T_GRND_R: ", numpy.min(CLM_Initial_File.variables["T_GRND_R"][:]), numpy.max(CLM_Initial_File.variables["T_GRND_R"][:])
            print "T_LAKE: ", numpy.min(CLM_Initial_File.variables["T_LAKE"][:,0]), numpy.max(CLM_Initial_File.variables['T_LAKE'][:,0])
            print "T_VEG: ", numpy.min(CLM_Initial_File.variables['T_VEG'][:]), numpy.max(CLM_Initial_File.variables['T_VEG'][:])
            print "TH2OSFC: ", numpy.min(CLM_Initial_File.variables['TH2OSFC'][:]), numpy.max(CLM_Initial_File.variables['TH2OSFC'][:])
            print "T_SOISNO: ", numpy.min(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+2]), numpy.max(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+2])
        #T_GRND_Original = CLM_Initial_File.variables["T_GRND"][:]
        #print CLM_Initial_File.variables["T_GRND"][1500],CLM_Initial_File.variables["T_SOISNO"][1500,5]
        if Def_ParFor:
            #print CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_GRND"][:],CLM_Ground_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue
            #print numpy.shape(CLM_Initial_File.variables["T_GRND"][:]),numpy.shape(ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_GRND"][:],CLM_Ground_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue))
            CLM_Initial_File.variables["T_GRND"][:] = ParFor.ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_GRND"][:]),CLM_Ground_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            CLM_Initial_File.variables["T_GRND_U"][:] = CLM_Initial_File.variables["T_GRND"][:]
            CLM_Initial_File.variables["T_GRND_R"][:] = CLM_Initial_File.variables["T_GRND"][:]
            CLM_Initial_File.variables["TH2OSFC"][:] = CLM_Initial_File.variables["T_GRND"][:]
            CLM_Initial_File.variables["T_VEG"][:] = ParFor.ParFor_Common(CLM_NA,pft_len,Row_Numbers,Col_Numbers,pfts1d_ixy,pfts1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_VEG"][:]),CLM_Vegetation_Temperature,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            #T_GRND_Post = CLM_Initial_File.variables["T_GRND"][:]
            T_SOISNO_Temp = numpy.asarray(CLM_Initial_File.variables["T_SOISNO"][:,:])
            T_LAKE_Temp = numpy.asarray(CLM_Initial_File.variables["T_LAKE"][:,:])
            
            for Soil_Layer_Index in range(Soil_Layer_Num):
                #print numpy.shape(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Index)])
                #print numpy.shape(ParFor_Ratio(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Index)],T_GRND_Post,T_GRND_Original,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue))
                #CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index] = ParFor_Ratio(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index],T_GRND_Post,T_GRND_Original,Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                T_SOISNO_Temp[:,Snow_Layer_Num+Soil_Layer_Index] = ParFor.ParFor_Common(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,numpy.asarray(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index][:]),CLM_Soil_Temperature[Soil_Layer_Index,:,:],Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                if Soil_Layer_Index <= 9:
                    T_LAKE_Temp[:,Soil_Layer_Index] = T_SOISNO_Temp[:,Snow_Layer_Num+Soil_Layer_Index]
            CLM_Initial_File.variables["T_SOISNO"][:,:] = T_SOISNO_Temp
            CLM_Initial_File.variables["T_LAKE"][:,:] = T_LAKE_Temp
            del T_SOISNO_Temp,T_LAKE_Temp
            
            if Feedback_Assim:# and (string.atoi(Stop_Month) >= 5) and (string.atoi(Stop_Month) <= 9):
                print "------------------------------------------- Update the Soil Moisture Data"
                if Def_Print >= 2:
                    print "Soil Liquid Water: ", numpy.min(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA]), numpy.max(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA])
                
                
                H2OSOI_LIQ_Temp = ParFor.ParFor_Soil_Moisture(CLM_NA,column_len,Row_Numbers,Col_Numbers,cols1d_ixy,cols1d_jxy,Soil_Layer_Index_DA,\
                                        Snow_Layer_Num,Soil_Layer_Num-Snow_Layer_Num,CLM_Soil_Moisture,CLM_Soil_Layer_Thickness,Density_of_liquid_water,Density_of_ice,\
                                        Teta_Saturated[Soil_Layer_Index_DA, :, :],cols1d_ityplun,CLM_Soil_Temperature,Freezing_temperature_of_fresh_water,\
                                        numpy.asarray(CLM_Initial_File.variables["H2OSOI_LIQ"][::]),Mask_Index_Sub,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
                
                # Only update the first 7 layers to avoid unreasonable subsurface flow
                CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num - 5)] = H2OSOI_LIQ_Temp[:,0:(Soil_Layer_Num - 5)]
                del H2OSOI_LIQ_Temp
                if Def_Print >= 2:
                    print "Soil Liquid Water: ", numpy.min(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA]), numpy.max(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA])
                        
        else:
            
            T_GRND_Temp = numpy.asarray(CLM_Initial_File.variables["T_GRND"][:])
            T_GRND_U_Temp = numpy.asarray(CLM_Initial_File.variables["T_GRND_U"][:])
            T_GRND_R_Temp = numpy.asarray(CLM_Initial_File.variables["T_GRND_R"][:])
            TH2OSFC_Temp = numpy.asarray(CLM_Initial_File.variables["TH2OSFC"][:])
            T_SOISNO_Temp = numpy.asarray(CLM_Initial_File.variables["T_SOISNO"][:,:])
            T_LAKE_Temp = numpy.asarray(CLM_Initial_File.variables["T_LAKE"][:,:])
            
            for Column_Index in range(column_len):
                #Grid_Index = cols1d_gi[Column_Index] - 1
                Row_Index = Row_Numbers - cols1d_jxy[Column_Index]
                Col_Index = cols1d_ixy[Column_Index] - 1
                
                if not Mask_Index_Sub[Row_Index,Col_Index]:
                    if CLM_Ground_Temperature[Row_Index, Col_Index] != CLM_NA:
                        [Column_Index] = CLM_Ground_Temperature[Row_Index, Col_Index]
                        T_GRND_Temp[Column_Index] = CLM_Ground_Temperature[Row_Index, Col_Index]
                        T_GRND_U_Temp[Column_Index] = CLM_Ground_Temperature[Row_Index, Col_Index]
                        T_GRND_R_Temp[Column_Index] = CLM_Ground_Temperature[Row_Index, Col_Index]
                        TH2OSFC_Temp[Column_Index] = CLM_Ground_Temperature[Row_Index, Col_Index]
                        #CLM_Initial_File.variables["T_SOISNO"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = CLM_Initial_File.variables["T_SOISNO"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (CLM_Initial_File.variables["T_GRND"][Column_Index] / T_GRND_Original[Column_Index])
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            T_SOISNO_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Temperature[Soil_Layer_Index, Row_Index, Col_Index]
                        #print CLM_Initial_File.variables["T_SOISNO"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)]
                        if Soil_Layer_Index <= 9:
                            T_LAKE_Temp[:,Soil_Layer_Index] = CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+Soil_Layer_Index]
            
            CLM_Initial_File.variables["T_GRND"][:] = T_GRND_Temp
            CLM_Initial_File.variables["T_GRND_U"][:] = T_GRND_U_Temp
            CLM_Initial_File.variables["T_GRND_R"][:] = T_GRND_R_Temp
            CLM_Initial_File.variables["TH2OSFC"][:] = TH2OSFC_Temp
            CLM_Initial_File.variables["T_SOISNO"][:,:] = T_SOISNO_Temp
            CLM_Initial_File.variables["T_LAKE"][:,:] = T_LAKE_Temp
            del T_GRND_Temp,T_GRND_U_Temp,T_GRND_R_Temp,TH2OSFC_Temp,T_SOISNO_Temp,T_LAKE_Temp
            
            T_Veg_Temp = numpy.asarray(CLM_Initial_File.variables["T_VEG"][:])
            for Pft_Index in range(pft_len):
                Row_Index = Row_Numbers - pfts1d_jxy[Pft_Index]
                Col_Index = pfts1d_ixy[Pft_Index] - 1
                if not Mask_Index_Sub[Row_Index,Col_Index]:
                    if CLM_Ground_Temperature[Row_Index, Col_Index] != CLM_NA:
                        T_Veg_Temp[Pft_Index] = CLM_Vegetation_Temperature[Row_Index,Col_Index]
            CLM_Initial_File.variables["T_VEG"][:] = T_Veg_Temp
            del T_Veg_Temp
            
            if Feedback_Assim:# and (string.atoi(Stop_Month) >= 5) and (string.atoi(Stop_Month) <= 9):
                print "------------------------------------------- Update the Soil Moisture Data"
                if Def_Print >= 2:
                    print "Soil Liquid Water: ", numpy.min(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA]), numpy.max(CLM_Initial_File.variables["H2OSOI_LIQ"][:,Snow_Layer_Num+Soil_Layer_Index_DA])
                
                H2OSOI_LIQ_Temp = numpy.asarray(CLM_Initial_File.variables["H2OSOI_LIQ"][:,:])
                
                for Column_Index in range(column_len):
                    Row_Index = Row_Numbers - cols1d_jxy[Column_Index]
                    Col_Index = cols1d_ixy[Column_Index] - 1
                    cols1d_ityplun_index = cols1d_ityplun[Column_Index]
                    
                    if cols1d_ityplun_index != 3:
                    
                        if not Mask_Index_Sub[Row_Index,Col_Index]:
                            if cols1d_ityplun_index == 1 or cols1d_ityplun_index == 5 or cols1d_ityplun_index == 8:
                                
                                # Transform the CLM Soil Moisture Unit (m3/m3) to (kg/m2) after data assimilation
                                #for Soil_Layer_Index in range(Soil_Layer_Index_DA,Soil_Layer_Index_DA+1,1):
                                for Soil_Layer_Index in range(Soil_Layer_Num - 5):  # Only update the first 7 layers to avoid unreasonable subsurface flow
                                    if CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] != NAvalue and CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] < CLM_NA:
                                        #H2OSOI_LIQ_Original = CLM_Initial_File.variables["H2OSOI_LIQ"][Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                        #H2OSOI_ICE_Original = CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                        if CLM_Soil_Temperature[Soil_Layer_Index,Row_Index, Col_Index] <= Freezing_temperature_of_fresh_water:
                                            
                                            H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water)
                                            #print H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index]
                                            # Update the Soil Ice
                                            #if CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] + CLM_Soil_Ice[Soil_Layer_Index,Row_Index,Col_Index] > Teta_Saturated[Soil_Layer_Index,Row_Index,Col_Index]:
                                            #    CLM_Soil_Ice[Snow_Layer_Num+Soil_Layer_Index,Row_Index,Col_Index] = max([Teta_Saturated[Soil_Layer_Index,Row_Index,Col_Index] - CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index],0.0])
                                            #    CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Ice[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_ice)
        
                                            #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / H2OSOI_LIQ_Temp[Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                            #CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                        else:
                                            H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water)
                                            #CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = 0.0
                                            #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] * (H2OSOI_LIQ_Original / H2OSOI_LIQ_Temp[Column_Index,(Snow_Layer_Num+Soil_Layer_Index_DA)])
                                            #for Soil_Layer_Index in range(Soil_Layer_Num - 5):
                                            #    if CLM_Soil_Temperature[Soil_Layer_Index,Row_Index, Col_Index] > Freezing_temperature_of_fresh_water:
                                            #        CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,Snow_Layer_Num+Soil_Layer_Index] = 0.0
                            else:
                                pass
        #                            if H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] > (CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index] * Density_of_liquid_water) or H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index] < 0.0:
        #                                print Soil_Layer_Index,H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num+Soil_Layer_Index],CLM_Soil_Moisture[Soil_Layer_Index,Row_Index,Col_Index],CLM_Soil_Layer_Thickness[Soil_Layer_Index,Row_Index,Col_Index], Density_of_liquid_water 
                                #H2OSOI_LIQ_Temp[Column_Index,Snow_Layer_Num:(Snow_Layer_Num+Soil_Layer_Num)] = 0.2 * (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
                                #H2OSOI_LIQ_Temp[Column_Index,range(Soil_Layer_Num)] = CLM_Soil_Moisture[range(Soil_Layer_Num),Row_Index,Col_Index] * (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
                                #print CLM_Initial_File.variables["H2OSOI_ICE"][Column_Index,:]
                                #print CLM_Soil_Moisture[range(Soil_Layer_Num),Row_Index,Col_Index], (CLM_Soil_Layer_Thickness[range(Soil_Layer_Num),Row_Index,Col_Index] * 1000.0) 
                CLM_Initial_File.variables["H2OSOI_LIQ"][:,:] = H2OSOI_LIQ_Temp
                del H2OSOI_LIQ_Temp

        #CLM_Initial_File.variables["T_GRND"][:] = 313.15
        #CLM_Initial_File.variables["T_VEG"][:] = 313.15
        #CLM_Initial_File.variables["T_SOISNO"][:,5] = 313.15
#        CLM_Initial_File.variables["T_SOISNO"][:,6] = 313.15
        #print CLM_Initial_File.variables["T_GRND"][1500],CLM_Initial_File.variables["T_SOISNO"][1500,5]
        if Def_Print >= 2:
            print "T_GRND: ", numpy.min(CLM_Initial_File.variables["T_GRND"][:]), numpy.max(CLM_Initial_File.variables["T_GRND"][:])
            print "T_GRND_U: ", numpy.min(CLM_Initial_File.variables["T_GRND_U"][:]), numpy.max(CLM_Initial_File.variables["T_GRND_U"][:])
            print "T_GRND_R: ", numpy.min(CLM_Initial_File.variables["T_GRND_R"][:]), numpy.max(CLM_Initial_File.variables["T_GRND_R"][:])
            print "T_LAKE: ", numpy.min(CLM_Initial_File.variables["T_LAKE"][:,0]), numpy.max(CLM_Initial_File.variables['T_LAKE'][:,0])
            print "T_VEG: ", numpy.min(CLM_Initial_File.variables['T_VEG'][:]), numpy.max(CLM_Initial_File.variables['T_VEG'][:])
            print "TH2OSFC: ", numpy.min(CLM_Initial_File.variables['TH2OSFC'][:]), numpy.max(CLM_Initial_File.variables['TH2OSFC'][:])
            print "T_SOISNO: ", numpy.min(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+2]), numpy.max(CLM_Initial_File.variables["T_SOISNO"][:,Snow_Layer_Num+2])
    
    CLM_Initial_File.sync()
    CLM_Initial_File.close()
    
    del Analysis_Grid_Array, CLM_Soil_Moisture_Ensemble_Mat, CLM_Soil_Temperature_Ensemble_Mat, CLM_Vegetation_Temperature_Ensemble_Mat
    del CLM_Ground_Temperature_Ensemble_Mat, CLM_2m_Air_Temperature_Ensemble_Mat, CLM_Snow_Depth_Ensemble_Mat, CLM_Snow_Water_Ensemble_Mat
    del CLM_INT_SNOW_Ensemble_Mat, CLM_FH2OSFC_Ensemble_Mat, CLM_Soil_Layer_Thickness, CLM_Soil_Moisture_Ensemble_Mat_Copy
    del CLM_Soil_Moisture, CLM_Soil_Temperature, CLM_Vegetation_Temperature, CLM_Ground_Temperature, CLM_Albedo,
    del CLM_Snow_Depth, CLM_Snow_Cover_Fraction, CLM_Snow_Water, CLM_INT_SNOW, CLM_FH2OSFC, CLM_Canopy_Water, CLM_Emissivity, CLM_2m_Air_Temperature, CLM_LAI
    del CLM_HTOP, CLM_HBOT, CLM_Water_Storage, CLM_Water_Table_Depth, CLM_Irrigation_Rate, CLM_Planting_Date, CLM_Harvest_Date
    del CLM_Snow_Density, CLM_Snow_Density_Whole, CLM_Soil_Moisture_Ratio_Ensemble_Mat, CLM_Soil_Moisture_Ratio
    del cols1d_ixy,cols1d_jxy,cols1d_ityplun,pfts1d_ixy,pfts1d_jxy,pfts1d_itypveg,pfts1d_ci,pfts1d_ityplun
    del Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point, STD_ELEV
        
    gc.collect()
    del gc.garbage[:]
    
    
    return 0

def Parameter_Space_Function(Model_Driver, Def_Print, Def_PP, active_nodes_server, job_server_node_array, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, Def_Debug, Def_Region, Def_First_Run, Def_Par_Optimized, \
                             Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Dim_ParFlow_Par, ParFlow_Layer_Num, Start_Month, maxpft, Soil_Texture_Layer_Opt_Num, Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict, \
                    Parameter_Optimization, Row_Numbers_String, Col_Numbers_String, Soil_Sand_Clay_Sum, Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  PFT_Dominant_Index, topo_slope,\
                    fsurdat_name, fpftcon_name, NC_FileName_Assimilation_2_Constant, DasPy_Path,
                    DAS_Data_Path, DAS_Output_Path, Region_Name, Datetime_Start, Datetime_Initial, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, DAS_Depends_Path, Def_ParFor, omp_get_num_procs_ParFor, r):
    
#    numpy.savez("Parameter_Space_Function.npz",Def_Debug, Def_Region, Def_First_Run, Def_Par_Optimized, Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Texture_Layer_Opt_Num, Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict, \
#                    Parameter_Optimization, Parameter_Soil_Space_Single, Parameter_Veg_Space_Single, Row_Numbers_String, Col_Numbers_String, Soil_Sand_Clay_Sum, Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  PFT_Par_Sens_Array, \
#                    DAS_Data_Path, Datetime_Start, Datetime_Initial, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, DAS_Depends_Path, Def_ParFor, omp_get_num_procs_ParFor)
        
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    Bulk_Density = NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Top_Region'][:,:]
    NC_File_Out_Assimilation_2_Constant.close()
    
    # 0-min 1-max
    Parameter_Range_Soil = numpy.zeros((2, Dim_Soil_Par),dtype=numpy.float32)
    Parameter_Range_Veg = numpy.zeros((2,Dim_Veg_Par),dtype=numpy.float32)
    Parameter_Range_PFT = numpy.zeros((2,Dim_PFT_Par),dtype=numpy.float32)
    Parameter_Range_Hard = numpy.zeros((2, Dim_Hard_Par),dtype=numpy.float32)
    
    #Parameter_Range_Hard[:,0] = numpy.asarray([1.0e-6,1e-2])
    #Parameter_Range_Hard[:,0] = numpy.asarray([-6.0,-2.0])  # log10 transform
    Parameter_Range_Hard[:,0] = numpy.asarray([-2.25,2.25])   # log10 factor - qdrai
    Parameter_Range_Hard[:,1] = numpy.asarray([-1.25,1.25])    # log2 fdrai factor
    Parameter_Range_Hard[:,2] = numpy.asarray([0.75,2.25])    # vcmax factor (CLM overestimate LE in the tropical region, underestimate LE in the mid-latitude and high-latitude)
    
    Parameter_Range_PFT[:,0] = numpy.asarray([0.75,2.25])
    Parameter_Range_PFT[:,1] = numpy.asarray([0.75,2.25])
    Parameter_Range_PFT[:,2] = numpy.asarray([0.75,2.25])
    
    if Def_Region == -2:
        Parameter_Range_PFT[:,0] = numpy.asarray([0.75,2.75])
        Parameter_Range_PFT[:,1] = numpy.asarray([0.75,2.75])
        Parameter_Range_PFT[:,2] = numpy.asarray([0.75,2.75])
    
    Parameter_Range_Soil[:,0] = numpy.asarray([6.0,90.0])
    Parameter_Range_Soil[:,1] = numpy.asarray([3.0,80.0])
    Parameter_Range_Soil[:,2] = numpy.asarray([1.0,130.0])
    Parameter_Range_Soil[:,3] = numpy.asarray([0.01,0.907])
    Parameter_Range_Soil[:,4] = numpy.asarray([1,20])
    
    
    Parameter_Range_Veg[:,0] = numpy.asarray([4.0,13.0])
    Parameter_Range_Veg[:,1] = numpy.asarray([1.0,5.0])
    Parameter_Range_Veg[:,2] = numpy.asarray([0.6,0.8])
    Parameter_Range_Veg[:,3] = numpy.asarray([0.04,0.15])
    Parameter_Range_Veg[:,4] = numpy.asarray([-450000.0,-200000.0])
    Parameter_Range_Veg[:,5] = numpy.asarray([-90000.0,-30000.0])
    Parameter_Range_Veg[:,6] = numpy.asarray([0.3,0.5])
    Parameter_Range_Veg[:,7] = numpy.asarray([0.05,0.15])
    Parameter_Range_Veg[:,8] = numpy.asarray([0.05,0.4])
    Parameter_Range_Veg[:,9] = numpy.asarray([0.025,0.075])
    Parameter_Range_Veg[:,10] = numpy.asarray([0.3,0.6])
    Parameter_Range_Veg[:,11] = numpy.asarray([0.1,0.35])
    Parameter_Range_Veg[:,12] = numpy.asarray([0.005,0.3])
    Parameter_Range_Veg[:,13] = numpy.asarray([0.005,0.15])
    Parameter_Range_Veg[:,14] = numpy.asarray([0.5,5.0])
    
    MONTHLY_LAI_Empirical = numpy.asarray([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                        [8.76, 9.16, 9.827, 10.093, 10.36, 10.76, 10.493, 10.227, 10.093, 9.827, 9.16, 8.76],
                                        [8.76, 9.16, 9.827, 10.093, 10.36, 10.76, 10.493, 10.227, 10.093, 9.827, 9.16, 8.76],
                                        [8.76, 9.16, 9.827, 10.093, 10.36, 10.76, 10.493, 10.227, 10.093, 9.827, 9.16, 8.76],
                                        [5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117],
                                        [5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117, 5.117],
                                        [0.52, 0.52, 0.867, 2.107, 4.507, 6.773, 7.173, 6.507, 5.04, 2.173, 0.867, 0.52],
                                        [0.52, 0.52, 0.867, 2.107, 4.507, 6.773, 7.173, 6.507, 5.04, 2.173, 0.867, 0.52],
                                        [0.52, 0.52, 0.867, 2.107, 4.507, 6.773, 7.173, 6.507, 5.04, 2.173, 0.867, 0.52],
                                        [0.580555, 0.6290065, 0.628558, 0.628546, 0.919255, 1.7685454, 2.5506969, 2.5535975, 1.7286418, 0.9703975, 0.726358, 0.6290065],
                                        [0.3999679, 0.4043968, 0.3138257, 0.2232945, 0.2498679, 0.3300675, 0.4323964, 0.7999234, 1.1668827, 0.7977234, 0.5038257, 0.4043968],
                                        [0.3999679, 0.4043968, 0.3138257, 0.2232945, 0.2498679, 0.3300675, 0.4323964, 0.7999234, 1.1668827, 0.7977234, 0.5038257, 0.4043968],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893],
                                        [0.782, 0.893, 1.004, 1.116, 1.782, 3.671, 4.782, 4.227, 2.004, 1.227, 1.004, 0.893]])
    
    Par_Index_Increment_Hard_Par = numpy.zeros((2, Dim_Hard_Par),dtype=numpy.float32)
    Par_Index_Increment_Hard_Par[:,0] = numpy.asarray([-0.25,0.25])
    Par_Index_Increment_Hard_Par[:,1] = numpy.asarray([-0.25,0.25])
    Par_Index_Increment_Hard_Par[:,2] = numpy.asarray([-0.25,0.25])
    
    Par_Index_Increment_Soil_Par = numpy.zeros((2, Dim_Soil_Par),dtype=numpy.float32)
    Par_Index_Increment_Soil_Par[:,0] = numpy.asarray([-10.0,10.0])
    Par_Index_Increment_Soil_Par[:,1] = numpy.asarray([-10.0,10.0])
    Par_Index_Increment_Soil_Par[:,2] = numpy.asarray([-10.0,10.0])
    Par_Index_Increment_Soil_Par[:,3] = numpy.asarray([-0.1,0.1])
    Par_Index_Increment_Soil_Par[:,4] = numpy.asarray([-2.0,2.0])
    
    Parameter_Range_Veg[:,10] = numpy.asarray([0.3,0.6])
    Parameter_Range_Veg[:,11] = numpy.asarray([0.1,0.35])
    Parameter_Range_Veg[:,12] = numpy.asarray([0.005,0.3])
    Parameter_Range_Veg[:,13] = numpy.asarray([0.005,0.15])
    Parameter_Range_Veg[:,14] = numpy.asarray([0.5,1.0])
    
    Par_Index_Increment_Veg_Par = numpy.zeros((2, Dim_Veg_Par),dtype=numpy.float32)
    Par_Index_Increment_Veg_Par[:,0] = numpy.asarray([-0.5,0.5]) * 2.0
    Par_Index_Increment_Veg_Par[:,1] = numpy.asarray([-0.2,0.2]) * 2.0
    Par_Index_Increment_Veg_Par[:,2] = numpy.asarray([-0.025,0.025]) * 2.0
    Par_Index_Increment_Veg_Par[:,3] = numpy.asarray([-0.005,0.005]) * 2.0
    Par_Index_Increment_Veg_Par[:,4] = numpy.asarray([-200.0,200.0]) * 2.0
    Par_Index_Increment_Veg_Par[:,5] = numpy.asarray([-20.0,20.0]) * 2.0
    Par_Index_Increment_Veg_Par[:,6] = numpy.asarray([-0.02,0.02]) * 2.0
    Par_Index_Increment_Veg_Par[:,7] = numpy.asarray([-0.005,0.005]) * 2.0
    Par_Index_Increment_Veg_Par[:,8] = numpy.asarray([-0.005,0.005]) * 2.0
    Par_Index_Increment_Veg_Par[:,9] = numpy.asarray([-0.002,0.002]) * 2.0
    Par_Index_Increment_Veg_Par[:,10] = numpy.asarray([-0.025,0.025]) * 2.0
    Par_Index_Increment_Veg_Par[:,11] = numpy.asarray([-0.005,0.005]) * 2.0
    Par_Index_Increment_Veg_Par[:,12] = numpy.asarray([-0.002,0.002]) * 2.0
    Par_Index_Increment_Veg_Par[:,13] = numpy.asarray([-0.002,0.002]) * 2.0
    Par_Index_Increment_Veg_Par[:,14] = numpy.asarray([-0.025,0.025]) * 2.0
    
    Par_Index_Increment_PFT_Par = numpy.zeros((2, Dim_PFT_Par),dtype=numpy.float32)
    Par_Index_Increment_PFT_Par[:,0] = numpy.asarray([-0.25,0.25])
    Par_Index_Increment_PFT_Par[:,1] = numpy.asarray([-0.25,0.25])
    Par_Index_Increment_PFT_Par[:,2] = numpy.asarray([-0.25,0.25])
    
    print "*******************************Parameter Perturbation Preparation*******************************"
    Parameter_ParFlow_Space_Single = numpy.zeros((Dim_ParFlow_Par, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Parameter_Soil_Space_Single = numpy.zeros((Dim_Soil_Par, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Parameter_Veg_Space_Single = numpy.zeros((Dim_Veg_Par, maxpft), dtype=numpy.float32)
    Parameter_PFT_Space_Single = numpy.zeros((Dim_PFT_Par, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Parameter_Hard_Space_Single = numpy.zeros((Dim_Hard_Par, Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    mksurfdata_NC_FileName_In = DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name
    pft_physiology_file_name = DAS_Data_Path + 'SysModel/CLM/inputdata/lnd/clm2/pftdata/'+fpftcon_name
    print("*****************************Read the surfdata.nc")
    mksurfdata_NC_File_In = netCDF4.Dataset(mksurfdata_NC_FileName_In, 'r')
    #print(mksurfdata_NC_File_In)
    
    #print(mksurfdata_NC_File_In)
    Sand_Ratio = numpy.zeros((len(mksurfdata_NC_File_In.dimensions['nlevsoi']),Row_Numbers,Col_Numbers))
    Clay_Ratio = numpy.zeros((len(mksurfdata_NC_File_In.dimensions['nlevsoi']),Row_Numbers,Col_Numbers))
    Organic_Ratio = numpy.zeros((8,Row_Numbers,Col_Numbers),dtype=numpy.float32)
    
    Sand_Ratio[0,:,:] = 1.0
    Clay_Ratio[0,:,:] = 1.0
    Organic_Ratio[0,:,:] = 1.0
    
    
    for i in range(1,len(mksurfdata_NC_File_In.dimensions['nlevsoi']),1):
        Sand_Ratio[i,:,:] = numpy.asarray(mksurfdata_NC_File_In.variables["PCT_SAND"][i,:,:])/numpy.asarray(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:])
        numexpr_a = numpy.asarray(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:])
        numexpr_b = 0.0
        numexpr_c = numpy.where(numexpr_a == numexpr_b)
        #print numexpr_c
        Sand_Ratio[i,:,:][numexpr_c] = 1.0
        Clay_Ratio[i,:,:] = numpy.asarray(mksurfdata_NC_File_In.variables["PCT_CLAY"][i,:,:])/numpy.asarray(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:])
        numexpr_a = numpy.asarray(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:])
        numexpr_b = 0.0
        numexpr_c = numpy.where(numexpr_a == numexpr_b)
        Clay_Ratio[i,:,:][numexpr_c] = 1.0
    for i in range(1,8,1):
        Organic_Ratio[i,:,:] = numpy.asarray(mksurfdata_NC_File_In.variables["ORGANIC"][i,:,:]) / numpy.asarray(mksurfdata_NC_File_In.variables["ORGANIC"][0,:,:])
        numexpr_a = numpy.asarray(mksurfdata_NC_File_In.variables["ORGANIC"][0,:,:])
        numexpr_b = 0.0
        numexpr_c = numpy.where(numexpr_a == numexpr_b)
        Organic_Ratio[i,:,:][numexpr_c] = 1.0            
    
    #print numpy.min(Sand_Ratio),numpy.max(Sand_Ratio)
    #os.abort()
    #####################
    r.assign("variance",100)
    r.assign("Ensemble_Number",3)                    
    r.assign("seed_value",123456)
    r.assign("Correlation_Scale",numpy.sqrt(numpy.square(Row_Numbers)+numpy.square(Col_Numbers))/2)
    r('set.seed(seed_value)')
    r('grf_data <- grf(n = Row_Numbers*Col_Numbers, grid = "reg", nx = Col_Numbers, ny = Row_Numbers, xlims = c(1, Col_Numbers), ylims = c(1, Row_Numbers), nsim = Ensemble_Number, cov.model = "exponential", \
        cov.pars = c(variance, Correlation_Scale), kappa = 0.5, nugget = 0.0, mean = 0, RF=TRUE)')
    #print r('dim(grf_data$data)')
    #print numpy.shape(numpy.array(r['grf_data$data']))
    GaussRF_Array = numpy.abs(numpy.asarray(r['grf_data$data']))
    GaussRF_Array = numpy.round(GaussRF_Array)
    numexpr_a = GaussRF_Array
    numexpr_b = 0
    numexpr_c = numpy.where(numexpr_a < numexpr_b)
    GaussRF_Array[numexpr_c] = 0
    numexpr_a = GaussRF_Array
    numexpr_b = 5
    numexpr_c = numpy.where(numexpr_a > numexpr_b)
    GaussRF_Array[numexpr_c] = 5
    
    Parameter_Soil_Space_Single[0,::] = numpy.flipud(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:])
    Parameter_Soil_Space_Single[1,::] = numpy.flipud(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:])
    Parameter_Soil_Space_Single[2,::] = numpy.flipud(mksurfdata_NC_File_In.variables["ORGANIC"][0,:,:])
    Parameter_Soil_Space_Single[Dim_Soil_Par-2,::] = numpy.flipud(mksurfdata_NC_File_In.variables['FMAX'][::])
    Parameter_Soil_Space_Single[Dim_Soil_Par-1,::] = numpy.flipud(mksurfdata_NC_File_In.variables["SOIL_COLOR"][::])
    
    Soil_Sand_Clay_Sum[0,::] = numpy.flipud(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:]) + numpy.flipud(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:])
    
    Parameter_PFT_Space_Single[0,::] = 1.0 # factor
    Parameter_PFT_Space_Single[1,::] = 1.0 # factor
    Parameter_PFT_Space_Single[2,::] = 1.0 # factor
    
    mksurfdata_NC_File_In.close()
    
    print("*****************************Read the pft-physiology.nc file")
    pft_physiology_file = netCDF4.Dataset(pft_physiology_file_name, 'r')
    pft_len = len(pft_physiology_file.dimensions['pft'])
    
    if Def_Debug:    # For Parameter Optimization and Synthetic study
        Parameter_Veg_Space_Single[0,:] = pft_physiology_file.variables["roota_par"][0:maxpft]# * 2.0
    else:
        Parameter_Veg_Space_Single[0,:] = pft_physiology_file.variables["roota_par"][0:maxpft]
        
    Parameter_Veg_Space_Single[1,:] = pft_physiology_file.variables["rootb_par"][0:maxpft]
    Parameter_Veg_Space_Single[2,:] = pft_physiology_file.variables["displar"][0:maxpft]
    
    if Def_Debug:    # For Parameter Optimization and Synthetic study
        Parameter_Veg_Space_Single[3,:] = pft_physiology_file.variables["z0mr"][0:maxpft]# * 2.0
    else:
        Parameter_Veg_Space_Single[3,:] = pft_physiology_file.variables["z0mr"][0:maxpft]
        
    Parameter_Veg_Space_Single[4,:] = pft_physiology_file.variables["smpsc"][0:maxpft]
    Parameter_Veg_Space_Single[5,:] = pft_physiology_file.variables["smpso"][0:maxpft]
    Parameter_Veg_Space_Single[6,:] = pft_physiology_file.variables["rholnir"][0:maxpft]
    Parameter_Veg_Space_Single[7,:] = pft_physiology_file.variables["rholvis"][0:maxpft]
    Parameter_Veg_Space_Single[8,:] = pft_physiology_file.variables["taulnir"][0:maxpft]
    Parameter_Veg_Space_Single[9,:] = pft_physiology_file.variables["taulvis"][0:maxpft]
    Parameter_Veg_Space_Single[10,:] = pft_physiology_file.variables["rhosnir"][0:maxpft]
    Parameter_Veg_Space_Single[11,:] = pft_physiology_file.variables["rhosvis"][0:maxpft]
    Parameter_Veg_Space_Single[12,:] = pft_physiology_file.variables["tausnir"][0:maxpft]
    Parameter_Veg_Space_Single[13,:] = pft_physiology_file.variables["tausvis"][0:maxpft]
    Parameter_Veg_Space_Single[14,:] = pft_physiology_file.variables["flnr"][0:maxpft]
    Parameter_Veg_Space_Single[15,:] = pft_physiology_file.variables["slatop"][0:maxpft]
    pft_physiology_file.close()


    print "Prepare Parameter_Space_Single.nc"
    if os.path.exists(NC_FileName_Parameter_Space_Single):
        os.remove(NC_FileName_Parameter_Space_Single)
            
    if Def_Print:
            print 'Write NetCDF File:',NC_FileName_Parameter_Space_Single
            
    NC_File_Parameter_Space_Single = netCDF4.Dataset(NC_FileName_Parameter_Space_Single, 'w', diskless=True, persist=True, format='NETCDF4')
    # Dim the dimensions of NetCDF
    NC_File_Parameter_Space_Single.createDimension('lon', Col_Numbers)
    NC_File_Parameter_Space_Single.createDimension('lat', Row_Numbers)
    NC_File_Parameter_Space_Single.createDimension('Dim_Soil_Par', Dim_Soil_Par)
    NC_File_Parameter_Space_Single.createDimension('Dim_Veg_Par', Dim_Veg_Par)
    NC_File_Parameter_Space_Single.createDimension('Dim_PFT_Par', Dim_PFT_Par)
    NC_File_Parameter_Space_Single.createDimension('Dim_Hard_Par', Dim_Hard_Par)
    NC_File_Parameter_Space_Single.createDimension('Dim_ParFlow_Par', Dim_ParFlow_Par)  # Perm Alpha N SRes SSat
    NC_File_Parameter_Space_Single.createDimension('maxpft', maxpft)
    NC_File_Parameter_Space_Single.createDimension('GaussRF_Array_Row', Row_Numbers*Col_Numbers)
    NC_File_Parameter_Space_Single.createDimension('GaussRF_Array_Col', 3)
    NC_File_Parameter_Space_Single.createDimension('Texture_Layer_Num', 10)
    NC_File_Parameter_Space_Single.createDimension('Organic_Layer_Num', 8)
    NC_File_Parameter_Space_Single.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
    NC_File_Parameter_Space_Single.createDimension('Month', 12)
    
    NC_File_Parameter_Space_Single.createVariable('Parameter_ParFlow_Space_Single','f4',('Dim_ParFlow_Par','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Parameter_ParFlow_Space_Single'][:,:,:] = Parameter_ParFlow_Space_Single
    
    NC_File_Parameter_Space_Single.createVariable('Parameter_Soil_Space_Single','f4',('Dim_Soil_Par','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][:,:,:] = Parameter_Soil_Space_Single
    
    NC_File_Parameter_Space_Single.createVariable('Parameter_Veg_Space_Single','f4',('Dim_Veg_Par','maxpft',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Parameter_Veg_Space_Single'][:,:] = Parameter_Veg_Space_Single
    
    NC_File_Parameter_Space_Single.createVariable('Parameter_PFT_Space_Single','f4',('Dim_PFT_Par','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Parameter_PFT_Space_Single'][:,:,:] = Parameter_PFT_Space_Single
    
    NC_File_Parameter_Space_Single.createVariable('Parameter_Hard_Space_Single','f4',('Dim_Hard_Par','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Parameter_Hard_Space_Single'][:,:,:] = Parameter_Hard_Space_Single
    
    NC_File_Parameter_Space_Single.createVariable('GaussRF_Array','f4',('GaussRF_Array_Row','GaussRF_Array_Col',),zlib=True)
    NC_File_Parameter_Space_Single.variables['GaussRF_Array'][:,:] = GaussRF_Array
    
    NC_File_Parameter_Space_Single.createVariable('Sand_Ratio','f4',('Texture_Layer_Num','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Sand_Ratio'][:,:,:] = Sand_Ratio
    
    NC_File_Parameter_Space_Single.createVariable('Clay_Ratio','f4',('Texture_Layer_Num','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Clay_Ratio'][:,:,:] = Clay_Ratio
    
    NC_File_Parameter_Space_Single.createVariable('Organic_Ratio','f4',('Organic_Layer_Num','lat','lon',),zlib=True)
    NC_File_Parameter_Space_Single.variables['Organic_Ratio'][:,:,:] = Organic_Ratio
    
    NC_File_Parameter_Space_Single.createVariable('MONTHLY_LAI_Empirical','f4',('maxpft','Month',),zlib=True)
    NC_File_Parameter_Space_Single.variables['MONTHLY_LAI_Empirical'][:,:] = MONTHLY_LAI_Empirical[0:maxpft,:]
    
    
    NC_File_Parameter_Space_Single.sync()            
    NC_File_Parameter_Space_Single.close()
    
    del Sand_Ratio, Clay_Ratio, Organic_Ratio, GaussRF_Array   
    print "################## Soil Parameters"
    '''
    Soil organic matter tends to increase as the clay content increases. 
    This increase depends on two mechanisms. 
    First, bonds between the surface of clay particles and organic matter retard the decomposition process. 
    Second, soils with higher clay content increase the potential for aggregate formation.
    '''
    Par_Index_Increment_Soil = numpy.zeros((Ensemble_Number,Dim_Soil_Par),dtype=numpy.float32)
    for Par_Index in range(Dim_Soil_Par):
        numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str((Par_Index))))
        Par_Index_Increment_Soil[:,Par_Index] = numpy.random.uniform(Par_Index_Increment_Soil_Par[0,Par_Index],Par_Index_Increment_Soil_Par[1,Par_Index],size=Ensemble_Number)
    
    
    # change the sign of increment   
    Par_Index_Increment_Soil_Sign = numpy.zeros((Ensemble_Number,Dim_Soil_Par),dtype=numpy.float32)
    Par_Index_Increment_Soil_Sign[::] = 1
    
    Par_Index_Increment_Soil_Sign[:,0][numpy.where(Par_Index_Increment_Soil[:,0] < 0)] = -1
     
    Par_Index_Increment_Soil_Sign[:,1] = -1 * Par_Index_Increment_Soil_Sign[:,0]
                 
    Par_Index_Increment_Soil_Sign[:,2] = -1 * Par_Index_Increment_Soil_Sign[:,0]
     
    Par_Index_Increment_Soil[:,0] = numpy.abs(Par_Index_Increment_Soil[:,0]) * Par_Index_Increment_Soil_Sign[:,0]
    Par_Index_Increment_Soil[:,1] = numpy.abs(Par_Index_Increment_Soil[:,1]) * Par_Index_Increment_Soil_Sign[:,1]
    Par_Index_Increment_Soil[:,2] = numpy.abs(Par_Index_Increment_Soil[:,2]) * Par_Index_Increment_Soil_Sign[:,2]
    
    Par_Index_Increment_PFT = numpy.zeros((Ensemble_Number,Dim_PFT_Par),dtype=numpy.float32)
    for Par_Index in range(Dim_PFT_Par):
        numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str(Par_Index)+str(Dim_PFT_Par)))
        Par_Index_Increment_PFT[:,Par_Index] = numpy.random.uniform(Par_Index_Increment_PFT_Par[0,Par_Index],Par_Index_Increment_PFT_Par[1,Par_Index],size=Ensemble_Number)
    
    # change the sign of increment   
    Par_Index_Increment_PFT_Sign = numpy.zeros((Ensemble_Number,Dim_PFT_Par),dtype=numpy.float32)
    Par_Index_Increment_PFT_Sign[::] = 1
    
    Par_Index_Increment_PFT_Sign[:,0][numpy.where(Par_Index_Increment_PFT[:,0] < 0)] = -1
    
    Par_Index_Increment_PFT_Sign[:,1] = 1 * Par_Index_Increment_PFT_Sign[:,0]
    
    Par_Index_Increment_PFT_Sign[:,2] = 1 * Par_Index_Increment_PFT_Sign[:,0]
    
    Par_Index_Increment_PFT[:,0] = numpy.abs(Par_Index_Increment_PFT[:,0]) * Par_Index_Increment_PFT_Sign[:,0]
    Par_Index_Increment_PFT[:,1] = numpy.abs(Par_Index_Increment_PFT[:,1]) * Par_Index_Increment_PFT_Sign[:,1]
    Par_Index_Increment_PFT[:,2] = numpy.abs(Par_Index_Increment_PFT[:,2]) * Par_Index_Increment_PFT_Sign[:,2]
    
    Par_Index_Increment_PFT_Predict = numpy.zeros((Ensemble_Number_Predict,Dim_PFT_Par),dtype=numpy.float32)
    for Par_Index in range(Dim_PFT_Par):
        numpy.random.seed(seed=string.atoi(str((Datetime_Start - Datetime_Initial).days)+str((Par_Index))+str(Dim_PFT_Par)))
        Par_Index_Increment_PFT_Predict[:,Par_Index] = numpy.random.uniform(Par_Index_Increment_PFT_Par[0,Par_Index],Par_Index_Increment_PFT_Par[1,Par_Index],size=Ensemble_Number_Predict)
    
    
    ###########################
    Par_Soil_Uniform_STD = numpy.zeros(Dim_Soil_Par,dtype=numpy.float32)
    for Par_Index in range(Dim_Soil_Par):
        Par_Soil_Uniform_STD[Par_Index] = numpy.sqrt(numpy.square(numpy.abs(Par_Index_Increment_Soil_Par[0,Par_Index])*2.0)/12.0)
    
    Par_Veg_Uniform_STD = numpy.zeros(Dim_Veg_Par,dtype=numpy.float32)
    for Par_Index in range(Dim_Veg_Par):
        Par_Veg_Uniform_STD[Par_Index] = numpy.sqrt(numpy.square(numpy.abs(Par_Index_Increment_Veg_Par[0,Par_Index])*2.0)/12.0)
    
    Par_PFT_Uniform_STD = numpy.zeros(Dim_PFT_Par,dtype=numpy.float32)
    for Par_Index in range(Dim_PFT_Par):
        Par_PFT_Uniform_STD[Par_Index] = numpy.sqrt(numpy.square(numpy.abs(Par_Index_Increment_PFT_Par[0,Par_Index])*2.0)/12.0)
    
    Par_Hard_Uniform_STD = numpy.zeros(Dim_Hard_Par,dtype=numpy.float32)
    for Par_Index in range(Dim_Hard_Par):
        Par_Hard_Uniform_STD[Par_Index] = numpy.sqrt(numpy.square(numpy.abs(Par_Index_Increment_Hard_Par[0,Par_Index])*2.0)/12.0)
    
    if Def_First_Run == 1 and Ensemble_Number > 1:
        
        print "*************************** Begin Parameter Generation"
        NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r+')
        
        
        print "Perturb the soil peroperties with spatial correlated noise"
        # change the sign of increment   
        Par_Index_Increment_Sign = numpy.zeros((Ensemble_Number,Dim_Soil_Par,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Par_Index_Increment_Sign[:,:,:,:] = 1
         
        Par_Index_Increment_Noise = numpy.zeros((Ensemble_Number,Dim_Soil_Par,Row_Numbers,Col_Numbers),dtype=numpy.float32)
         
        for Par_Index in range(Dim_Soil_Par):
            #r.assign("variance",numpy.square(numpy.abs(Par_Index_Increment_Soil_Par[0,Par_Index])/2.0))
            r.assign("variance",numpy.square(numpy.abs(Par_Index_Increment_Soil_Par[0,Par_Index])*2.0)/12.0)
            r.assign("Ensemble_Number",Ensemble_Number)                    
            r.assign("seed_value",str(Par_Index+1))
            r.assign("Correlation_Scale",numpy.sqrt(numpy.square(Row_Numbers)+numpy.square(Col_Numbers))/2)
            r('set.seed(seed_value)')
            r('grf_data <- grf(n = Row_Numbers*Col_Numbers, grid = "reg", nx = Col_Numbers, ny = Row_Numbers, xlims = c(1, Col_Numbers), ylims = c(1, Row_Numbers), nsim = Ensemble_Number, cov.model = "exponential", \
                cov.pars = c(variance, Correlation_Scale), \
                kappa = 0.5, nugget = 0.0, mean = 0, RF=TRUE)')
            #print r('dim(grf_data$data)')
            #print numpy.shape(numpy.array(r['grf_data$data']))
            GaussRF_Array = numpy.asarray(r['grf_data$data'])
                         
            for Ens_Index in range(Ensemble_Number):
                Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:] = numpy.reshape(GaussRF_Array[:,Ens_Index],(Row_Numbers,Col_Numbers))
                #print numpy.min(Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:]),numpy.max(Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:])
                Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:][numpy.where(Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:] > Par_Index_Increment_Soil_Par[1,Par_Index])] = Par_Index_Increment_Soil_Par[1,Par_Index]
                Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:][numpy.where(Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:] < Par_Index_Increment_Soil_Par[0,Par_Index])] = Par_Index_Increment_Soil_Par[0,Par_Index]                    
             
        for Ens_Index in range(Ensemble_Number):               
            # change the sign of increment, increase sand, decrease clay and organic
             
            Par_Index_Increment_Sign[:,0,:,:][numpy.where(Par_Index_Increment_Noise[:,0,:,:] < 0)] = -1            
            Par_Index_Increment_Sign[:,1,:,:] = -1 * Par_Index_Increment_Sign[:,0,:,:]            
            Par_Index_Increment_Sign[:,2,:,:] = -1 * Par_Index_Increment_Sign[:,0,:,:]
             
            Par_Index_Increment_Noise[:,0,:,:] = numpy.abs(Par_Index_Increment_Noise[:,0,:,:]) * Par_Index_Increment_Sign[:,0,:,:]
            Par_Index_Increment_Noise[:,1,:,:] = numpy.abs(Par_Index_Increment_Noise[:,1,:,:]) * Par_Index_Increment_Sign[:,1,:,:]
            Par_Index_Increment_Noise[:,2,:,:] = numpy.abs(Par_Index_Increment_Noise[:,2,:,:]) * Par_Index_Increment_Sign[:,2,:,:]
        
        Parameter_Soil_Space_Ensemble = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:])
        
        for Par_Index in range(Dim_Soil_Par):
            if (Par_Index <= 2) or ((Par_Index > 2) and (Soil_Par_Sens_Array[0][Par_Index] or Soil_Par_Sens_Array[1][Par_Index])):
                for Ens_Index in range(Ensemble_Number):
                    Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = Parameter_Soil_Space_Single[Par_Index,::] + Par_Index_Increment_Noise[Ens_Index,Par_Index,:,:]
                    
                    numexpr_a = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]
                    numexpr_b = Parameter_Range_Soil[0,Par_Index]
                    numexpr_c = numpy.where(numexpr_a < numexpr_b)
                    #print numexpr_b,numexpr_c
                    numexpr_a[numexpr_c] = Parameter_Range_Soil[0,Par_Index]
                    Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = numexpr_a
                    #print numpy.min(Parameter_Soil_Space_Ensemble[:,Par_Index,:,:])
                    numexpr_a = Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:]
                    numexpr_b = Parameter_Range_Soil[1,Par_Index]
                    numexpr_c = numpy.where(numexpr_a > numexpr_b)
                    numexpr_a[numexpr_c] = Parameter_Range_Soil[1,Par_Index]
                    Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = numexpr_a
            else:
                for Ens_Index in range(Ensemble_Number):
                    Parameter_Soil_Space_Ensemble[Ens_Index,Par_Index,:,:] = Parameter_Soil_Space_Single[Par_Index,::]
        
        NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:] = Parameter_Soil_Space_Ensemble
        del Parameter_Soil_Space_Ensemble
        
        del Par_Index_Increment_Sign, Par_Index_Increment_Noise
        
        print "###################### Check Texture 1"  
        if Def_ParFor:
            #print Dim_Soil_Par, Ensemble_Number, Row_Numbers, Col_Numbers, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Soil_Space_Ensemble, DAS_Depends_Path, omp_get_num_procs_ParFor
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:] = ParFor.ParFor_Texture_Check(Dim_Soil_Par, Ensemble_Number, Row_Numbers, Col_Numbers, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:], DAS_Depends_Path, omp_get_num_procs_ParFor)
        else:
            Parameter_Soil_Space_Ensemble = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:])
            
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
    #                        if Texture_Sum < Soil_Sand_Clay_Sum[Soil_Layer_Index_Sub,Row_Index,Col_Index]:
    #                            Ratio = Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index] / (Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index]+Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index])
    #                            Diff = Soil_Sand_Clay_Sum[Soil_Layer_Index_Sub,Row_Index,Col_Index] - Texture_Sum
    #                            Diff_Part1 = Ratio*Diff
    #                            Diff_Part2 = Diff - Diff_Part1
    #                            Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub,Row_Index,Col_Index] += Diff_Part1
    #                            Parameter_Soil_Space_Ensemble[Ens_Index,Soil_Layer_Index_Sub+Soil_Texture_Layer_Opt_Num,Row_Index,Col_Index] += Diff_Part2
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:] = Parameter_Soil_Space_Ensemble
            del Parameter_Soil_Space_Ensemble
            
        NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:] = numpy.asarray(numpy.round(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:]),dtype=numpy.integer)
        numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:])
        numexpr_b = 1
        numexpr_c = numpy.where(numexpr_a < numexpr_b)
        numexpr_a[numexpr_c] = 1
        NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:] = numexpr_a
        numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:])
        numexpr_b = 20
        numexpr_c = numpy.where(numexpr_a > numexpr_b)
        numexpr_a[numexpr_c] = 20
        NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,Dim_Soil_Par-1,:,:] = numexpr_a
        numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,2*Soil_Texture_Layer_Opt_Num:3*Soil_Texture_Layer_Opt_Num,:,:])
        numexpr_b = 130.0
        numexpr_c = numpy.where(numexpr_a > numexpr_b)
        numexpr_a[numexpr_c] = 130.0
        NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,2*Soil_Texture_Layer_Opt_Num:3*Soil_Texture_Layer_Opt_Num,:,:] = numexpr_a
        
        
        if (numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) >= 1):
            Parameter_PFT_Space_Ensemble_Temp = numpy.zeros_like(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:],dtype=numpy.float32)
            
            for Par_Index in range(Dim_PFT_Par):
                
                if PFT_Par_Sens_Array[0][Par_Index] or PFT_Par_Sens_Array[1][Par_Index]:
                    for Ens_Index in range(Ensemble_Number):
                        Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:] = Parameter_PFT_Space_Single[Par_Index,::] + Par_Index_Increment_PFT[Ens_Index,Par_Index]
                        #Parameter_PFT_Space_Ensemble_Temp[:,Par_Index,:,:] = numpy.random.uniform(Parameter_Range_PFT[0,Par_Index], Parameter_Range_PFT[1,Par_Index])
                        #print "before", numpy.min(Parameter_PFT_Space_Ensemble_Temp[:,Par_Index,:,:])
                        numexpr_a = numpy.asarray(Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:])
                        numexpr_b = Parameter_Range_PFT[0,Par_Index]
                        numexpr_c = numpy.where(numexpr_a < numexpr_b)
                        #print numexpr_b,numexpr_c
                        numexpr_a[numexpr_c] = Parameter_Range_PFT[0,Par_Index]
                        Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:] = numexpr_a
                        #print numpy.min(Parameter_PFT_Space_Ensemble_Temp[:,Par_Index,:,:])
                        numexpr_a = numpy.asarray(Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:])
                        numexpr_b = Parameter_Range_PFT[1,Par_Index]
                        numexpr_c = numpy.where(numexpr_a > numexpr_b)
                        numexpr_a[numexpr_c] = Parameter_Range_PFT[1,Par_Index]
                        Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:] = numexpr_a
                        if Parameter_Optimization == 1:
                            for Ens_Index in range(Ensemble_Number_Predict):
                                NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:] = Parameter_PFT_Space_Single[Par_Index,::] + Par_Index_Increment_PFT_Predict[Ens_Index,Par_Index]
                                #NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][:,Par_Index,:,:] = numpy.random.uniform(Parameter_Range_PFT[0,Par_Index], Parameter_Range_PFT[1,Par_Index])
                                #print "before", numpy.min(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][:,Par_Index,:,:])
                                numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:])
                                numexpr_b = Parameter_Range_PFT[0,Par_Index]
                                numexpr_c = numpy.where(numexpr_a < numexpr_b)
                                #print numexpr_b,numexpr_c
                                numexpr_a[numexpr_c] = Parameter_Range_PFT[0,Par_Index]
                                NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:] = numexpr_a
                                #print numpy.min(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][:,Par_Index,:,:])
                                numexpr_a = numpy.asarray(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:])
                                numexpr_b = Parameter_Range_PFT[1,Par_Index]
                                numexpr_c = numpy.where(numexpr_a > numexpr_b)
                                numexpr_a[numexpr_c] = Parameter_Range_PFT[1,Par_Index]
                                NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:] = numexpr_a
                else:
                    for Ens_Index in range(Ensemble_Number):
                        Parameter_PFT_Space_Ensemble_Temp[Ens_Index,Par_Index,:,:] = Parameter_PFT_Space_Single[Par_Index,::]
                    if Parameter_Optimization == 1:
                        for Ens_Index in range(Ensemble_Number_Predict):
                            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble_Predict'][Ens_Index,Par_Index,:,:] = Parameter_PFT_Space_Single[Par_Index,::]
                
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:] = Parameter_PFT_Space_Ensemble_Temp
            del Parameter_PFT_Space_Ensemble_Temp
                
        NC_File_Out_Assimilation_2_Parameter.sync()
        NC_File_Out_Assimilation_2_Parameter.close()
    
    numexpr_a = []
    numexpr_b = []
    numexpr_c = []
    
    del Parameter_ParFlow_Space_Single, Parameter_Soil_Space_Single, Parameter_Veg_Space_Single, Parameter_PFT_Space_Single
    
    gc.collect()
    del gc.garbage[:]
    
    print "*************************** Finish Parameter Generation"
    return Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Par_Index_Increment_Soil_Par, Par_Soil_Uniform_STD, Par_Veg_Uniform_STD, Par_PFT_Uniform_STD, Par_Hard_Uniform_STD


def Run_CMEM(Ens_Index, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, \
             NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, \
             Forcing_File_Path_Array, history_file_name, CMEM_Work_Path_Array, mksrf_edgen, mksrf_edges, mksrf_edgew, mksrf_edgee, \
             Grid_Resolution_GEO, Row_Numbers, Col_Numbers, Datetime_Start, Datetime_Stop, Datetime_Stop_Init,
             Def_Region, DAS_Data_Path, Region_Name, Constant_File_Name, Soil_Layer_Num, NAvalue, CLM_NA, DasPy_Path, \
             Def_Print, Density_of_liquid_water, Run_Dir_Array, Row_Numbers_String, Col_Numbers_String, Station_XY_Index,
             Tair_Time, Soil_Texture_Layer_Opt_Num, SensorType_Name, Mask_Index, Station_XY, DAS_Output_Path):
    
    Call_CMEM = imp.load_source("Call_CMEM",DasPy_Path+"ObsModel/CMEM/Call_CMEM.py")
    
    NC_FileName_CMEM_Par = DAS_Output_Path+"Analysis/"+Region_Name+"/CMEM_Par.nc"
    print 'Read NetCDF File:',NC_FileName_CMEM_Par
    NC_File_CMEM_Par = netCDF4.Dataset(NC_FileName_CMEM_Par, 'r')    
    ECOCVL_Mat = NC_File_CMEM_Par.variables['ECOCVL_Mat'][:,:]
    ECOCVH_Mat = NC_File_CMEM_Par.variables['ECOCVH_Mat'][:,:]
    ECOTVL_Mat = NC_File_CMEM_Par.variables['ECOTVL_Mat'][:,:]
    ECOTVH_Mat = NC_File_CMEM_Par.variables['ECOTVH_Mat'][:,:]
    ECOWAT_Mat = NC_File_CMEM_Par.variables['ECOWAT_Mat'][:,:]
    NC_File_CMEM_Par.close()
    
    # For evaluation parameter estimation only in CLM
    #NC_FileName_Assimilation_2_Parameter = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_First.nc"
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    DEM_Data = NC_File_Out_Assimilation_2_Constant.variables['DEM_Data'][:,:]
    NC_File_Out_Assimilation_2_Constant.close()
    
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,Ens_Index]
    #CLM_Soil_Ice_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Snow_Depth_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Snow_Water_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,Ens_Index]
    #CLM_ROOTFR_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_LAI = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][Ens_Index,12,:,:]
    NC_File_Out_Assimilation_2_Initial.close()
    
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')        
    CLM_2m_Air_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    NC_File_Out_Assimilation_2_Diagnostic.close()
    
    #NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
    #Parameter_Soil_Space_Ensemble = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][Ens_Index,:,:,:]
    #NC_File_Out_Assimilation_2_Parameter.close()

#     NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
#     CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,Ens_Index]
#     CLM_Soil_Temperature_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:],axis=3)
#     #CLM_Soil_Ice_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_Ensemble_Mat'][:,:,:,Ens_Index]
#     CLM_Vegetation_Temperature_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:],axis=2)
#     CLM_Ground_Temperature_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:],axis=2)
#     CLM_Snow_Depth_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:],axis=2)
#     CLM_Snow_Water_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:],axis=2)
#     #CLM_ROOTFR_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,Ens_Index]
#     CLM_LAI = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:,12,:,:],axis=0)
#     NC_File_Out_Assimilation_2_Initial.close()
#     
    #NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r+')        
    #CLM_2m_Air_Temperature_Ensemble_Mat = numpy.mean(NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:],axis=2)
    #NC_File_Out_Assimilation_2_Diagnostic.close()
    
    # Mean values are better for parameter estimation
    NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
    Parameter_Soil_Space_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:],axis=0)
    NC_File_Out_Assimilation_2_Parameter.close()
    
    Forcing_File_Path_CMEM = Forcing_File_Path_Array[Ens_Index]
    CMEM_Work_Path = CMEM_Work_Path_Array[Ens_Index]
    
    CLM_Soil_Moisture = CLM_Soil_Moisture_Ensemble_Mat
    CLM_Soil_Temperature = CLM_Soil_Temperature_Ensemble_Mat
    CLM_Ground_Temperature = CLM_Ground_Temperature_Ensemble_Mat
    CLM_Vegetation_Temperature = CLM_Vegetation_Temperature_Ensemble_Mat
    CLM_2m_Air_Temperature = CLM_2m_Air_Temperature_Ensemble_Mat
    CLM_Snow_Water = CLM_Snow_Water_Ensemble_Mat
    CLM_Snow_Depth = CLM_Snow_Depth_Ensemble_Mat
    
    Sand_Mat = Parameter_Soil_Space_Ensemble[0,::]
    Clay_Mat = Parameter_Soil_Space_Ensemble[0+Soil_Texture_Layer_Opt_Num,::]
    
    if Def_Print >= 2:
        print "numpy.max(CLM_Soil_Moisture),numpy.min(CLM_Soil_Moisture)",numpy.max(CLM_Soil_Moisture[:,~Mask_Index[0,::]]),numpy.min(CLM_Soil_Moisture[:,~Mask_Index[0,::]])
        print "numpy.max(CLM_Soil_Temperature),numpy.min(CLM_Soil_Temperature)",numpy.max(CLM_Soil_Temperature[:,~Mask_Index[0,::]]),numpy.min(CLM_Soil_Temperature[:,~Mask_Index[0,::]])
        print "numpy.max(CLM_Ground_Temperature),numpy.min(CLM_Ground_Temperature)",numpy.max(CLM_Ground_Temperature[~Mask_Index[0,::]]),numpy.min(CLM_Ground_Temperature[~Mask_Index[0,::]])
        print "numpy.max(CLM_Vegetation_Temperature),numpy.min(CLM_Vegetation_Temperature)",numpy.max(CLM_Vegetation_Temperature[~Mask_Index[0,::]]),numpy.min(CLM_Vegetation_Temperature[~Mask_Index[0,::]])
        print "numpy.max(CLM_Snow_Water),numpy.min(CLM_Snow_Water)",numpy.max(CLM_Snow_Water[~Mask_Index[0,::]]),numpy.min(CLM_Snow_Water[~Mask_Index[0,::]])
        print "numpy.max(CLM_Snow_Depth),numpy.min(CLM_Snow_Depth)",numpy.max(CLM_Snow_Depth[~Mask_Index[0,::]]),numpy.min(CLM_Snow_Depth[~Mask_Index[0,::]])
        print "numpy.max(Clay_Mat),numpy.min(Clay_Mat)",numpy.max(Clay_Mat[~Mask_Index[0,::]]),numpy.min(Clay_Mat[~Mask_Index[0,::]])
        print "numpy.max(Sand_Mat),numpy.min(Sand_Mat)",numpy.max(Sand_Mat[~Mask_Index[0,::]]),numpy.min(Sand_Mat[~Mask_Index[0,::]])
        
    Longitude = numpy.linspace(mksrf_edgew+Grid_Resolution_GEO[0]/2.0,mksrf_edgee-Grid_Resolution_GEO[0]/2.0,Col_Numbers)
    Latitude = numpy.linspace(mksrf_edges+Grid_Resolution_GEO[1]/2.0,mksrf_edgen-Grid_Resolution_GEO[1]/2.0,Row_Numbers)
        
    #NA_Index = CLM_Soil_Moisture[numpy.where(CLM_Soil_Moisture == NAvalue)]
    #RSN_Mat = CLM_Snow_Water / CLM_Snow_Depth       # snow density (kg/m3)
    RSN_Mat = 100.0 * numpy.ones_like((CLM_Snow_Water))      # snow density (kg/m3)
    SD_Mat = CLM_Snow_Water / 1000.0        # snow depth in water equiv (m)
    # soil layer thickness (default is TESSEL: 0.07, 0.21, 0.72)
    STL_Mat = CLM_Soil_Temperature[:, :, :]    # Soil temperature level (K)
    SWVL_Mat = CLM_Soil_Moisture[:, :, :]     # volumetric soil moist level (m3/m3)
    TSKIN_Mat = CLM_Ground_Temperature[::]  # Skin temperature (K)
    Tair_Mat = CLM_Vegetation_Temperature[::]
    ECOLAIL_Mat = CLM_LAI[::]
    
    Z_Mat = DEM_Data / 1000.0           # geopotential at surface in km
    Frequency = 1.4
    View_Angle = 0.0
    if SensorType_Name == 'SMOS':
        Frequency = 1.4
        View_Angle = 50.0
    
    print "CMEM_Work_Path",CMEM_Work_Path,SensorType_Name,Frequency,View_Angle
    TBH, TBV, EFFECTIVE_TEMP = Call_CMEM.Call_CMEM(Def_Print, CMEM_Work_Path, Row_Numbers, Col_Numbers, Latitude, Longitude, Tair_Mat, Tair_Time, Clay_Mat, Sand_Mat, ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat, ECOLAIL_Mat,
                                       RSN_Mat, SD_Mat, STL_Mat, SWVL_Mat, TSKIN_Mat, Z_Mat, Frequency, View_Angle)
    
    TBH_Flipud = numpy.flipud(TBH)
    TBH_Flipud[Mask_Index[0,::]] = -9999.0
    #print numpy.where(TBH_Flipud == -10000)
    #Model_Value[0] = TBH_Flipud
    #print "numpy.max(TBH_Flipud),numpy.min(TBH_Flipud)",numpy.max(TBH_Flipud),numpy.min(TBH_Flipud)
    
#     print "TBH_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]",TBH_Flipud[Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     
#     print "CLM_Soil_Moisture[:,Station_XY_Index[0][1],Station_XY_Index[0][0]]",CLM_Soil_Moisture[:,Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     
#     print "CLM_Soil_Temperature[:,Station_XY_Index[0][1],Station_XY_Index[0][0]]",CLM_Soil_Temperature[:,Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     
#     print "CLM_Ground_Temperature[Station_XY_Index[0][1],Station_XY_Index[0][0]]",CLM_Ground_Temperature[Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     #print numpy.shape(Tair_Mat),Tair_Mat
#     print "Tair_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]",Tair_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]
# #            
#     print "Clay_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]",Clay_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     print "Sand_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]",Sand_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     print "ECOCVL_Mat",ECOCVL_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]    # fraction of low veg
#     print "ECOCVH_Mat",ECOCVH_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]    # fraction of high veg
#     print "ECOTVL_Mat",ECOTVL_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]      # low veg type 
#     print "ECOTVH_Mat",ECOTVH_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]      # high veg type
#     print "ECOWAT_Mat",ECOWAT_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]    # water fraction
#     print "ECOLAIL_Mat",ECOLAIL_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]   # lai of low veg for each pixel
#     print "RSN_Mat",RSN_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]      # snow density (kg/m3)
#     print "SD_Mat",SD_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]
#     print "Z_Mat",Z_Mat[Station_XY_Index[0][1],Station_XY_Index[0][0]]
                
    if Def_Print >= 1:
        print "numpy.max(TBH_Flipud),numpy.min(TBH_Flipud)",numpy.max(TBH_Flipud[~Mask_Index[0,::]]),numpy.min(TBH_Flipud[~Mask_Index[0,::]])
        if Def_Print >= 3:
            for Staton_Index in range(numpy.size(Station_XY)/2):
                print "Sand_Mat",Sand_Mat[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]],"Clay_Mat",Clay_Mat[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]]
                print "Tair_Mat",Tair_Mat[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]],"ECOLAIL_Mat",ECOLAIL_Mat[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]]
                print "STL_Mat",STL_Mat[1,Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]],"STL_Mat",STL_Mat[2,Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]]
                print "SWVL_Mat",SWVL_Mat[1,Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]],"SWVL_Mat",SWVL_Mat[2,Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]]
                print "TSKIN_Mat",TSKIN_Mat[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]],"TBH_Flipud",TBH_Flipud[Station_XY_Index[Staton_Index][1],Station_XY_Index[Staton_Index][0]]
                print "\n"
    if Def_Print >= 3:
        print "numpy.where(TBH_Flipud<50.0)"
        print "CLM_Soil_Moisture,CLM_Soil_Moisture)",CLM_Soil_Moisture[:,numpy.where(TBH_Flipud<50.0)],CLM_Soil_Moisture[:,numpy.where(TBH_Flipud<50.0)]
        print "CLM_Soil_Temperature,CLM_Soil_Temperature)",CLM_Soil_Temperature[:,numpy.where(TBH_Flipud<50.0)],CLM_Soil_Temperature[:,numpy.where(TBH_Flipud<50.0)]
        print "CLM_Ground_Temperature,CLM_Ground_Temperature)",CLM_Ground_Temperature[numpy.where(TBH_Flipud<50.0)],CLM_Ground_Temperature[numpy.where(TBH_Flipud<50.0)]
        print "CLM_Vegetation_Temperature,CLM_Vegetation_Temperature)",CLM_Vegetation_Temperature[numpy.where(TBH_Flipud<50.0)],CLM_Vegetation_Temperature[numpy.where(TBH_Flipud<50.0)]
        print "CLM_Snow_Water,CLM_Snow_Water)",CLM_Snow_Water[numpy.where(TBH_Flipud<50.0)],CLM_Snow_Water[numpy.where(TBH_Flipud<50.0)]
        print "CLM_Snow_Depth,CLM_Snow_Depth)",CLM_Snow_Depth[numpy.where(TBH_Flipud<50.0)],CLM_Snow_Depth[numpy.where(TBH_Flipud<50.0)]
        print "Clay_Mat,Clay_Mat)",Clay_Mat[numpy.where(TBH_Flipud<50.0)],Clay_Mat[numpy.where(TBH_Flipud<50.0)]
        print "Sand_Mat,Sand_Mat)",Sand_Mat[numpy.where(TBH_Flipud<50.0)],Sand_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOCVL_Mat,ECOCVL_Mat)",ECOCVL_Mat[numpy.where(TBH_Flipud<50.0)],ECOCVL_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOCVH_Mat,ECOCVH_Mat)",ECOCVH_Mat[numpy.where(TBH_Flipud<50.0)],ECOCVH_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOTVL_Mat,ECOTVL_Mat)",ECOTVL_Mat[numpy.where(TBH_Flipud<50.0)],ECOTVL_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOTVH_Mat,ECOTVH_Mat)",ECOTVH_Mat[numpy.where(TBH_Flipud<50.0)],ECOTVH_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOWAT_Mat,ECOWAT_Mat)",ECOWAT_Mat[numpy.where(TBH_Flipud<50.0)],ECOWAT_Mat[numpy.where(TBH_Flipud<50.0)]
        print "ECOLAIL_Mat,ECOLAIL_Mat)",ECOLAIL_Mat[numpy.where(TBH_Flipud<50.0)],ECOLAIL_Mat[numpy.where(TBH_Flipud<50.0)]
        
        
        print "numpy.where(TBH_Flipud>350.0)"
        print "CLM_Soil_Moisture),CLM_Soil_Moisture)",CLM_Soil_Moisture[:,numpy.where(TBH_Flipud>350.0)],CLM_Soil_Moisture[:,numpy.where(TBH_Flipud>350.0)]
        print "CLM_Soil_Temperature),CLM_Soil_Temperature)",CLM_Soil_Temperature[:,numpy.where(TBH_Flipud>350.0)],CLM_Soil_Temperature[:,numpy.where(TBH_Flipud>350.0)]
        print "CLM_Ground_Temperature),CLM_Ground_Temperature)",CLM_Ground_Temperature[numpy.where(TBH_Flipud>350.0)],CLM_Ground_Temperature[numpy.where(TBH_Flipud>350.0)]
        print "CLM_Vegetation_Temperature),CLM_Vegetation_Temperature)",CLM_Vegetation_Temperature[numpy.where(TBH_Flipud>350.0)],CLM_Vegetation_Temperature[numpy.where(TBH_Flipud>350.0)]
        print "CLM_Snow_Water),CLM_Snow_Water)",CLM_Snow_Water[numpy.where(TBH_Flipud>350.0)],CLM_Snow_Water[numpy.where(TBH_Flipud>350.0)]
        print "CLM_Snow_Depth),CLM_Snow_Depth)",CLM_Snow_Depth[numpy.where(TBH_Flipud>350.0)],CLM_Snow_Depth[numpy.where(TBH_Flipud>350.0)]
        print "Clay_Mat),Clay_Mat)",Clay_Mat[numpy.where(TBH_Flipud>350.0)],Clay_Mat[numpy.where(TBH_Flipud>350.0)]
        print "Sand_Mat),Sand_Mat)",Sand_Mat[numpy.where(TBH_Flipud>350.0)],Sand_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOCVL_Mat),ECOCVL_Mat)",ECOCVL_Mat[numpy.where(TBH_Flipud>350.0)],ECOCVL_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOCVH_Mat),ECOCVH_Mat)",ECOCVH_Mat[numpy.where(TBH_Flipud>350.0)],ECOCVH_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOTVL_Mat),ECOTVL_Mat)",ECOTVL_Mat[numpy.where(TBH_Flipud>350.0)],ECOTVL_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOTVH_Mat),ECOTVH_Mat)",ECOTVH_Mat[numpy.where(TBH_Flipud>350.0)],ECOTVH_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOWAT_Mat),ECOWAT_Mat)",ECOWAT_Mat[numpy.where(TBH_Flipud>350.0)],ECOWAT_Mat[numpy.where(TBH_Flipud>350.0)]
        print "ECOLAIL_Mat),ECOLAIL_Mat)",ECOLAIL_Mat[numpy.where(TBH_Flipud>350.0)],ECOLAIL_Mat[numpy.where(TBH_Flipud>350.0)]
        
            
    os.chdir(DasPy_Path)
    
    del CLM_Soil_Moisture_Ensemble_Mat, CLM_Soil_Temperature_Ensemble_Mat, CLM_Vegetation_Temperature_Ensemble_Mat, CLM_Ground_Temperature_Ensemble_Mat
    del CLM_2m_Air_Temperature_Ensemble_Mat, CLM_Snow_Depth_Ensemble_Mat, CLM_Snow_Water_Ensemble_Mat, Parameter_Soil_Space_Ensemble, CLM_LAI, TBH, TBV, EFFECTIVE_TEMP
    del ECOCVL_Mat,ECOCVH_Mat,ECOTVL_Mat,ECOTVH_Mat,ECOWAT_Mat,DEM_Data
    del CLM_Soil_Moisture,CLM_Soil_Temperature,CLM_Ground_Temperature,CLM_Vegetation_Temperature,CLM_2m_Air_Temperature
    del CLM_Snow_Water,CLM_Snow_Depth,Sand_Mat,Clay_Mat, RSN_Mat,SD_Mat,STL_Mat,SWVL_Mat,TSKIN_Mat,Tair_Mat,ECOLAIL_Mat,Z_Mat
    del Call_CMEM
    
    gc.collect()
    del gc.garbage[:]
    
    return TBH_Flipud


def Read_CMEM_Par(Def_Region, DAS_Data_Path, DAS_Output_Path, Region_Name, Row_Numbers, Col_Numbers, \
                  Datetime_Stop, maxpft, PCT_PFT_High, PCT_PFT_Low, PCT_PFT_WATER, NC_FileName_Assimilation_2_Constant):
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    
    NC_FileName_CMEM_Par = DAS_Output_Path+"Analysis/"+Region_Name+"/CMEM_Par.nc"
    if os.path.exists(NC_FileName_CMEM_Par):
        os.remove(NC_FileName_CMEM_Par)
    
    ECOCVL_Mat = PCT_PFT_Low    # fraction of low veg
    ECOCVH_Mat = PCT_PFT_High   # fraction of high veg
    
    #No vegetation: 0
    #High vegetation: 1 Decidious forests; 2 Coniferous forests; 3 Rain forests
    #Low vegetation: 4 C3 Grasslands; 5 C4 Grasslands; 6 C3 Crops; 7 C4 Crops

    ECOTVL_Mat = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)    # low veg type
    ECOTVH_Mat = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)    # high veg type
    for Row_Index in range(Row_Numbers):
        for Col_Index in range(Col_Numbers):
            Type_Index = numpy.argmax(NC_File_Out_Assimilation_2_Constant.variables['PCT_PFT'][:, Row_Index, Col_Index])
            Type_Index_High = numpy.argmax(NC_File_Out_Assimilation_2_Constant.variables['PCT_PFT'][1:9, Row_Index, Col_Index])
            Type_Index_Low = numpy.argmax(NC_File_Out_Assimilation_2_Constant.variables['PCT_PFT'][9:maxpft, Row_Index, Col_Index])

            if (Type_Index_High == 0 or Type_Index_High == 1 or Type_Index_High == 2):
                ECOTVH_Mat[Row_Index, Col_Index] = 2
            elif (Type_Index_High == 3 or Type_Index_High == 4):
                ECOTVH_Mat[Row_Index, Col_Index] = 3
            elif (Type_Index_High == 5 or Type_Index_High == 6 or Type_Index_High == 7):
                ECOTVH_Mat[Row_Index, Col_Index] = 1
            
            if (Type_Index_Low == 0 or Type_Index_Low == 1 or Type_Index_Low == 2 or Type_Index_Low == 3 or Type_Index_Low == 4 or Type_Index_Low == 5):
                ECOTVL_Mat[Row_Index, Col_Index] = 4
            elif (Type_Index_Low == 6 or Type_Index_Low == 7):
                ECOTVL_Mat[Row_Index, Col_Index] = 5
    
    ECOWAT_Mat = PCT_PFT_WATER
    
    print 'Write NetCDF File:',NC_FileName_CMEM_Par
    NC_File_CMEM_Par = netCDF4.Dataset(NC_FileName_CMEM_Par, 'w', diskless=True, persist=True, format='NETCDF4')
    # Dim the dimensions of NetCDF
    NC_File_CMEM_Par.createDimension('lon', Col_Numbers)
    NC_File_CMEM_Par.createDimension('lat', Row_Numbers)
    
    NC_File_CMEM_Par.createVariable('ECOCVL_Mat','f4',('lat','lon',),zlib=True)
    NC_File_CMEM_Par.variables['ECOCVL_Mat'][:,:] = ECOCVL_Mat
    NC_File_CMEM_Par.createVariable('ECOCVH_Mat','f4',('lat','lon',),zlib=True)
    NC_File_CMEM_Par.variables['ECOCVH_Mat'][:,:] = ECOCVH_Mat
    NC_File_CMEM_Par.createVariable('ECOTVL_Mat','f4',('lat','lon',),zlib=True)
    NC_File_CMEM_Par.variables['ECOTVL_Mat'][:,:] = ECOTVL_Mat
    NC_File_CMEM_Par.createVariable('ECOTVH_Mat','f4',('lat','lon',),zlib=True)
    NC_File_CMEM_Par.variables['ECOTVH_Mat'][:,:] = ECOTVH_Mat
    NC_File_CMEM_Par.createVariable('ECOWAT_Mat','f4',('lat','lon',),zlib=True)
    NC_File_CMEM_Par.variables['ECOWAT_Mat'][:,:] = ECOWAT_Mat
    
    NC_File_CMEM_Par.close()
    
    NC_File_Out_Assimilation_2_Constant.close()
    
    del ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat
    
    gc.collect()
    del gc.garbage[:]
    
    return 0


def Run_COSMOS(Ens_Index, DAS_Depends_Path, Def_Region, Def_PP, Def_Print, N0, nlyr, Soil_Layer_Num, Grid_Resolution_CEA, DasPy_Path, \
               Row_Numbers, Col_Numbers, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, omp_get_num_procs_ParFor, *vartuple):
    
    COSMOS = imp.load_source("COSMOS",DasPy_Path+"ObsModel/COSMOS/COSMOS.py")
    COSMIC_Py = imp.load_source("COSMIC_Py",DasPy_Path+"ObsModel/COSMOS/COSMIC_Py.py")
    window = imp.load_source("window",DasPy_Path+"Algorithm/window.py")
    COSMIC = imp.load_dynamic("COSMIC",DasPy_Path+"ObsModel/COSMOS/COSMIC.so")
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    CLM_Soil_Layer_Thickness_Cumsum = NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness_Cumsum'][:,:,:]
    DEM_Data = NC_File_Out_Assimilation_2_Constant.variables['DEM_Data'][:,:]
    Land_Mask_Data = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][:,:]
    Bulk_Density_Top_Region = NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Top_Region'][:,:]
    Bulk_Density_Sub_Region = NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Sub_Region'][:,:]
    NC_File_Out_Assimilation_2_Constant.close()
    
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')          
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    CLM_Soil_Moisture = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,Ens_Index]
    CLM_Air_Pressure = NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Air_Pressure_Ensemble_Mat'][:,:,Ens_Index]
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Diagnostic.close()
    
    if Def_Region == 3 or Def_Region == 8: # for rur, we use the field parameters
        NC_FileName_Observation = "/lustre/jwork/jicg41/jicg4128/DAS_Data/Analysis/Rur/Observation_0.nc"
        NC_File_Observation = netCDF4.Dataset(NC_FileName_Observation, 'r')
        bd = NC_File_Observation.variables['Observation_Misc'][0,:,:]
        lw = NC_File_Observation.variables['Observation_Misc'][1,:,:]
        Ncosmic = NC_File_Observation.variables['Observation_Misc'][2,:,:]
        NC_File_Observation.close()
    
    Air_Pressure = CLM_Air_Pressure
    Observation_Frequency = 1 # Hours
    
    #Sensor_Layer, Sensor_Depth, Sensor_Soil_Moisture, Neutron_COSMOS = COSMOS(N0, nlyr, Def_Print, Observation_Frequency, DEM_Data, CLM_Soil_Moisture_Smoothed, CLM_Soil_Layer_Thickness, Air_Pressure, Bulk_Density_Top_Region, Bulk_Density_Sub_Region, Row_Numbers, Col_Numbers, Soil_Layer_Num)
    
    Mask = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.bool)
    Mask[::] = True
    numexpr_a = CLM_Soil_Moisture[0,:,:]
    numexpr_b = -9999.0
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    Mask[numexpr_c] = False
    Mask_Col = Mask.flatten()
    Mask_Index = numpy.where(Mask_Col == True)
    Mask_Size = numpy.size(Mask_Index)
    
    CLM_Soil_Layer_Thickness_Cumsum = CLM_Soil_Layer_Thickness_Cumsum * 100 # from meters to cm
    
    soil_moisture = numpy.zeros((Mask_Size,Soil_Layer_Num-5),dtype=numpy.float32)
    layerz = numpy.zeros((Mask_Size,Soil_Layer_Num-5),dtype=numpy.float32)
    
    for Soil_Layer_Index in range(Soil_Layer_Num-5):
        soil_moisture[:,Soil_Layer_Index] = CLM_Soil_Moisture[Soil_Layer_Index,:,:].flatten()[Mask_Index]
        layerz[:,Soil_Layer_Index] = CLM_Soil_Layer_Thickness_Cumsum[Soil_Layer_Index,:,:].flatten()[Mask_Index]
    
    if Def_Region == 3 or Def_Region == 8:
        bd = bd.flatten()[Mask_Index]
        n = Ncosmic.flatten()[Mask_Index]
        lattwat = lw.flatten()[Mask_Index]
    else:
        bd = Bulk_Density_Top_Region.flatten()[Mask_Index]
        n = numpy.zeros(Mask_Size,dtype=numpy.float32)
        lattwat = numpy.zeros(Mask_Size,dtype=numpy.float32)
        
    n[:] = N0
    lattwat[:] = 0.03
    nthreads = omp_get_num_procs_ParFor
    
    if Def_Print:
        print "numpy.min(n),numpy.max(n),numpy.min(nlyr),numpy.max(nlyr),numpy.min(soil_moisture),numpy.max(soil_moisture)"
        print numpy.min(n),numpy.max(n),numpy.min(nlyr),numpy.max(nlyr),numpy.min(soil_moisture),numpy.max(soil_moisture)
        print "numpy.min(layerz),numpy.max(layerz),numpy.min(bd),numpy.max(bd),numpy.min(lattwat),numpy.max(lattwat),nthreads"
        print numpy.min(layerz),numpy.max(layerz),numpy.min(bd),numpy.max(bd),numpy.min(lattwat),numpy.max(lattwat),nthreads
        
    Neutron_COSMOS, Sensor_Soil_Moisture, Sensor_Depth = COSMIC_Py.COSMIC_Py(COSMIC,n,nlyr,soil_moisture,layerz,bd,lattwat,nthreads)
    
    if Def_Print:
        print numpy.min(Neutron_COSMOS),numpy.max(Neutron_COSMOS),numpy.min(Sensor_Soil_Moisture),numpy.max(Sensor_Soil_Moisture)
        
    Neutron_Matirx = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    Neutron_Matirx_Col = Neutron_Matirx.flatten()
    Neutron_Matirx_Col[Mask_Index] = Neutron_COSMOS
    Neutron_Matirx = numpy.reshape(Neutron_Matirx_Col,(Row_Numbers, Col_Numbers))
    
    Neutron_Matirx_Smoothed = numpy.zeros_like(Neutron_Matirx)
    
    Data = numpy.copy(Neutron_Matirx[:, :])
    invalid = numpy.zeros_like(Data,dtype=numpy.bool)
    invalid[::] = False
    numexpr_a = Data
    numexpr_b = -9999.0
    numexpr_c = numpy.where(numexpr_a == numexpr_b)
    invalid[numexpr_c] = True
    ind = scipy.ndimage.morphology.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    Data = Data[tuple(ind)]
    
    Neutron_Matirx_Smoothed[:, :] = Neutron_Matirx[:, :]
    
    del CLM_Soil_Moisture, CLM_Air_Pressure, CLM_Soil_Layer_Thickness_Cumsum
    del soil_moisture, layerz, bd, n, lattwat
    del DEM_Data,Land_Mask_Data,Bulk_Density_Top_Region,Bulk_Density_Sub_Region
    del Data,invalid,ind,Air_Pressure,Mask,Mask_Col,Mask_Index,Neutron_COSMOS, Sensor_Soil_Moisture, Sensor_Depth
    del COSMIC_Py,window,COSMIC
    numexpr_a = []
    numexpr_b = []
    numexpr_c = []
    
    gc.collect()
    del gc.garbage[:]
    
    return Neutron_Matirx_Smoothed


def Run_TSM(Ens_Index, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor, NAvalue, Variable_Assimilation_Flag, Variable_List,\
            DasPy_Path, Def_Print, SensorVariable, Observation_View_Zenith_Angle, NC_FileName_Assimilation_2_Constant, \
            NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias):
    
    Clumping_Index = imp.load_source("Clumping_Index",DasPy_Path+"ObsModel/LST/Clumping_Index.py")
    
    NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    CLM_Vegetation_Temperature = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_Ground_Temperature = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,Ens_Index]
    CLM_LAI = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][Ens_Index,Variable_List.index('LAI'),:,:]
    #CLM_LAI = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:,Variable_List.index('LAI'),:,:],axis=0)
    Teta_Residual = NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][:,:,:]
    Teta_Saturated = NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][:,:,:]
    Teta_Field_Capacity = NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][:,:,:]
    Teta_Wilting_Point = NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][:,:,:]
    Land_Mask_Data = NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][::]
    PFT_Dominant_Index = NC_File_Out_Assimilation_2_Constant.variables['PFT_Dominant_Index'][::]
    #Canopy_Height_Width_Ratio = NC_File_Out_Assimilation_2_Constant.variables['Canopy_Height_Width_Ratio'][::]
    NC_File_Out_Assimilation_2_Initial.close()
    NC_File_Out_Assimilation_2_Constant.close()
    
    if Def_Print:
        print "----------------------------- Generate Fused Surface Temperature"
    
    #if Variable_Assimilation_Flag[Variable_List.index("Sensible_Heat")] == 1:
    #    SensorVariable_Index = SensorVariable.index('Sensible_Heat')
    if Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")] == 1:
        SensorVariable_Index = SensorVariable.index('Surface_Temperature')
    
    # Estimate of bulk surface radiometric temperature based on Kustas, W. and M. Anderson (2009)
    #Fractional_Vegetation_Cover = 1.0 - numpy.exp(-0.5*MONTHLY_LAI[string.atoi(Stop_Month)-1,:,:]/numpy.cos(Observation_View_Zenith_Angle / 180.0 * numpy.pi))
    if Def_Print:
        print "SensorVariable_Index",SensorVariable_Index
    
    Azimuth_Angle = numpy.zeros_like(PFT_Dominant_Index,dtype=numpy.float32)
    Azimuth_Angle[:,:] = 45.0
    Canopy_Height_Width_Ratio = numpy.ones_like(PFT_Dominant_Index,dtype=numpy.float32)
    Omega = Clumping_Index.Clumping_Index(CLM_LAI, Observation_View_Zenith_Angle[SensorVariable_Index,:,:],Azimuth_Angle,Canopy_Height_Width_Ratio)
    
    # Correct the Forest Omega
    Omega[numpy.where(PFT_Dominant_Index <= 8)] = 1.0
    
    Fractional_Vegetation_Cover = 1.0 - numpy.exp(-0.5*Omega*CLM_LAI/numpy.cos(Observation_View_Zenith_Angle[SensorVariable_Index,:,:] / 180.0 * numpy.pi))
    #print numpy.max(MONTHLY_LAI[string.atoi(Stop_Month)-1,:,:]),numpy.min(MONTHLY_LAI[string.atoi(Stop_Month)-1,:,:])
    if Def_Print:
        print "Observation_View_Zenith_Angle",numpy.max(Observation_View_Zenith_Angle[SensorVariable_Index,:,:]),numpy.min(Observation_View_Zenith_Angle[SensorVariable_Index,:,:])
        print "Omega",numpy.max(Omega),numpy.min(Omega)
        print "CLM_LAI",numpy.max(CLM_LAI[numpy.where(Land_Mask_Data!=NAvalue)]),numpy.min(CLM_LAI[numpy.where(Land_Mask_Data!=NAvalue)])
        print "Fractional_Vegetation_Cover",numpy.max(Fractional_Vegetation_Cover[numpy.where(Land_Mask_Data!=NAvalue)]),numpy.min(Fractional_Vegetation_Cover[numpy.where(Land_Mask_Data!=NAvalue)])
        if Def_Print >= 2:
            print "Observation_View_Zenith_Angle",Observation_View_Zenith_Angle,numpy.cos(Observation_View_Zenith_Angle / 180.0 * numpy.pi)
            print "Fractional_Vegetation_Cover",Fractional_Vegetation_Cover
                
    TSM_Temperature = numpy.power(Fractional_Vegetation_Cover * numpy.power(CLM_Vegetation_Temperature, 4) + (1.0 - Fractional_Vegetation_Cover) * numpy.power(CLM_Ground_Temperature, 4), 0.25)
    
    NA_Flag = False
    #print CLM_Ground_Temperature
    Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,TSM_Temperature,NA_Flag, 'Soil_Temperature', Variable_Assimilation_Flag, Variable_List, 
                    Teta_Residual[0,:,:], Teta_Saturated[0,:,:], Teta_Field_Capacity[0,:,:], Teta_Wilting_Point[0,:,:], NAvalue)
    
    if Def_Print:
        print "TSM_Temperature",numpy.max(TSM_Temperature[numpy.where(Land_Mask_Data!=NAvalue)]),numpy.min(TSM_Temperature[numpy.where(Land_Mask_Data!=NAvalue)])
    
    del Canopy_Height_Width_Ratio, PFT_Dominant_Index, Azimuth_Angle
    del CLM_Vegetation_Temperature, CLM_Ground_Temperature, CLM_LAI, Omega, Fractional_Vegetation_Cover, Land_Mask_Data, Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point
    
    gc.collect()
    del gc.garbage[:]
    
    return TSM_Temperature


def Run_CLM(Model_Driver, Def_PP, Do_DA_Flag, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized,   Def_Region, Def_Initial, \
            Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, 
            Run_Dir, Run_Dir_Array, Model_Path, CLM_Flag, domain_name, rdirc_name, fatmgrid_name, \
            fatmlndfrc_name, fndepdat_name, fsurdat_name, fsnowoptics_name, fsnowaging_name, fglcmask_name, finidat_name, flndtopo_name, \
            fpftcon_name, popd_streams_name, light_streams_name, N_Steps, Row_Numbers_String, Col_Numbers_String, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, \
            Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, DasPy_Path, DAS_Data_Path, Forcing_File_Path, Forcing_File_Path_Array, 
            dtime, ntasks_CLM, rootpe_CLM, nthreads_CLM, Ensemble_Number, num_processors,
            COUP_OAS_PFL, CESM_Init_Flag, fcomm, fcomm_null, fcomm_rank):
    
    if CLM_Flag == "3D":
        
        #=============================== datm_atm_in parameters
        domain_file_path = DAS_Data_Path + "SysModel/CLM/tools/"
        
        domain_file_lnd_path = DAS_Data_Path + "SysModel/CLM/tools/"
        
        align_year = Start_Year
        first_year = Start_Year
        last_year = Stop_Year
        
        #=============================== drv_in parameters
        case_name = repr(Region_Name)
        hostname = repr(socket.gethostname())
        orb_iyear_ad = Start_Year
        start_type = "'startup'"
        # startup, continue, branch, hybrid
        # startup - arbitrary initialization determined by components (default)
        # hybrid - initialization occurs from the restart/initial files of a previous reference case, the start date can be changed with respect to reference case
        # branch - initialization occurs from restart files of a previous reference case, cannot change start date with respect to reference case
        
        username = repr(getpass.getuser())
        if Def_SpinUp:
            atm_cpl_dt = dtime
            lnd_cpl_dt = dtime
            ocn_cpl_dt = 86400
            ice_cpl_dt = dtime
            glc_cpl_dt = 86400
            rof_cpl_dt = 10800
            wav_cpl_dt = dtime
            
        else:
            atm_cpl_dt = dtime
            lnd_cpl_dt = dtime
            ocn_cpl_dt = 86400
            ice_cpl_dt = dtime
            glc_cpl_dt = 86400
            rof_cpl_dt = 10800
            wav_cpl_dt = dtime
        end_restart = '.true.'
        restart_option = "'end'"
        start_tod = str((Datetime_Start - Datetime_Start_Init).seconds)           # Start time of day of run (seconds).
        start_ymd = Start_Year + Start_Month + Start_Day    # Start date of run (yyyymmddformat)
        stop_ymd = Stop_Year + Stop_Month + Stop_Day     # Stop date(YYYYMMDD)
        if (Datetime_Stop - Datetime_Stop_Init).days > 0:
            stop_tod = str((Datetime_Stop - Datetime_Stop_Init).days * 86400)            # Stop time of day (seconds)
        elif (Datetime_Stop - Datetime_Stop_Init).seconds > 0:
            stop_tod = str((Datetime_Stop - Datetime_Stop_Init).seconds)            # Stop time of day (seconds)
        else:
            stop_tod = '0'
        
        stop_tod_string = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
        
        #=============================== lnd_in parameters
        
        rtm_nsteps = 1          # Frequency in number of time steps at which RTM is called.
        aero_file_path = DAS_Data_Path + "SysModel/CLM/inputdata/atm/cam/chem/trop_mozart_aero/aero/"
        aero_file_name = "aerosoldep_monthly_2000_mean_1.9x2.5_c090421.nc"
        megan_factors_file_path = DAS_Data_Path + "SysModel/CLM/inputdata/atm/cam/chem/trop_mozart/emis/"
        megan_factors_file_name = "megan21_emis_factors_c20120313.nc"
        frivinp_rtm_path = DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"
        frivinp_rtm_name = "rdirc_0.5x0.5_simyr2000_slpmxvl_c120717.nc"
        seq_maps_file_name = "seq_maps.rc"
        
        wrtdia = ".false."      # TRUE if want diagnostic of global radiative temperature written to CLM log file 
        
        # Change the Output Samples in each Histroy File
        if Def_Par_Sensitivity or Def_Par_Correlation:
            hist_mfilt = N_Steps  # The primary history file will have hist_mfilt time samples on every tape
        else:
            hist_mfilt = 1  # The primary history file will have hist_mfilt time samples on every tape
        
        hist_nhtfrq = 1 # Per file history write frequency in steps if positive, (0=monthly), negative means hours Default: 0,-24,-24,-24,-24,-24
        
        #hist_nhtfrq = 0
        hist_crtinic = "'NONE'"   # Frequency with which initial datasets will be generated. 
                                    # Valid values are 'YEARLY','MONTHLY','DAILY','6-HOURLY' or 'NONE'.
        
        hist_dov2xy = '.true.' # Per tape spatial averaging flag. 
                                # If set to true, produces grid-average history fields on output tape. If set to false,one-dimensional fields are produced.
        #if Def_Region == -1:    # For Drip Irrigation Area, We Need the Column Based Output
        #    hist_dov2xy = '.false.'
        hist_ndens = 2          # Per file history file density (i.e. output precision) (1=double precision, 2=single precision)
        hist_type1d_pertape = "'','GRID'"  # Per tape one dimensional output type. Only used if one dimensional output is selected for the given tape ( via the setting of HIST DOV2XY). Valid values are
                                        # 'GRID','LAND','COLS','PFTS'.
        hist_empty_htapes = '.true.'    # If set to true, all the history tapes are empty by default. Only variables explicitly listed by the user will be output.
        # List of fields to include on the respective history tape. See tables 9-18 for the list of default fields on the primary history tape. Namelist specfication can take one of two forms.
        # The user may simply specify the name of the field to be included on the history tape (in which case the default time averaging for that field will be used).
        # For example,HIST FINCL2='TV',will add the field TV to the second history tape with whatever default time averaging was specified for TV.
        # Alternatively, the user may specify the field name, followed by a ":" followed by the time averaging flag desired (valid flags are 'I' for instantaneous, 'A' for average, 'M'
        # for minimum, and 'X' for maximum). For example, HIST FINCL2 = 'TV:I" will add the field TV with instantaneous output to the second history tape
        
        hist_avgflag_pertape = "'I'"  # A Average, over the output interval.
                                    #    I Instantaneous, output the value at the output interval.
                                    #    X Maximum, over the output interval.
                                    #    M Minimum, over the output interval.
        
        hist_fincl1 = "'TSA','PSurf','ONSET_FDD','RAIN','SNOW','QDRIP','H2OSFC','FH2OSFC','Rnet','FSA','PSurf','QSOIL','QVEGE','QVEGT','QIRRIG','HTOP','HBOT','SOILC','NPP','NEE','GPP','TV','TG','TG_U','TG_R',\
        'TSOI','TLAKE','FSH_V','FSH_G','Qh','Qle','TLAI','TSAI','H2OSNO','H2OCAN','H2OSOI', 'SOILICE','SNOWDP','FSNO','INT_SNOW','QINFL','FGR','FGEV','ZWT','WA','TWS','TSA','TBOT','QOVER','QRGWL','QSNWCPLIQ','QDRAI','QH2OSFC'"
        
        # List of fields to exclude from the respective history tape. The field name must appear in the Master Field List of the default history tape (the primary tape)
        hist_fexcl1 = "''"
        
        ################################# Call the CLM Model #################################
        Call_CLM_3D(Def_First_Run, Def_CESM_Multi_Instance, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir, Run_Dir_Array, Ensemble_Number, num_processors, DasPy_Path, Model_Path, Model_Driver, Def_SpinUp, Def_PP, Def_Print, align_year, first_year, last_year, \
                    domain_file_path, Forcing_File_Path, Forcing_File_Path_Array, domain_file_lnd_path, domain_name, rdirc_name, aero_file_path, aero_file_name, megan_factors_file_path, megan_factors_file_name,frivinp_rtm_path,frivinp_rtm_name, \
                    case_name, hostname, orb_iyear_ad, start_type, username, atm_cpl_dt, lnd_cpl_dt, ocn_cpl_dt, ice_cpl_dt, glc_cpl_dt, rof_cpl_dt, wav_cpl_dt, end_restart, restart_option, start_tod, start_ymd, stop_tod, stop_ymd, ntasks_CLM, rootpe_CLM, nthreads_CLM, \
                    dtime, rtm_nsteps, fatmgrid_name, fatmlndfrc_name, fglcmask_name, finidat_name, flndtopo_name, fndepdat_name, fpftcon_name, fsnowaging_name, fsnowoptics_name, fsurdat_name, popd_streams_name, light_streams_name, \
                    wrtdia, hist_nhtfrq, hist_mfilt, hist_crtinic, hist_dov2xy, hist_ndens, hist_type1d_pertape, hist_empty_htapes, hist_avgflag_pertape, hist_fincl1, hist_fexcl1,\
                    Region_Name, Stop_Year, Stop_Month, Stop_Day, stop_tod_string, seq_maps_file_name, Row_Numbers_String, Col_Numbers_String, DAS_Data_Path,
                    COUP_OAS_PFL, CESM_Init_Flag, fcomm, fcomm_null, fcomm_rank)
    
    gc.collect()
    del gc.garbage[:]
    
    return 0




def Check_Outliers(DasPy_Path, Land_Mask_Data, Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,Analysis_Grid,NA_Flag, SensorVariable,                         
                        Variable_Assimilation_Flag, Variable_List, Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point, NAvalue):
    
#    numpy.savez("Check_Outliers.npz",Def_ParFor,DAS_Depends_Path,omp_get_num_procs_ParFor,Analysis_Grid,NA_Flag, SensorVariable, Soil_Moisture_DA_Flag, Surface_Temperature_DA_Flag, Vegetation_Temperature_DA_Flag, Albedo_DA_Flag, Snow_Depth_DA_Flag, Snow_Cover_Fraction_DA_Flag,Canopy_Water_DA_Flag,
#                        Snow_Water_Equivalent_DA_Flag, Crop_Planting_DA_Flag, Crop_Harvest_DA_Flag, LAI_DA_Flag, Latent_Heat_DA_Flag, Latent_Heat_Daily_DA_Flag, Sensible_Heat_DA_Flag, Water_Storage_DA_Flag, Water_Table_Depth_DA_Flag, Emissivity_DA_Flag,
#                        Teta_Residual, Teta_Saturated, Teta_Field_Capacity, Teta_Wilting_Point, NAvalue)
    
    ParFor = imp.load_source("ParFor",DasPy_Path+"ParFor.py")
    
    '''
        NA_Flag means whether we assign the unreasonable grid as NAvalue or real value
    '''
    Soil_Moisture_Minimum = numpy.zeros_like(Analysis_Grid)
    Soil_Moisture_Maximum = numpy.zeros_like(Analysis_Grid)
    Soil_Moisture_Minimum[::] = Teta_Residual
    Soil_Moisture_Maximum[::] = Teta_Saturated
    #Soil_Moisture_Maximum[::] = (Teta_Field_Capacity + 0.7*(Teta_Saturated-Teta_Field_Capacity)) # sat will results in high ground water recharge
    
    Ground_Temperature_Minimum = numpy.zeros_like(Analysis_Grid)
    Ground_Temperature_Maximum = numpy.zeros_like(Analysis_Grid)
    Ground_Temperature_Minimum[::] = 273.15 - 60.0
    Ground_Temperature_Maximum[::] = 273.15 + 60.0
    Vegetation_Temperature_Minimum = numpy.zeros_like(Analysis_Grid)
    Vegetation_Temperature_Maximum = numpy.zeros_like(Analysis_Grid)
    Vegetation_Temperature_Minimum[::] = 273.15 - 60.0
    Vegetation_Temperature_Maximum[::] = 273.15 + 60.0
    
    
    
    Row_Numbers,Col_Numbers = numpy.shape(Analysis_Grid)
    # Check the Analysis Boundary
    if not NA_Flag:
        if SensorVariable == "Soil_Moisture":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Soil_Moisture_Minimum,Soil_Moisture_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Soil_Moisture_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = Soil_Moisture_Minimum[Row_Index, Col_Index]
                            elif Analysis_Grid[Row_Index, Col_Index] > Soil_Moisture_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = Soil_Moisture_Maximum[Row_Index, Col_Index]
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            
        elif SensorVariable == "Surface_Temperature" or SensorVariable == "Ground_Temperature" or SensorVariable == "Soil_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Ground_Temperature_Minimum,Ground_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Ground_Temperature_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = Ground_Temperature_Minimum[Row_Index, Col_Index]
                            elif Analysis_Grid[Row_Index, Col_Index] > Ground_Temperature_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = Ground_Temperature_Maximum[Row_Index, Col_Index]
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            
        elif SensorVariable == "Vegetation_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Vegetation_Temperature_Minimum,Vegetation_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index,Col_Index] < Vegetation_Temperature_Minimum[Row_Index,Col_Index]:
                                Analysis_Grid[Row_Index,Col_Index] = Vegetation_Temperature_Minimum[Row_Index,Col_Index]
                            elif Analysis_Grid[Row_Index,Col_Index] >  Vegetation_Temperature_Maximum[Row_Index,Col_Index]:
                                Analysis_Grid[Row_Index,Col_Index] = Vegetation_Temperature_Maximum[Row_Index,Col_Index]
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
        
        
        elif SensorVariable == "2m_Air_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers_NA(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Ground_Temperature_Minimum,Ground_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Ground_Temperature_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            elif Analysis_Grid[Row_Index, Col_Index] > Ground_Temperature_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
    else:
        
        if SensorVariable == "Soil_Moisture":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers_NA(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Soil_Moisture_Minimum,Soil_Moisture_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Soil_Moisture_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            elif Analysis_Grid[Row_Index, Col_Index] > Soil_Moisture_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
                
        elif SensorVariable == "Surface_Temperature" or SensorVariable == "Ground_Temperature" or SensorVariable == "Soil_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers_NA(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Ground_Temperature_Minimum,Ground_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Ground_Temperature_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            elif Analysis_Grid[Row_Index, Col_Index] > Ground_Temperature_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
        elif SensorVariable == "Vegetation_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers_NA(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Vegetation_Temperature_Minimum,Vegetation_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index,Col_Index] < Vegetation_Temperature_Minimum[Row_Index,Col_Index]:
                                Analysis_Grid[Row_Index,Col_Index] = NAvalue
                            elif Analysis_Grid[Row_Index,Col_Index] >  Vegetation_Temperature_Maximum[Row_Index,Col_Index]:
                                Analysis_Grid[Row_Index,Col_Index] = NAvalue
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
        
        
        elif SensorVariable == "2m_Air_Temperature":
            if Def_ParFor:
                ParFor.ParFor_Check_Outliers_NA(Land_Mask_Data, Analysis_Grid,Row_Numbers,Col_Numbers,Ground_Temperature_Minimum,Ground_Temperature_Maximum,DAS_Depends_Path,omp_get_num_procs_ParFor,NAvalue)
            else:
                for Row_Index in range(Row_Numbers):
                    for Col_Index in range(Col_Numbers):
                        if Land_Mask_Data[Row_Index, Col_Index] != NAvalue:
                            if Analysis_Grid[Row_Index, Col_Index] < Ground_Temperature_Minimum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                            elif Analysis_Grid[Row_Index, Col_Index] > Ground_Temperature_Maximum[Row_Index, Col_Index]:
                                Analysis_Grid[Row_Index, Col_Index] = NAvalue
                        else:
                            Analysis_Grid[Row_Index, Col_Index] = NAvalue
    
    del Soil_Moisture_Minimum,Soil_Moisture_Maximum,Ground_Temperature_Minimum,Ground_Temperature_Maximum,Vegetation_Temperature_Minimum,Vegetation_Temperature_Maximum
    
    
    gc.collect()
    del gc.garbage[:]
    
    return 0
