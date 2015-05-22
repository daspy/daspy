'''
This Program Use the Remote Sensed 'Soil Moisture', 'Canopy Water', 'Ground Surface Temperature', 'Albedo', 'Snow Depth', 'Snow Cover Fraction', 'Leaf Area Index', 'Water Storage Variations', 'Ground Table Depth', 'Run Off' data to Update the Community Land Model.
Excution Flows:
1. Define the Time; 2. Read the Observation Time Records File; 3. Determine when to Do Assimilation on which Variable;
4. Run CLM to the Observation Time; 5. Read the CLM Initial File and Observation File;
6. Analyze the CLM Output and Observation Data; 7. Call LETKF to Geft the Analysis;
8. Update CLM Initial File; 9. Go To 1.
m
The Observation Data are Contained in the 'DAS_Data/Observation' Directory. The 'Observation_Time.txt' File Contains the Observation Information Line by Line.

The multi-scale observation will be processed by the discrete wavelets transformation (DWT).

The Observation variance is estimated by the geostatistics or image noise method.

We only use one ensemble during the simulation, after that we perturb the model state to get the ensembles for assimilation.

1=vis, 2=nir

'''
from mpi4py import MPI
import multiprocessing, shutil
from DAS_Assim_Common import *
from DAS_Assim import *
from DAS_Misc import *
from DAS_Driver_Common import *
from DAS_Utilities import *

# Data Assimilation, Parameter Estimation, Bias Estimation

def DAS_Driver(mpi4py_comm, mpi4py_null, mpi4py_rank,  mpi4py_size, mpi4py_name, Model_Driver,Do_DA_Flag, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Texture_Layer_Opt_Num, Observation_Box, LAI_Year_String, MODIS_LAI_Data_ID,\
             Num_of_Days_Monthly, Start_Year, Start_Month, Start_Day, Start_Hour, Start_Minute, End_Year, End_Month, End_Day, End_Hour, End_Minute, Datetime_Start, Datetime_Start_Init, \
             Datetime_End, Datetime_End_Init, Datetime_Initial, UTC_Zone, CLM_NA, NAvalue, Assim_Algorithm_Name, Station_XY, Station_XY_Index, dtime,\
             NSLOTS, Feedback_Assim, Parameter_Optimization, Parameter_Regularization, Def_CDF_Matching, Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, Post_Inflation_Alpha, Def_Snow_Effects, N0, nlyr,\
             Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,
             Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
             Observation_Time_File_Path, Def_CESM_Multi_Instance, Constant_File_Name_Header, finidat_initial_CLM, finidat_initial_PFCLM, Def_PP, DAS_Fortran_Lib, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, \
             Def_ParFor, Def_Region, Def_Initial, Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, CLM_Flag,  Def_ReBEL, Def_Localization, \
             Num_Local_Obs_State, Num_Local_Obs_Par,  Num_Local_Obs_Bias, eps, msw_infl, Def_Multiresolution, Def_Write_Initial, Ensemble_Number, Ensemble_Number_Predict,  Call_Gstat_Flag, Write_DA_File_Flag, Use_Mask_Flag, Def_Figure_Output,\
             Forcing_File_Path_Home, Soil_Layer_Num, Snow_Layer_Num, ParFlow_Layer_Num, maxpft, numrad, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water, Plot_Analysis, Def_Debug, Initial_Perturbation, \
             Weather_Forecast_Days, PicHeight, PicWidth, RegionName, Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, Grid_Resolution_CEA_String, xllcenter, yllcenter, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution, Proj_String, \
             Grid_Resolution_CEA, Grid_Resolution_GEO, mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, \
             Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Region_Name, Run_Dir_Home, Model_Path, Hydraulic_File_Name, Mask_File, Observation_Path, DAS_Data_Path, DasPy_Path, DAS_Output_Path, DAS_Depends_Path, 
             octave, r, plt, cm, colors, inset_axes, fm, legend):
        
    #print DasPy_Path+"ObsModel/COSMOS/COSMIC_Py.py"
    #print DasPy_Path+"ObsModel/COSMOS/COSMIC_Py.py"
    COSMIC_Py = imp.load_source("COSMIC_Py",DasPy_Path+"ObsModel/COSMOS/COSMIC_Py.py")
    memory_profiler = imp.load_source("memory_profiler",DasPy_Path+"Utilities/memory_profiler.py") 
    COSMIC = imp.load_dynamic("COSMIC",DasPy_Path+"ObsModel/COSMOS/COSMIC.so")
    
    num_processors = multiprocessing.cpu_count()
    
    start = time.time()
    
    UTM_Zone = int(round(mksrf_edgee/15.0))    # The Difference Between the UTM Time and the Local Time
    
    diskless_flag = True
    persist_flag = True
    
    ############################################ PFCLM ####################################################
    COUP_OAS_PFL = False        # whether to run coupled ParFlow or not
    if Model_Driver == "PFCLM":
        COUP_OAS_PFL = True
    CESM_Init_Flag = 1          # whether is the first time to run ccsm_driver
    ############################################# PFCLM ###################################################
    
    gelmna_threshold = 1.0  # Threshold to Stop Parameter Optimization (No Stop)
    #gelmna_threshold = 1.08  # Threshold to Stop Parameter Optimization
    
    if mpi4py_rank == 0:
        if Def_Region == -1:
            forcing_file_name = Forcing_File_Path_Home +"_Ens1/"+  Start_Year + '-' + Start_Month + '-' + Start_Day + '.nc'
            if not os.path.exists(forcing_file_name):
                print "Forcing file",forcing_file_name,"does not exists!!!!!!!!"
                print "Please Change the Start Date and Time."
                os.abort()
        else:
            forcing_file_name = Forcing_File_Path_Home +"/"+  Start_Year + '-' + Start_Month + '-' + Start_Day + '.nc'
            if not os.path.exists(forcing_file_name):
                print "Forcing file",forcing_file_name,"does not exists!!!!!!!!"
                print "Please Change the Start Date and Time."
                if Def_SpinUp != 1: # If we do multi years spinup, we could use one year foring to simulate multi years
                    os.abort()
        
        PP_Port = 23335 + int(numpy.random.uniform(1000*Def_Region,2000*Def_Region))
        active_nodes_server = []
    else:
        forcing_file_name = None
        PP_Port = None
        active_nodes_server = None
        
    
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()    
        PP_Port = mpi4py_comm.bcast(PP_Port)
        active_nodes_server = mpi4py_comm.bcast(active_nodes_server)
        forcing_file_name = mpi4py_comm.bcast(forcing_file_name)
    
    if mpi4py_rank == 0:
        print "mpi4py_rank",mpi4py_rank,"PP_Port",PP_Port
    
    job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
    if mpi4py_rank == 0:
        restart_pp_server = (len(job_server_node_array) < 1)
        print "-------------- restart_pp_server",restart_pp_server
    else:
        restart_pp_server = None
    
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
        restart_pp_server = mpi4py_comm.bcast(restart_pp_server)
    
    while Def_PP and restart_pp_server and Ensemble_Number > 1:
        job_server_node_array = Stop_ppserver(Def_PP, DAS_Depends_Path, job_server_node_array, NSLOTS, DasPy_Path, active_nodes_server, PP_Servers_Per_Node)
        job_server_node_array, active_nodes_server, PROCS_PER_NODE, PP_Port, PP_Servers_Per_Node = Start_ppserver(mpi4py_comm, mpi4py_rank, mpi4py_name, DAS_Output_Path, Ensemble_Number, DAS_Depends_Path, active_nodes_server, Def_Region, NSLOTS, Def_Print, DasPy_Path, Def_PP, Def_CESM_Multi_Instance, PP_Port)
        
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
    
    if mpi4py_rank == 0:
        print "job_server_node_array",job_server_node_array
        print "active_nodes_server",active_nodes_server
    
    if mpi4py_rank == 0:
        print "mpi4py_comm,mpi4py_rank,mpi4py_size,mpi4py_name",mpi4py_comm,mpi4py_rank,mpi4py_size,mpi4py_name
        
    # MPI Split into Ensemble_Number Groups
    mpi4py_comm_split = []
    mpipy_comm_decomposition = 1
    if Def_PP == 2:
        if mpi4py_rank == 0:
            mpipy_comm_decomposition = mpi4py_size/Ensemble_Number
        else:
            mpipy_comm_decomposition = None
        mpipy_comm_decomposition = mpi4py_comm.bcast(mpipy_comm_decomposition)
        
        if Ensemble_Number > 1:
            color = mpi4py_rank/mpipy_comm_decomposition
            key = mpi4py_rank
        else:
            color = 0
            key = 0
            
        mpi4py_comm_split = mpi4py_comm.Split(color, key)
        
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
            
    if mpi4py_rank == 0:
        print "mpipy_comm_decomposition,mpi4py_comm_split",mpipy_comm_decomposition,mpi4py_comm_split
        if Def_PP == 2:
            print "mpi4py_comm_split.py2f()",mpi4py_comm_split.py2f(),"mpi4py_comm_split.Get_size()",mpi4py_comm_split.Get_size()
    
    if Def_PP == 2:
        ntasks_CLM[:] = mpi4py_comm_split.Get_size()
    
    ###################################################################
    if Def_Region >= 9:
        Grid_Resolution_GEO_Global = 360.0/43200.0
        Resolution_Name = "1KM"
    else:
        Grid_Resolution_GEO_Global = 360.0/8640.0
        Resolution_Name = "5KM"
    
    # Path for DA
    
    if mpi4py_rank == 0:
        
        # Decrease the float bit
        MODEL_CEA_X = MODEL_CEA_X - MODEL_X_Left
        MODEL_CEA_Y = MODEL_CEA_Y - MODEL_Y_Lower
        MODEL_X_Left = numpy.min(MODEL_CEA_X)
        MODEL_X_Right = numpy.max(MODEL_CEA_X)
        MODEL_Y_Lower = numpy.min(MODEL_CEA_Y)
        MODEL_Y_Upper = numpy.max(MODEL_CEA_Y)
        MODEL_X_Right = MODEL_X_Right - MODEL_X_Left
        MODEL_X_Left = MODEL_X_Left - MODEL_X_Left    
        MODEL_Y_Upper = MODEL_Y_Upper - MODEL_Y_Lower
        MODEL_Y_Lower = MODEL_Y_Lower - MODEL_Y_Lower
        
        r.assign('MODEL_X_Left', MODEL_X_Left)
        r.assign('MODEL_X_Right', MODEL_X_Right)
        r.assign('MODEL_Y_Lower', MODEL_Y_Lower)
        r.assign('MODEL_Y_Upper', MODEL_Y_Upper)
        print "MODEL_X_Left,MODEL_X_Right,MODEL_Y_Lower,MODEL_Y_Upper"
        print MODEL_X_Left,MODEL_X_Right,MODEL_Y_Lower,MODEL_Y_Upper
        
        Run_Dir_Multi_Instance = []
        Run_Dir_Array = []
        Forcing_File_Path_Array = []
        Forcing_File_Path_Array_Par = []
    
        if Ensemble_Number == 1:
            if not os.path.exists(Run_Dir_Home):
                os.makedirs(Run_Dir_Home)
            Run_Dir_Multi_Instance = Run_Dir_Home+"/"
            Run_Dir_Array.append(Run_Dir_Home+"/")
            Forcing_File_Path_Array.append(Forcing_File_Path_Home)
            Forcing_File_Path_Array_Par.append(Forcing_File_Path_Home)
        else:
            Run_Dir_Multi_Instance = Run_Dir_Home+"_Ens/"
            for Ens_Index in range(Ensemble_Number):
                Run_Dir_Temp = Run_Dir_Home+"_Ens"+str(Ens_Index+1)+"/"
                if Def_CESM_Multi_Instance == 1:
                    Run_Dir_Temp = Run_Dir_Multi_Instance
                Forcing_File_Path_Temp = Forcing_File_Path_Home+"_Ens"+str(Ens_Index+1)+"/"
                if not os.path.exists(Run_Dir_Temp):
                    os.makedirs(Run_Dir_Temp)
                Run_Dir_Array.append(Run_Dir_Temp)
                Forcing_File_Path_Array.append(Forcing_File_Path_Temp)
                Forcing_File_Path_Array_Par.append(Forcing_File_Path_Home)
        
        if not os.path.exists(Run_Dir_Multi_Instance):
            os.makedirs(Run_Dir_Multi_Instance)
    
    else:
        Run_Dir_Multi_Instance = None
        Run_Dir_Array = None
        Forcing_File_Path_Array = None
        Forcing_File_Path_Array_Par = None
    
    if Def_PP == 2:
        Run_Dir_Multi_Instance = mpi4py_comm.bcast(Run_Dir_Multi_Instance)
        Run_Dir_Array = mpi4py_comm.bcast(Run_Dir_Array)
        Forcing_File_Path_Array = mpi4py_comm.bcast(Forcing_File_Path_Array)
        Forcing_File_Path_Array_Par = mpi4py_comm.bcast(Forcing_File_Path_Array_Par)
    
    ######################################################### NC Files
    NC_FileName_Assimilation_2_Constant = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Constant.nc"
    NC_FileName_Assimilation_2_Diagnostic = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Diagnostic.nc"
    
    NC_FileName_Assimilation_2_Initial = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Initial.nc"
    NC_FileName_Assimilation_2_Initial_Copy = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Initial_Copy.nc"
    
    NC_FileName_Assimilation_2_Bias = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Bias.nc"
    NC_FileName_Assimilation_2_Bias_Copy = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Bias_Copy.nc"
    
    NC_FileName_Assimilation_2_Bias_Monthly = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Bias_Monthly.nc"
    NC_FileName_Assimilation_2_Bias_Monthly_Copy = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Bias_Monthly_Copy.nc"
    
    NC_FileName_Estimated_Bias = DAS_Output_Path+"Analysis/"+Region_Name+"/Estimated_Bias.nc"
    
    NC_FileName_Assimilation_2_Parameter = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter.nc"
    NC_FileName_Assimilation_2_Parameter_Copy = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_Copy.nc"
    NC_FileName_Assimilation_2_Parameter_Obs_Dim = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_Obs_Dim.nc"
    
    NC_FileName_Assimilation_2_Parameter_Monthly = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_Monthly.nc"
    NC_FileName_Assimilation_2_Parameter_Monthly_Copy = DAS_Output_Path+"Analysis/"+Region_Name+"/Assimilation_2_Parameter_Monthly_Copy.nc"
    
    NC_FileName_Optimized_Parameter = DAS_Output_Path+"Analysis/"+Region_Name+"/Optimized_Parameter.nc"
    NC_FileName_Soil_Moisture_Difference = DAS_Output_Path+"Analysis/"+Region_Name+"/Soil_Moisture_Difference.nc"
    
    NC_FileName_Parameter_Space_Single = DAS_Output_Path+"Analysis/"+Region_Name+"/Parameter_Space_Single.nc"
    DAS_File_Name_List = [NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic,
                           NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy,
                           NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy,
                           NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy,
                           NC_FileName_Estimated_Bias,
                           NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy, NC_FileName_Assimilation_2_Parameter_Obs_Dim,
                           NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy,
                           NC_FileName_Optimized_Parameter, NC_FileName_Soil_Moisture_Difference, NC_FileName_Parameter_Space_Single]
    
    ########################################################################################
    
    Soil_Thickness = numpy.asarray([0.01751282, 0.02757897, 0.04547003, 0.07496741, 0.1236004, 0.2037826, 0.3359806, 0.5539384, 0.91329, 1.50576, 2.48258, 4.093082, 6.748351, 11.12615, 13.85115])
        
    Variable_List = ["Soil_Moisture","Surface_Temperature","Vegetation_Temperature","Canopy_Water","Albedo_BSA_Band_vis","Albedo_BSA_Band_nir","Albedo_WSA_Band_vis",
                         "Albedo_WSA_Band_nir","Emissivity","Snow_Depth","Snow_Cover_Fraction","Snow_Water_Equivalent","LAI","Sensible_Heat"]
    
    Dim_CLM_State = len(Variable_List)
    Variable_Assimilation_Flag = numpy.zeros(Dim_CLM_State,dtype=numpy.float32)
    
    Initial_Perturbation_SM_Flag = numpy.array([0 for i in range(15)])   # whether to perturb the initial data
    Initial_Perturbation_ST_Flag = numpy.array([0 for i in range(15)])   # whether to perturb the initial data
        
    ################################################################################ CLM Input File Names
    finidat_initial_CLM_Copy = finidat_initial_CLM
    fndepdat_name = "fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc"
    fatmgrid_name = "griddata_"+Row_Numbers_String+"x"+Col_Numbers_String+".nc"
    fatmlndfrc_name = "domain.lnd_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    #fatmlndfrc_name = "fracdata_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    fsurdat_name = "surfdata_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    fglcmask_name = "glcmaskdata_"+Row_Numbers_String+"x"+Col_Numbers_String+"_Gland20km.nc"
    flndtopo_name = "topodata_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    fsnowoptics_name = "snicar_optics_5bnd_c090915.nc"
    fsnowaging_name = "snicar_drdt_bst_fit_60_c070416.nc"
    fpftcon_name = "pft-physiology.c130503.nc"
    domain_name = "domain.lnd_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    rdirc_name = "rdirc_0.5x0.5_simyr2000_c101124.nc"
    popd_streams_name = "clmforc.Li_2012_hdm_0.5x0.5_AVHRR_simyr1850-2010_c130401.nc"
    light_streams_name = "clmforc.Li_2012_climo1995-2011.T62.lnfm_c130327.nc"
    
    CLM_File_Name_List = [fndepdat_name, fatmgrid_name, fatmlndfrc_name, fsurdat_name, fglcmask_name, flndtopo_name,
                          fsnowoptics_name, fsnowaging_name, fpftcon_name, domain_name, rdirc_name, popd_streams_name, light_streams_name]
    
    
    ##################DA Parameters
    Dim_Observation_Quantity = 4
    
    # Irrigation Parameters
    irrig_nsteps_per_day = 1
    PFT_Num = 1
    PFT_Type_Index = 4
    irrig_nsteps_per_day = 3600.0 / dtime * Irrigation_Hours * PFT_Num
    
    N_Steps = []
    
    SensorType = []
    SensorVariable = []
    SensorQuantity = []
    Variable_ID = []
    QC_ID = []
    SensorResolution = []
    Observation_File_Name = []
    Soil_Layer_Index_DA = 0
    
    Def_First_Run_RTM = 1   # whether is RTM called at the first time
    
    if mpi4py_rank == 0:
        ##### Some Index Variables
        column_len = []
        pft_len = []
        finidat_name_string = []
        if Do_DA_Flag:
            
            finidat_name_string = Run_Dir_Home+"_Ens" + str(1) +"/"+ finidat_initial_CLM
            
            print '============================= Open the Model Initial File and Read the Index Data ==========================================='
            #------------------------------------------- Read the CLM Initial File
            print "Open Initial File:", finidat_name_string
            try:
                CLM_Initial_File = netCDF4.Dataset(finidat_name_string, 'r')
                column_len = len(CLM_Initial_File.dimensions['column'])
                pft_len = len(CLM_Initial_File.dimensions['pft'])
                CLM_Initial_File.close()
            except:
                print finidat_name_string,"not exists!!!!!!!!!!!!!!!!!!!!!"
                os.abort()
    else:
        finidat_name_string = None
        column_len = None
        pft_len = None
    
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
        finidat_name_string = mpi4py_comm.bcast(finidat_name_string)
        column_len = mpi4py_comm.bcast(column_len)
        pft_len = mpi4py_comm.bcast(pft_len)
        
    if mpi4py_rank == 0:
        ##### Some Index Variables
        
        
        #############################################################
        
        if not os.path.exists(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name):
            os.makedirs(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name)
    
        ####################################################################################
        
        Land_Mask_Data = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Teta_Residual = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Teta_Saturated = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Teta_Field_Capacity = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Teta_Wilting_Point = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        sucsat = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        bsw = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        watdry = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        watopt = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        watfc = numpy.ones((Soil_Layer_Num,Row_Numbers,Col_Numbers),dtype=numpy.float32)
        
        smpso = []
        smpsc = []
        
        Sand_Top_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Clay_Top_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Organic_Top_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Bulk_Density_Top_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Sand_Sub_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Clay_Sub_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Organic_Sub_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        Bulk_Density_Sub_Region = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        
        DEM_Data = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        
        PCT_PFT = numpy.zeros((maxpft, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        STD_ELEV = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)    # standard deviation of the elevation within a grid cell
        topo_slope = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        PCT_LAKE = numpy.ones((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        
        Irrigation_Grid_Flag = numpy.zeros((Row_Numbers,Col_Numbers),dtype=numpy.bool)
        
        micro_sigma = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        # Calculate the Vegetation Fraction
        PCT_Veg = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        PCT_PFT_High = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        PCT_PFT_Low = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        PCT_PFT_WATER = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        PFT_Dominant_Index = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Crop_Sum = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Bare_Grid_Index = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        print "Open Hydraulic_File_Name",Hydraulic_File_Name
        try:
            Hydraulic_File = netCDF4.Dataset(Hydraulic_File_Name, 'r')
            Teta_Residual = Hydraulic_File.variables['RES'][:,:,:]
            Teta_Saturated = Hydraulic_File.variables['SAT'][:,:,:]
            Teta_Field_Capacity = Hydraulic_File.variables['FC'][:,:,:]
            Teta_Wilting_Point = Hydraulic_File.variables['WP'][:,:,:]
            sucsat = Hydraulic_File.variables["sucsat"][:,:,:]
            bsw = Hydraulic_File.variables["bsw"][:,:,:]
            watdry = Hydraulic_File.variables["watdry"][:,:,:]
            watopt = Hydraulic_File.variables["watopt"][:,:,:]
            watfc = Hydraulic_File.variables["watfc"][:,:,:]
            Hydraulic_File.close()
        except:
            print "Open Hydraulic_File_Name",Hydraulic_File_Name,"Failed!!"
            os.abort()
        
        # Make sure the minimum soil moisture larger than 0.05, because we assume the maximum bias is 0.05
        Teta_Residual[numpy.where(Teta_Residual < 0.05)] = 0.05
        
        pft_physiology_file_name = DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name
        pft_physiology_file = netCDF4.Dataset(pft_physiology_file_name, "r")
        smpso = pft_physiology_file.variables["smpso"][:]
        smpsc = pft_physiology_file.variables["smpsc"][:]
        #smpso = numpy.zeros(len(pft_physiology_file.dimensions['pft']))
        #smpso[:] = -50000.0
        #print z0mr,PFT_Height_Top[1:17]
        pft_physiology_file.close()
        
        # Create Plot Folder
        if not os.path.exists(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name):
            os.makedirs(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name)
        
        if Use_Mask_Flag:
            # find the model grids which are need to be assiilated.
            # Read the Mask Grid
            print "Read the Land Water Mask Data"
            mkdatadomain_NC_FileName_In = DAS_Data_Path+"SysModel/CLM/tools/" + fatmlndfrc_name
            print "mkdatadomain_NC_FileName_In",mkdatadomain_NC_FileName_In
            mkdatadomain_NC_File_In = netCDF4.Dataset(mkdatadomain_NC_FileName_In, 'r+')
            Land_Mask_Data = numpy.flipud(mkdatadomain_NC_File_In.variables['mask'][:,:])            
            mkdatadomain_NC_File_In.close()
        
        Land_Mask_Data[numpy.where(Land_Mask_Data == 0.0)] = NAvalue
        
        Data = numpy.ma.masked_where(Land_Mask_Data == NAvalue, Land_Mask_Data)
        
        fig1 = plt.figure(figsize=(15, 10), dpi=80)
        ax = fig1.add_subplot(1, 1, 1)
        im1 = ax.imshow(Data, cmap=cm.jet)
        ax.set_title('Land_Mask_Data')
        plt.grid(True)
        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Land_Mask_Data.png")
        plt.close('all')
        
        del Data
        
        ####################################################################
        Corner_Row_Index = 0
        Corner_Col_Index = 0
        
        if Def_Region == -2:
            Grid_Resolution_GEO_Global = 360.0/432000.0*5.0
            Resolution_Name = "500m"
        elif Def_Region == -1:
            Grid_Resolution_GEO_Global = 360.0/432000.0
            Resolution_Name = "100m"
        elif Def_Region <= 8:
            Grid_Resolution_GEO_Global = 360.0/43200.0
            Resolution_Name = "1KM"
        elif Def_Region >= 9:
            Grid_Resolution_GEO_Global = 360.0/8640.0
            Resolution_Name = "5KM"
        
        Sand_Top_Region, Clay_Top_Region, Organic_Top_Region, Bulk_Density_Top_Region, Sand_Sub_Region, Clay_Sub_Region, Organic_Sub_Region, Bulk_Density_Sub_Region \
        = Read_Soil_Texture(Def_Region, DAS_Data_Path, Resolution_Name, Region_Name, Row_Numbers, Col_Numbers, Corner_Row_Index, Corner_Col_Index)
        
        print "*******************************Read CLM mksurfdata"
        mksurfdata_NC_FileName_In = DAS_Data_Path+"SysModel/CLM/tools/" + fsurdat_name
        print "mksurfdata_NC_FileName_In",mksurfdata_NC_FileName_In
        print "************************Open*******************",mksurfdata_NC_FileName_In
        mksurfdata_NC_File_In = netCDF4.Dataset(mksurfdata_NC_FileName_In, 'r')
        #print mksurfdata_NC_File_In.variables
        STD_ELEV = numpy.flipud(mksurfdata_NC_File_In.variables['STD_ELEV'][::])    # standard deviation of the elevation within a grid cell
        topo_slope = numpy.flipud(mksurfdata_NC_File_In.variables['SLOPE'][::])    # mean topographic slop
        DEM_Data = numpy.flipud(mksurfdata_NC_File_In.variables['TOPO'][::])    # mean elevation on land
        # check for near zero slopes, set minimum value
        topo_slope[numpy.where(topo_slope < 0.2)] = 0.2
        for Pft_index in range(maxpft):
            #print numpy.shape(mksurfdata_NC_File_In.variables['PCT_PFT'])
            PCT_PFT[Pft_index,::] = numpy.flipud(mksurfdata_NC_File_In.variables['PCT_PFT'][Pft_index,:,:])
        PCT_LAKE = numpy.flipud(mksurfdata_NC_File_In.variables['PCT_LAKE'][::])
        mksurfdata_NC_File_In.close()
        
        #################################
        
        ###########################3
        # microtopographic parameter, units are meters
        minslope=0.05
        slopemax=0.4
        maxslope=(slopemax - minslope)/(slopemax)
    
        # try smooth function of slope
        slopebeta=3.0
        slopemax=0.4
        slope0=slopemax**(-1.0/slopebeta)
        micro_sigma = (topo_slope + slope0)**(-slopebeta)
                
        # Calculate the Vegetation Fraction
        PCT_Veg = numpy.sum(PCT_PFT[1:maxpft,:,:],axis=0)
        PCT_PFT_High = numpy.sum(PCT_PFT[1:9,:,:],axis=0) / 100.0
        PCT_PFT_Low = numpy.sum(PCT_PFT[9:maxpft,:,:],axis=0) / 100.0
        PCT_PFT_WATER = PCT_LAKE / 100.0
        PFT_Dominant_Index = numpy.argmax(PCT_PFT,axis=0)
        numpy.savetxt("PFT_Dominant_Index_"+Region_Name+".txt",PFT_Dominant_Index)
        w,h = plt.figaspect(float(Row_Numbers)/Col_Numbers)
        fig1 = plt.figure(figsize=(w,h))
        ax1 = fig1.add_subplot(1,1,1)
        im1 = ax1.imshow(PFT_Dominant_Index, cmap=cm.jet, interpolation='bilinear')
        plt.colorbar(im1)
        if Def_Figure_Output:
            plt.savefig("DataBase/PFT_Dominant_Index_"+Region_Name+".png")
        #plt.show()
        
        
        Bare_Grid_Index = numpy.where(PFT_Dominant_Index == 0)
        #Bare_Grid_Index = numpy.where(PCT_PFT[0,:,:] == 100)
        print "numpy.size(Bare_Grid_Index)",numpy.size(Bare_Grid_Index)
        for Soil_Layer_Index in range(Soil_Layer_Num):
            watopt[Soil_Layer_Index,:,:][Bare_Grid_Index] = watfc[Soil_Layer_Index,:,:][Bare_Grid_Index]
            watdry[Soil_Layer_Index,:,:][Bare_Grid_Index] = Teta_Residual[Soil_Layer_Index,:,:][Bare_Grid_Index]
            #watopt[Soil_Layer_Index,:,:] = Teta_Saturated[Soil_Layer_Index,:,:]
            #watdry[Soil_Layer_Index,:,:] = Teta_Residual[Soil_Layer_Index,:,:]
        
        print "************************"
        print "numpy.max(watopt),numpy.min(watopt),numpy.max(watdry),numpy.min(watdry)"
        print numpy.max(watopt),numpy.min(watopt),numpy.max(watdry),numpy.min(watdry)
        #    w,h = plt.figaspect(float(Row_Numbers)/Col_Numbers)
        #    fig1 = plt.figure(figsize=(w,h))
        #    ax1 = fig1.add_subplot(1,1,1)
        #    im1 = ax1.imshow(MONTHLY_LAI[6,:,:], cmap=cm.jet, interpolation='bilinear')
        #    plt.colorbar(im1)
        #    if Def_Figure_Output:
        #        plt.savefig("SysModel/CLM/Surfdata_Figure/"+Region_Name+"_PFT/"+"MONTHLY_LAI.png")
        #    plt.show()
        
        
        ##############################################################################
        # COSMOS Circle Mask
        Mask_X_COSMOS = MODEL_CEA_X
        Mask_Y_COSMOS = MODEL_CEA_Y
        
        COSMOS_Circle_Plot = numpy.zeros((Row_Numbers,Col_Numbers),dtype=numpy.float32)
        
        COSMOS_Circle_Array = []
        COSMOS_Circle_Index_Array = []
        COSMOS_Circle_Num_Array = []
        #for Station_Index in range(numpy.size(Station_XY)/2):
        for Station_Index in range(numpy.size(Station_XY)/2-12):
            COSMOS_Circle = numpy.zeros((Row_Numbers,Col_Numbers),dtype=numpy.bool)
            COSMOS_Circle[::] = False
        
            print "Station_"+str(Station_Index+1),Station_XY[Station_Index][0],Station_XY[Station_Index][1]
            r.assign('X_Coordiates',Station_XY[Station_Index][0])
            r.assign('Y_Coordiates',Station_XY[Station_Index][1])                    
            r('xy <- cbind(X_Coordiates,Y_Coordiates)')
            print r['xy']
            
            print "========================== GEO to CEA"
            r('xy_cea <- project(xy,"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +ellps=WGS84 +no_defs")')
            print 'x,y',r['xy_cea'][0][0],r['xy_cea'][0][1]
            ii = r['xy_cea'][0][0]
            jj = r['xy_cea'][0][1]
            dist = numpy.sqrt(abs(ii - Mask_X_COSMOS) ** 2 + abs(jj - Mask_Y_COSMOS) ** 2)
            COSMOS_Circle_Index = numpy.where(dist <= 300)
            COSMOS_Circle[COSMOS_Circle_Index] = True
            COSMOS_Circle_Array.append(COSMOS_Circle)
            COSMOS_Circle_Index_Array.append(COSMOS_Circle_Index)
            #print COSMOS_Circle_Index,numpy.size(COSMOS_Circle_Index)
            COSMOS_Circle_Num = numpy.zeros_like(COSMOS_Circle_Index)
            COSMOS_Circle_Num = numpy.floor(dist[COSMOS_Circle_Index] / Grid_Resolution_CEA)
            #print COSMOS_Circle_Num
            COSMOS_Circle_Num_Array.append(COSMOS_Circle_Num)
        
            COSMOS_Circle_Plot[COSMOS_Circle] = 1.0
        
        if Plot_Analysis:
            fig1 = plt.figure(figsize=(15, 10), dpi=80)
            ax = fig1.add_subplot(1, 1, 1)
            im1 = ax.imshow(COSMOS_Circle_Plot, cmap=cm.jet)
            ax.set_title('COSMOS_Circle_Plot')
            plt.grid(True)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/COSMOS_Circle_Plot.png")
            plt.close('all')
        
        ##########################################################
        LONGXY_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        LATIXY_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        
        longitudes = numpy.linspace(mksrf_edgew+Grid_Resolution_GEO[0]/2.0,mksrf_edgee-Grid_Resolution_GEO[0]/2.0,Col_Numbers)
        latitudes = numpy.linspace(mksrf_edges+Grid_Resolution_GEO[1]/2.0,mksrf_edgen-Grid_Resolution_GEO[1]/2.0,Row_Numbers)
        LONGXY_Row = longitudes
        LATIXY_Col = latitudes
        for row in range(Row_Numbers):
            LONGXY_Mat[row,:] = LONGXY_Row
        for col in range(Col_Numbers):
            LATIXY_Mat[:,col] = LATIXY_Col
        #print LATIXY_Col
                    
        Mask_Index = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers), dtype=numpy.bool)
        Mask_Index[:,:,:] = False
        
        #------------------------------------------- Data Assimilation Flags
        print "Soil Moisture Products: SMAP(10km), AMSR-E(25km), SMOS(40km), ASCAT(12.5km,25km), MODIS(1km), ASAR(120m), PALSAR(60m)"
        print "Soil Temperature Products: MODIS Terra and Aqua(1km)"        
        
        ############################################################## For Bias Estimation #######################################
        Bias_Remove_Start_Time_Array = ['' for i in range(Dim_CLM_State)]
        
        # Flag to check wehter the Observation Bias of each observation type and each observation ensemble has been perturbed
        Observation_Bias_Initialization_Flag = numpy.zeros((Dim_CLM_State,Dim_Observation_Quantity,Ensemble_Number),dtype=numpy.float32)
        
        Model_Bias_Optimized = numpy.zeros((Ensemble_Number, Dim_CLM_State, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Observation_Bias_Optimized = numpy.zeros((Ensemble_Number, Dim_CLM_State, Dim_Observation_Quantity, numpy.size(Station_XY)/2), dtype=numpy.float32)
        
        # Bias Estimation Range and Standard Deviation Defination
        Model_Bias_Range = numpy.zeros((Dim_CLM_State,2),dtype=numpy.float32)
        Observation_Bias_Range = numpy.zeros((Dim_CLM_State,Dim_Observation_Quantity,2),dtype=numpy.float32)
        
        Model_Bias_Range_STD = numpy.zeros((Dim_CLM_State,2),dtype=numpy.float32)
        Observation_Bias_Range_STD = numpy.zeros((Dim_CLM_State,Dim_Observation_Quantity,2),dtype=numpy.float32)
        Model_Bias_STD = numpy.zeros(Dim_CLM_State,dtype=numpy.float32)
        Observation_Bias_STD = numpy.zeros((Dim_CLM_State,Dim_Observation_Quantity),dtype=numpy.float32)
        
        # Model State Ensemble Inflation STD
        Model_State_Inflation_Range = numpy.zeros((Dim_CLM_State,2),dtype=numpy.float32)
        Model_State_Inflation_Range_STD = numpy.zeros(Dim_CLM_State,dtype=numpy.float32)
        
        ########### Simulate Model Error
        Additive_Noise_SM_Par = numpy.zeros((10,11),dtype=numpy.float32)
        Additive_Noise_SM_Par[::] = numpy.array([[1.00E-3, 1.0, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6, 0.4, 0.4, 0.4],
                                            [7.00E-4, 0.7, 1.0, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6, 0.4, 0.4],
                                            [5.00E-4, 0.7, 0.7, 1.0, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6, 0.4],
                                            [3.00E-4, 0.6, 0.7, 0.7, 1.0, 0.7, 0.7, 0.6, 0.6, 0.6, 0.6],
                                            [2.00E-5, 0.6, 0.6, 0.7, 0.7, 1.0, 0.7, 0.7, 0.6, 0.6, 0.6],
                                            [2.00E-5, 0.6, 0.6, 0.6, 0.7, 0.7, 1.0, 0.7, 0.7, 0.6, 0.6],
                                            [2.00E-5, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 1.0, 0.7, 0.7, 0.6],
                                            [1.50E-6, 0.4, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 1.0, 0.7, 0.7],
                                            [1.50E-6, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 1.0, 0.7],
                                            [5.00E-8, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 1.0]])
        #print numpy.shape(Additive_Noise_SM_Par)
        
        Additive_Noise_SM = numpy.zeros((Ensemble_Number,Soil_Layer_Num-5),dtype=numpy.float32)
        Additive_Noise_ST = numpy.zeros((Ensemble_Number,2),dtype=numpy.float32)
        
        Irrigation_Grid_Flag_Array = []
        
        cols1d_ixy = numpy.zeros(column_len, dtype=numpy.integer)
        cols1d_jxy = numpy.zeros(column_len, dtype=numpy.integer)
        cols1d_ityplun = numpy.zeros(column_len, dtype=numpy.integer)
        pfts1d_ixy = numpy.zeros(pft_len, dtype=numpy.integer)
        pfts1d_jxy = numpy.zeros(pft_len, dtype=numpy.integer)
        pfts1d_itypveg = numpy.zeros(pft_len, dtype=numpy.integer)
        pfts1d_ci = numpy.zeros(pft_len, dtype=numpy.integer)
        pfts1d_ityplun = numpy.zeros(pft_len, dtype=numpy.integer)
        
        ##### Some Index Variables
        if Do_DA_Flag:
            
            finidat_name_string = Run_Dir_Home+"_Ens" + str(1) +"/"+ finidat_initial_CLM
            
            print '============================= Open the Model Initial File and Read the Index Data ==========================================='
            #------------------------------------------- Read the CLM Initial File
            print "Open Initial File:", finidat_name_string
            try:
                CLM_Initial_File = netCDF4.Dataset(finidat_name_string, 'r')
                cols1d_ixy[:] = CLM_Initial_File.variables['cols1d_ixy'][:]
                cols1d_jxy[:] = CLM_Initial_File.variables['cols1d_jxy'][:]
                cols1d_ityplun[:] = CLM_Initial_File.variables['cols1d_ityplun'][:]
                
                #numpy.savetxt('cols1d_ixy',cols1d_ixy)
                #numpy.savetxt('cols1d_jxy',cols1d_jxy)
                pfts1d_ixy[:] = CLM_Initial_File.variables['pfts1d_ixy'][:]
                pfts1d_jxy[:] = CLM_Initial_File.variables['pfts1d_jxy'][:]
                pfts1d_itypveg[:] = CLM_Initial_File.variables['pfts1d_itypveg'][:]
                pfts1d_ci[:] = CLM_Initial_File.variables['pfts1d_ci'][:]
                pfts1d_ityplun[:] = CLM_Initial_File.variables['pfts1d_ityplun'][:]
                
                CLM_Initial_File.close()
            
            except:
                print finidat_name_string,"not exists!!!!!!!!!!!!!!!!!!!!!"
                os.abort()
        
        Analysis_Variable_Name = ['' for i in range(Dim_CLM_State)]
        
        Soil_Sand_Clay_Sum = numpy.zeros((Soil_Texture_Layer_Opt_Num, Row_Numbers, Col_Numbers), dtype=numpy.float32)
        
        
        print "################ Go to CLM"
        Parameter_Soil_Optimized = numpy.zeros((Ensemble_Number, Dim_Soil_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Parameter_PFT_Optimized = numpy.zeros((Ensemble_Number, Dim_PFT_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Parameter_Hard_Optimized = numpy.zeros((Ensemble_Number, Dim_Hard_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Parameter_Soil_PSRF = numpy.zeros((Dim_Soil_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Parameter_PFT_PSRF = numpy.zeros((Dim_PFT_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        Parameter_Hard_PSRF = numpy.zeros((Dim_Hard_Par, numpy.size(Station_XY)/2), dtype=numpy.float32)
        
        Parameter_Optimization_First_Flag = True
        
        
        Mean_Index_Prop_Grid_Array_Sys = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Model_Variance = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask = numpy.zeros((Dim_CLM_State, 3, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        Mask_X = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask_Y = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask_X[::] = MODEL_CEA_X
        Mask_Y[::] = MODEL_CEA_Y
        
        ###################################### CMEM Matrix
        
        
        Clay_Fraction = []
        Sand_Fraction = []
        Soil_Density = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
                                
        CMEM_Work_Path_Array = []
        Clay_Mat = []
        Sand_Mat = []
        ECOCVL_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOCVH_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVL_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVH_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOWAT_Mat = numpy.zeros((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        
        # Folder to Save Ensemble Mean
        Mean_Dir = Run_Dir_Home+"_Ens_Mean"
        if not os.path.exists(Mean_Dir):
            os.makedirs(Mean_Dir)
                
        if not os.path.exists(DAS_Output_Path+"Analysis/"+Region_Name):
            os.makedirs(DAS_Output_Path+"Analysis/"+Region_Name)
        
        for Block_Index in range(Sub_Block_Ratio_Row*Sub_Block_Ratio_Col):
            if not os.path.exists(DAS_Output_Path+"Analysis/"+Region_Name+"/Block_"+str(Block_Index+1)):
                os.makedirs(DAS_Output_Path+"Analysis/"+Region_Name+"/Block_"+str(Block_Index+1))
        
        
        NC_File_In = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name, 'r')
        File_Format1 = NC_File_In.file_format
        NC_File_In.close()
        
        NC_File_In = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/tools/"+fatmlndfrc_name, 'r')
        File_Format2 = NC_File_In.file_format
        NC_File_In.close()
        
        NC_File_In = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name, 'r')
        File_Format3 = NC_File_In.file_format
        NC_File_In.close()
        
        NC_File_In = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"+rdirc_name, 'r')
        File_Format4 = NC_File_In.file_format
        NC_File_In.close()
        
        NC_File_In = netCDF4.Dataset(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/ndepdata/"+fndepdat_name, 'r')
        File_Format5 = NC_File_In.file_format
        NC_File_In.close()
        
        if Def_First_Run == 1 and ((File_Format1 == 'NETCDF3_CLASSIC') or (File_Format1 == 'NETCDF3_64BIT') and \
                              (File_Format2 == 'NETCDF3_CLASSIC') or (File_Format2 == 'NETCDF3_64BIT') and \
                              (File_Format3 == 'NETCDF3_CLASSIC') or (File_Format3 == 'NETCDF3_64BIT') and \
                              (File_Format4 == 'NETCDF3_CLASSIC') or (File_Format4 == 'NETCDF3_64BIT') and \
                              (File_Format5 == 'NETCDF3_CLASSIC') or (File_Format5 == 'NETCDF3_64BIT') ):
            
            print "Convert netCDF3 input to netCDF4 for CLM"
            
            subprocess.call(DAS_Depends_Path+"bin/nccopy -k 3 "+DAS_Data_Path + "SysModel/CLM/tools/"+fatmlndfrc_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fatmlndfrc_name,shell=True)
            os.remove(DAS_Data_Path + "SysModel/CLM/tools/"+fatmlndfrc_name)
            subprocess.call(DAS_Depends_Path+"bin/nccopy -d 4 "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fatmlndfrc_name+" "+DAS_Data_Path + "SysModel/CLM/tools/"+fatmlndfrc_name,shell=True)
            
            subprocess.call(DAS_Depends_Path+"bin/nccopy -k 3 "+DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fsurdat_name,shell=True)
            os.remove(DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name)
            subprocess.call(DAS_Depends_Path+"bin/nccopy -d 4 "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fsurdat_name+" "+DAS_Data_Path + "SysModel/CLM/tools/"+fsurdat_name,shell=True)
            
            subprocess.call(DAS_Depends_Path+"bin/nccopy -k 3 "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fpftcon_name,shell=True)
            os.remove(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name)
            subprocess.call(DAS_Depends_Path+"bin/nccopy -d 4 "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fpftcon_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/pftdata/"+fpftcon_name,shell=True)
            
            subprocess.call(DAS_Depends_Path+"bin/nccopy -k 3 "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"+rdirc_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/"+rdirc_name,shell=True)
            os.remove(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"+rdirc_name)
            subprocess.call(DAS_Depends_Path+"bin/nccopy -d 4 "+DAS_Data_Path + "SysModel/CLM/inputdata/"+rdirc_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/rtmdata/"+rdirc_name,shell=True)
            
            subprocess.call(DAS_Depends_Path+"bin/nccopy -k 3 "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/ndepdata/"+fndepdat_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fndepdat_name,shell=True)
            os.remove(DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/ndepdata/"+fndepdat_name)
            subprocess.call(DAS_Depends_Path+"bin/nccopy -d 4 "+DAS_Data_Path + "SysModel/CLM/inputdata/"+fndepdat_name+" "+DAS_Data_Path + "SysModel/CLM/inputdata/lnd/clm2/ndepdata/"+fndepdat_name,shell=True)        
        
        if (Def_First_Run == -1) and Ensemble_Number > 1:
            copyLargeFile(NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Initial)
            copyLargeFile(NC_FileName_Assimilation_2_Parameter_Copy, NC_FileName_Assimilation_2_Parameter)
            copyLargeFile(NC_FileName_Assimilation_2_Parameter_Monthly_Copy, NC_FileName_Assimilation_2_Parameter_Monthly)
            copyLargeFile(NC_FileName_Assimilation_2_Bias_Copy, NC_FileName_Assimilation_2_Bias)
            copyLargeFile(NC_FileName_Assimilation_2_Bias_Monthly_Copy, NC_FileName_Assimilation_2_Bias_Monthly)
            
        
        if Def_First_Run == 1:
            
            print "**************** Prepare Initial netCDF file"
            if os.path.exists(NC_FileName_Assimilation_2_Constant):
                os.remove(NC_FileName_Assimilation_2_Constant)
            
            print 'Write NetCDF File:',NC_FileName_Assimilation_2_Constant
            NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_Assimilation_2_Constant.createDimension('lon', Col_Numbers)
            NC_File_Out_Assimilation_2_Constant.createDimension('lat', Row_Numbers)
            NC_File_Out_Assimilation_2_Constant.createDimension('Soil_Layer_Num', Soil_Layer_Num)
            NC_File_Out_Assimilation_2_Constant.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
            NC_File_Out_Assimilation_2_Constant.createDimension('Ensemble_Number', Ensemble_Number)
            NC_File_Out_Assimilation_2_Constant.createDimension('Dim_CLM_State', Dim_CLM_State)
            NC_File_Out_Assimilation_2_Constant.createDimension('maxpft', maxpft)
            
            NC_File_Out_Assimilation_2_Constant.createVariable('Land_Mask_Data','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Land_Mask_Data'][:,:] = Land_Mask_Data
            NC_File_Out_Assimilation_2_Constant.createVariable('PCT_LAKE','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['PCT_LAKE'][:,:] = PCT_LAKE
            NC_File_Out_Assimilation_2_Constant.createVariable('PCT_Veg','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['PCT_Veg'][:,:] = PCT_Veg
            NC_File_Out_Assimilation_2_Constant.createVariable('PCT_PFT','f4',('maxpft','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['PCT_PFT'][:,:] = PCT_PFT
            #print "numpy.mean(numpy.sum(PCT_PFT[:,:,:],axis=0))",numpy.mean(numpy.sum(PCT_PFT[:,:,:],axis=0))
            NC_File_Out_Assimilation_2_Constant.createVariable('STD_ELEV','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['STD_ELEV'][:,:] = STD_ELEV
            NC_File_Out_Assimilation_2_Constant.createVariable('DEM_Data','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['DEM_Data'][:,:] = DEM_Data
            NC_File_Out_Assimilation_2_Constant.createVariable('Bulk_Density_Top_Region','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Top_Region'][:,:] = Bulk_Density_Top_Region
            NC_File_Out_Assimilation_2_Constant.createVariable('Bulk_Density_Sub_Region','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Bulk_Density_Sub_Region'][:,:] = Bulk_Density_Sub_Region
                    
            NC_File_Out_Assimilation_2_Constant.createVariable('Teta_Residual','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Teta_Residual'][:,:,:] = Teta_Residual
            
            NC_File_Out_Assimilation_2_Constant.createVariable('Teta_Saturated','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Teta_Saturated'][:,:,:] = Teta_Saturated
            
            NC_File_Out_Assimilation_2_Constant.createVariable('Teta_Field_Capacity','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Teta_Field_Capacity'][:,:,:] = Teta_Field_Capacity
            
            NC_File_Out_Assimilation_2_Constant.createVariable('Teta_Wilting_Point','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['Teta_Wilting_Point'][:,:,:] = Teta_Wilting_Point
            
            NC_File_Out_Assimilation_2_Constant.createVariable('watopt','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['watopt'][:,:,:] = watopt
            
            NC_File_Out_Assimilation_2_Constant.createVariable('watdry','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['watdry'][:,:,:] = watdry
            
            NC_File_Out_Assimilation_2_Constant.createVariable('watfc','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['watfc'][:,:,:] = watfc
            
            NC_File_Out_Assimilation_2_Constant.createVariable('PFT_Dominant_Index','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.variables['PFT_Dominant_Index'][:,:] = PFT_Dominant_Index
            
            del Land_Mask_Data,PCT_LAKE,PCT_Veg,PCT_PFT,STD_ELEV,DEM_Data
            del Bulk_Density_Top_Region,Bulk_Density_Sub_Region,Teta_Residual,Teta_Saturated,Teta_Field_Capacity,Teta_Wilting_Point,watopt,watdry,watfc
            
            NC_File_Out_Assimilation_2_Constant.createVariable('CLM_Soil_Layer_Thickness','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.createVariable('CLM_Soil_Layer_Thickness_Cumsum','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.createVariable('Soil_Layer_Thickness_Ratio_Moisture','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Constant.createVariable('Soil_Layer_Thickness_Ratio_Temperature','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            
            
            NC_File_Out_Assimilation_2_Constant.sync()
            NC_File_Out_Assimilation_2_Constant.close()
        
            NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r+', format='NETCDF4')        
            # Meters
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][0, :, :] = Soil_Thickness[0]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][1, :, :] = Soil_Thickness[1]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][2, :, :] = Soil_Thickness[2]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][3, :, :] = Soil_Thickness[3]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][4, :, :] = Soil_Thickness[4]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][5, :, :] = Soil_Thickness[5]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][6, :, :] = Soil_Thickness[6]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][7, :, :] = Soil_Thickness[7]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][8, :, :] = Soil_Thickness[8]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][9, :, :] = Soil_Thickness[9]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][10, :, :] = Soil_Thickness[10]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][11, :, :] = Soil_Thickness[11]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][12, :, :] = Soil_Thickness[12]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][13, :, :] = Soil_Thickness[13]
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][14, :, :] = Soil_Thickness[14]
            
            NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness_Cumsum'][:,:,:] = numpy.cumsum(numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][:,:,:]), axis=0)
            
            for Soil_Layer_Index in range(Soil_Layer_Num):
                #NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Moisture'][Soil_Layer_Index, :, :] = numpy.exp(-1.0*NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness_Cumsum'][Soil_Layer_Index,:,:])
                NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Moisture'][Soil_Layer_Index, :, :] = numpy.exp(-1.0*NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][Soil_Layer_Index,:,:])
                
            for Soil_Layer_Index in range(Soil_Layer_Num):
                #NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][Soil_Layer_Index, :, :] = numpy.exp(-1.0*numpy.log(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness_Cumsum'][Soil_Layer_Index,:,:]))
                NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][Soil_Layer_Index, :, :] = numpy.exp(-1.0*numpy.log(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][Soil_Layer_Index,:,:]))
            
            Ratio_Temp = NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][0, :, :]
            for Soil_Layer_Index in range(Soil_Layer_Num):
                NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][Soil_Layer_Index, :, :] = \
                numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['Soil_Layer_Thickness_Ratio_Temperature'][Soil_Layer_Index, :, :]) / Ratio_Temp
            
            NC_File_Out_Assimilation_2_Constant.sync()
            NC_File_Out_Assimilation_2_Constant.close()
                    
            print "**************** Prepare Initial netCDF file"
            if os.path.exists(NC_FileName_Assimilation_2_Diagnostic):
                os.remove(NC_FileName_Assimilation_2_Diagnostic)
            
            print 'Write NetCDF File:',NC_FileName_Assimilation_2_Diagnostic
            NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('lon', Col_Numbers)
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('lat', Row_Numbers)
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('Soil_Layer_Num', Soil_Layer_Num)
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('Ensemble_Number', Ensemble_Number)
            NC_File_Out_Assimilation_2_Diagnostic.createDimension('Dim_CLM_State', Dim_CLM_State)
                        
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Initial_SM_Noise','f4',('Ensemble_Number','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_SM_Noise'][:,:,:] = 0.0
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Initial_ST_Noise','f4',('Ensemble_Number','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.variables['Initial_ST_Noise'][:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Mask_Index','i4',('Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('CLM_Soil_Temperature_Ratio_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('CLM_Soil_Moisture_Ratio_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat'][:,:,:] = 1.0
            NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Temperature_Ratio_Ensemble_Mat'][:,:,:] = 1.0
            NC_File_Out_Assimilation_2_Diagnostic.variables['CLM_Soil_Moisture_Ratio_Ensemble_Mat_MultiScale'][:,:,:] = 1.0
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Analysis_Grid_Array','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Innovation_State','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Increments_State','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('Observation','f4',('Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('CLM_2m_Air_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.createVariable('CLM_Air_Pressure_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Diagnostic.sync()
            NC_File_Out_Assimilation_2_Diagnostic.close()
            
            print "**************** Prepare Initial netCDF file"
            if os.path.exists(NC_FileName_Assimilation_2_Initial):
                os.remove(NC_FileName_Assimilation_2_Initial)
            if os.path.exists(NC_FileName_Assimilation_2_Initial_Copy):
                os.remove(NC_FileName_Assimilation_2_Initial_Copy)
            
            print 'Write NetCDF File:',NC_FileName_Assimilation_2_Initial
            NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_Assimilation_2_Initial.createDimension('lon', Col_Numbers)
            NC_File_Out_Assimilation_2_Initial.createDimension('lat', Row_Numbers)
            NC_File_Out_Assimilation_2_Initial.createDimension('Soil_Layer_Num', Soil_Layer_Num)
            NC_File_Out_Assimilation_2_Initial.createDimension('Ensemble_Number', Ensemble_Number)
            NC_File_Out_Assimilation_2_Initial.createDimension('Dim_CLM_State', Dim_CLM_State)
            
            NC_File_Out_Assimilation_2_Initial.createVariable('Prop_Grid_Array_Sys_Parallel','f4',('Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('Prop_Grid_Array_H_Trans_Paralle','f4',('Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Moisture_Ensemble_Mat_Parallel','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Temperature_Ensemble_Mat_Parallel','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Initial.createVariable('Prop_Grid_Array_Sys','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('Prop_Grid_Array_H_Trans','f4',('Ensemble_Number','Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Moisture_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Temperature_Ensemble_Mat','f4',('Soil_Layer_Num','lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Vegetation_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Ground_Temperature_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Snow_Depth_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Snow_Water_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_INT_SNOW_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_FH2OSFC_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            #onset freezing degree days counters
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Onset_Freezing_Degree_Days_Counter_Ensemble_Mat','f4',('lat','lon','Ensemble_Number',),zlib=True,least_significant_digit=None)
            
            NC_File_Out_Assimilation_2_Initial.createVariable('Prop_Grid_Array_Sys_parm_infl','f4',('Dim_CLM_State','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Moisture_parm_infl','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Soil_Temperature_parm_infl','f4',('Soil_Layer_Num','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Vegetation_Temperature_parm_infl','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Ground_Temperature_parm_infl','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Initial.createVariable('CLM_Surface_Temperature_parm_infl','f4',('lat','lon',),zlib=True,least_significant_digit=None)
            
            
            NC_File_Out_Assimilation_2_Initial.sync()
            NC_File_Out_Assimilation_2_Initial.close()
        
            NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r+', format='NETCDF4')
            NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys_parm_infl'][:,:,:] = numpy.abs(msw_infl[0])     
            NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_parm_infl'][:,:,:] = numpy.abs(msw_infl[0])
            NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_parm_infl'][:,:,:] = numpy.abs(msw_infl[0])
            NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_parm_infl'][:,:] = numpy.abs(msw_infl[0])
            NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_parm_infl'][:,:] = numpy.abs(msw_infl[0])
            NC_File_Out_Assimilation_2_Initial.variables['CLM_Surface_Temperature_parm_infl'][:,:] = numpy.abs(msw_infl[0])
                
            NC_File_Out_Assimilation_2_Initial.sync()
            NC_File_Out_Assimilation_2_Initial.close()
            
            
            print "**************** Prepare Parameter netCDF file"
            if os.path.exists(NC_FileName_Assimilation_2_Parameter):
                os.remove(NC_FileName_Assimilation_2_Parameter)
            if os.path.exists(NC_FileName_Assimilation_2_Parameter_Copy):
                os.remove(NC_FileName_Assimilation_2_Parameter_Copy)
            if os.path.exists(NC_FileName_Assimilation_2_Parameter_Obs_Dim):
                os.remove(NC_FileName_Assimilation_2_Parameter_Obs_Dim)
            
            if os.path.exists(NC_FileName_Assimilation_2_Parameter_Monthly):
                os.remove(NC_FileName_Assimilation_2_Parameter_Monthly)
            if os.path.exists(NC_FileName_Assimilation_2_Parameter_Monthly_Copy):
                os.remove(NC_FileName_Assimilation_2_Parameter_Monthly_Copy)
                
            print 'Write NetCDF File:',NC_FileName_Assimilation_2_Parameter
            NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_Assimilation_2_Parameter.createDimension('lon', Col_Numbers)
            NC_File_Out_Assimilation_2_Parameter.createDimension('lat', Row_Numbers)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Soil_Layer_Num', Soil_Layer_Num)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Ensemble_Number', Ensemble_Number)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Ensemble_Number_Predict', Ensemble_Number_Predict)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Dim_CLM_State', Dim_CLM_State)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Dim_Soil_Par', Dim_Soil_Par)
            NC_File_Out_Assimilation_2_Parameter.createDimension('Dim_PFT_Par', Dim_PFT_Par)
            NC_File_Out_Assimilation_2_Parameter.createDimension('maxpft', maxpft)
            
            NC_File_Out_Assimilation_2_Parameter.createVariable('Parameter_Soil_Space_Ensemble','f4',('Ensemble_Number','Dim_Soil_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Parameter.createVariable('Parameter_Soil_Space_parm_infl','f4',('Dim_Soil_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_parm_infl'][:,:,:] = numpy.abs(msw_infl[1])
        
            NC_File_Out_Assimilation_2_Parameter.createVariable('Parameter_PFT_Space_Ensemble','f4',('Ensemble_Number','Dim_PFT_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,:,:] = 0.0
            NC_File_Out_Assimilation_2_Parameter.createVariable('Parameter_PFT_Space_parm_infl','f4',('Dim_PFT_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_parm_infl'][:,:,:] = numpy.abs(msw_infl[1])
            
            NC_File_Out_Assimilation_2_Parameter.sync()            
            NC_File_Out_Assimilation_2_Parameter.close()
            
            print 'Write NetCDF File:',NC_FileName_Assimilation_2_Parameter_Obs_Dim
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter_Obs_Dim, 'w', diskless=True, persist=True, format='NETCDF4')
            # Dim the dimensions of NetCDF
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('lon', Col_Numbers)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('lat', Row_Numbers)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Soil_Layer_Num', Soil_Layer_Num)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('ParFlow_Layer_Num', ParFlow_Layer_Num)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Ensemble_Number', Ensemble_Number)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Ensemble_Number_Predict', Ensemble_Number_Predict)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Dim_CLM_State', Dim_CLM_State)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Dim_Soil_Par', Dim_Soil_Par)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Dim_Veg_Par', Dim_Veg_Par)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Dim_PFT_Par', Dim_PFT_Par)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('Dim_Hard_Par', Dim_Hard_Par)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createDimension('maxpft', maxpft)
            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createVariable('Parameter_Soil_Space_Ensemble_Obs_Dim','f4',('Ensemble_Number','Dim_Soil_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Soil_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createVariable('Parameter_Veg_Space_Ensemble_Obs_Dim','f4',('Ensemble_Number','Dim_Veg_Par','maxpft',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Veg_Space_Ensemble_Obs_Dim'][:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createVariable('Parameter_Hard_Space_Ensemble_Obs_Dim','f4',('Ensemble_Number','Dim_Hard_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Hard_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createVariable('Parameter_Veg_Space_Ensemble_Matrix_Obs_Dim','f4',('Ensemble_Number','Dim_Veg_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_Veg_Space_Ensemble_Matrix_Obs_Dim'][:,:,:,:] = 0.0
            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.createVariable('Parameter_PFT_Space_Ensemble_Obs_Dim','f4',('Ensemble_Number','Dim_PFT_Par','lat','lon',),zlib=True,least_significant_digit=None)
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.variables['Parameter_PFT_Space_Ensemble_Obs_Dim'][:,:,:,:] = 0.0
            
                
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.sync()            
            NC_File_Out_Assimilation_2_Parameter_Obs_Dim.close()
            
                
            if Parameter_Optimization:
                if os.path.exists(NC_FileName_Optimized_Parameter):
                    os.remove(NC_FileName_Optimized_Parameter)
                print 'Write NetCDF File:',NC_FileName_Optimized_Parameter
                NC_File_Out_Optimized_Parameter = netCDF4.Dataset(NC_FileName_Optimized_Parameter, 'w', diskless=True, persist=True, format='NETCDF4')
                # Dim the dimensions of NetCDF
                NC_File_Out_Optimized_Parameter.createDimension('lon', Col_Numbers)
                NC_File_Out_Optimized_Parameter.createDimension('lat', Row_Numbers)
                NC_File_Out_Optimized_Parameter.createDimension('Soil_Layer_Num', Soil_Layer_Num)
                NC_File_Out_Optimized_Parameter.createDimension('Ensemble_Number', Ensemble_Number)
                NC_File_Out_Optimized_Parameter.createDimension('Dim_Soil_Par', Dim_Soil_Par)
                NC_File_Out_Optimized_Parameter.createDimension('Dim_PFT_Par', Dim_PFT_Par)
                NC_File_Out_Optimized_Parameter.createDimension('maxpft', maxpft)
                NC_File_Out_Optimized_Parameter.createDimension('Station_Dim',numpy.size(Station_XY)/2)
                NC_File_Out_Optimized_Parameter.createDimension('time_soil', None)
                NC_File_Out_Optimized_Parameter.createDimension('time_pft', None)
                
                NC_File_Out_Optimized_Parameter.createVariable('Parameter_Soil_Optimized','f4',('time_soil', 'Ensemble_Number','Dim_Soil_Par','Station_Dim',),zlib=True)
                NC_File_Out_Optimized_Parameter.createVariable('Parameter_PFT_Optimized','f4',('time_pft', 'Ensemble_Number','Dim_PFT_Par','Station_Dim',),zlib=True)
                
                NC_File_Out_Optimized_Parameter.sync()            
                NC_File_Out_Optimized_Parameter.close()
        
        
        Dim_ParFlow_Par = 1
        # Get the Perturbed Parameter Space    
        if Ensemble_Number > 1 or Def_Par_Optimized:
            Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Par_Index_Increment_Soil_Par, Par_Soil_Uniform_STD, Par_Veg_Uniform_STD, Par_PFT_Uniform_STD, Par_Hard_Uniform_STD = \
            Parameter_Space_Function(Model_Driver, Def_Print, Def_PP,active_nodes_server, job_server_node_array, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single, Def_Debug, Def_Region, Def_First_Run, Def_Par_Optimized, 
                                     Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Dim_ParFlow_Par, ParFlow_Layer_Num, Start_Month, maxpft, Soil_Texture_Layer_Opt_Num, Row_Numbers, Col_Numbers, Ensemble_Number,  Ensemble_Number_Predict, \
                                    Parameter_Optimization, Row_Numbers_String, Col_Numbers_String, Soil_Sand_Clay_Sum, Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  PFT_Dominant_Index, topo_slope,\
                                    fsurdat_name, fpftcon_name, NC_FileName_Assimilation_2_Constant, DasPy_Path,
                                    DAS_Data_Path, DAS_Output_Path, Region_Name, Datetime_Start, Datetime_Initial, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, 
                                    DAS_Depends_Path, Def_ParFor, omp_get_num_procs_ParFor, r)
        
        if Def_First_Run == 1:
            Optimized_Parameter_Index = numpy.zeros(4,dtype=numpy.integer)
            Bias_Record_Index = numpy.zeros(2,dtype=numpy.integer)
            Soil_Moisture_Diff_Index = 0
        else:
            Optimized_Parameter_Index = numpy.zeros(4,dtype=numpy.integer)
            Bias_Record_Index = numpy.zeros(2,dtype=numpy.integer)
            
            if Parameter_Optimization:
                NC_File_Out_Optimized_Parameter = netCDF4.Dataset(NC_FileName_Optimized_Parameter, 'r')
                Optimized_Parameter_Index[0] = len(NC_File_Out_Optimized_Parameter.dimensions['time_soil']) - 1
                Optimized_Parameter_Index[2] = len(NC_File_Out_Optimized_Parameter.dimensions['time_pft']) - 1
                NC_File_Out_Optimized_Parameter.close()
                
            if (numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1):
                NC_File_Out_Estimated_Bias = netCDF4.Dataset(NC_FileName_Estimated_Bias, 'r')
                Bias_Record_Index[0] = len(NC_File_Out_Estimated_Bias.dimensions['time']) - 1
                Bias_Record_Index[1] = len(NC_File_Out_Estimated_Bias.dimensions['time']) - 1
                NC_File_Out_Estimated_Bias.close()
                
            NC_File_Out_Soil_Moisture_Difference = netCDF4.Dataset(NC_FileName_Soil_Moisture_Difference, 'r')
            if len(NC_File_Out_Soil_Moisture_Difference.dimensions['time']) >= 1:
                Soil_Moisture_Diff_Index = len(NC_File_Out_Soil_Moisture_Difference.dimensions['time']) - 1
            NC_File_Out_Soil_Moisture_Difference.close()
        
        NC_File_Out_Assimilation_2_Constant = netCDF4.Dataset(NC_FileName_Assimilation_2_Constant, 'r')
        CLM_Soil_Layer_Thickness = numpy.asarray(NC_File_Out_Assimilation_2_Constant.variables['CLM_Soil_Layer_Thickness'][:,:,:])
        NC_File_Out_Assimilation_2_Constant.close()
        
        Analysis_Grid = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Localization_Map_Mask = numpy.zeros((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        ObsModel_Mat_Masked = numpy.zeros((Row_Numbers, Col_Numbers),dtype=numpy.float32)
    
    else:
        LONGXY_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        LATIXY_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        
        Mean_Index_Prop_Grid_Array_Sys = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Model_Variance = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask = numpy.empty((Dim_CLM_State, 3, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask_Index = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers), dtype=numpy.bool)
        Mask_X = numpy.empty((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Mask_Y = numpy.empty((Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        Mean_Dir = None
        Bias_Remove_Start_Time_Array = None
        Observation_Bias_Initialization_Flag = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, Ensemble_Number), dtype=numpy.float32)
        Observation_Bias_Optimized = numpy.empty((Ensemble_Number, Dim_CLM_State, Dim_Observation_Quantity, numpy.size(Station_XY) / 2), dtype=numpy.float32)
        Model_Bias_Range = numpy.empty((Dim_CLM_State, 2), dtype=numpy.float32)
        Observation_Bias_Range = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, 2), dtype=numpy.float32)
        Model_Bias_Range_STD = numpy.empty((Dim_CLM_State, 2), dtype=numpy.float32)
        Observation_Bias_Range_STD = numpy.empty((Dim_CLM_State, Dim_Observation_Quantity, 2), dtype=numpy.float32)
        Model_Bias_STD = numpy.empty(Dim_CLM_State, dtype=numpy.float32)
        Observation_Bias_STD = numpy.empty((Dim_CLM_State,Dim_Observation_Quantity), dtype=numpy.float32)
        Model_State_Inflation_Range = numpy.empty((Dim_CLM_State,2), dtype=numpy.float32)
        Model_State_Inflation_Range_STD = numpy.empty(Dim_CLM_State, dtype=numpy.float32)
        Additive_Noise_SM_Par = numpy.empty((10,11), dtype=numpy.float32)
        Additive_Noise_SM = numpy.empty((Ensemble_Number, Soil_Layer_Num - 5), dtype=numpy.float32)
        Additive_Noise_ST = numpy.empty((Ensemble_Number, 2), dtype=numpy.float32)
        Irrigation_Grid_Flag_Array = None
        cols1d_ixy = numpy.empty(column_len, dtype=numpy.integer)
        cols1d_jxy = numpy.empty(column_len, dtype=numpy.integer)
        cols1d_ityplun = numpy.empty(column_len, dtype=numpy.integer)
        pfts1d_ixy = numpy.empty(pft_len, dtype=numpy.integer)
        pfts1d_jxy = numpy.empty(pft_len, dtype=numpy.integer)
        pfts1d_itypveg = numpy.empty(pft_len, dtype=numpy.integer)
        pfts1d_ci = numpy.empty(pft_len, dtype=numpy.integer)
        pfts1d_ityplun = numpy.empty(pft_len, dtype=numpy.integer)
        
        PCT_PFT_High = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        PCT_PFT_Low = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        PCT_PFT_WATER = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOCVL_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOCVH_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVL_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOTVH_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        ECOWAT_Mat = numpy.empty((Row_Numbers, Col_Numbers), dtype=numpy.float32)
        Soil_Density = numpy.empty((Row_Numbers, Col_Numbers),dtype=numpy.float32)
                                
        CMEM_Work_Path_Array = None
        
        Analysis_Variable_Name = None
        Soil_Sand_Clay_Sum = numpy.empty((Soil_Texture_Layer_Opt_Num, Row_Numbers, Col_Numbers), dtype=numpy.float32)
        
        Parameter_Range_Soil = numpy.empty((2, Dim_Soil_Par),dtype=numpy.float32)
        Parameter_Range_Veg = numpy.empty((2,Dim_Veg_Par),dtype=numpy.float32)
        Parameter_Range_PFT = numpy.empty((2,Dim_PFT_Par),dtype=numpy.float32)
        Parameter_Range_Hard = numpy.empty((2, Dim_Hard_Par),dtype=numpy.float32)
        Par_Index_Increment_Soil_Par = numpy.empty((2, Dim_Soil_Par),dtype=numpy.float32)
        Par_Soil_Uniform_STD = numpy.empty(Dim_Soil_Par,dtype=numpy.float32)
        Par_Veg_Uniform_STD = numpy.empty(Dim_Veg_Par,dtype=numpy.float32)
        Par_PFT_Uniform_STD = numpy.empty(Dim_PFT_Par,dtype=numpy.float32)
        Par_Hard_Uniform_STD = numpy.empty(Dim_Hard_Par,dtype=numpy.float32)
        
        Optimized_Parameter_Index = numpy.empty(4,dtype=numpy.integer)
        Bias_Record_Index = numpy.empty(2,dtype=numpy.integer)
        Soil_Moisture_Diff_Index = None
        
        COSMOS_Circle_Array = None
        COSMOS_Circle_Index_Array = None
        COSMOS_Circle_Num_Array = None
        
        CLM_Soil_Layer_Thickness = numpy.empty((Soil_Layer_Num, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        
        Analysis_Grid = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        Localization_Map_Mask = numpy.empty((Dim_CLM_State, Row_Numbers, Col_Numbers),dtype=numpy.float32)
        ObsModel_Mat_Masked = numpy.empty((Row_Numbers, Col_Numbers),dtype=numpy.float32)
            
    Two_Step_Bias_Estimation_Flag = 0   # whether the bias estimation has been done
    Two_Step_Bias_Estimation_Active = False # whether it is in the process of bias estimation step
    Def_First_Run_Bias = 1
    
    if Def_PP == 2:
        Mean_Dir = mpi4py_comm.bcast(Mean_Dir)
        mpi4py_comm.bcast(Bias_Remove_Start_Time_Array)
        mpi4py_comm.Bcast([Mean_Index_Prop_Grid_Array_Sys,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_Variance,MPI.FLOAT])
        mpi4py_comm.Bcast([Mask,MPI.FLOAT])
        mpi4py_comm.Bcast([Mask_Index,MPI.BOOL])
        
        mpi4py_comm.Bcast([Mask_X,MPI.FLOAT])
        mpi4py_comm.Bcast([Mask_Y,MPI.FLOAT])
        mpi4py_comm.Bcast([LONGXY_Mat,MPI.FLOAT])
        mpi4py_comm.Bcast([LATIXY_Mat,MPI.FLOAT])
        
        mpi4py_comm.Bcast([Observation_Bias_Initialization_Flag,MPI.FLOAT])
        mpi4py_comm.Bcast([Observation_Bias_Optimized,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_Bias_Range,MPI.FLOAT])
        mpi4py_comm.Bcast([Observation_Bias_Range,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_Bias_Range_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Observation_Bias_Range_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_Bias_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Observation_Bias_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_State_Inflation_Range,MPI.FLOAT])
        mpi4py_comm.Bcast([Model_State_Inflation_Range_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Additive_Noise_SM_Par,MPI.FLOAT])
        mpi4py_comm.Bcast([Additive_Noise_SM,MPI.FLOAT])
        mpi4py_comm.Bcast([Additive_Noise_ST,MPI.FLOAT])
        
        Irrigation_Grid_Flag_Array = mpi4py_comm.bcast(Irrigation_Grid_Flag_Array)
        mpi4py_comm.Bcast([cols1d_ixy,MPI.INT])
        mpi4py_comm.Bcast([cols1d_jxy,MPI.INT])
        mpi4py_comm.Bcast([cols1d_ityplun,MPI.INT])
        mpi4py_comm.Bcast([pfts1d_ixy,MPI.INT])
        mpi4py_comm.Bcast([pfts1d_jxy,MPI.INT])
        mpi4py_comm.Bcast([pfts1d_itypveg,MPI.INT])
        mpi4py_comm.Bcast([pfts1d_ci,MPI.INT])
        mpi4py_comm.Bcast([pfts1d_ityplun,MPI.INT])
#         
        mpi4py_comm.Bcast([PCT_PFT_High,MPI.FLOAT])
        mpi4py_comm.Bcast([PCT_PFT_Low,MPI.FLOAT])
        mpi4py_comm.Bcast([PCT_PFT_WATER,MPI.FLOAT])
        mpi4py_comm.Bcast([ECOCVL_Mat,MPI.FLOAT])
        mpi4py_comm.Bcast([ECOCVH_Mat,MPI.FLOAT])
        mpi4py_comm.Bcast([ECOTVL_Mat,MPI.FLOAT])
        mpi4py_comm.Bcast([ECOTVH_Mat,MPI.FLOAT])
        mpi4py_comm.Bcast([ECOWAT_Mat,MPI.FLOAT])
         
        mpi4py_comm.Bcast([Soil_Density,MPI.FLOAT])
        CMEM_Work_Path_Array = mpi4py_comm.bcast(CMEM_Work_Path_Array) 
         
        Analysis_Variable_Name = mpi4py_comm.bcast(Analysis_Variable_Name)
        mpi4py_comm.Bcast([Soil_Sand_Clay_Sum,MPI.FLOAT])
        mpi4py_comm.Bcast([Parameter_Range_Soil,MPI.FLOAT])
        mpi4py_comm.Bcast([Parameter_Range_Veg,MPI.FLOAT])
        mpi4py_comm.Bcast([Parameter_Range_PFT,MPI.FLOAT])
        mpi4py_comm.Bcast([Parameter_Range_Hard,MPI.FLOAT])
        mpi4py_comm.Bcast([Par_Index_Increment_Soil_Par,MPI.FLOAT])
        mpi4py_comm.Bcast([Par_Soil_Uniform_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Par_Veg_Uniform_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Par_PFT_Uniform_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Par_Hard_Uniform_STD,MPI.FLOAT])
        mpi4py_comm.Bcast([Optimized_Parameter_Index,MPI.INT])
        mpi4py_comm.Bcast([Bias_Record_Index,MPI.INT])
        Soil_Moisture_Diff_Index = mpi4py_comm.bcast(Soil_Moisture_Diff_Index)
        mpi4py_comm.Bcast([CLM_Soil_Layer_Thickness,MPI.FLOAT])
        mpi4py_comm.Bcast([Analysis_Grid,MPI.FLOAT])
        mpi4py_comm.Bcast([Localization_Map_Mask,MPI.FLOAT])
        mpi4py_comm.Bcast([ObsModel_Mat_Masked,MPI.FLOAT])
         
        COSMOS_Circle_Array = mpi4py_comm.bcast(COSMOS_Circle_Array)
        COSMOS_Circle_Index_Array = mpi4py_comm.bcast(COSMOS_Circle_Index_Array)
        COSMOS_Circle_Num_Array = mpi4py_comm.bcast(COSMOS_Circle_Num_Array)
        
    if Def_PP == 2:
        mpi4py_comm.barrier()
        mpi4py_comm.Barrier()
    
    #================================================= Do Data Assimilation
    if Do_DA_Flag:
        if mpi4py_rank == 0:
            #---------------------------------------------- Read Observation Time -----------------------------
            Observation_Time_File = open(Observation_Time_File_Path + '/Observation_Time.txt', 'r')
            Observation_Time_File_Header = Observation_Time_File.readline()
            print "Observation_Time_File_Header",Observation_Time_File_Header
            Observation_Time_Lines = Observation_Time_File.readlines()
        else:
            Observation_Time_Lines = None
        
        if Def_PP == 2:
            Observation_Time_Lines = mpi4py_comm.bcast(Observation_Time_Lines)
        
        #print len(Observation_Time_Lines)
        Observation_Index = 0
        Forward = True
        
        while Observation_Index < len(Observation_Time_Lines):
            
            if mpi4py_rank == 0:
                #print Observation_Index,len(Observation_Time_Lines)
                # Find the Same Time Observation
                Variable_Assimilation_Flag = numpy.zeros(Dim_CLM_State)
                
                while Forward:
                    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    print "+++++++++++++++ Find Synchronous Observation ++++++++++++++++++++++++++++++++++++"
                    print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                    Observation_Time_Line = Observation_Time_Lines[Observation_Index]
                    Observation_Time_Line_Split = string.split(Observation_Time_Line)
                    print Observation_Time_Line_Split
                    
                    # Model Stop Time
                    Stop_Year = Observation_Time_Line_Split[6].zfill(4)
                    Stop_Month = Observation_Time_Line_Split[7].zfill(2)
                    Stop_Day = Observation_Time_Line_Split[8].zfill(2)
                    Stop_Hour = Observation_Time_Line_Split[9].zfill(2)
                    Stop_Minute = Observation_Time_Line_Split[10].zfill(2)
                    Datetime_Stop_First = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month), string.atoi(Stop_Day), string.atoi(Stop_Hour), 00)
                    if string.atoi(Stop_Minute) >= 60 or string.atoi(Stop_Minute) < 0:
                        sys.exit('The Observation is Wrong in ' + Stop_Year + ' ' + Stop_Month + ' ' + Stop_Day + ' ' + Stop_Hour + ' ' + Stop_Minute)
                    elif string.atoi(Stop_Minute) <= 30:
                        Datetime_Stop = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month), string.atoi(Stop_Day), string.atoi(Stop_Hour), 00)
                    else:
                        if (string.atoi(Stop_Hour) + 1 < 24):
                            Datetime_Stop = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month), string.atoi(Stop_Day), string.atoi(Stop_Hour) + 1, 00)
                        elif (string.atoi(Stop_Hour) + 1 == 24) and (string.atoi(Stop_Day) + 1 <= Num_of_Days_Monthly[string.atoi(Stop_Month) - 1]):
                            Datetime_Stop = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month), string.atoi(Stop_Day) + 1, 00, 00)
                        elif (string.atoi(Stop_Hour) + 1 == 24) and (string.atoi(Stop_Day) + 1 > Num_of_Days_Monthly[string.atoi(Stop_Month) - 1])  and (string.atoi(Stop_Month) + 1 <= 12):
                            Datetime_Stop = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month) + 1, 01, 00, 00)
                        elif (string.atoi(Stop_Hour) + 1 == 24) and (string.atoi(Stop_Day) + 1 > Num_of_Days_Monthly[string.atoi(Stop_Month) - 1])  and (string.atoi(Stop_Month) + 1 > 12):
                            Datetime_Stop = datetime.datetime(string.atoi(Stop_Year) + 1, 01, 01, 00, 00)
                    #Datetime_Stop_Init = datetime.datetime(Datetime_Stop.year-1,12,31,23,00)
                    Datetime_Stop_Init = datetime.datetime(Datetime_Stop.year, Datetime_Stop.month, Datetime_Stop.day, 00, 00)
                    
                    print "Datetime_Start, Datetime_Stop, Datetime_Stop_Init, (Datetime_Stop - Datetime_Stop_Init).seconds"  
                    print Datetime_Start, Datetime_Stop, Datetime_Stop_Init, (Datetime_Stop - Datetime_Stop_Init).seconds   
                    
                    # Because there is a calendar bug in cesm1.1.1, so we skip 02-29 for data assimilation
                    if (calendar.isleap(string.atoi(Stop_Year)) and (datetime.datetime(string.atoi(Stop_Year),string.atoi(Stop_Month),string.atoi(Stop_Day)) == datetime.datetime(string.atoi(Stop_Year),2,29))) or (Datetime_Stop <= Datetime_Start):
                        Observation_Index = Observation_Index + 1
                        if Observation_Index >= len(Observation_Time_Lines):
                            Forward = False
                            break
                        else:
                            continue
                    
                    print "Observation_Index",Observation_Index,"len(Observation_Time_Lines)",len(Observation_Time_Lines)
                    #sys.exit()
                    # Judge the Sensor Type
                    SensorType.append(Observation_Time_Line_Split[0]) # MODIS AMSR-E SMOS ASCAT ASAR
                    SensorVariable.append(Observation_Time_Line_Split[1])
                    SensorQuantity.append(Observation_Time_Line_Split[2])
                    Variable_ID.append(Observation_Time_Line_Split[3])
                    QC_ID.append(Observation_Time_Line_Split[4])
                    SensorResolution.append(string.atof(Observation_Time_Line_Split[5]))
                    Observation_File_Name.append(Observation_Time_Line_Split[11])
                    
                    SensorVariable_Temp = Observation_Time_Line_Split[1]
                    DA_Flag = 0 # if DA_Flag=1, Do Data Assimilation
                    
                    print SensorVariable_Temp,"Variable_List.index(SensorVariable_Temp)",Variable_List.index(SensorVariable_Temp)
                    
                    if SensorVariable_Temp == 'Albedo':
                        Variable_Assimilation_Flag[4] = 1
                        Variable_Assimilation_Flag[5] = 1
                        Variable_Assimilation_Flag[6] = 1
                        Variable_Assimilation_Flag[7] = 1
                    elif Variable_List.index(SensorVariable_Temp) >= 0:
                        Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Temp)] = 1
                    else:
                        print "Wrong SensorVariable_Temp Specification",SensorVariable_Temp
                    
                    # Recore the last bias remove time of each observation type
                    if Def_First_Run_Bias or (Bias_Remove_Start_Time_Array[Variable_List.index(SensorVariable_Temp)] == ''):
                        Bias_Remove_Start_Time_Array[Variable_List.index(SensorVariable_Temp)] = Datetime_Start
                    
                    print ""
                    print "-----Bias_Remove_Start_Time_Array",Bias_Remove_Start_Time_Array
                    print ""
                    
                    if Def_Print:
                        print "------------Variable_Assimilation_Flag",Variable_Assimilation_Flag
                        
                    if numpy.sum(Variable_Assimilation_Flag) > 0:
                        DA_Flag = 1                    
                    
                    if Def_Print:
                        print "len(SensorType),len(SensorVariable),len(SensorQuantity),len(Variable_ID)",len(SensorType),len(SensorVariable),len(SensorQuantity),len(Variable_ID)                
                    
                    if Observation_Index + 1 >= len(Observation_Time_Lines):
                        Forward = False
                        break
                    else:
                        Observation_Index = Observation_Index + 1
                        Observation_Time_Line = Observation_Time_Lines[Observation_Index]
                        Observation_Time_Line_Split = string.split(Observation_Time_Line)
                        #print Observation_Time_Line_Split
                        Stop_Year = Observation_Time_Line_Split[6].zfill(4)
                        Stop_Month = Observation_Time_Line_Split[7].zfill(2)
                        Stop_Day = Observation_Time_Line_Split[8].zfill(2)
                        Stop_Hour = Observation_Time_Line_Split[9].zfill(2)
                        Stop_Minute = Observation_Time_Line_Split[10].zfill(2)
                        Datetime_Stop_Second = datetime.datetime(string.atoi(Stop_Year), string.atoi(Stop_Month), string.atoi(Stop_Day), string.atoi(Stop_Hour), 00)
                        #print Datetime_Stop_First,Datetime_Stop_Second
                        if Datetime_Stop_First == Datetime_Stop_Second:
                            Forward = True
                        else:
                            Observation_Index = Observation_Index - 1
                            Observation_Time_Line = Observation_Time_Lines[Observation_Index]
                            Observation_Time_Line_Split = string.split(Observation_Time_Line)
                            #print Observation_Time_Line_Split
                            Stop_Year = Observation_Time_Line_Split[6].zfill(4)
                            Stop_Month = Observation_Time_Line_Split[7].zfill(2)
                            Stop_Day = Observation_Time_Line_Split[8].zfill(2)
                            Stop_Hour = Observation_Time_Line_Split[9].zfill(2)
                            Stop_Minute = Observation_Time_Line_Split[10].zfill(2)
                            break   
            
            else:
                Variable_Assimilation_Flag = None
                Stop_Year = None
                Stop_Month = None
                Stop_Day = None
                Stop_Hour = None
                Stop_Minute = None
                Datetime_Stop = None
                Datetime_Stop_Init = None
                Observation_Index = None
                SensorType = None
                SensorVariable = None
                SensorQuantity = None
                Variable_ID = None
                QC_ID = None
                SensorResolution = None
                Observation_File_Name = None
                
                SensorVariable_Temp = None
                DA_Flag = None
            
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
                Variable_Assimilation_Flag = mpi4py_comm.bcast(Variable_Assimilation_Flag)
                Stop_Year = mpi4py_comm.bcast(Stop_Year)
                Stop_Month = mpi4py_comm.bcast(Stop_Month)
                Stop_Day = mpi4py_comm.bcast(Stop_Day)
                Stop_Hour = mpi4py_comm.bcast(Stop_Hour)
                Stop_Minute = mpi4py_comm.bcast(Stop_Minute)
                Datetime_Stop = mpi4py_comm.bcast(Datetime_Stop)
                Datetime_Stop_Init = mpi4py_comm.bcast(Datetime_Stop_Init)
                Observation_Index = mpi4py_comm.bcast(Observation_Index)
                SensorType = mpi4py_comm.bcast(SensorType)
                SensorVariable = mpi4py_comm.bcast(SensorVariable)
                SensorQuantity = mpi4py_comm.bcast(SensorQuantity)
                Variable_ID = mpi4py_comm.bcast(Variable_ID)
                QC_ID = mpi4py_comm.bcast(QC_ID)
                SensorResolution = mpi4py_comm.bcast(SensorResolution)
                Observation_File_Name = mpi4py_comm.bcast(Observation_File_Name)
            
            Dim_Obs_Type = len(SensorType)
            if mpi4py_rank == 0:
                print "******************************************* Dim_Obs_Type",Dim_Obs_Type
            
            Month_String = "-" + Stop_Month.zfill(2)
            Day_String = "-" + Stop_Day.zfill(2)
            Hour_String = "-" + Stop_Hour.zfill(2)
            DateString_Plot=Stop_Year+ Month_String+Day_String+Hour_String
            
            if mpi4py_rank == 0:
                print "#######################Initial file is",finidat_initial_CLM
            stop_tod_string = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
            history_file_name = Region_Name + '.clm2.h0.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
            finidat_name = Region_Name + '.clm2.r.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
            
            if mpi4py_rank == 0:
                print "*************************************************Prepare Observation Matrix***************************************************************"
                Observation_Matrix = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_Variance = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_Latitude = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_Longitude = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_View_Zenith_Angle = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_View_Time = numpy.zeros((len(SensorType), Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_NLons = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_NLats = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_X_Left = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_X_Right = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_Y_Lower = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_Y_Upper = numpy.zeros(len(SensorType),dtype=numpy.float32)
                Observation_Misc = numpy.zeros((len(SensorType), 10, Row_Numbers, Col_Numbers),dtype=numpy.float32)
                Observation_Corelation_Par = numpy.zeros((len(SensorType), 5, 2),dtype=numpy.float32)   
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
                
            if Def_PP == 2:
                mpi4py_comm.Bcast([Observation_Matrix,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Variance,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Latitude,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Longitude,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_View_Zenith_Angle,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_View_Time,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_NLons,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_NLats,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_X_Left,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_X_Right,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Y_Lower,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Y_Upper,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Misc,MPI.FLOAT])
                mpi4py_comm.Bcast([Observation_Corelation_Par,MPI.FLOAT])
            
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
                
            if Datetime_Stop > Datetime_Start and Datetime_Stop <= Datetime_End and Dim_Obs_Type > 0:
                if mpi4py_rank == 0:
                    print "Datetime_Start, Datetime_Stop",Datetime_Start, Datetime_Stop   
                
                os.chdir(DasPy_Path)
                
                if Datetime_Start != Datetime_Stop: # If There are severl Observation at the Same Time, CLM Only Need to Be Run Once.
                    if mpi4py_rank == 0:
                        print "*************************************************Start Online Forcing Perturbation***************************************************************"
#                         Forcing_Perturbation_Online(Def_PP,Ensemble_Number,DasPy_Path,Forcing_File_Path_Home,Row_Numbers,Col_Numbers,Grid_Resolution_GEO,
#                                                 mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, LATIXY_Mat, LONGXY_Mat, Forcepert_ntrmdt,
#                                                 Datetime_Start, Datetime_Stop, active_nodes_server, job_server_node_array, Def_Print)
        
                if (numpy.size(numpy.where(Bias_Estimation_Option_Model == 1)) >= 1) or (numpy.size(numpy.where(Bias_Estimation_Option_Obs == 1)) >= 1):  # Model/Observation Bias Estimation
                    if Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")] or Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")] or Variable_Assimilation_Flag[Variable_List.index("Vegetation_Temperature")] or Variable_Assimilation_Flag[Variable_List.index("Latent_Heat")]\
                     or Variable_Assimilation_Flag[Variable_List.index("Latent_Heat_Daily")] or Variable_Assimilation_Flag[Variable_List.index("Sensible_Heat")]:
                        Two_Step_Bias_Estimation_Active = True  #If it is true, read the ensemble mean in the read_history function, use the bias ensembles to generate the model ensembles
                        # In theory, only the parameter estimation or bias estimation can be activated, both can not be activated at the same time
                
                if mpi4py_rank == 0:
                    print "===================================== Call DAS_Driver_Common"
                Observation_Matrix, Observation_Longitude, Observation_Latitude, Observation_Variance, Observation_NLons, Observation_NLats, \
                Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Misc, Observation_View_Zenith_Angle, Observation_View_Time, Observation_Corelation_Par, \
                Analysis_Variable_Name, Constant_File_Name, Def_Par_Optimized, Soil_Layer_Index_DA,\
                Mask, Mask_Index, Model_Variance, Def_First_Run_RTM, ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat, \
                Mean_Index_Prop_Grid_Array_Sys, Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Observation_Bias_Initialization_Flag, job_server_node_array, active_nodes_server = \
                DAS_Driver_Common(mpi4py_comm, mpi4py_null, mpi4py_rank, mpi4py_name, mpi4py_comm_split, mpipy_comm_decomposition, Model_Driver, PDAF_Assim_Framework, PDAF_Filter_Type, Def_PP, Def_CESM_Multi_Instance, Def_Par_Optimized, Def_Par_Sensitivity, Def_Par_Correlation, Do_DA_Flag, CLM_NA, NAvalue, finidat_initial_CLM, finidat_initial_CLM_Copy, Def_ParFor, Def_Region, 
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
                                Initial_Perturbation_ST_Flag, Def_First_Run, Def_SpinUp, active_nodes_server, job_server_node_array, PP_Servers_Per_Node,
                                Forcing_File_Path_Home, SensorType, SensorVariable, SensorQuantity, SensorResolution, Variable_ID, QC_ID, Observation_File_Name, Dim_Obs_Type,
                                Write_DA_File_Flag, Use_Mask_Flag, Mask_File, Def_ReBEL, Def_Localization, Call_Gstat_Flag, Plot_Analysis, Num_Local_Obs_State,
                                Observation_Path,  Grid_Resolution_CEA, Grid_Resolution_GEO, mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen,
                                LATIXY_Mat, LONGXY_Mat, MODEL_CEA_X, MODEL_CEA_Y, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower,MODEL_Y_Upper, LAI_Year_String, Month_String, 
                                DAS_Output_Path, Dim_CLM_State, Parameter_Optimization, Def_Debug, PCT_PFT_High, PCT_PFT_Low, PCT_PFT_WATER, Soil_Layer_Index_DA,
                                ECOCVL_Mat, ECOCVH_Mat, ECOTVL_Mat, ECOTVH_Mat, ECOWAT_Mat, column_len, pft_len, cols1d_jxy, cols1d_ixy, pfts1d_jxy, pfts1d_ixy, 
                                Constant_File_Name_Header, finidat_name_string, Feedback_Assim, numrad,  
                                Density_of_ice, omp_get_num_procs_ParFor, Model_Variance, Additive_Noise_SM_Par, Additive_Noise_SM, Additive_Noise_ST,
                                Variable_List, Variable_Assimilation_Flag, Analysis_Variable_Name, Mask, Mask_X, Mask_Y, N0, nlyr, Station_XY, Station_XY_Index, Def_First_Run_RTM, Soil_Density, CMEM_Work_Path_Array,
                                DAS_File_Name_List, COUP_OAS_PFL, CESM_Init_Flag,
                                MODIS_LAI_Data_ID, Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, PFT_Par_Sens_Array, UTC_Zone, plt, cm, colors, DateString_Plot, octave, r, COSMIC_Py, [], memory_profiler, COSMIC, Observation_Time_File_Path)
                
                if mpi4py_rank == 0:
                    ###################################################################
                    Weather_Forecast_Days_State = copy.copy(Weather_Forecast_Days)
                    Weather_Forecast_Days_Par = 0
                    
                    
                    # Record how many iterations we have for each parameter
                    Soil_Par_Accum_Dim = 0  
                    Veg_Par_Accum_Dim = 0
                    PFT_Par_Accum_Dim = 0
                    Hard_Par_Accum_Dim = 0
                    
                    
                    for Observation_Matrix_Index in range(Dim_Obs_Type):
                        print "Read",str(Observation_Matrix_Index+1)+"th","Observation Matrix!"
                        SensorType_Sub = SensorType[Observation_Matrix_Index]
                        SensorVariable_Sub = SensorVariable[Observation_Matrix_Index]
                        SensorQuantity_Sub = SensorQuantity[Observation_Matrix_Index]
                        SensorResolution_Sub = SensorResolution[Observation_Matrix_Index]
                        Variable_ID_Sub = Variable_ID[Observation_Matrix_Index]
                        QC_ID_Sub = QC_ID[Observation_Matrix_Index]
                        
                        print "SensorType_Sub,SensorVariable_Sub,SensorQuantity_Sub,SensorResolution_Sub,Variable_ID_Sub,QC_ID_Sub"
                        print SensorType_Sub,SensorVariable_Sub,SensorQuantity_Sub,SensorResolution_Sub,Variable_ID_Sub,QC_ID_Sub
                        
                        if SensorVariable_Sub == "Irrigation_Scheduling":
                            print "Skip to next observation because of irrigation"
                            continue
                        
                        if SensorVariable_Sub != "Albedo":
                            Prop_Grid_Array_Sys_Index = Variable_List.index(SensorVariable_Sub)
                        else:
                            Prop_Grid_Array_Sys_Index = Variable_List.index(Variable_ID_Sub)
                        
                        print "SensorVariable_Sub,Prop_Grid_Array_Sys_Index",SensorVariable_Sub,Prop_Grid_Array_Sys_Index
                        print ""
                        
                        if SensorQuantity_Sub == "K":
                            if SensorVariable_Sub == "Soil_Moisture":
                                SensorQuantity_Index = 0
                            if SensorVariable_Sub == "Surface_Temperature" or SensorVariable_Sub == "Vegetation_Temperature":
                                SensorQuantity_Index = 3
                        if SensorQuantity_Sub == "DB":
                            SensorQuantity_Index = 1
                        if SensorQuantity_Sub == "Neutron_Count":
                            SensorQuantity_Index = 2
                        if SensorQuantity_Sub == "m3/m3":
                            SensorQuantity_Index = 3
                        
                        print "SensorQuantity_Sub,SensorQuantity_Index",SensorQuantity_Sub,SensorQuantity_Index
                        print ""
                        
                        print "Mask the Observation where the Model State is invalid"
                        Observation_Matrix[Observation_Matrix_Index,:,:][Mask_Index[Prop_Grid_Array_Sys_Index,::]] = -9999.0
                        
                        
                        if SensorType_Sub == "COSMOS" or SensorType_Sub == "InSitu":
                            Def_Localization_Original = Def_Localization
                            Def_Localization = 1.0
                            
                        if Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")]:
                            Soil_Par_Sens = Soil_Par_Sens_Array[0]
                            Veg_Par_Sens = Veg_Par_Sens_Array[0]
                            PFT_Par_Sens = PFT_Par_Sens_Array[0]
                            Hard_Par_Sens = Hard_Par_Sens_Array[0]
                        elif Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")]:
                            Soil_Par_Sens = Soil_Par_Sens_Array[1]
                            Veg_Par_Sens = Veg_Par_Sens_Array[1]
                            PFT_Par_Sens = PFT_Par_Sens_Array[1]
                            Hard_Par_Sens = Hard_Par_Sens_Array[1]
                        else:
                            Soil_Par_Sens = numpy.array([False, False, False, False, False],dtype=numpy.bool)
                            Veg_Par_Sens = numpy.array([False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],dtype=numpy.bool)
                            PFT_Par_Sens = numpy.array([False, False, False],dtype=numpy.bool)
                            Hard_Par_Sens = numpy.array([False, False],dtype=numpy.bool)
                            
                        Soil_Par_Sens_Dim = numpy.size(numpy.where(Soil_Par_Sens == True))
                        Veg_Par_Sens_Dim = numpy.size(numpy.where(Veg_Par_Sens == True))
                        PFT_Par_Sens_Dim = numpy.size(numpy.where(PFT_Par_Sens == True))
                        Hard_Par_Sens_Dim = numpy.size(numpy.where(Hard_Par_Sens == True))
                        
                        print "++++++++++++++++++++++ Soil_Par_Sens_Dim",Soil_Par_Sens_Dim,"+++++++++++++++++++++++ PFT_Par_Sens_Dim",PFT_Par_Sens_Dim
                                
                        Par_Soil_Uniform_STD_Sub = Par_Soil_Uniform_STD[Soil_Par_Sens]
                        Par_PFT_Uniform_STD_Sub = Par_PFT_Uniform_STD[PFT_Par_Sens]
                        
                        print "Par_Soil_Uniform_STD_Sub",Par_Soil_Uniform_STD_Sub,"Par_PFT_Uniform_STD_Sub",Par_PFT_Uniform_STD_Sub
                        ##################################################################################################################   
                        
                        NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
                        Mask_Index_Vector = ~Mask_Index[Prop_Grid_Array_Sys_Index,:,:].flatten()
                        Model_Variance_Sub = Model_Variance[Prop_Grid_Array_Sys_Index, :, :]
                        Prop_Grid_Array_Sys_Sub = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :]
                        Prop_Grid_Array_H_Trans_Sub = NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_H_Trans'][:, Prop_Grid_Array_Sys_Index, :, :]
                        
                        Mask_Sub = Mask[Prop_Grid_Array_Sys_Index, :, :, :]
                        Mask_Index_Sub = Mask_Index[Prop_Grid_Array_Sys_Index, :, :]
                                            
                        E0_SysModel_Mask = numpy.zeros((numpy.size(numpy.where(Mask_Index_Vector==True)),Ensemble_Number),dtype=numpy.float32)
                        E0_ObsModel_Mask = numpy.zeros((numpy.size(numpy.where(Mask_Index_Vector==True)),Ensemble_Number),dtype=numpy.float32)
                        for Ens_Index in range(Ensemble_Number):
                            E0_SysModel_Mask[:,Ens_Index] = Prop_Grid_Array_Sys_Sub[Ens_Index, :,:][~Mask_Index_Sub].flatten()
                            E0_ObsModel_Mask[:,Ens_Index] = Prop_Grid_Array_H_Trans_Sub[Ens_Index, :,:][~Mask_Index_Sub].flatten()
                        SysModel_Variance_Value = numpy.var(E0_SysModel_Mask,axis=1)
                        Mean_SysModel_Variance_Value = numpy.mean(SysModel_Variance_Value)
                        ObsModel_Variance_Value = numpy.var(E0_ObsModel_Mask,axis=1)
                        Mean_ObsModel_Variance_Value = numpy.mean(ObsModel_Variance_Value)
                        
                        for Soil_Layer_Index in range(Soil_Layer_Num):
                            Mean_SysModel_Variance_Value_Layer = numpy.mean(numpy.var(E0_SysModel_Mask,axis=1))
                            
                            if SensorVariable_Sub == "Soil_Moisture":
                                if numpy.sqrt(Mean_SysModel_Variance_Value_Layer) < 0.02:
                                    Initial_Perturbation_SM_Flag[Soil_Layer_Index] = 1
                                else:
                                    Initial_Perturbation_SM_Flag[Soil_Layer_Index] = 0
                            #print Initial_Perturbation_SM_Flag[0]
                            elif SensorVariable_Sub == "Surface_Temperature":
                                if numpy.sqrt(Mean_SysModel_Variance_Value_Layer) < 1.0:
                                    Initial_Perturbation_ST_Flag[Soil_Layer_Index] = 1
                                else:
                                    Initial_Perturbation_ST_Flag[Soil_Layer_Index] = 0
                                        
                        if Write_DA_File_Flag:
                            numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_SysModel_Mask.txt", E0_SysModel_Mask)
                            numpy.savetxt(DasPy_Path+"Analysis/DAS_Temp/E0_ObsModel_Mask.txt", E0_ObsModel_Mask)
                            
                        if Def_Print:
                            print "******************************************************************************"
                            print "SysModel Mean Variance is:", Mean_SysModel_Variance_Value
                            print "Max SysModel Variance is:", numpy.max(SysModel_Variance_Value),"Median SysModel Variance is:", numpy.median(SysModel_Variance_Value),"Min SysModel Variance is:", numpy.min(SysModel_Variance_Value)
                            print "ObsModel Mean Variance is:", Mean_ObsModel_Variance_Value
                            print "Max ObsModel Variance is:", numpy.max(ObsModel_Variance_Value),"Median ObsModel Variance is:", numpy.median(ObsModel_Variance_Value),"Min ObsModel Variance is:", numpy.min(ObsModel_Variance_Value)
                            print "******************************************************************************"
                        
                        SysModel_Mat = numpy.zeros((Row_Numbers,Col_Numbers))
                        SysModel_Mat_Col = SysModel_Mat.flatten()
                        SysModel_Mat_Col[Mask_Index_Vector] = numpy.mean(E0_SysModel_Mask,axis=1)
                        SysModel_Mat = numpy.reshape(SysModel_Mat_Col, (Row_Numbers, -1))
                        SysModel_Mat = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index,:,:], SysModel_Mat)
                        
                        ObsModel_Variance_Mat = numpy.zeros((Row_Numbers,Col_Numbers))
                        ObsModel_Variance_Mat_Col = ObsModel_Variance_Mat.flatten()
                        ObsModel_Variance_Mat_Col[Mask_Index_Vector] = ObsModel_Variance_Value
                        ObsModel_Variance_Mat = numpy.reshape(ObsModel_Variance_Mat_Col, (Row_Numbers, -1))
                        ObsModel_Variance_Mat = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index,:,:], ObsModel_Variance_Mat)
                            
                        ObsModel_Mat = numpy.zeros((Row_Numbers,Col_Numbers))
                        ObsModel_Mat_Col = ObsModel_Mat.flatten()
                        ObsModel_Mat_Col[Mask_Index_Vector] = numpy.mean(E0_ObsModel_Mask,axis=1)
                        ObsModel_Mat = numpy.reshape(ObsModel_Mat_Col, (Row_Numbers, -1))
                        ObsModel_Mat = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index,:,:], ObsModel_Mat)
                            
                        if Plot_Analysis:
                            import matplotlib.pyplot as plt
                            import matplotlib.cm as cm
                            import matplotlib.colors as colors
                            from mpl_toolkits.axes_grid.inset_locator import inset_axes
                            
                            w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)
                            
                            SysModel_Mat_Masked = numpy.ma.masked_values(SysModel_Mat, NAvalue)
                            ObsModel_Mat_Masked = numpy.ma.masked_values(ObsModel_Mat, NAvalue)
                            Observation_Matrix_Masked = numpy.ma.masked_values(Observation_Matrix[Observation_Matrix_Index,:,:],NAvalue)
                            ObsModel_Variance_Mat_Masked = numpy.ma.masked_values(ObsModel_Variance_Mat,NAvalue)
                            
                            Variable_Min = numpy.zeros(Dim_Obs_Type)
                            Variable_Max = numpy.zeros(Dim_Obs_Type)
                            
                            Variable_Min[Observation_Matrix_Index] = numpy.min(Observation_Matrix_Masked)
                            Variable_Max[Observation_Matrix_Index] = numpy.max(Observation_Matrix_Masked)
                            print "Variable_Min_Obs",Variable_Min[Observation_Matrix_Index],"Variable_Max_Obs",Variable_Max[Observation_Matrix_Index]
        
                            
                            Variable_Min = numpy.min(ObsModel_Variance_Mat_Masked)
                            Variable_Max = numpy.max(ObsModel_Variance_Mat_Masked)
                                
                            if Variable_Min != Variable_Max and (not numpy.isnan(Variable_Min)) and (not numpy.isnan(Variable_Max)):
                                
                                fig1 = plt.figure(figsize=(w*2, h*2), dpi=80)
                                fig1.suptitle(DateString_Plot, fontsize=16)
                                
                                Variable_Min = numpy.min(SysModel_Mat_Masked)
                                Variable_Max = numpy.max(SysModel_Mat_Masked)
                                ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
                                color_boun_list = []
                                color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
                                for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                                    color_bound[0] += color_bound[2]
                                    color_boun_list.append(color_bound[0])
                                
                                ax = fig1.add_subplot(2, 2, 1)
                                im1 = ax.imshow(SysModel_Mat_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                                plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                                ax.set_title('SysModel_Value')
                                plt.grid(True)
                                
                                Variable_Min = numpy.min(ObsModel_Mat_Masked)
                                Variable_Max = numpy.max(ObsModel_Mat_Masked)
                                ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
                                color_boun_list = []
                                color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
                                for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                                    color_bound[0] += color_bound[2]
                                    color_boun_list.append(color_bound[0])
                                
                                ax = fig1.add_subplot(2, 2, 2)
                                im1 = ax.imshow(ObsModel_Mat_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                                plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                                ax.set_title('ObsModel_Value')
                                plt.grid(True)
                                
                                ax = fig1.add_subplot(2, 2, 3)
                                im1 = ax.imshow(Observation_Matrix_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                                plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                                ax.set_title('Observation_Matrix')
                                plt.grid(True)
                                
                                Variable_Min = numpy.min(ObsModel_Variance_Mat_Masked)
                                Variable_Max = numpy.max(ObsModel_Variance_Mat_Masked)
                            
                                print "Variable_Max,Variable_Min",Variable_Max,Variable_Min
                                ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
                                color_boun_list = []
                                color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
                                for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                                    color_bound[0] += color_bound[2]
                                    color_boun_list.append(color_bound[0])
                                
                                ax = fig1.add_subplot(2, 2, 4)
                                im1 = ax.imshow(ObsModel_Variance_Mat_Masked, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                                plt.colorbar(im1, ticks=ticks, orientation='horizontal')
                                ax.set_title('ObsModel_Variance_Value')
                                plt.grid(True)
                                
                                plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/SysModel_ObsModel_Observation_"+str(Observation_Matrix_Index)+".png")
                                plt.close('all')
                                #os.abort()
                                
                        if Def_Print >= 2:
                            print "numpy.mean(Prop_Grid_Array_Sys[:, Prop_Grid_Array_Sys_Index, :, :],axis=0)",numpy.shape(numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0))
                        Model_State = numpy.mean(NC_File_Out_Assimilation_2_Initial.variables['Prop_Grid_Array_Sys'][:, Prop_Grid_Array_Sys_Index, :, :],axis=0)[~Mask_Index[Prop_Grid_Array_Sys_Index,::]].flatten()
                        
                        NC_File_Out_Assimilation_2_Initial.close()
                        
                        print "******************** Prepare the Input For Block_Assim"
                        NC_FileName_Block_Assim_Common = DAS_Output_Path+"Analysis/"+Region_Name+"/Block_Assim_Common.nc"
                        if os.path.exists(NC_FileName_Block_Assim_Common):
                            os.remove(NC_FileName_Block_Assim_Common)
                            
                        if Def_Print:
                            print 'Write NetCDF File:',NC_FileName_Block_Assim_Common
                            
                        NC_File_Block_Assim_Common = netCDF4.Dataset(NC_FileName_Block_Assim_Common, 'w', diskless=True, persist=True, format='NETCDF4')
                        # Dim the dimensions of NetCDF
                        NC_File_Block_Assim_Common.createDimension('Ensemble_Number', Ensemble_Number)
                        NC_File_Block_Assim_Common.createDimension('lon', Col_Numbers)
                        NC_File_Block_Assim_Common.createDimension('lat', Row_Numbers)
                        NC_File_Block_Assim_Common.createDimension('Dim_CLM_State', Dim_CLM_State)
                        NC_File_Block_Assim_Common.createDimension('Dim_Obs_Type', Dim_Obs_Type)
                        NC_File_Block_Assim_Common.createDimension('Mask_Dim', 3)
                        
                        NC_File_Block_Assim_Common.createVariable('Mask_Sub','f4',('Mask_Dim','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Mask_Sub'][:,:,:] = Mask_Sub
                        
                        NC_File_Block_Assim_Common.createVariable('Mask_Index_Sub','i4',('lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Mask_Index_Sub'][:,:] = Mask_Index_Sub
                        
                        NC_File_Block_Assim_Common.createVariable('Model_Variance','f4',('Dim_CLM_State','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Model_Variance'][:,:,:] = Model_Variance
                        
                        NC_File_Block_Assim_Common.createVariable('Observation_Variance','f4',('Dim_Obs_Type','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Observation_Variance'][:,:,:] = Observation_Variance
                        
                        NC_File_Block_Assim_Common.createVariable('Observation_Latitude','f4',('Dim_Obs_Type','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Observation_Latitude'][:,:,:] = Observation_Latitude
                        
                        NC_File_Block_Assim_Common.createVariable('Observation_Longitude','f4',('Dim_Obs_Type','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Observation_Longitude'][:,:,:] = Observation_Longitude
                        
                        NC_File_Block_Assim_Common.createVariable('Observation_Matrix','f4',('Dim_Obs_Type','lat','lon',),zlib=True)
                        NC_File_Block_Assim_Common.variables['Observation_Matrix'][:,:,:] = Observation_Matrix
                        
                        NC_File_Block_Assim_Common.sync()      
                        NC_File_Block_Assim_Common.close()
        
                                       
                        if Parameter_Optimization:                   
                             
                            if Parameter_Optimization_First_Flag and Def_First_Run:
                                NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
                                for Station_Index in range(numpy.size(Station_XY)/2):
                                    Parameter_Soil_Optimized[:,:,Station_Index] = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                                    Parameter_PFT_Optimized[:,:,Station_Index] = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                                
                                NC_File_Out_Optimized_Parameter = netCDF4.Dataset(NC_FileName_Optimized_Parameter, 'a')
                                NC_File_Out_Optimized_Parameter.variables['Parameter_Soil_Optimized'][Optimized_Parameter_Index[0],:,:,:] = Parameter_Soil_Optimized
                                NC_File_Out_Optimized_Parameter.variables['Parameter_PFT_Optimized'][Optimized_Parameter_Index[2],:,:,:] = Parameter_PFT_Optimized
                                NC_File_Out_Optimized_Parameter.sync()
                                NC_File_Out_Optimized_Parameter.close()
                                
                                NC_File_Out_Assimilation_2_Parameter.close()
                                
                            Parameter_Optimization_First_Flag = False
                             
                            print "****************************Optimizing the Parameters When Soil_Moisture_DA_Flag or Surface_Temperature_DA_Flag or Vegetation_Temperature_DA_Flag is True!!***********8"
                             
                            if Variable_Assimilation_Flag[Variable_List.index("Soil_Moisture")]:
                                Soil_Par_Sens = Soil_Par_Sens_Array[0]
                                Veg_Par_Sens = Veg_Par_Sens_Array[0]
                                PFT_Par_Sens = PFT_Par_Sens_Array[0]
                                Hard_Par_Sens = Hard_Par_Sens_Array[0]
                            elif Variable_Assimilation_Flag[Variable_List.index("Surface_Temperature")]:
                                Soil_Par_Sens = Soil_Par_Sens_Array[1]
                                Veg_Par_Sens = Veg_Par_Sens_Array[1]
                                PFT_Par_Sens = PFT_Par_Sens_Array[1]
                                Hard_Par_Sens = Hard_Par_Sens_Array[1]
                            else:
                                Soil_Par_Sens = numpy.array([False, False, False, False, False],dtype=numpy.bool)
                                Veg_Par_Sens = numpy.array([False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],dtype=numpy.bool)
                                PFT_Par_Sens = numpy.array([False, False, False],dtype=numpy.bool)
                                Hard_Par_Sens = numpy.array([False, False],dtype=numpy.bool)
    
                            Soil_Par_Sens_Dim = numpy.size(numpy.where(Soil_Par_Sens == True))
                            Veg_Par_Sens_Dim = numpy.size(numpy.where(Veg_Par_Sens == True))
                            PFT_Par_Sens_Dim = numpy.size(numpy.where(PFT_Par_Sens == True))
                            Hard_Par_Sens_Dim = numpy.size(numpy.where(Hard_Par_Sens == True))
                             
                            if Soil_Par_Sens_Dim >= 1 or PFT_Par_Sens_Dim >= 1:                                
                                
                                print "*************************************************Start Parameter Optimization***************************************************************"
                                
                                Normal_Score_Trans_Par = 0  # No Normal Scaore Transformation for Parameter Estimation
                                
                                Def_Par_Optimized, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, job_server_node_array, active_nodes_server, Optimized_Parameter_Index = \
                                Parameter_Update(mpi4py_comm, mpi4py_rank, mpi4py_name, gelmna_threshold, Optimized_Parameter_Index, Model_Driver, NSLOTS, Def_PP, Def_First_Run, Def_Print, Feedback_Assim, Def_Par_Optimized, Parameter_Optimization, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, [], Par_PFT_Uniform_STD_Sub, [], Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, \
                                                 SensorQuantity_Sub, SensorType_Sub, SensorVariable_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub, Variable_List, maxpft, \
                                              Row_Numbers, Col_Numbers, Ensemble_Number, Ensemble_Number_Predict, Dim_Obs_Type, Observation_Matrix, Observation_Longitude, Observation_Latitude, job_server_node_array, active_nodes_server, ntasks_CLM, \
                                              Mask, Mask_Index, NAvalue, COSMOS_Circle_Array, COSMOS_Circle_Index_Array, COSMOS_Circle_Num_Array, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Par_Index_Increment_Soil_Par, DasPy_Path, \
                                              Variable_Assimilation_Flag, DAS_Depends_Path, Def_ParFor, omp_get_num_procs_ParFor, Def_CDF_Matching, Normal_Score_Trans_Par, PDAF_Assim_Framework, PDAF_Filter_Type, PP_Servers_Per_Node, Def_CESM_Multi_Instance, PP_Port, \
                                              Plot_Analysis, Soil_Layer_Index_DA, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, Post_Inflation_Alpha[1], \
                                              Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Datetime_Start, Datetime_Initial, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, Write_DA_File_Flag, r, Observation_Box, Def_Region, Dim_CLM_State, Num_Local_Obs_Par, Model_Variance, DateString_Plot,
                                            Def_Multiresolution, Def_ReBEL, Def_Localization, Assim_Algorithm_Name, eps[1], msw_infl[1], Region_Name, Call_Gstat_Flag, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,Proj_String, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution,
                                            dtime, Irrigation_Hours, column_len, Weather_Forecast_Days, Datetime_End, Hydraulic_File_Name, fpftcon_name, Run_Dir_Array, \
                                            Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Dim_Observation_Quantity, 
                                            Snow_Layer_Num, Def_Write_Initial, cols1d_ixy, cols1d_jxy, cols1d_ityplun, pfts1d_ityplun, Freezing_temperature_of_fresh_water, Density_of_ice, N0, nlyr,
                                            Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,
                                            Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,
                                            diskless_flag, persist_flag, Irrig_Scheduling, Run_Dir_Home, Start_Month, Stop_Year, Stop_Month, Stop_Day, Stop_Hour, UTC_Zone, finidat_name, Density_of_liquid_water, Irrigation_Grid_Flag_Array,
                                            mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Datetime_Stop, Datetime_Stop_Init, CLM_NA,
                                            Observation_Variance, Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Corelation_Par, octave, Station_XY, Station_XY_Index, Soil_Layer_Num, Analysis_Variable_Name,
                                            Analysis_Grid, Localization_Map_Mask, ObsModel_Mat, ObsModel_Variance_Mat, Mask_Sub, Mask_Index_Sub, Mask_Index_Vector, Observation_Matrix_Index, Prop_Grid_Array_Sys_Index, Model_State,
                                            SensorQuantity_Index, E0_ObsModel_Mask, Soil_Par_Sens, Veg_Par_Sens, PFT_Par_Sens, Hard_Par_Sens, Soil_Par_Accum_Dim, Veg_Par_Accum_Dim, PFT_Par_Accum_Dim, Hard_Par_Accum_Dim, ParFlow_Layer_Num,
                                            Forcing_File_Path_Home, Observation_Path, DAS_Data_Path, Grid_Resolution_CEA, Grid_Resolution_GEO, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy,  \
                                            NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Bias_Copy, NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy,
                                            NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Parameter_Copy, NC_FileName_Assimilation_2_Parameter_Obs_Dim, NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy, NC_FileName_Parameter_Space_Single, DAS_Output_Path, \
                                            COSMIC_Py, [], memory_profiler, COSMIC, Observation_Time_File_Path)
                                 
                                           
                                NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
                                for Station_Index in range(numpy.size(Station_XY)/2):
                                    if (numpy.size(numpy.where(numpy.asarray(Soil_Par_Sens_Array) == True)) >= 1):
                                        Parameter_Soil_Optimized[:,:,Station_Index] = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,:,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                                    if (numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) >= 1):
                                        Parameter_PFT_Optimized[:,:,Station_Index] = NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,:,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                                    
                                NC_File_Out_Optimized_Parameter = netCDF4.Dataset(NC_FileName_Optimized_Parameter, 'a')
                                if (numpy.size(numpy.where(numpy.asarray(Soil_Par_Sens_Array) == True)) >= 1):
                                    NC_File_Out_Optimized_Parameter.variables['Parameter_Soil_Optimized'][Optimized_Parameter_Index[0],:,:,:] = Parameter_Soil_Optimized
                                if (numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) >= 1):
                                    NC_File_Out_Optimized_Parameter.variables['Parameter_PFT_Optimized'][Optimized_Parameter_Index[2],:,:,:] = Parameter_PFT_Optimized
                                                               
                                NC_File_Out_Optimized_Parameter.sync()
                                NC_File_Out_Optimized_Parameter.close()
                                 
                                NC_File_Out_Assimilation_2_Parameter.close()
                                 
                                #print "Parameter_Soil_Space_Ensemble[:,13,Station_XY_Index[0][1],Station_XY_Index[0][0]]",Parameter_Soil_Space_Ensemble[:,13,Station_XY_Index[0][1],Station_XY_Index[0][0]]
                                #print "Parameter_Soil_Optimized_Array",Parameter_Soil_Optimized_Array[0][:,13,0],Parameter_Soil_Optimized_Array[-1][:,13,0]
                                 
                                if Plot_Analysis:           
                                    print "################################################# Plot Parameter Results"
                                    Plot_Parameters(Def_Print, fm, legend, plt, cm, colors, r, Def_Region, DasPy_Path, Region_Name, Row_Numbers, Col_Numbers, DAS_Data_Path, Row_Numbers_String, Col_Numbers_String, Dim_Soil_Par, Dim_Veg_Par, Start_Month, DateString_Plot, Variable_Assimilation_Flag, Mask_Index, PFT_Dominant_Index,
                                                    Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Ensemble_Number, Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Station_XY, Station_XY_Index, NC_FileName_Assimilation_2_Parameter, NC_FileName_Optimized_Parameter, NC_FileName_Parameter_Space_Single)               
                                                         
                        
                        #os.abort()
                        print "*************************************************Start Data Assimilation***************************************************************"
                        print ""
                        Analysis_Grid, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, job_server_node_array, active_nodes_server = \
                        Assimilation_Update(mpi4py_comm, mpi4py_rank, mpi4py_name, Model_Driver, NSLOTS, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, Irrig_Scheduling, Irrigation_Hours,  Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Model_Path, CLM_Flag, Def_PP, job_server_node_array, active_nodes_server,
                                              Start_Year,Start_Month,Start_Day,Stop_Year,Stop_Month,Stop_Day,Stop_Hour, UTC_Zone, Datetime_Start,Datetime_Start_Init,Datetime_Stop,Datetime_Stop_Init, Datetime_End, Datetime_Initial, Weather_Forecast_Days, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water, N0, nlyr,
                                              DAS_Data_Path, DAS_Depends_Path, DasPy_Path, Def_Multiresolution, Def_ReBEL, Def_Localization, Num_Local_Obs_State, eps[0], msw_infl[0], Plot_Analysis, Def_Figure_Output, DateString_Plot,
                                              Def_Write_Initial,DA_Flag, Write_DA_File_Flag, Mask, Mask_Index, COSMOS_Circle_Array, COSMOS_Circle_Index_Array, COSMOS_Circle_Num_Array, Call_Gstat_Flag, mksrf_edgee, mksrf_edges, mksrf_edgew, mksrf_edgen, Station_XY_Index, Station_XY, Observation_Box,
                                              Variable_Assimilation_Flag, Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, Model_Variance, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, PP_Servers_Per_Node, Def_CESM_Multi_Instance, PP_Port,
                                              Z_Resolution, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper,Proj_String, MODEL_CEA_X, MODEL_CEA_Y, Hydraulic_File_Name, Assim_Algorithm_Name, Low_Ratio_Par, High_Ratio_Par, Post_Inflation_Alpha[0], irrig_nsteps_per_day, PFT_Num, 
                                              Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset, fpftcon_name, Crop_Sum, 
                                              Model_State_Inflation_Range, Model_State_Inflation_Range_STD, Model_Bias_Range, Observation_Bias_Range, Model_Bias_Range_STD, Observation_Bias_Range_STD, Model_Bias_STD, Observation_Bias_STD, Dim_Observation_Quantity, 
                                              Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array, finidat_name,
                                              Ensemble_Number, Ensemble_Number_Predict,  Soil_Layer_Num, Snow_Layer_Num, maxpft, Forcing_File_Path_Home, dtime, Observation_Path, Dim_CLM_State, Dim_Obs_Type, CLM_NA, NAvalue, Variable_List, ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor,
                                              Grid_Resolution_CEA, Grid_Resolution_GEO, SensorQuantity_Sub, SensorType_Sub, SensorVariable_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub,Analysis_Variable_Name, Soil_Layer_Index_DA,
                                              Observation_Matrix, Observation_Variance, Observation_Latitude, Observation_Longitude, Observation_NLons, Observation_NLats, Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Corelation_Par,
                                              octave, r, Def_CDF_Matching, numrad, cols1d_ixy, cols1d_jxy, pfts1d_ixy, pfts1d_jxy, cols1d_ityplun, pfts1d_ityplun, column_len, pft_len, pfts1d_itypveg, pfts1d_ci,
                                              diskless_flag, persist_flag, Forcing_File_Path_Home, Forcing_File_Path_Array, history_file_name, Constant_File_Name, Run_Dir_Array, Feedback_Assim, 
                                              Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Texture_Layer_Opt_Num, Soil_Sand_Clay_Sum, Parameter_Range_Soil, Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Parameter_Regularization, Par_Soil_Uniform_STD_Sub, [], Par_PFT_Uniform_STD_Sub, [], \
                                              Analysis_Grid, Localization_Map_Mask, ObsModel_Mat, ObsModel_Variance_Mat, Prop_Grid_Array_Sys_Index, Observation_Matrix_Index, Mask_Sub, Mask_Index_Sub, 
                                              SensorQuantity_Index, Soil_Par_Sens_Dim, Veg_Par_Sens_Dim, PFT_Par_Sens_Dim, Hard_Par_Sens_Dim, ParFlow_Layer_Num,
                                              NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial,  NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Assimilation_2_Initial_Copy, NC_FileName_Assimilation_2_Bias_Copy, 
                                              NC_FileName_Assimilation_2_Bias_Monthly, NC_FileName_Assimilation_2_Bias_Monthly_Copy, NC_FileName_Assimilation_2_Parameter_Monthly, NC_FileName_Assimilation_2_Parameter_Monthly_Copy, 
                                              NC_FileName_Parameter_Space_Single, DAS_Output_Path, COSMIC_Py, [], memory_profiler, COSMIC, finidat_name_string, Observation_Time_File_Path)
                       
                        if Plot_Analysis:
                            print "###################################### Plot the Updated Model States"
                                
                            Plot_States(octave, fm, legend, plt, cm, colors, Def_Region, Region_Name, Plot_Analysis, DasPy_Path, Soil_Layer_Num, Ensemble_Number, Row_Numbers, Col_Numbers, Dim_Obs_Type, Observation_Matrix, NAvalue, Def_Print,
                                        Prop_Grid_Array_Sys_Index, Observation_Matrix_Index, SensorType_Sub, SensorVariable_Sub, SensorQuantity_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub, Variable_List,
                                        Mean_Index_Prop_Grid_Array_Sys, Mask_Index, Analysis_Grid, Analysis_Variable_Name, DateString_Plot, Variable_Assimilation_Flag,
                                        NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy, ObsModel_Mat_Masked, Observation_File_Name)
                        
                        
                        ##########################################                                        
                        
                        if SensorType_Sub == "COSMOS" or SensorType_Sub == "InSitu":
                            Def_Localization = Def_Localization_Original
                    
                # MPI
                if Def_PP == 2:
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
    
                if Parameter_Optimization:
                    Def_Par_Optimized = 1
                                                                   
                #raw_input("Press Enter to continue...")
                
                Def_First_Run = 0
                Def_First_Run_Bias = 0
                Two_Step_Bias_Estimation_Flag = 0
                
                # don't init cesm mpi again
                CESM_Init_Flag = 0
                
            else:
                if mpi4py_rank == 0:
                    print "************************ Datetime_Stop > Datetime_End, Skip to Simulation *****************"
                break
            
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
    
            if mpi4py_rank == 0:
                if (Datetime_Stop > Datetime_Start):
                
                    print "############################### Delete the Histroy Ensembles to Save Space in Daily Step"
                    # Folder to Save Ensemble Mean
                    Mean_Dir = Run_Dir_Home+"_Ens_Mean"
                    if not os.path.exists(Mean_Dir):
                        os.makedirs(Mean_Dir)
                    
                    #Datetime_Start_Mean = datetime.datetime(Datetime_Start.year,Datetime_Start.month,Datetime_Start.day,0,0)
                    #Datetime_Stop_Mean = datetime.datetime((Datetime_Stop-datetime.timedelta(days=1)).year,(Datetime_Stop-datetime.timedelta(days=1)).month,(Datetime_Stop-datetime.timedelta(days=1)).day,23,0)
                    
                    Datetime_Start_Mean = Datetime_Start
                    Datetime_Stop_Mean = Datetime_Stop
                    
                    stop_tod_string_final = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
                    
                       
                print "-------------------- Remove the old Initial Files"   # Because of the Get_Ens_Mean cannot do it
                Datetime_Start_Temp = Datetime_Start-datetime.timedelta(hours=1)
                Datetime_Start_Temp_Init_Mean = datetime.datetime(Datetime_Start_Temp.year, Datetime_Start_Temp.month, Datetime_Start_Temp.day, 00, 00)
                stop_tod_string = str((Datetime_Start_Temp - Datetime_Start_Temp_Init_Mean).seconds).zfill(5)
                restart_file_name_last_step = Region_Name + '.clm2.r.' + str(Datetime_Start_Temp.year) + '-' + str(Datetime_Start_Temp.month).zfill(2) + '-' + str(Datetime_Start_Temp.day).zfill(2) + '-' + stop_tod_string + '.nc'
                Command_String = "rm -irdf "+Run_Dir_Home+"*/"+restart_file_name_last_step
                print Command_String
                if (restart_file_name_last_step != finidat_initial_CLM_Copy):
                    subprocess.call(Command_String,shell=True)
                
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
    
            #print Parameter_Range_PFT
            # After the Assimilation, We should use the new initial file
            Def_Initial = 1
            finidat_initial_CLM = finidat_name
            Datetime_Start = Datetime_Stop + datetime.timedelta(hours=1)
            Datetime_Start_Init = datetime.datetime(Datetime_Start.year,Datetime_Start.month,Datetime_Start.day,00,00)
            Start_Year = str(Datetime_Start.year).zfill(4)
            Start_Month = str(Datetime_Start.month).zfill(2)
            Start_Day = str(Datetime_Start.day).zfill(2)
            Start_Hour = str(Datetime_Start.hour).zfill(2)
            
            Constant_File_Name_Header = Region_Name + ".clm2.h0."+Start_Year+"-"+Start_Month+"-"+Start_Day+"-"+str((Datetime_Start - Datetime_Start_Init).seconds).zfill(5)+".nc"
            if mpi4py_rank == 0:
                print "============================================================================================================================"
                print "Constant_File_Name_Header",Constant_File_Name_Header
                print "============================================================================================================================"
            
            # We only read the column and pft index during the first run
            Def_Read_Index = 0
            
            Def_First_Run = 0
            CESM_Init_Flag = 0
#    
            SensorType = []
            SensorVariable = []
            SensorQuantity = []
            Variable_ID = []
            QC_ID = []
            SensorResolution = []
            Observation_File_Name = []
            
            #------------------------------------------- Data Assimilation Flags
            
            Variable_Assimilation_Flag = numpy.zeros(Dim_CLM_State)
            
            if mpi4py_rank == 0:
                del Observation_Matrix, Observation_Longitude, Observation_Latitude, Observation_Variance, Observation_NLons, Observation_NLats
                del Observation_X_Left, Observation_X_Right, Observation_Y_Lower, Observation_Y_Upper, Observation_Misc, Observation_View_Zenith_Angle, Observation_View_Time, Observation_Corelation_Par
                        
                plt.close('all')
                gc.collect()
                del gc.garbage[:]
                
                r('gc(TRUE)')
                
                mem_usage = memory_profiler.memory_usage(proc=-1, interval=.1, timeout=None)
                print "mem_usage",mem_usage
            
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
    
            Observation_Index = Observation_Index + 1
            
            end = time.time()
            if mpi4py_rank == 0:
                print 'Time Is: ', (end - start) / 3600.0, 'Hours'
        
        if mpi4py_rank == 0:
            Observation_Time_File.close()
        
        if Datetime_Start > Datetime_End:
            Datetime_Stop = Datetime_Start
            Stop_Year = Start_Year
            Stop_Month = Start_Month
            Stop_Day = Start_Day
            Stop_Hour = Start_Hour
                    
        if Datetime_Start < Datetime_End:
            Stop_Year = End_Year
            Stop_Month = End_Month
            Stop_Day = End_Day
            Stop_Hour = End_Hour
            Datetime_Stop = Datetime_End
            Datetime_Stop_Init = Datetime_End_Init
            
            if mpi4py_rank == 0:
                print "**************** Drive CLM from",Datetime_Start,"to",Datetime_Stop,"after data assimilation"
            
            Do_DA_Flag = 0
            
            if mpi4py_rank == 0:
            
                
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
                            job_server_node = job_server_node_array[numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]
                            #job_server_node = job_server_node_array[Node_Index]
                            if Ens_Index > Ensemble_Number - 1:
                                break
                            job_server_node_results.append(job_server_node.submit(Prepare_Model_Operator, args=(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                    Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                    CLM_File_Name_List, Parameter_Range_Soil,
                                                  Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, 
                                                  DAS_Data_Path, DasPy_Path, DAS_Output_Path, Forcing_File_Path_Array, dtime, Variable_Assimilation_Flag, Variable_List,\
                                                  Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, 
                                                  DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Variable_Assimilation_Flag[Variable_List.index("Irrigation_Scheduling")],\
                                                  omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                  Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                                                  NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single,),
                                          depfuncs=(Run_CLM, Call_CLM_3D, Write_datm_atm_in, Write_datm_streams_txt, Write_presaero_stream_txt, Write_lnd_in, Write_rof_in, Write_Config_Files, Write_drv_in, Write_seq_maps),
                                          modules=("numpy", "netCDF4",  "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "signal", "gc", "imp", "getpass", "calendar", "glob","scipy.stats","scipy.signal",'scipy.weave',), group='Prepare_Model_Operator'))
                            
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
                        Prepare_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                    Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                    CLM_File_Name_List, Parameter_Range_Soil,
                                                  Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, 
                                                  DAS_Data_Path, DasPy_Path, DAS_Output_Path, Forcing_File_Path_Array, dtime, Variable_Assimilation_Flag, Variable_List,\
                                                  Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, 
                                                  DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Variable_Assimilation_Flag[Variable_List.index("Irrigation_Scheduling")],\
                                                  omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                  Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, \
                                                  NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Assimilation_2_Parameter, NC_FileName_Parameter_Space_Single)    
            
            if Def_PP == 2:
                mpi4py_comm.barrier()
                mpi4py_comm.Barrier()
            
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
                            job_server_node = job_server_node_array[numpy.min([Job_Num_Per_Node_Index+Node_Index*len(job_server_node_array)/len(active_nodes_server),len(job_server_node_array)-1])]
                            #job_server_node = job_server_node_array[Node_Index]
                            if Ens_Index > Ensemble_Number - 1:
                                break
                                                                                                            
                            job_server_node_results.append(job_server_node.submit(Call_Model_Operator, args=(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                                                                                 Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                                                                                 CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                                                                                                              Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Variable_Assimilation_Flag[Variable_List.index("Irrigation_Scheduling")],\
                                                                                                              Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                                                                                              Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                                                                                              NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null),
                                          depfuncs=(Run_CLM, Call_CLM_3D, Write_datm_atm_in, Write_datm_streams_txt, Write_presaero_stream_txt, Write_lnd_in, Write_rof_in, Write_Config_Files, Write_drv_in, Write_seq_maps),
                                          modules=("numpy", "netCDF4",  "sys", "os", "re", "unittest", "time", "datetime", "shutil", "fnmatch", "subprocess", "string", "socket", "signal", "gc", "imp", "getpass", "calendar", "glob","scipy.stats","scipy.signal",'scipy.weave',), group='Call_Model_Operator'))
                            
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
                                              Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Variable_Assimilation_Flag[Variable_List.index("Irrigation_Scheduling")],\
                                              Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                              Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                              NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single,COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null)    
                        
                    mpi4py_comm.barrier()
                    mpi4py_comm.Barrier()
                
                else:
                    print "*************************************************** Run DAS Sequentially"
                    for Ens_Index in range(Ensemble_Number):
                        Call_Model_Operator(Ens_Index, Model_Driver, Def_CESM_Multi_Instance, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized, Do_DA_Flag, Def_Debug, CLM_NA, NAvalue, finidat_initial_CLM, Def_ParFor, Def_Region, Def_Initial, \
                                                Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, Region_Name, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir_Array, Model_Path, CLM_Flag, num_processors,
                                                CLM_File_Name_List, Start_Year, Start_Month, Start_Day, Stop_Year, Stop_Month, Stop_Day, Datetime_Start, Datetime_Start_Init, Datetime_Stop, Datetime_Stop_Init, Datetime_End, Datetime_Initial, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Array, dtime,\
                                              Def_PP, N_Steps, Ensemble_Number, Ensemble_Number_Predict,  Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, DAS_Depends_Path, maxpft, ntasks_CLM, rootpe_CLM, nthreads_CLM, Weather_Forecast_Days, Variable_Assimilation_Flag[Variable_List.index("Irrigation_Scheduling")],\
                                              Low_Ratio_Par, High_Ratio_Par, Soil_Texture_Layer_Opt_Num, Def_Snow_Effects, PFT_Par_Sens_Array,\
                                              Soil_Thickness, Soil_Layer_Num, Snow_Layer_Num, Density_of_liquid_water, Initial_Perturbation, Initial_Perturbation_SM_Flag, Initial_Perturbation_ST_Flag, NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, 
                                              NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Bias, NC_FileName_Parameter_Space_Single, COUP_OAS_PFL, CESM_Init_Flag, mpi4py_comm_split, mpi4py_null)     
        
        
        
        
        #Datetime_Start_Mean = datetime.datetime(Datetime_Start.year,Datetime_Start.month,Datetime_Start.day,0,0)
        stop_tod_string_final = str((Datetime_Stop - Datetime_Stop_Init).seconds).zfill(5)
        
        Datetime_Start_Mean = Datetime_Start
        Datetime_Stop_Mean = Datetime_Stop
            
        
        end = time.time()
        if mpi4py_rank == 0:
            print 'Time Is: ', (end - start) / 3600.0, 'Hours'
    
    