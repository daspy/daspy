'''
    CLM and WRF Coupled System
'''

Def_PP = 0 # (0: Serial, 1: ParallelPython 2: MPI4Py)
mpi4py_comm = []
mpi4py_null = []
mpi4py_rank = 0
mpi4py_size = 0
mpi4py_name = []
if Def_PP == 2:
    from mpi4py import MPI
    try:
        import dill
        MPI.pickle.dumps = dill.dumps
        MPI.pickle.loads = dill.loads
    except:
        pass
    
    mpi4py_rank = MPI.COMM_WORLD.Get_rank()
    mpi4py_size = MPI.COMM_WORLD.Get_size()
    mpi4py_comm = MPI.COMM_WORLD
    mpi4py_null = MPI.COMM_NULL
    mpi4py_name = MPI.Get_processor_name()
    
#print "mpi4py_comm,mpi4py_rank,mpi4py_size,mpi4py_name",mpi4py_comm,mpi4py_rank,mpi4py_size,mpi4py_name

import os, sys, time, datetime, math, gc, subprocess, glob, string, shutil, copy, imp
import warnings, multiprocessing, socket, getpass, pickle, ctypes, platform
import numpy, scipy, netCDF4
sys.path.append('Utilities/Soil')
sys.path.append('Utilities')
sys.path.append('Algorithm')
sys.path.append('Algorithm/DAS')
sys.path.append('ForcingData')

Def_Figure_Output = 1
if Def_Figure_Output:
    # Generate images without having a window popup
    import matplotlib
    matplotlib.use('Agg')
    # Generate images without having a window popup
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from pylab import legend
import matplotlib.font_manager as fm

from DAS_Initialize import *
from DAS_Driver import *

#print os.getenv("PYTHONPATH")
start = time.time()
# Enable automatic garbage collection
gc.enable()

if mpi4py_rank == 0:

    print """Usage: python DAS.py [ncpus_main] - such as python DAS.py Ensember_number
        [ncpus] - the number of workers to run in parallel, 
        if omitted it will be set to the number of processors in the system
        
        Babaohe: python DAS.py 10
    """


Def_Region = 3    
Model_Driver = "CLM_45" 
PicHeight, PicWidth, RegionName, Row_Numbers, Col_Numbers, Grid_Resolution_CEA, Grid_Resolution_GEO, \
mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, Region_Name, Run_Dir_Home, DAS_Output_Path, Hydraulic_File_Name, \
Mask_File, Observation_Path, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Home, DAS_Depends_Path, geog_data_path, \
WRF_WPS_Path, WRF_WRF_Path, Station_XY, Station_XY_Index, r, octave, Row_Numbers_String, Col_Numbers_String, Grid_Resolution_CEA_String,\
xllcenter, yllcenter, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution, Proj_String, UTC_Zone = \
DAS_Initialize(Model_Driver, Def_Region, mpi4py_rank)

if Def_PP == 2:
    mpi4py_comm.barrier()
    mpi4py_comm.Barrier()

#print "mpi4py_rank",mpi4py_rank

Observation_Time_File_Path = DasPy_Path + "Examples/Rur/Only_LST_Par_LAI"

Def_CESM_Multi_Instance = 0 # for future
if mpi4py_rank == 0:
    print "*********************************************** Common Configuration"
Def_Run_DAS_Model    = 1     # Do Data Assimilation
Def_Run_Model        = 0     # for future
Def_Run_WRF         = 0     # for future
Def_Irrigation_Opt  = 0     # for future

Def_Snow_Effects    = 0     # for future

Feedback_Assim      = 0     # Whether to use LST update SM or use SM to update LST
Parameter_Optimization = 2  # Define whether to call the parameter optimization module (0: No 1: SODA 2: Augmentation)
Parameter_Regularization = 1.0  # for future
Def_Par_Sensitivity = 0     # for future
Def_Par_Correlation = 0     # for future
Def_Par_Optimized  = 1     # Define whether to use the optimized parameters

Def_First_Run       = 1  # 0 for restart run, 1 for first run, -1 for recover run if 0 fails. Define whether it is the first run
                        # It controls the copy and perturbation of surface data
Ensemble_Number         = 2    # Run CLM in Ensemble
Ensemble_Number_Predict = 100  # for future

Normal_Score_Trans      = 0 # for future
PDAF_Assim_Framework    = 0 # for future
PDAF_Filter_Type        = 5 # for future

if mpi4py_rank == 0:
    print "**********************************************************CLM******************************************************************"

Def_SpinUp      = 0  	# for future
Def_Print       = 1   # (0: No printed information) (1: short information) (2: medium output statistics for Debug) (3: full output statistics for Debug)
Plot_Analysis   = 1   # whether to plot the results (1: Plot few, 2: Plot more)
Def_Debug       = 0     # for future
Write_DA_File_Flag   = 0  # Define whether to write the assimilation files
Initial_Perturbation = 0    # Whether to perturb the initial state before starting CLM model to prevent filter divergence
Def_Localization    = 1  # Whether to use the Observation Horizontal Correlation in the general data fusion algorithm
                        # Or in LETKF, if 0, the observation variance will not be divided by the correlation coefficients
Observation_Box     = numpy.min([int(numpy.sqrt(Row_Numbers**2+Col_Numbers**2) / 2.0), 10])
# The Local Observation Window Increment for the Large Area, How many boundary grid cells should be considered
# Number of Observation used in State Local Analysis (16 is the best one for soil moisture)
#["Soil_Moisture","Surface_Temperature","Vegetation_Temperature","Canopy_Water","Albedo_BSA_Band_vis","Albedo_BSA_Band_nir","Albedo_WSA_Band_vis",
#                         "Albedo_WSA_Band_nir","Emissivity","Snow_Depth","Snow_Cover_Fraction","Snow_Water_Equivalent","LAI","Crop_Planting_Date","Crop_Harvest_Date",
#                         "Water_Storage","Water_Table","Irrigation_Rate","Irrigation_Scheduling"]
Num_Local_Obs_State = numpy.asarray([9, 9])
# Number of Observation used in Parameter Local Analysis (5 is better for soil moisture)
Num_Local_Obs_Par   = numpy.asarray([9, 9])
# Number of Observation used in Bias Local Analysis (5 is better for soil moisture)
Num_Local_Obs_Bias   = numpy.asarray([9, 9])
eps           = numpy.asarray([0.1, 0.1, 0.1])   # Threshold to Select Correlated Observations (State, Parameter and Bias)
msw_infl             = numpy.asarray([1.01, 1.01, 1.01])  #  (State, Parameter and Bias)inflation mode switch  #  < 0 : adaptive inflation    #  > 0 : fixed inflation value
                            # During dry period, the soil moisture needs the inflation

if mpi4py_rank == 0:
    print "-------------------------The Observation Box is",Observation_Box

Post_Inflation_Alpha = numpy.asarray([1.0, 1.0, 1.0])  # State, Parameter and Bias Inflation Alpha
if mpi4py_rank == 0:
    print "---------------------- Define the COSMOS Max Counting Rate"
N0 = 1132
nlyr = 300
Call_Gstat_Flag = 0  # for future

Def_Multiresolution = 1  # 1: DWT; 2: Contourlet
Def_ReBEL           = 1  # Whether to 1: use the Bayesian Filtering  0: use Optimal Interpolation(OI) algorithm 2: Direct Insertion (DI)
Def_Write_Initial   = 1  # whether to Output the Assimilation Results to the CLM Initial Files
Independent_Obs     = 1   # If Independent_Obsrevations_Flag = 0: Means the Observations are dependdent; Independent_Obsrevations_Flag = 1: Means the Observations are Independent

Def_CDF_Matching       = 0  # for future
Bias_Estimation_Option_Model = numpy.asarray([0, 0])     # for future
Bias_Estimation_Option_Obs = numpy.asarray([0, 0])     # for future

Low_Ratio_Par = 0.8         # Lower limit to perturb the parameters
High_Ratio_Par = 1.2        # Upper limit to perturb the parameters
Low_Ratio_Par_Uniform = 0.2  # for future
High_Ratio_Par_Uniform = 1.8  # for future

if mpi4py_rank == 0:
    print "************************************** Irrigation Configuration"
Irrig_Scheduling      = 0 # for future
Irrigation_Hours      = 0 # for future
Weather_Forecast_Days = 0 # for future
if Def_Region == -1:
    Irrig_Scheduling  = 2    # for future
    Irrigation_Hours      = 2  # for future
    Weather_Forecast_Days = 1    # for future
    
if mpi4py_rank == 0:
    print "************************************** Datetime Configuration"
# Model Start Time
Start_Year      = '2012'
Start_Month     = '01'
Start_Day       = '01'
Start_Hour      = '00'
Start_Minute    = '00'
# Model Start Time
Datetime_Start = datetime.datetime(string.atoi(Start_Year), string.atoi(Start_Month), string.atoi(Start_Day), string.atoi(Start_Hour), string.atoi(Start_Minute))
Datetime_Start_Init = datetime.datetime(string.atoi(Start_Year), string.atoi(Start_Month), string.atoi(Start_Day), 00, 00)

# Model End Time
End_Year        = '2012'
End_Month       = '01'
End_Day         = '03'
End_Hour        = '23'
End_Minute      = '00'
Datetime_End = datetime.datetime(string.atoi(End_Year), string.atoi(End_Month), string.atoi(End_Day), string.atoi(End_Hour), string.atoi(End_Minute))
Datetime_End_Init = datetime.datetime(string.atoi(End_Year), string.atoi(End_Month), string.atoi(End_Day), 00, 00)

LAI_Year_String = '2012'    # for future
MODIS_LAI_Data_ID = 'Lai_1km'   # for future

if Ensemble_Number == 1:
    Def_PP = 0
DAS_Fortran_Lib = []
if Def_PP and Plot_Analysis >= 2:
    Plot_Analysis   = 1
Def_ParFor      = 1  # Whether to Use Weave and OpenMP to Speed Up the For Loop
Def_Initial     = 1  # Define whether to use the provided initial files, 1: use 0: not use
Use_Mask_Flag   = 1  # Define whether to use the watershed boundary

#          0     1      2          3        4       5      6        7       8          9        10    11        12        13            14
lftype = ['ekf', 'ukf', 'cdkf', 'srukf', 'srcdkf', 'pf', 'gspf', 'sppf', 'gmsppf', 'letkf', 'enkf', '4dletkf', 'aenkf', 'letkf_ddsm','letkoi']
Assim_Algorithm_Name = lftype[9]
if mpi4py_rank == 0:
    print "Assimilation Algorithm is",Assim_Algorithm_Name

#*************************** Soil Parameter
#'PCT_SAND', 'PCT_CLAY', 'ORGANIC'
Soil_Par_Sens_Array = ['' for i in range(2)]
#*************************** Vegetation Parameter
Veg_Par_Sens_Array = ['' for i in range(2)]
#***************************** PFT Monthly Parameters
# 'LAI' 'SAI' 'HTOP'
PFT_Par_Sens_Array = ['' for i in range(2)]
#***************************** Hard Coded Parameters
Hard_Par_Sens_Array = ['' for i in range(2)]


#SensorVariable_Sub == "Soil_Moisture":
Soil_Par_Sens_Array[0] = numpy.array([True, True, True, False, False],dtype=numpy.bool)
Veg_Par_Sens_Array[0] = numpy.array([False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],dtype=numpy.bool)# for future
PFT_Par_Sens_Array[0] = numpy.array([False, False, False],dtype=numpy.bool)
Hard_Par_Sens_Array[0] = numpy.array([False, False, False, False, False],dtype=numpy.bool)# for future

#SensorVariable_Sub == "Surface_Temperature":
Soil_Par_Sens_Array[1] = numpy.array([False, False, False, False, False],dtype=numpy.bool)
Veg_Par_Sens_Array[1] = numpy.array([False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],dtype=numpy.bool)# for future
PFT_Par_Sens_Array[1] = numpy.array([True, True, False],dtype=numpy.bool)
Hard_Par_Sens_Array[1] = numpy.array([False, False, False, False, False],dtype=numpy.bool)# for future

NAvalue = -9999.0
CLM_NA = 1e36
Dim_Soil_Par = numpy.size(Soil_Par_Sens_Array[0])
Dim_Veg_Par = numpy.size(Veg_Par_Sens_Array[0])# for future
Dim_PFT_Par = numpy.size(PFT_Par_Sens_Array[0])
Dim_Hard_Par = numpy.size(Hard_Par_Sens_Array[0]) # for future

Soil_Texture_Layer_Opt_Num = 1

dtime = 3600    # Model time step (seconds)

Def_wget        = 1
Def_geogrid     = 1
Def_ungrib      = 1
Def_metgrid     = 1
Def_WRF         = 1
Def_Read_wrfout = 0
OMP_NUM_THREADS_WRF = '1'

Num_of_Days_Monthly, Datetime_Initial, NSLOTS, Constant_File_Name_Header, finidat_initial_CLM, finidat_initial_PFCLM, \
Soil_Layer_Num, Snow_Layer_Num, ParFlow_Layer_Num, maxpft, numrad, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water,\
ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor, Model_Path, CLM_Flag,\
Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,\
Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array = \
DAS_Config(mpi4py_rank, Model_Driver, Start_Year, Start_Month, Start_Day, Start_Hour, Start_Minute, Datetime_Start, \
                Def_CESM_Multi_Instance, Ensemble_Number, Region_Name, Def_Initial, Run_Dir_Home,\
                Datetime_Start_Init, Def_ParFor, DAS_Data_Path, Def_Region,\
                Def_PP, Row_Numbers, Col_Numbers, Def_Print, DAS_Depends_Path, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type)

if mpi4py_rank == 0:
    end = time.time()
    print 'Time Is: ', (end - start), 'Seconds'

if mpi4py_rank == 0:
    print "============================================================================================================================"
    print "=====================================================Call DAS CLM=========================================================="
    print "============================================================================================================================"

if Def_Run_DAS_Model:
    Do_DA_Flag = 1  # If Do_DA_Flag = 0: Run the CLM Only; If Do_DA_Flag = 1: Run the CLM and Do Data Assimilation
    DAS_Driver(mpi4py_comm, mpi4py_null, mpi4py_rank,  mpi4py_size, mpi4py_name, Model_Driver,Do_DA_Flag, Def_Par_Sensitivity, Def_Par_Correlation, Def_Par_Optimized,   Dim_Soil_Par, Dim_Veg_Par, Dim_PFT_Par, Dim_Hard_Par, Soil_Texture_Layer_Opt_Num, Observation_Box, LAI_Year_String, MODIS_LAI_Data_ID,\
             Num_of_Days_Monthly, Start_Year, Start_Month, Start_Day, Start_Hour, Start_Minute, End_Year, End_Month, End_Day, End_Hour, End_Minute, Datetime_Start, Datetime_Start_Init, \
             Datetime_End, Datetime_End_Init, Datetime_Initial, UTC_Zone, CLM_NA, NAvalue, Assim_Algorithm_Name, Station_XY, Station_XY_Index, dtime, \
             NSLOTS, Feedback_Assim, Parameter_Optimization, Parameter_Regularization, Def_CDF_Matching, Bias_Estimation_Option_Model, Bias_Estimation_Option_Obs, Post_Inflation_Alpha, Def_Snow_Effects, N0, nlyr,\
             Sub_Block_Ratio_Row, Sub_Block_Ratio_Col, Sub_Block_Index_Row_Mat_Vector, Sub_Block_Index_Col_Mat_Vector, Row_Offset, Col_Offset,\
             Row_Numbers_SubBlock_Array, Col_Numbers_SubBlock_Array, Sub_Block_Row_Start_Array, Sub_Block_Row_End_Array, Sub_Block_Col_Start_Array, Sub_Block_Col_End_Array,\
             Observation_Time_File_Path, Def_CESM_Multi_Instance, Constant_File_Name_Header, finidat_initial_CLM, finidat_initial_PFCLM, Def_PP, DAS_Fortran_Lib, Normal_Score_Trans, PDAF_Assim_Framework, PDAF_Filter_Type, Def_ParFor, Def_Region, Def_Initial, Irrig_Scheduling, Irrigation_Hours,  Def_SpinUp, Def_First_Run, Def_Print, CLM_Flag,  Def_ReBEL, Def_Localization, \
             Num_Local_Obs_State, Num_Local_Obs_Par,  Num_Local_Obs_Bias, eps, msw_infl, Def_Multiresolution, Def_Write_Initial, Ensemble_Number, Ensemble_Number_Predict,  Call_Gstat_Flag, Write_DA_File_Flag, Use_Mask_Flag, Def_Figure_Output,\
             Forcing_File_Path_Home, Soil_Layer_Num, Snow_Layer_Num, ParFlow_Layer_Num, maxpft, numrad, Density_of_liquid_water, Density_of_ice, Freezing_temperature_of_fresh_water, Plot_Analysis, Def_Debug, Initial_Perturbation, \
             Weather_Forecast_Days, PicHeight, PicWidth, RegionName, Row_Numbers, Col_Numbers, Row_Numbers_String, Col_Numbers_String, Grid_Resolution_CEA_String, xllcenter, yllcenter, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution, Proj_String, \
             Grid_Resolution_CEA, Grid_Resolution_GEO, mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, ntasks_CLM, rootpe_CLM, nthreads_CLM, omp_get_num_procs_ParFor, Low_Ratio_Par, High_Ratio_Par, Low_Ratio_Par_Uniform, High_Ratio_Par_Uniform, \
             Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Region_Name, Run_Dir_Home, Model_Path, Hydraulic_File_Name, Mask_File, Observation_Path, DAS_Data_Path, DasPy_Path, DAS_Output_Path, DAS_Depends_Path, octave, r, plt, cm, colors, inset_axes, fm, legend)

if mpi4py_rank == 0:
    print "============================================================================================================================"
    print "=======================================================Finishing============================================================"
    print "============================================================================================================================"
