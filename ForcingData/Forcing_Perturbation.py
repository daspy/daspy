import os, sys, time, datetime, random, math, gc, subprocess, glob, string, shutil, warnings, commands
import numpy, scipy, scipy.io.netcdf, netCDF4, pp
from scipy import ogrid, sin, mgrid, ndimage, array, signal, stats
from numpy import min,max


Def_Figure_Output = 1
if Def_Figure_Output:
    # Generate images without having a window popup
    import matplotlib
    matplotlib.use('Agg')
    # Generate images without having a window popup

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.font_manager as fm
from pylab import legend

sys.path.append('../../')
sys.path.append('../../ForcingData')
sys.path.append('../../Utilities')
from DAS_Initialize import *

'''
Precipitation forcing was perturbed using multiplicative noise with a log-normal distribution (mean=1,std=0.25); 
short-wave radiative forcing was perturbed using multiplicative noise with normal distribution [N(mean=1, std=0.7)]; 
air temperature forcing and long-wave forcing data were perturbed using additive noises with normal distribution [N(mean=0, std=2.5 K) and N(mean=0, std=10 W.m-2) respectively]. 
The above perturbations are independent. 
The precipitation perturbation multiplication factor was limited between 0 and 4
where the actual precipitation value was further prevented to exceed the true precipitation value with (+-)5mm/hour 
in ensemble generation.
 The short-wave perturbation multiplication factor was limited between 0.2 and 1.8.
Temperature and long-wave radiation perturbations were limited to (+-)4 times their respective standard deviations.
'''
# Enable automatic garbage collection
gc.enable()

# Indicate which year will be processed!
Year_String = '2012'

Def_Region = -2    #-3(LW) -2 (HiWATER_M) -1 (AGADAPT) 0(Hulugou),1(Babaohe),2(ZWA),3(Rur),4(Heihe),5(Shiyang),6(Haihe),7(Tahe),8(China),9(Global)
Model_Driver = "CLM_45" # CLM_40 CLM_45 CLM_CN PFCLM TSMP 0: cesm_sp 1: cesm_cn, 2: cesm_ad_spinup, 3: cesm_exit_spinup, 4: cesm_cn_spinup 5: cesm_irrig
PicHeight, PicWidth, RegionName, Row_Numbers, Col_Numbers, Grid_Resolution_CEA, Grid_Resolution_GEO, \
mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, Region_Name, Run_Dir_Home, DAS_Output_Path, Hydraulic_File_Name, \
Mask_File, Observation_Path, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Home, DAS_Depends_Path, geog_data_path, \
WRF_WPS_Path, WRF_WRF_Path, Station_XY, Station_XY_Index, r, octave, Row_Numbers_String, Col_Numbers_String, Grid_Resolution_CEA_String,\
xllcenter, yllcenter, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution, Proj_String, UTC_Zone = \
DAS_Initialize(Model_Driver, Def_Region)

#print "mpi4py_rank",mpi4py_rank

#--------------------------------- UTC Time ----------------------------------
Start_Year = 2012
Start_Month = 1
Start_Day = 1
Start_Hour = 0
Start_Minute = 0
    
Stop_Year = 2012
Stop_Month = 8
Stop_Day = 31
Stop_Hour = 23
Stop_Minute = 0

NTime_Out_Daily = 24

Ensemble_Number = 20

NC_Output_Path=DAS_Data_Path+"ForcingData/Bilinear_1km_1hour_"+Region_Name+"/CLM"
#NC_Output_Path="/disk4/hanxujun/DAS_Data/ForcingData/Micromet_1km_1hour_"+Region_Name+"/CLM"

#################################### Define the Perturbation Methods #############################################

#std min max
TBOT_Perturbation_Par = numpy.array([5.0, -5.0, 5.0])
SHUM_Perturbation_Par = numpy.array([0.0, 0.8, 1.2])
PRECTmms_Perturbation_Par = numpy.array([0.0, 0.8, 1.2])
FSDS_Perturbation_Par = numpy.array([10.0, 0.9, 1.1])
PSRF_Perturbation_Par = numpy.array([10.0, -100.0, 100.0])
FLDS_Perturbation_Par = numpy.array([10.0, 0.9, 1.1])
WIND_Speed_Perturbation_Par = numpy.array([1.0, 0.8, 1.2])
##################################################################################################################

Datetime_Start  = datetime.datetime(string.atoi(Year_String),12,31,23,0)
Datetime_Initial  = datetime.datetime(string.atoi(Year_String),1,1,0,0)
N_Steps = ((Datetime_Start - Datetime_Initial).days + 1) * 24
print "There are",N_Steps,"N_Steps"

TBOT_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)
SHUM_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)
PRECTmms_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)
FSDS_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)
FLDS_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)
WIND_Ensemble = numpy.zeros((N_Steps,Ensemble_Number),dtype=numpy.float32)

PRECTmms_Ensemble_Multiplier = numpy.zeros(N_Steps,dtype=numpy.float32) 

TBOT_Series = numpy.zeros(N_Steps,dtype=numpy.float32)
TBOT_Series_Filter = numpy.zeros(N_Steps,dtype=numpy.float32)
SHUM_Series = numpy.zeros(N_Steps,dtype=numpy.float32)
PRECTmms_Series = numpy.zeros(N_Steps,dtype=numpy.float32)
FSDS_Series = numpy.zeros(N_Steps,dtype=numpy.float32)
FLDS_Series = numpy.zeros(N_Steps,dtype=numpy.float32)
WIND_Series = numpy.zeros(N_Steps,dtype=numpy.float32)

Stop_Year_Loop = Stop_Year + 1

for Year_Index in range(Start_Year,Stop_Year_Loop):
    
    Stop_Month_Interpolation = Stop_Month
    if Stop_Year > Year_Index:
        Stop_Month_Interpolation = 12
    elif Stop_Month <= Start_Month and Stop_Year > Year_Index:
        Stop_Month_Interpolation = Stop_Month
    
    # Model Start Time
    Datetime_Start = datetime.datetime(Year_Index,Start_Month,Start_Day,Start_Hour,Start_Minute)
    
    # Model Stop Time
    Datetime_Stop = datetime.datetime(Year_Index,Stop_Month_Interpolation,Stop_Day,Stop_Hour,Stop_Minute)
    
    # Model Initial Time
    Initial_Year = Year_Index
    Initial_Month = 1
    Initial_Day = 1
    Initial_Hour = 0
    Initial_Minute = 0
    Datetime_Initial = datetime.datetime(Initial_Year,Initial_Month,Initial_Day,Initial_Hour,Initial_Minute)
    
    #If you want to simulate one month: time range should be 1.1 2.1
    
    Time_Interval_Hourly = datetime.timedelta(hours=1)
    Time_Interval_Daily = datetime.timedelta(hours=24)
    
    Start_Day_Index = (Datetime_Start - Datetime_Initial).days
    End_Day_Index = (Datetime_Stop - Datetime_Initial).days
    
    Number_of_Days = End_Day_Index - Start_Day_Index
    print '\n'
    print 'There are', Number_of_Days, 'Days Need to Be Generated from the',str(Start_Day_Index+1)+'th','Day!'
    print '\n'
    
    # Get the Number of Days in Each Day
    Num_of_Days_Monthly = numpy.zeros(12)
    Num_of_Days_Monthly = [31,28,31,30,31,30,31,31,30,31,30,31]
    Num_of_Days_Monthly[1] = (datetime.datetime(Initial_Year,3,1,0)-datetime.datetime(Initial_Year,2,1,0)).days
    print Num_of_Days_Monthly
        
    start = time.time()
    
    
    #================================= Generate Met Files
    for Month_Index in range(Start_Month,Stop_Month_Interpolation+1):
        print '\n'
        print 'The ',Year_Index,str(Month_Index)+'th','Month is Being Processing!'
        print '\n'
        
        Num_of_Records_Daily = 24    # Number of the Records in Each Day
        
        Num_Days = Num_of_Days_Monthly[Month_Index-1] - Start_Day + 1
        
        NTime_Out_Monthly   = int(24 * Num_Days)
        NTime_Out_Daily = 24
        
    #    Num_Days = 1    # For Testing
    #    NTime_Out_Monthly = 2
    #    NTime_Out_Daily = 1
        
        Year = numpy.zeros(NTime_Out_Monthly,dtype=numpy.integer)
        Month = numpy.zeros(NTime_Out_Monthly,dtype=numpy.integer)
        Day = numpy.zeros(NTime_Out_Monthly,dtype=numpy.integer)
        Hour = numpy.zeros(NTime_Out_Monthly)
    
        for Time_Step in range(NTime_Out_Monthly):
            Year[Time_Step]    = Datetime_Start.year
            Month[Time_Step]   = Datetime_Start.month
            Day[Time_Step]     = Datetime_Start.day
            Hour[Time_Step]    = Datetime_Start.hour
            Datetime_Start = Datetime_Start + Time_Interval_Hourly
            
        
        for Day_Index in range(Start_Day,Start_Day+Num_Days):
                                
            if Month_Index < 10:
                Month_String = "-0" + str(Month_Index)
            elif Month_Index >= 10 and Month_Index <= 12:
                Month_String = "-" + str(Month_Index)
        
            if Day_Index < 10:
                Day_String = "-0" + str(Day_Index)
            elif Day_Index >= 10:
                Day_String = "-" + str(Day_Index)
        
        
            FileName = NC_Output_Path + "/" + str(Year_Index) + Month_String + Day_String + ".nc"
            
            print "Open the Original Forcing Data Files",FileName
            Data_File = netCDF4.Dataset(FileName, 'r')
            TBOT_Whole = Data_File.variables['TBOT'][:,:,:]
            SHUM_Whole = Data_File.variables['SHUM'][:,:,:]
            PRECTmms_Whole = Data_File.variables['PRECTmms'][:,:,:]
            FSDS_Whole = Data_File.variables['FSDS'][:,:,:]
            PSRF_Whole = Data_File.variables['PSRF'][:,:,:]
            FLDS_Whole = Data_File.variables['FLDS'][:,:,:]
            WIND_Whole = Data_File.variables['WIND'][:,:,:]
            Data_File.close()
            
            
            Random_TBOT_ARMA = numpy.random.normal(loc=0.0,scale=2.5,size=(Ensemble_Number,NTime_Out_Daily,Row_Numbers,Col_Numbers))
            Random_FSDS_ARMA = numpy.random.normal(loc=1.0,scale=0.7,size=(Ensemble_Number,NTime_Out_Daily,Row_Numbers,Col_Numbers))
            Random_FLDS_ARMA = numpy.random.normal(loc=0.0,scale=10.0,size=(Ensemble_Number,NTime_Out_Daily,Row_Numbers,Col_Numbers))
            Random_Prec_ARMA = numpy.random.lognormal(mean=1.0,sigma=0.25,size=(Ensemble_Number,NTime_Out_Daily,Row_Numbers,Col_Numbers))
            
            for Time_Step in range(NTime_Out_Daily):
                
                for Ens_Index in range(Ensemble_Number):

                    Random_TBOT_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_TBOT_ARMA[Ens_Index,Time_Step,:,:] < -5.0)] = -5.0
                    Random_TBOT_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_TBOT_ARMA[Ens_Index,Time_Step,:,:] >  5.0)] =  5.0
                    Random_FSDS_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_FSDS_ARMA[Ens_Index,Time_Step,:,:] <  0.2)] =  0.2
                    Random_FSDS_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_FSDS_ARMA[Ens_Index,Time_Step,:,:] >  1.8)] =  1.8
                    Random_FLDS_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_FLDS_ARMA[Ens_Index,Time_Step,:,:] < -40.0)] = -40.0
                    Random_FLDS_ARMA[Ens_Index,Time_Step,:,:][numpy.where(Random_FLDS_ARMA[Ens_Index,Time_Step,:,:] >  40.0)] =  40.0
                
        #        fig1 = plt.figure()
        #        ax1 = fig1.add_subplot(4, 1, 1)
        #        ax1.plot(Random_TBOT_ARMA[:,Time_Step,0,0])
        #        ax1.set_title('TBOT')
        #        ax1 = fig1.add_subplot(4, 1, 2)
        #        ax1.plot(Random_FSDS_ARMA[:,Time_Step,0,0])
        #        ax1.set_title('FSDS')
        #        ax1 = fig1.add_subplot(4, 1, 3)
        #        ax1.plot(Random_FLDS_ARMA[:,Time_Step,0,0])
        #        ax1.set_title('FLDS')
        #        ax1 = fig1.add_subplot(4, 1, 4)
        #        ax1.plot(Random_Prec_ARMA[:,Time_Step,0,0])
        #        ax1.set_title('Prec')
        #        plt.show()
        #    
        #    fig1 = plt.figure()
        #    ax1 = fig1.add_subplot(4, 1, 1)
        #    ax1.plot(Random_TBOT_ARMA[Ens_Index,:,0,0])
        #    ax1.set_title('TBOT')
        #    ax1 = fig1.add_subplot(4, 1, 2)
        #    ax1.plot(Random_FSDS_ARMA[Ens_Index,:,0,0])
        #    ax1.set_title('FSDS')
        #    ax1 = fig1.add_subplot(4, 1, 3)
        #    ax1.plot(Random_FLDS_ARMA[Ens_Index,:,0,0])
        #    ax1.set_title('FLDS')
        #    ax1 = fig1.add_subplot(4, 1, 4)
        #    ax1.plot(Random_Prec_ARMA[Ens_Index,:,0,0])
        #    ax1.set_title('Prec')
        #    plt.show()
            
            for Ens_Index in range(Ensemble_Number):
                Ens_Ouput_Path = NC_Output_Path+"_Ens"+str(Ens_Index+1)
                if not os.path.exists(Ens_Ouput_Path):
                    os.makedirs(Ens_Ouput_Path)
                
                CLM_NC_FileName_Out = Ens_Ouput_Path + "/" + str(Year_Index) + Month_String + Day_String + ".nc"
                
                print 'Write CLM NetCDF File:',CLM_NC_FileName_Out
                
                if os.path.exists(CLM_NC_FileName_Out):
                    os.remove(CLM_NC_FileName_Out)
                
                CLM_NC_File_Out = netCDF4.Dataset(CLM_NC_FileName_Out, 'w', diskless=True, persist=True, format='NETCDF4')
                # Dim the dimensions of netCDF: 4 dimensions
                CLM_NC_File_Out.createDimension('scalar', 1)
                CLM_NC_File_Out.createDimension('lon', Col_Numbers)
                CLM_NC_File_Out.createDimension('lat', Row_Numbers)
                CLM_NC_File_Out.createDimension('time', NTime_Out_Daily)
                
                EDGEE = CLM_NC_File_Out.createVariable('EDGEE','f4',('scalar',),zlib=True)
                EDGEE.long_name = 'eastern edge in atmospheric data'
                EDGEE.units = 'degrees east'        
                EDGEE.mode = "time-invariant"
                
                EDGEN = CLM_NC_File_Out.createVariable('EDGEN','f4',('scalar',),zlib=True)
                EDGEN.long_name = 'northern edge in atmospheric data'
                EDGEN.units = 'degrees north'        
                EDGEN.mode = "time-invariant"
                
                EDGES = CLM_NC_File_Out.createVariable('EDGES','f4',('scalar',),zlib=True)
                EDGES.long_name = 'southern edge in atmospheric data'
                EDGES.units = 'degrees north'        
                EDGES.mode = "time-invariant"
                
                EDGEW = CLM_NC_File_Out.createVariable('EDGEW','f4',('scalar',),zlib=True)
                EDGEW.long_name = 'western edge in atmospheric data'
                EDGEW.units = 'degrees east'        
                EDGEW.mode = "time-invariant"
                
                LATIXY = CLM_NC_File_Out.createVariable('LATIXY','f4',('lat','lon',),zlib=True)
                LATIXY.long_name = 'latitude'
                LATIXY.units = 'degrees north'        
                LATIXY.mode = "time-invariant"
                
                LONGXY = CLM_NC_File_Out.createVariable('LONGXY','f4',('lat','lon',),zlib=True)
                LONGXY.long_name = 'longitude'
                LONGXY.units = 'degrees east'
                LONGXY.mode = "time-invariant"
                            
                TBOT = CLM_NC_File_Out.createVariable('TBOT','f4',('time','lat','lon',),zlib=True)
                TBOT.long_name = 'temperature at the lowest atm level (TBOT)'
                TBOT.units = 'K'        
                TBOT.mode = "time-dependent"
                
                SHUM = CLM_NC_File_Out.createVariable('SHUM','f4',('time','lat','lon',),zlib=True)
                SHUM.long_name = 'specific humidity at the lowest atm level (SHUM)'
                SHUM.units = 'kg/kg'        
                SHUM.mode = "time-dependent"
                
                PRECTmms = CLM_NC_File_Out.createVariable('PRECTmms','f4',('time','lat','lon',),zlib=True)
                PRECTmms.long_name = 'precipitation (PRECT)'
                PRECTmms.units = 'mm/s'        
                PRECTmms.mode = "time-dependent"
                
                FSDS = CLM_NC_File_Out.createVariable('FSDS','f4',('time','lat','lon',),zlib=True)
                FSDS.long_name = 'incident solar (FSDS)'
                FSDS.units = 'W/m2'        
                FSDS.mode = "time-dependent"
                
                PSRF = CLM_NC_File_Out.createVariable('PSRF','f4',('time','lat','lon',),zlib=True)
                PSRF.long_name = 'surface pressure at the lowest atm level (PSRF)'
                PSRF.units = 'Pa'
                PSRF.mode = "time-dependent"
            
                FLDS = CLM_NC_File_Out.createVariable('FLDS','f4',('time','lat','lon',),zlib=True)
                FLDS.long_name = 'incident longwave radiation (FLDS)'
                FLDS.units = 'W/m2'        
                FLDS.mode = "time-dependent"
                        
                WIND = CLM_NC_File_Out.createVariable('WIND','f4',('time','lat','lon',),zlib=True)
                WIND.long_name = 'wind at the lowest atm level (WIND)'
                WIND.units = 'm/s'        
                WIND.mode = "time-dependent"
                
                times = CLM_NC_File_Out.createVariable('time','f4',('time',),zlib=True)
                times.long_name = "observation time"
                times.units = "hours since "+str(Year_Index)+Month_String+Day_String+" "+"00:00:00"
                times.calendar = "gregorian"
                
                CLM_NC_File_Out.title = "GLDAS atmospheric data at 1-hr resolution"
                CLM_NC_File_Out.conventions = "CF-1.0"
                #print CLM_NC_File_Out.case_title
                
                # Assign the Values
                longitudes = numpy.linspace(mksrf_edgew+Grid_Resolution_GEO[0]/2.0,mksrf_edgee-Grid_Resolution_GEO[0]/2.0,Col_Numbers)
                latitudes = numpy.linspace(mksrf_edges+Grid_Resolution_GEO[1]/2.0,mksrf_edgen-Grid_Resolution_GEO[1]/2.0,Row_Numbers)
                times[:] = numpy.linspace(0,23,NTime_Out_Daily)
                
                #print numpy.linspace(0,23,NTime_Out_Daily),times[:]
                
                EDGEW[:] = mksrf_edgew
                EDGEE[:] = mksrf_edgee
                EDGES[:] = mksrf_edges
                EDGEN[:] = mksrf_edgen
                LONGXY_Row = longitudes
                LATIXY_Col = latitudes
                for row in range(Row_Numbers):
                    LONGXY[row,:] = LONGXY_Row
                for col in range(Col_Numbers):
                    LATIXY[:,col] = LATIXY_Col
                #print LATIXY_Col
                
                
                for Time_Step in range(NTime_Out_Daily):
                    
                    #numpy.random.seed(123456789 + (Time_Step + (Day_Index-1)*NTime_Out_Daily))
                    TBOT[Time_Step,:,:] = TBOT_Whole[Time_Step,::] + Random_TBOT_ARMA[Ens_Index,Time_Step,:,:]
                    PRECTmms_Temp = PRECTmms_Whole[Time_Step,::] * Random_Prec_ARMA[Ens_Index,Time_Step,:,:]
                    
                    SHUM[Time_Step,:,:] = SHUM_Whole[Time_Step,::]
                    
                    Prec_Copy = numpy.copy(PRECTmms_Whole[Time_Step,:,:])
                    
                    Prec_Threshold = 5.0 / 3600.0
                    if numpy.size(numpy.where((Prec_Copy-PRECTmms_Temp) > Prec_Threshold)) > 0:
                        PRECTmms_Temp[numpy.where((Prec_Copy-PRECTmms_Temp) > Prec_Threshold)] = Prec_Copy[numpy.where((Prec_Copy-PRECTmms_Temp) > Prec_Threshold)] - Prec_Threshold
                    if numpy.size(numpy.where((Prec_Copy-PRECTmms_Temp) < -1.0*Prec_Threshold)) > 0:
                        PRECTmms_Temp[numpy.where((Prec_Copy-PRECTmms_Temp) < -1.0*Prec_Threshold)] = Prec_Copy[numpy.where((Prec_Copy-PRECTmms_Temp) < -1.0*Prec_Threshold)] + Prec_Threshold
                    PRECTmms[Time_Step,:,:] = PRECTmms_Temp
                    #print numpy.where((Prec_Copy-PRECTmms_Temp) < -5.0)
                    
                    FSDS_Temp = FSDS_Whole[Time_Step,::] * Random_FSDS_ARMA[Ens_Index,Time_Step,:,:]
                    FSDS_Temp[numpy.where(FSDS_Temp > 1361.0)] = 1361.0
                    FSDS[Time_Step,:,:] = FSDS_Temp
                    
                    PSRF[Time_Step,:,:] = PSRF_Whole[Time_Step,::]
                    #numpy.random.seed(123456789 + (Time_Step + (Day_Index-1)*NTime_Out_Daily))
                    FLDS[Time_Step,:,:] = FLDS_Whole[Time_Step,::] + Random_FLDS_ARMA[Ens_Index,Time_Step,:,:]
                    
                    #numpy.random.seed(123456789 + (Time_Step + (Day_Index-1)*NTime_Out_Daily))
                    WIND[Time_Step,:,:] = WIND_Whole[Time_Step,::]
                    #print numpy.max(TBOT[Time_Step,:,:]),numpy.min(TBOT[Time_Step,:,:])
                    
            #        for Row_Index in range(Row_Numbers):
            #            for Col_Index in range(Col_Numbers):
            #                TBOT[:,Row_Index,Col_Index] = savitzky_golay(TBOT[:,Row_Index,Col_Index], window_size=31, order=4)
            #                SHUM[:,Row_Index,Col_Index] = savitzky_golay(SHUM[:,Row_Index,Col_Index], window_size=31, order=4)
            #                FSDS[:,Row_Index,Col_Index] = savitzky_golay(FSDS[:,Row_Index,Col_Index], window_size=31, order=4)
            #                FLDS[:,Row_Index,Col_Index] = savitzky_golay(FLDS[:,Row_Index,Col_Index], window_size=31, order=4)
                        
            #                fig1 = plt.figure()
            #                ax1 = fig1.add_subplot(1, 1, 1)
            #                x=numpy.arange(NTime_Out_Daily)
            #                ax1.plot(x,FSDS[:,Row_Index,Col_Index])
            #                ax1.plot(x,savitzky_golay(FSDS[:,Row_Index,Col_Index], window_size=31, order=4))
            #                ax1.set_title('TBOT')
            #                plt.show()
                
                for Time_Step in range(NTime_Out_Daily):
                    TBOT_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = TBOT_Whole[Time_Step,::][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    SHUM_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = SHUM_Whole[Time_Step,::][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    PRECTmms_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = PRECTmms_Whole[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    FSDS_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = FSDS_Whole[Time_Step,::][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    FLDS_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = FLDS_Whole[Time_Step,::][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    WIND_Series[Time_Step + (Day_Index-1)*NTime_Out_Daily] = WIND_Whole[Time_Step,::][Station_XY_Index[0][1],Station_XY_Index[0][0]]
            
                    TBOT_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = TBOT[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    SHUM_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = SHUM[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    PRECTmms_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = PRECTmms[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    FSDS_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = FSDS[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    FLDS_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = FLDS[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                    WIND_Ensemble[(Time_Step + (Day_Index-1)*NTime_Out_Daily),Ens_Index] = WIND[Time_Step,:,:][Station_XY_Index[0][1],Station_XY_Index[0][0]]
                
                CLM_NC_File_Out.sync()
                CLM_NC_File_Out.close()
                    
                for Time_Step in range(NTime_Out_Daily):
                    PRECTmms_Ensemble_Multiplier[(Time_Step + (Day_Index-1)*NTime_Out_Daily)] = numpy.mean(Random_Prec_ARMA[:,Time_Step])
                  
        #    fig1 = plt.figure()
        #    ax1 = fig1.add_subplot(1, 1, 1)
        #    x=numpy.arange(NTime_Out_Daily)
        #    for Ens_Index in range(Ensemble_Number):
        #        #ax1.plot(x,TBOT_Ensemble[0:24,Ens_Index])
        #        ax1.plot(x,savitzky_golay(TBOT_Ensemble[0:24,Ens_Index], window_size=31, order=4))
        #    ax1.set_title('TBOT')
        #    plt.show()
            time.sleep(3)
            del Random_TBOT_ARMA, Random_FSDS_ARMA, Random_FLDS_ARMA, Random_Prec_ARMA
            gc.collect()
            del gc.garbage[:]
            
            
        Start_Day = 1
    Start_Month = 1
            
                    
            
numpy.savetxt("../../Analysis/DAS_Temp/TBOT_Ensemble.txt",TBOT_Ensemble)
numpy.savetxt("../../Analysis/DAS_Temp/SHUM_Ensemble.txt",SHUM_Ensemble)
numpy.savetxt("../../Analysis/DAS_Temp/PRECTmms_Ensemble.txt",PRECTmms_Ensemble)
numpy.savetxt("../../Analysis/DAS_Temp/FSDS_Ensemble.txt",FSDS_Ensemble)
numpy.savetxt("../../Analysis/DAS_Temp/FLDS_Ensemble.txt",FLDS_Ensemble)
numpy.savetxt("../../Analysis/DAS_Temp/TBOT_Ensemble.txt",WIND_Ensemble)

numpy.savetxt("../../ForcingData/Ensembles/"+Region_Name+ "_PRECTmms_Ensemble_Multiplier_"+Year_String+".txt",PRECTmms_Ensemble_Multiplier)

x=numpy.arange(N_Steps)
fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
for Ens_Index in range(Ensemble_Number):
    ax1.plot(x,TBOT_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
ax1.set_title('TBOT')
plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_TBOT_Ensembles_"+Year_String+".png")

#fig2 = plt.figure()
#ax1 = fig2.add_subplot(1, 1, 1)
#for Ens_Index in range(Ensemble_Number):
#    ax1.plot(x,SHUM_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
#ax1.set_title('SHUM')
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_SHUM_Ensembles_"+Year_String+".png")

fig3 = plt.figure()
ax1 = fig3.add_subplot(1, 1, 1)
for Ens_Index in range(Ensemble_Number):
    ax1.plot(x,PRECTmms_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
ax1.set_title('PRECTmms')
plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_PRECTmms_Ensembles_"+Year_String+".png")

fig4 = plt.figure()
ax1 = fig4.add_subplot(1, 1, 1)
for Ens_Index in range(Ensemble_Number):
    ax1.plot(x,FSDS_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
ax1.set_title('FSDS')
plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_FSDS_Ensembles_"+Year_String+".png")

#fig5 = plt.figure()
#ax1 = fig5.add_subplot(1, 1, 1)
#for Ens_Index in range(Ensemble_Number):
#    ax1.plot(x,WIND_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
#ax1.set_title('WIND')
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_WIND_Ensembles_"+Year_String+".png")

fig6 = plt.figure()
ax1 = fig6.add_subplot(1, 1, 1)
for Ens_Index in range(Ensemble_Number):
    ax1.plot(x,FLDS_Ensemble[:,Ens_Index],color='k', linewidth=1 ,alpha=4.0/float(Ensemble_Number))
ax1.set_title('FLDS')
plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_FLDS_Ensembles_"+Year_String+".png")

#
#Index = numpy.where(TBOT_Series>200)
#
#x=numpy.arange(N_Steps)
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x[Index],TBOT_Series[Index])
#ax1.plot(x[Index],TBOT_Series_Filter[Index])
#ax1.plot(x[Index],numpy.mean(TBOT_Ensemble[:,:],axis=1)[Index])
#ax1.set_title('TBOT')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Filter', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_TBOT_Mean_Filter_"+Year_String+".png")
#
#x=numpy.arange(N_Steps)
#fig1 = plt.figure()
#ax1 = fig1.add_subplot(1, 1, 1)
#ax1.plot(x,TBOT_Series)
#ax1.plot(x,numpy.mean(TBOT_Ensemble[:,:],axis=1))
#ax1.set_title('TBOT')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Filter', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_TBOT_Mean_"+Year_String+".png")
#
#fig2 = plt.figure()
#ax1 = fig2.add_subplot(1, 1, 1)
#ax1.plot(x,SHUM_Series)
#ax1.plot(x,numpy.mean(SHUM_Ensemble[:,:],axis=1))
#ax1.set_title('SHUM')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_SHUM_Mean_"+Year_String+".png")
#
#fig3 = plt.figure()
#ax1 = fig3.add_subplot(1, 1, 1)
#ax1.plot(x,PRECTmms_Series)
#ax1.plot(x,numpy.mean(PRECTmms_Ensemble[:,:],axis=1))
#ax1.set_title('PRECTmms')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_PRECTmms_Mean_"+Year_String+".png")
#
#fig4 = plt.figure()
#ax1 = fig4.add_subplot(1, 1, 1)
#ax1.plot(x,FSDS_Series)
#ax1.plot(x,numpy.mean(FSDS_Ensemble[:,:],axis=1))
#ax1.set_title('FSDS')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_FSDS_Mean_"+Year_String+".png")
#
#fig5 = plt.figure()
#ax1 = fig5.add_subplot(1, 1, 1)
#ax1.plot(x,WIND_Series)
#ax1.plot(x,numpy.mean(WIND_Ensemble[:,:],axis=1))
#ax1.set_title('WIND')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_WIND_Mean_"+Year_String+".png")
#
#fig6 = plt.figure()
#ax1 = fig6.add_subplot(1, 1, 1)
#ax1.plot(x,FLDS_Series)
#ax1.plot(x,numpy.mean(FLDS_Ensemble[:,:],axis=1))
#ax1.set_title('FLDS')
#prop = fm.FontProperties(size=16) 
#legend(('True', 'Mean'), 'upper right', shadow=True,prop=prop)
#plt.savefig("../../ForcingData/Ensembles/"+Region_Name+ "_FLDS_Mean_"+Year_String+".png")

plt.show()

    
    