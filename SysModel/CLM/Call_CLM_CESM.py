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
1. Han, X., Li, X., He, G., Kumbhar, P., Montzka, C., Kollet, S., Miyoshi, T., Rosolem, R., Zhang, Y., Vereecken, H., and Franssen, H. J. H.: DasPy 1.0 &ndash; the Open Source Multivariate Land Data Assimilation Framework in combination with the Community Land Model 4.5, Geosci. Model Dev. Discuss., 8, 7395-7444, 2015.
2. Han, X., Franssen, H. J. H., Rosolem, R., Jin, R., Li, X., and Vereecken, H.: Correction of systematic model forcing bias of CLM using assimilation of cosmic-ray Neutrons and land surface temperature: a study in the Heihe Catchment, China, Hydrology and Earth System Sciences, 19, 615-629, 2015a.
3. Han, X., Franssen, H. J. H., Montzka, C., and Vereecken, H.: Soil moisture and soil properties estimation in the Community Land Model with synthetic brightness temperature observations, Water Resour Res, 50, 6081-6105, 2014a.
4. Han, X., Franssen, H. J. H., Li, X., Zhang, Y. L., Montzka, C., and Vereecken, H.: Joint Assimilation of Surface Temperature and L-Band Microwave Brightness Temperature in Land Data Assimilation, Vadose Zone J, 12, 0, 2013.
'''

import os, sys, time, datetime, calendar, subprocess, string, signal, socket, imp
import numpy

def Write_seq_maps(seq_maps_file_name, DAS_Data_Path, Row_Numbers_String, Col_Numbers_String, Region_Name):
    seq_maps_file = open(seq_maps_file_name,'w')
    seq_maps_file.write("##################################################################\n")
    seq_maps_file.write("#\n")
    seq_maps_file.write("# seq_maps.rc\n")
    seq_maps_file.write("#\n")
    seq_maps_file.write("# This is a resource file which lists the names of mapping\n")
    seq_maps_file.write("# weight files to use in a sequential CCSM run (mapname).\n")
    seq_maps_file.write("# You can also set when data is rearranged in the mapping (maptype).\n")
    seq_maps_file.write("#\n")
    seq_maps_file.write("# This file is read during the map_model2model_init calls.\n")
    seq_maps_file.write("#\n")
    seq_maps_file.write("# For maptype:  X = Rearrange the input so that the output\n")
    seq_maps_file.write("#                   is on the correct processor.\n")
    seq_maps_file.write("#               Y = Rearrange the output and sum partial outputs\n")
    seq_maps_file.write("#                   if necessary\n")
    seq_maps_file.write("#\n")
    seq_maps_file.write("# NOTE:  For bfb on different processor counts, set all maptypes to X.\n")
    seq_maps_file.write("################################################################## \n")
    seq_maps_file.write("atm2ice_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("atm2ice_fmaptype: 'X'\n")
    seq_maps_file.write("atm2ice_smapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_blin.nc")+"\n")
    seq_maps_file.write("atm2ice_smaptype: 'X'\n")
    seq_maps_file.write("atm2ice_vmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_patc.nc")+"\n")
    seq_maps_file.write("atm2ice_vmaptype: 'X'\n")
    seq_maps_file.write("atm2lnd_fmapname: 'idmap'\n")
    seq_maps_file.write("atm2lnd_fmaptype: 'X'\n")
    seq_maps_file.write("atm2lnd_smapname: 'idmap'\n")
    seq_maps_file.write("atm2lnd_smaptype: 'X'\n")
    seq_maps_file.write("atm2ocn_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("atm2ocn_fmaptype: 'X'\n")
    seq_maps_file.write("atm2ocn_smapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_blin.nc")+"\n")
    seq_maps_file.write("atm2ocn_smaptype: 'X'\n")
    seq_maps_file.write("atm2ocn_vmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_patc.nc")+"\n")
    seq_maps_file.write("atm2ocn_vmaptype: 'X'\n")
    seq_maps_file.write("atm2wav_smapname: 'idmap'\n")
    seq_maps_file.write("atm2wav_smaptype: 'Y'\n")
    seq_maps_file.write("ice2atm_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("ice2atm_fmaptype: 'Y'\n")
    seq_maps_file.write("ice2atm_smapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("ice2atm_smaptype: 'Y'\n")
    seq_maps_file.write("ice2wav_smapname: 'idmap'\n")
    seq_maps_file.write("ice2wav_smaptype: 'Y'\n")
    seq_maps_file.write("lnd2atm_fmapname: 'idmap'\n")
    seq_maps_file.write("lnd2atm_fmaptype: 'Y'\n")
    seq_maps_file.write("lnd2atm_smapname: 'idmap'\n")
    seq_maps_file.write("lnd2atm_smaptype: 'Y'\n")
    seq_maps_file.write("lnd2rof_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("lnd2rof_fmaptype: 'X'\n")
    seq_maps_file.write("ocn2atm_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("ocn2atm_fmaptype: 'Y'\n")
    seq_maps_file.write("ocn2atm_smapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("ocn2atm_smaptype: 'Y'\n")
    seq_maps_file.write("ocn2wav_smapname: 'idmap'\n")
    seq_maps_file.write("ocn2wav_smaptype: 'Y'\n")
    seq_maps_file.write("rof2lnd_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("rof2lnd_fmaptype: 'Y'\n")
    seq_maps_file.write("rof2lnd_smapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("rof2lnd_smaptype: 'Y'\n")
    seq_maps_file.write("rof2ocn_fmapname: "+repr(DAS_Data_Path+"/SysModel/CLM/tools/map_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_TO_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+"_aave.nc")+"\n")
    seq_maps_file.write("rof2ocn_fmaptype: 'Y'\n")
    seq_maps_file.write("rof2ocn_rmapname: ' '\n")
    seq_maps_file.write("rof2ocn_rmaptype: 'Y'\n")
    seq_maps_file.write("wav2ocn_smapname: 'idmap'\n")
    seq_maps_file.write("wav2ocn_smaptype: 'X'\n")
    seq_maps_file.write("/\n")
    seq_maps_file.close()

def Write_datm_atm_in(datm_atm_in_file_name, datm_streams_txt_file_name, presaero_stream_txt_file_name, domain_file_path,domain_name, rdirc_name,align_year,first_year,last_year):
    datm_atm_in_file = open(datm_atm_in_file_name,'w')
    datm_atm_in_file.write("&shr_strdata_nml\n")
    datm_atm_in_file.write(" dataMode = 'CLMNCEP'\n")
    datm_atm_in_file.write(" domainFile = "+repr(domain_file_path+domain_name)+ "\n")
    datm_atm_in_file.write(" dtlimit = 1.0e30,1.0e30\n")
    datm_atm_in_file.write(" fillalgo = 'copy','copy'\n")
    datm_atm_in_file.write(" fillmask = 'nomask','nomask'\n")
    datm_atm_in_file.write(" mapalgo  = 'bilinear','bilinear'\n")
    datm_atm_in_file.write(" mapmask  = 'nomask','nomask'\n")
    datm_atm_in_file.write(" streams  = "+repr(datm_streams_txt_file_name + " " +align_year+" "+first_year+" "+last_year+" ")+","+"\n")
    datm_atm_in_file.write("            "+repr(presaero_stream_txt_file_name +" 1 1 1")+"\n")
    datm_atm_in_file.write(" taxMode  = 'cycle','cycle'\n")
    datm_atm_in_file.write(" tintalgo = 'linear','linear'\n")
    datm_atm_in_file.write(" vectors  = 'null'\n")
    datm_atm_in_file.write("/\n")
    datm_atm_in_file.close()

def Write_datm_streams_txt(datm_streams_txt_file_name, Def_SpinUp, domain_file_path,domain_name, rdirc_name,forcing_file_path,start_ymd,stop_ymd):
    datm_streams_txt = open(datm_streams_txt_file_name,'w')
    datm_streams_txt.write("      <dataSource>\n")
    datm_streams_txt.write("         GENERIC\n")
    datm_streams_txt.write("      </dataSource>\n")
    datm_streams_txt.write("      <domainInfo>\n")
    datm_streams_txt.write("         <variableNames>\n")
    datm_streams_txt.write("            time    time\n")
    datm_streams_txt.write("            xc      lon\n")
    datm_streams_txt.write("            yc      lat\n")
    datm_streams_txt.write("            area    area\n")
    datm_streams_txt.write("            mask    mask\n")
    datm_streams_txt.write("         </variableNames>\n")
    datm_streams_txt.write("         <filePath>\n")
    datm_streams_txt.write("            "+domain_file_path+"\n")
    datm_streams_txt.write("         </filePath>\n")
    datm_streams_txt.write("         <fileNames>\n")
    datm_streams_txt.write("            "+domain_name+"\n")
    datm_streams_txt.write("         </fileNames>\n")
    datm_streams_txt.write("      </domainInfo>\n")
    datm_streams_txt.write("      <fieldInfo>\n")
    datm_streams_txt.write("         <variableNames>\n")
    datm_streams_txt.write("            PRECTmms precn\n")
    datm_streams_txt.write("            FSDS     swdn\n")
    datm_streams_txt.write("            FLDS     lwdn\n")
    datm_streams_txt.write("            TBOT     tbot\n")
    datm_streams_txt.write("            WIND     wind\n")
    datm_streams_txt.write("            SHUM     shum\n")
    datm_streams_txt.write("            PSRF     pbot\n")
    datm_streams_txt.write("            ZBOT     z\n")
    datm_streams_txt.write("         </variableNames>\n")
    datm_streams_txt.write("         <filePath>\n")
    datm_streams_txt.write("            "+forcing_file_path+"\n")
    datm_streams_txt.write("         </filePath>\n")
    datm_streams_txt.write("         <fileNames>\n")
    
    #print start_ymd,str.split(start_ymd),stop_ymd,str.split(stop_ymd)
    #print str.split(start_ymd)[0][0:4],str.split(start_ymd)[0][4:6],str.split(start_ymd)[0][6:8]
    
    Datetime_Start = datetime.datetime(string.atoi(str.split(start_ymd)[0][0:4]), string.atoi(str.split(start_ymd)[0][4:6]), string.atoi(str.split(start_ymd)[0][6:8]), 00, 00)
    Datetime_Stop = datetime.datetime(string.atoi(str.split(stop_ymd)[0][0:4]), string.atoi(str.split(stop_ymd)[0][4:6]), string.atoi(str.split(stop_ymd)[0][6:8]), 00, 00)

    #print Datetime_Start,Datetime_Stop,(Datetime_Stop - Datetime_Start).days
        
            
    Datetime_Stop_Temp = Datetime_Start
    while Datetime_Stop_Temp <= Datetime_Stop:
        
#            if calendar.isleap(Datetime_Stop_Temp.year):
#                if Datetime_Stop_Temp == datetime.datetime(Datetime_Stop_Temp.year,2,29):
#                    # Add the Delta Days
#                    data_delta = datetime.timedelta(days=1)
#                    Datetime_Stop_Temp = Datetime_Stop_Temp + data_delta
                
        if Datetime_Stop_Temp.day < 10:
            Day_String = '00'+str(Datetime_Stop_Temp.day)
        elif Datetime_Stop_Temp.day > 10 and Datetime_Stop_Temp.day < 100:
            Day_String = '0'+str(Datetime_Stop_Temp.day)
        else:
            Day_String = str(Datetime_Stop_Temp.day)
        
        # If the Month or Day is one-digit number, then add '0' in front of it.
        if len(str(Datetime_Stop_Temp.month)) == 1:
            Month = '0' + str(Datetime_Stop_Temp.month)
        else:
            Month = str(Datetime_Stop_Temp.month)
        if len(str(Datetime_Stop_Temp.day)) == 1:
            Day = '0' + str(Datetime_Stop_Temp.day)
        else:
            Day = str(Datetime_Stop_Temp.day)
        
        # MODIS Date
        #print Datetime_Stop_Temp.year,Datetime_Stop_Temp.month,Datetime_Stop_Temp.day
        #print str(Datetime_Stop_Temp.year),Month,Day
        OutputDate = str(Datetime_Stop_Temp.year) + '-' + Month + '-' + Day
        
        NC_FileName = str(Datetime_Stop_Temp.year) + '-' + Month + '-' + Day + ".nc"                
        datm_streams_txt.write("            "+NC_FileName+"\n")
            
        # Add the Delta Days
        data_delta = datetime.timedelta(days=1)
        Datetime_Stop_Temp = Datetime_Stop_Temp + data_delta
                         
                    
                
    datm_streams_txt.write("         </fileNames>\n")
    datm_streams_txt.write("         <offset>\n")
    datm_streams_txt.write("            0\n")
    datm_streams_txt.write("         </offset>\n")
    datm_streams_txt.write("      </fieldInfo>\n")
    datm_streams_txt.write("\n")
    datm_streams_txt.close()

def Write_presaero_stream_txt(presaero_stream_txt_file_name,aero_file_path,aero_file_name):
    presaero_stream_txt = open(presaero_stream_txt_file_name,'w')
    presaero_stream_txt.write("      <dataSource>\n")
    presaero_stream_txt.write("         GENERIC\n")
    presaero_stream_txt.write("      </dataSource>\n")
    presaero_stream_txt.write("      <domainInfo>\n")
    presaero_stream_txt.write("         <variableNames>\n")
    presaero_stream_txt.write("            time    time\n")
    presaero_stream_txt.write("            lon      lon\n")
    presaero_stream_txt.write("            lat      lat\n")
    presaero_stream_txt.write("            area    area\n")
    presaero_stream_txt.write("            mask    mask\n")
    presaero_stream_txt.write("         </variableNames>\n")
    presaero_stream_txt.write("         <filePath>\n")
    presaero_stream_txt.write("            "+aero_file_path+"\n")
    presaero_stream_txt.write("         </filePath>\n")
    presaero_stream_txt.write("         <fileNames>\n")
    presaero_stream_txt.write("            "+aero_file_name+"\n")
    presaero_stream_txt.write("         </fileNames>\n")
    presaero_stream_txt.write("      </domainInfo>\n")
    presaero_stream_txt.write("      <fieldInfo>\n")
    presaero_stream_txt.write("         <variableNames>\n")
    presaero_stream_txt.write("            BCDEPWET   bcphiwet\n")
    presaero_stream_txt.write("            BCPHODRY   bcphodry\n")
    presaero_stream_txt.write("            BCPHIDRY   bcphidry\n")
    presaero_stream_txt.write("            OCDEPWET   ocphiwet\n")
    presaero_stream_txt.write("            OCPHIDRY   ocphidry\n")
    presaero_stream_txt.write("            OCPHODRY   ocphodry\n")
    presaero_stream_txt.write("            DSTX01WD   dstwet1\n")
    presaero_stream_txt.write("            DSTX01DD   dstdry1\n")
    presaero_stream_txt.write("            DSTX02WD   dstwet2\n")
    presaero_stream_txt.write("            DSTX02DD   dstdry2\n")
    presaero_stream_txt.write("            DSTX03WD   dstwet3\n")
    presaero_stream_txt.write("            DSTX03DD   dstdry3\n")
    presaero_stream_txt.write("            DSTX04WD   dstdry4\n")
    presaero_stream_txt.write("            DSTX04DD   dstwet4\n")
    presaero_stream_txt.write("         </variableNames>\n")
    presaero_stream_txt.write("         <filePath>\n")
    presaero_stream_txt.write("            "+aero_file_path+"\n")
    presaero_stream_txt.write("         </filePath>\n")
    presaero_stream_txt.write("         <offset>\n")
    presaero_stream_txt.write("            0\n")
    presaero_stream_txt.write("         </offset>\n")
    presaero_stream_txt.write("         <fileNames>\n")
    presaero_stream_txt.write("            "+aero_file_name+"\n")
    presaero_stream_txt.write("         </fileNames>\n")
    presaero_stream_txt.write("         <offset>\n")
    presaero_stream_txt.write("            0\n")
    presaero_stream_txt.write("         </offset>\n")
    presaero_stream_txt.write("      </fieldInfo>\n")
    presaero_stream_txt.write("\n")
    presaero_stream_txt.close()

def Write_drv_in(Def_PP, Model_Driver, Def_CESM_Multi_Instance,Ensemble_Number,num_processors,case_name,hostname,orb_iyear,start_type,username,
                 atm_cpl_dt,lnd_cpl_dt,ocn_cpl_dt,ice_cpl_dt,glc_cpl_dt, rof_cpl_dt, wav_cpl_dt,
                 end_restart,restart_option,start_tod,start_ymd,stop_tod,stop_ymd,ntasks_CLM,rootpe_CLM,nthreads_CLM):
    drv_in_file = open("drv_in",'w')
    drv_in_file.write("&seq_cplflds_inparm\n")
    drv_in_file.write(" flds_co2_dmsa = .false.\n")
    drv_in_file.write(" flds_co2a = .false.\n")
    drv_in_file.write(" flds_co2b = .false.\n")
    drv_in_file.write(" flds_co2c = .false.\n")
    drv_in_file.write(" glc_nec = 0\n")
    drv_in_file.write("/\n")
    drv_in_file.write("&seq_cplflds_userspec\n")
    drv_in_file.write(" cplflds_custom = ''\n")
    drv_in_file.write("/\n")

    drv_in_file.write("&seq_infodata_inparm\n")
    drv_in_file.write(" aoflux_grid = 'ocn'\n")
    drv_in_file.write(" bfbflag = .false.\n")
    drv_in_file.write(" brnch_retain_casename = .false.\n")
    drv_in_file.write(" budget_ann = 1\n")
    drv_in_file.write(" budget_daily = 0\n")
    drv_in_file.write(" budget_inst = 0\n")
    drv_in_file.write(" budget_ltann = 1\n")
    drv_in_file.write(" budget_ltend = 0\n")
    drv_in_file.write(" budget_month = 1\n")
    drv_in_file.write(" case_desc = 'UNSET'\n")
    drv_in_file.write(" case_name = "+case_name+"\n")
    drv_in_file.write(" cpl_cdf64 = .true.\n")
    drv_in_file.write(" cpl_decomp = 0\n")
    drv_in_file.write(" do_budgets = .false.\n")
    drv_in_file.write(" do_histinit = .false.\n")
    drv_in_file.write(" drv_threading = .false.\n")
    drv_in_file.write(" eps_aarea = 9.0e-07\n")
    drv_in_file.write(" eps_agrid = 1.0e-12\n")
    drv_in_file.write(" eps_amask = 1.0e-13\n")
    drv_in_file.write(" eps_frac = 1.0e-02\n")
    drv_in_file.write(" eps_oarea = 1.0e-01\n")
    drv_in_file.write(" eps_ogrid = 1.0e-02\n")
    drv_in_file.write(" eps_omask = 1.0e-06\n")
    drv_in_file.write(" flux_albav = .false.\n")
    drv_in_file.write(" flux_epbal = 'off'\n")
    drv_in_file.write(" histaux_a2x = .false.\n")
    drv_in_file.write(" histaux_a2x24hr = .false.\n")
    drv_in_file.write(" histaux_a2x3hr = .false.\n")
    drv_in_file.write(" histaux_a2x3hrp = .false.\n")
    drv_in_file.write(" histaux_l2x = .false.\n")
    drv_in_file.write(" histaux_r2x = .false.\n")
    drv_in_file.write(" histaux_s2x1yr = .false.\n")
    drv_in_file.write(" hostname = "+hostname+"\n")
    drv_in_file.write(" info_debug = 1\n")
    drv_in_file.write(" mct_usealltoall = .false.\n")
    drv_in_file.write(" mct_usevector = .false.\n")
    drv_in_file.write(" model_version = 'cesm1_2_1'\n")
    drv_in_file.write(" ocean_tight_coupling = .false.\n")
    drv_in_file.write(" orb_iyear = "+orb_iyear+"\n")
    drv_in_file.write(" orb_iyear_align = "+orb_iyear+"\n")
    drv_in_file.write(" orb_mode = 'fixed_year'\n")
    drv_in_file.write(" run_barriers = .false.\n")
    drv_in_file.write(" samegrid_al = .true.\n")
    drv_in_file.write(" samegrid_ao = .false.\n")
    drv_in_file.write(" samegrid_aw = .false.\n")
    drv_in_file.write(" samegrid_ow = .false.\n")
    drv_in_file.write(" samegrid_ro = .false.\n")
    drv_in_file.write(" shr_map_dopole = .true.\n")
    drv_in_file.write(" start_type = "+start_type+"\n")
    drv_in_file.write(" tchkpt_dir = './timing/checkpoints'\n")
    drv_in_file.write(" timing_dir = './timing'\n")
    drv_in_file.write(" username = "+username+"\n")
    drv_in_file.write(" vect_map = 'cart3d'\n")
    drv_in_file.write("/\n")
    
    drv_in_file.write("&seq_timemgr_inparm\n")
    drv_in_file.write(" atm_cpl_dt = "+str(atm_cpl_dt)+"\n")
    drv_in_file.write(" calendar = 'GREGORIAN'\n")
    #drv_in_file.write(" calendar = 'NO_LEAP'\n")
    drv_in_file.write(" end_restart = "+end_restart+"\n")
    drv_in_file.write(" glc_cpl_dt = "+str(glc_cpl_dt)+"\n")
    drv_in_file.write(" histavg_n = -999\n")
    drv_in_file.write(" histavg_option = 'never'\n")
    drv_in_file.write(" histavg_ymd = -999\n")
    drv_in_file.write(" history_n = -999\n")
    drv_in_file.write(" history_option = 'never'\n")
    drv_in_file.write(" history_ymd = -999\n")
    drv_in_file.write(" ice_cpl_dt = "+str(ice_cpl_dt)+"\n")
    drv_in_file.write(" lnd_cpl_dt = "+str(lnd_cpl_dt)+"\n")
    drv_in_file.write(" ocn_cpl_dt = "+str(ocn_cpl_dt)+"\n")
    drv_in_file.write(" rof_cpl_dt = "+str(rof_cpl_dt)+"\n")
    drv_in_file.write(" start_tod = "+start_tod+"\n")
    drv_in_file.write(" start_ymd = "+start_ymd+"\n")
    drv_in_file.write(" stop_option = 'date'\n")
    drv_in_file.write(" stop_tod = "+stop_tod+"\n")
    drv_in_file.write(" stop_ymd = "+stop_ymd+"\n")
    drv_in_file.write(" end_restart = "+end_restart+"\n")
    drv_in_file.write(" tprof_n = -999\n")
    drv_in_file.write(" tprof_option = 'never'\n")
    drv_in_file.write(" tprof_ymd = -999\n")
    drv_in_file.write(" wav_cpl_dt = "+str(wav_cpl_dt)+"\n")
    drv_in_file.write("/\n")
    
    drv_in_file.write("&ccsm_pes\n")
    drv_in_file.write(" atm_layout = 'concurrent'\n")
    drv_in_file.write(" atm_ntasks = "+str(int(ntasks_CLM[0]))+"\n")
    drv_in_file.write(" atm_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" atm_pestride = 1\n")
    drv_in_file.write(" atm_rootpe = "+str(int(rootpe_CLM[0]))+"\n")
    drv_in_file.write(" lnd_layout = 'concurrent'\n")
    drv_in_file.write(" lnd_ntasks = "+str(int(ntasks_CLM[1]))+"\n")
    drv_in_file.write(" lnd_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" lnd_pestride = 1\n")
    drv_in_file.write(" lnd_rootpe = "+str(int(rootpe_CLM[1]))+"\n")
    drv_in_file.write(" cpl_ntasks = "+str(int(ntasks_CLM[2]))+"\n")
    drv_in_file.write(" cpl_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" cpl_pestride = 1\n")
    drv_in_file.write(" cpl_rootpe = "+str(int(rootpe_CLM[2]))+"\n")
    drv_in_file.write(" glc_layout = 'concurrent'\n")
    drv_in_file.write(" glc_ntasks = "+str(int(ntasks_CLM[3]))+"\n")
    drv_in_file.write(" glc_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" glc_pestride = 1\n")
    drv_in_file.write(" glc_rootpe = "+str(int(rootpe_CLM[3]))+"\n")
    drv_in_file.write(" ice_layout = 'concurrent'\n")
    drv_in_file.write(" ice_ntasks = "+str(int(ntasks_CLM[4]))+"\n")
    drv_in_file.write(" ice_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" ice_pestride = 1\n")
    drv_in_file.write(" ice_rootpe = "+str(int(rootpe_CLM[4]))+"\n")
    drv_in_file.write(" ocn_layout = 'concurrent'\n")
    drv_in_file.write(" ocn_ntasks = "+str(int(ntasks_CLM[5]))+"\n")
    drv_in_file.write(" ocn_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" ocn_pestride = 1\n")
    drv_in_file.write(" ocn_rootpe = "+str(int(rootpe_CLM[5]))+"\n")
    drv_in_file.write(" rof_layout = 'concurrent'\n")
    drv_in_file.write(" rof_ntasks = "+str(int(ntasks_CLM[6]))+"\n")
    drv_in_file.write(" rof_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" rof_pestride = 1\n")
    drv_in_file.write(" rof_rootpe = "+str(int(rootpe_CLM[6]))+"\n")
    drv_in_file.write(" wav_layout = 'concurrent'\n")
    drv_in_file.write(" wav_ntasks = "+str(int(ntasks_CLM[7]))+"\n")
    drv_in_file.write(" wav_nthreads = "+str(int(nthreads_CLM))+"\n")
    drv_in_file.write(" wav_pestride = 1\n")
    drv_in_file.write(" wav_rootpe = "+str(int(rootpe_CLM[7]))+"\n")
    drv_in_file.write("/\n")
    
    drv_in_file.write("&prof_inparm\n")
    drv_in_file.write(" profile_barrier = .false.\n")
    drv_in_file.write(" profile_depth_limit = 12\n")
    drv_in_file.write(" profile_detail_limit = 0\n")
    drv_in_file.write(" profile_disable = .false.\n")
    drv_in_file.write(" profile_global_stats = .false.\n")
    drv_in_file.write(" profile_single_file = .false.\n")
    if Def_PP == 2 or Def_CESM_Multi_Instance:
        drv_in_file.write(" profile_timer = 4\n")
    else:
        drv_in_file.write(" profile_timer = 1\n")
    drv_in_file.write("/\n")
    
    drv_in_file.write("&pio_default_inparm\n")
    drv_in_file.write(" pio_async_interface = .false.\n")
    if Def_PP or Def_CESM_Multi_Instance:
        drv_in_file.write(" pio_blocksize = -1\n")
        drv_in_file.write(" pio_buffer_size_limit = -1\n")
        drv_in_file.write(" pio_debug_level = 0\n")
        drv_in_file.write(" pio_numiotasks = 1\n")  # only 1 works for netcdf4c
        drv_in_file.write(" pio_root = 1\n")
        drv_in_file.write(" pio_stride = 1\n")  # only 1 works for netcdf4c
        drv_in_file.write(" pio_typename = 'netcdf4c'\n")
    else:
        drv_in_file.write(" pio_blocksize = -1\n")
        drv_in_file.write(" pio_buffer_size_limit = -1\n")
        drv_in_file.write(" pio_debug_level = 0\n")
        drv_in_file.write(" pio_numiotasks = -1\n")
        drv_in_file.write(" pio_root = 1\n")
        drv_in_file.write(" pio_stride = -1\n")
        drv_in_file.write(" pio_typename = 'netcdf'\n")
    drv_in_file.write("/\n")
    
    drv_in_file.close()

def Write_45_drv_flds_in(drv_flds_in_file_name, megan_factors_file_path, megan_factors_file_name):
    drv_flds_in_file = open(drv_flds_in_file_name,'w')
    drv_flds_in_file.write("&drydep_inparm\n")
    drv_flds_in_file.write("/\n")
    drv_flds_in_file.write("&megan_emis_nl")
    drv_flds_in_file.write(" megan_factors_file = "+repr(megan_factors_file_path+megan_factors_file_name)+ "\n")
    drv_flds_in_file.write(" megan_specifier = 'ISOP = isoprene', 'C10H16 = pinene_a + carene_3 + thujene_a', 'CH3OH = methanol', 'C2H5OH = ethanol',")
    drv_flds_in_file.write("         'CH2O = formaldehyde', 'CH3CHO = acetaldehyde', 'CH3COOH = acetic_acid', 'CH3COCH3 = acetone'")
    drv_flds_in_file.write("/")
    drv_flds_in_file.close()

def Write_rof_in(rof_in_file_name,frivinp_rtm_path,frivinp_rtm_name):
#    rof_in_file = open(rof_in_file_name,'w')
#    rof_in_file.write("&rtm_inparm\n")
#    rof_in_file.write(" finidat_rtm = ' '\n")
#    rof_in_file.write(" flood_mode = 'NULL'\n")
#    rof_in_file.write(" frivinp_rtm = "+repr(frivinp_rtm_path+frivinp_rtm_name)+ "\n")
#    rof_in_file.write(" ice_runoff = .true.\n")
#    rof_in_file.write(" rtm_tstep = 10800\n")
#    rof_in_file.write(" rtmhist_mfilt = 30\n")
#    rof_in_file.write(" rtmhist_ndens = 2\n")
#    rof_in_file.write(" rtmhist_nhtfrq = 0\n")
#    rof_in_file.write("/\n")
#    rof_in_file.close()
    rof_in_file = open(rof_in_file_name,'w')
    rof_in_file.write("&rtm_inparm\n")
    rof_in_file.write(" rtm_effvel = 'ACTIVE'\n")
    rof_in_file.write(" rtm_mode = 'NULL'\n")
    rof_in_file.write("/\n")
    rof_in_file.close()

def Write_lnd_in(Run_Dir, lnd_in_file_name,Model_Driver,dtime,rtm_nsteps,domain_file_lnd_path, domain_name, rdirc_name,fatmgrid_name,fatmlndfrc_name,fglcmask_name,finidat_name,\
                    flndtopo_name,fndepdat_name,fpftcon_name,frivinp_rtm_name,fsnowaging_name,fsnowoptics_name,fsurdat_name, popd_streams_name, light_streams_name,\
                 wrtdia,hist_nhtfrq,hist_mfilt,hist_crtinic,hist_dov2xy,hist_ndens,hist_type1d_pertape,hist_empty_htapes,hist_avgflag_pertape,hist_fincl1,hist_fexcl1,first_year):
    lnd_in_file = open(lnd_in_file_name,'w')
    lnd_in_file.write("&clm_inparm\n")
    lnd_in_file.write(" albice = 0.60,0.40\n")
    lnd_in_file.write(" co2_ppmv = 367.0\n")
    lnd_in_file.write(" co2_type = 'constant'\n")
    lnd_in_file.write(" create_crop_landunit = .false.\n")
    lnd_in_file.write(" dtime = "+str(dtime)+"\n")
    #lnd_in_file.write(" fatmgrid_name = "+fatmgrid_name+"\n")
    lnd_in_file.write(" fatmlndfrc = "+repr(domain_file_lnd_path+fatmlndfrc_name)+"\n")
    #lnd_in_file.write(" fglcmask_name  = "+fglcmask_name+"\n")
    if finidat_name == "":
        lnd_in_file.write(" finidat = ''\n")
    else:
        lnd_in_file.write(" finidat = "+repr(Run_Dir+finidat_name)+"\n")
    #lnd_in_file.write(" flndtopo = "+repr(Run_Dir+flndtopo_name)+"\n")
    lnd_in_file.write(" fpftcon = "+repr(Run_Dir+fpftcon_name)+"\n")
    lnd_in_file.write(" fsnowaging = "+repr(Run_Dir+fsnowaging_name)+"\n")
    lnd_in_file.write(" fsnowoptics = "+repr(Run_Dir+fsnowoptics_name)+"\n")
    lnd_in_file.write(" fsurdat = "+repr(Run_Dir+fsurdat_name)+"\n")
    lnd_in_file.write(" maxpatch_glcmec = 0\n")
    if Model_Driver == "CLM_45":
        lnd_in_file.write(" more_vertlayers = .false.\n")
    lnd_in_file.write(" nsegspc = 1\n")             # default is 20, but for parflow it should be 1, then the decomposition is right to the row number
    lnd_in_file.write(" hist_nhtfrq =  " + str(hist_nhtfrq) + "\n")
    lnd_in_file.write(" hist_mfilt =  " + str(hist_mfilt) + "\n")
    #lnd_in_file.write(" hist_crtinic   = " + hist_crtinic + "\n")
    lnd_in_file.write(" hist_dov2xy = "+hist_dov2xy+"\n")
    lnd_in_file.write(" hist_ndens = " + str(hist_ndens) + "\n")
    lnd_in_file.write(" hist_type1d_pertape = " + hist_type1d_pertape + "\n")
    lnd_in_file.write(" hist_empty_htapes = " + hist_empty_htapes + "\n")
    lnd_in_file.write(" hist_avgflag_pertape = " + hist_avgflag_pertape + "\n")
    lnd_in_file.write(" hist_fincl1 = " + hist_fincl1 + "\n")
    lnd_in_file.write(" hist_fexcl1 = " + hist_fexcl1 + "\n")
    lnd_in_file.write(" outnc_large_files = .true.\n")
    
    if Model_Driver == "CLM_BGC_SpinUp":
        lnd_in_file.write(" spinup_state = 1\n")
    else:
        if Model_Driver == "CLM_BGC":
            lnd_in_file.write(" spinup_state = 0\n")
    
    if Model_Driver == "CLM_CN":
        lnd_in_file.write(" suplnitro = 'PROG_CROP_ONLY'\n")
    lnd_in_file.write(" urban_hac = 'ON'\n")
    lnd_in_file.write(" urban_traffic = .false.\n")
    lnd_in_file.write("/\n")
    lnd_in_file.write("&ndepdyn_nml\n")
    if Model_Driver == "CLM_CN" or Model_Driver == "cesm_ad_spinup":        
        lnd_in_file.write(" ndepmapalgo = 'bilinear'\n")
        lnd_in_file.write(" stream_fldfilename_ndep  = "+repr(Run_Dir+fndepdat_name)+"\n")
        lnd_in_file.write(" stream_year_first_ndep  = 1850\n")
        lnd_in_file.write(" stream_year_last_ndep  = 2005\n")
 
    lnd_in_file.write("/\n")
    lnd_in_file.write("&popd_streams\n")
    if Model_Driver == "CLM_CN" or Model_Driver == "cesm_ad_spinup":  
        lnd_in_file.write(" popdensmapalgo = 'bilinear'\n")
        lnd_in_file.write(" stream_fldfilename_popdens  = "+repr(Run_Dir+popd_streams_name)+"\n")
        lnd_in_file.write(" stream_year_first_popdens  = 1850\n")
        lnd_in_file.write(" stream_year_last_popdens  = 2010\n")
    lnd_in_file.write("/\n")
    lnd_in_file.write("&light_streams\n")
    if Model_Driver == "CLM_CN" or Model_Driver == "cesm_ad_spinup":  
        lnd_in_file.write(" lightngmapalgo = 'bilinear'\n")
        lnd_in_file.write(" stream_fldfilename_lightng  = "+repr(Run_Dir+light_streams_name)+"\n")
        lnd_in_file.write(" stream_year_first_lightng  = 0001\n")
        lnd_in_file.write(" stream_year_last_lightng  = 0001\n")
    lnd_in_file.write("/\n")
    
    lnd_in_file.write("&clm_hydrology1_inparm\n")
    lnd_in_file.write(" oldfflag = 0\n")    #Use old snow cover fraction from Niu et al. 2007
    lnd_in_file.write("/\n")
    lnd_in_file.write("&clm_soilhydrology_inparm\n")
    lnd_in_file.write(" h2osfcflag = 0\n")  #If surface water is active or not
    lnd_in_file.write(" origflag = 0\n")    #Use original CLM4 soil hydraulic properties
    lnd_in_file.write("/\n")
    if Model_Driver == "CLM_CN" or Model_Driver == "cesm_ad_spinup":  
        lnd_in_file.write(" &ch4par_in\n")
        lnd_in_file.write(" fin_use_fsat = .true.\n")
        lnd_in_file.write(" /\n")

    lnd_in_file.write("#!--------------------------------------------------------------------------------------------------------------------------\n")
    lnd_in_file.write("#! lnd_in:: Comment:\n")
    lnd_in_file.write("#! This namelist was created using the following command-line:\n")
    lnd_in_file.write("#!     /lustrefs/lzhpc84/Library/cesm1_2_0/models/lnd/clm/bld/CLM build-namelist -infile /lustrefs/lzhpc84/Library/cesm1_2_0/scripts/sp_clm_ens_2/Buildconf/clmconf/cesm_namelist -csmdata /lustrefs/lzhpc84/DAS_Data/SysModel/CLM/inputdata -inputdata /lustrefs/lzhpc84/Library/cesm1_2_0/scripts/sp_clm_ens_2/Buildconf/clm.input_data_list -ignore_ic_year -namelist &clm_inparm start_ymd = 00010101  / -use_case 2000_control -res 1.9x2.5 -clm_start_type startup -clm_startfile I2000CN_f19_g16_c100503.clm2.r.0001-01-01-00000.nc -l_ncpl 48 -lnd_frac /lustrefs/lzhpc84/DAS_Data/SysModel/CLM/inputdata/share/domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc -glc_nec 0 -co2_ppmv 367.0 -co2_type constant -config /lustrefs/lzhpc84/Library/cesm1_2_0/scripts/sp_clm_ens_2/Buildconf/clmconf/config_cache.xml\n")
    lnd_in_file.write("#! For help on options use: /lustrefs/lzhpc84/Library/cesm1_2_0/models/lnd/clm/bld/CLM build-namelist -help\n")
    lnd_in_file.write("#!--------------------------------------------------------------------------------------------------------------------------\n")
    lnd_in_file.close()

def Write_Config_Files(datm_in_file_name,datm_atm_in_file_name,atm_modelio_file_name,cpl_modelio_file_name,glc_modelio_file_name,ice_modelio_file_name,lnd_modelio_file_name,ocn_modelio_file_name,rof_modelio_file_name,
                       wav_modelio_file_name, logfile_atm, logfile_cpl, logfile_glc, logfile_ice, logfile_lnd, logfile_ocn, logfile_rof, logfile_wav, Run_Dir):
    
    if not os.path.exists("timing/checkpoints"):
        os.makedirs("timing/checkpoints")
    
    datm_in = open(datm_in_file_name,'w')
    datm_in.write("&datm_nml\n")
    datm_in.write(" atm_in = "+repr(datm_atm_in_file_name)+"\n")
    datm_in.write(" decomp = '1d'\n")
    datm_in.write(" iradsw   = 1\n")
    datm_in.write(" presaero = .true.\n")
    datm_in.write(" restfilm = 'undefined'\n")
    datm_in.write(" restfils = 'undefined'\n")    
    datm_in.write("  /\n")
    datm_in.close()
    
    atm_modelio = open(atm_modelio_file_name,'w')
    atm_modelio.write('&modelio\n')
    atm_modelio.write(' diri = "."\n')
    atm_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    atm_modelio.write(' logfile = '+repr(logfile_atm)+'\n')
    atm_modelio.write('/\n')
    atm_modelio.write('&pio_inparm\n')
    atm_modelio.write(' pio_numiotasks = -99\n')
    atm_modelio.write(' pio_root = -99\n')
    atm_modelio.write(' pio_stride = -99\n')
    atm_modelio.write(' pio_typename = "nothing"\n')
    atm_modelio.write('/\n')
    atm_modelio.close()

    cpl_modelio = open(cpl_modelio_file_name,'w')
    cpl_modelio.write('&modelio\n')
    cpl_modelio.write(' diri = "."\n')
    cpl_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    cpl_modelio.write(' logfile = '+repr(logfile_cpl)+'\n')
    cpl_modelio.write('/\n')
    cpl_modelio.write('&pio_inparm\n')
    cpl_modelio.write(' pio_numiotasks = -99\n')
    cpl_modelio.write(' pio_root = -99\n')
    cpl_modelio.write(' pio_stride = -99\n')
    cpl_modelio.write(' pio_typename = "nothing"\n')
    cpl_modelio.write('/\n')
    cpl_modelio.close()

    glc_modelio = open(glc_modelio_file_name,'w')
    glc_modelio.write('&modelio\n')
    glc_modelio.write(' diri = "."\n')
    glc_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    glc_modelio.write(' logfile = '+repr(logfile_glc)+'\n')
    glc_modelio.write('/\n')
    glc_modelio.write('&pio_inparm\n')
    glc_modelio.write(' pio_numiotasks = -99\n')
    glc_modelio.write(' pio_root = -99\n')
    glc_modelio.write(' pio_stride = -99\n')
    glc_modelio.write(' pio_typename = "nothing"\n')
    glc_modelio.write('/\n')
    glc_modelio.close()
    
    ice_modelio = open(ice_modelio_file_name,'w')
    ice_modelio.write('&modelio\n')
    ice_modelio.write(' diri = "."\n')
    ice_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    ice_modelio.write(' logfile = '+repr(logfile_ice)+'\n')
    ice_modelio.write('/\n')
    ice_modelio.write('&pio_inparm\n')
    ice_modelio.write(' pio_numiotasks = -99\n')
    ice_modelio.write(' pio_root = -99\n')
    ice_modelio.write(' pio_stride = -99\n')
    ice_modelio.write(' pio_typename = "nothing"\n')
    ice_modelio.write('/\n')
    ice_modelio.close()

    lnd_modelio = open(lnd_modelio_file_name,'w')
    lnd_modelio.write('&modelio\n')
    lnd_modelio.write(' diri = "."\n')
    lnd_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    lnd_modelio.write(' logfile = '+repr(logfile_lnd)+'\n')
    lnd_modelio.write('/\n')
    lnd_modelio.write('&pio_inparm\n')
    lnd_modelio.write(' pio_numiotasks = -99\n')
    lnd_modelio.write(' pio_root = -99\n')
    lnd_modelio.write(' pio_stride = -99\n')
    lnd_modelio.write(' pio_typename = "nothing"\n')
    lnd_modelio.write('/\n')
    lnd_modelio.close()
    
    ocn_modelio = open(ocn_modelio_file_name,'w')
    ocn_modelio.write('&modelio\n')
    ocn_modelio.write(' diri = "."\n')
    ocn_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    ocn_modelio.write(' logfile = '+repr(logfile_ocn)+'\n')
    ocn_modelio.write('/\n')
    ocn_modelio.write('&pio_inparm\n')
    ocn_modelio.write(' pio_numiotasks = -99\n')
    ocn_modelio.write(' pio_root = 0\n')
    ocn_modelio.write(' pio_stride = -99\n')
    ocn_modelio.write(' pio_typename = "nothing"\n')
    ocn_modelio.write('/\n')
    ocn_modelio.close()

    rof_modelio = open(rof_modelio_file_name,'w')
    rof_modelio.write('&modelio\n')
    rof_modelio.write(' diri = "."\n')
    rof_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    rof_modelio.write(' logfile = '+repr(logfile_rof)+'\n')
    rof_modelio.write('/\n')
    rof_modelio.write('&pio_inparm\n')
    rof_modelio.write(' pio_numiotasks = -99\n')
    rof_modelio.write(' pio_root = -99\n')
    rof_modelio.write(' pio_stride = -99\n')
    rof_modelio.write(' pio_typename = "nothing"\n')
    rof_modelio.write('/\n')
    rof_modelio.close()
    
    wav_modelio = open(wav_modelio_file_name,'w')
    wav_modelio.write('&modelio\n')
    wav_modelio.write(' diri = "."\n')
    wav_modelio.write(' diro = '+repr(Run_Dir)+'\n')
    wav_modelio.write(' logfile = '+repr(logfile_wav)+'\n')
    wav_modelio.write('/\n')
    wav_modelio.write('&pio_inparm\n')
    wav_modelio.write(' pio_numiotasks = -99\n')
    wav_modelio.write(' pio_root = -99\n')
    wav_modelio.write(' pio_stride = -99\n')
    wav_modelio.write(' pio_typename = "nothing"\n')
    wav_modelio.write('/\n')
    wav_modelio.close()

def Call_CLM_3D(Def_First_Run,Def_CESM_Multi_Instance, Run_Dir_Home, Run_Dir_Multi_Instance, Run_Dir, Run_Dir_Array, Ensemble_Number, num_processors, DasPy_Path, Model_Path, Model_Driver, Def_SpinUp, Def_PP, Def_Print,align_year,first_year,last_year,\
                    domain_file_path,Forcing_File_Path, Forcing_File_Path_Array, domain_file_lnd_path, domain_name, rdirc_name, aero_file_path,aero_file_name, megan_factors_file_path, megan_factors_file_name,frivinp_rtm_path,frivinp_rtm_name, \
                    case_name,hostname,orb_iyear_ad,start_type,username,
                    atm_cpl_dt,lnd_cpl_dt,ocn_cpl_dt,ice_cpl_dt,glc_cpl_dt,rof_cpl_dt, wav_cpl_dt,
                    end_restart,restart_option,start_tod,start_ymd,stop_tod,stop_ymd,ntasks_CLM,rootpe_CLM,nthreads_CLM,\
                    dtime,rtm_nsteps,fatmgrid_name,fatmlndfrc_name,fglcmask_name, finidat_name, flndtopo_name,fndepdat_name,fpftcon_name,fsnowaging_name,fsnowoptics_name,fsurdat_name, popd_streams_name, light_streams_name,\
                    wrtdia,hist_nhtfrq,hist_mfilt,hist_crtinic,hist_dov2xy,hist_ndens,hist_type1d_pertape,hist_empty_htapes,hist_avgflag_pertape, hist_fincl1,hist_fexcl1,\
                    Region_Name, Stop_Year, Stop_Month, Stop_Day, stop_tod_string, seq_maps_file_name, Row_Numbers_String, Col_Numbers_String, DAS_Data_Path,
                    COUP_OAS_PFL, CESM_Init_Flag, fcomm, fcomm_null, fcomm_rank):
        
    if Def_Print:
        print "fcomm_rank",fcomm_rank
        print "Run_Dir",Run_Dir
    
    if Def_CESM_Multi_Instance == 1:
        history_file_name = Run_Dir_Multi_Instance + Region_Name + '.clm2_0001' + '.h0.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
        
        if Ensemble_Number == 1:
            os.chdir(Run_Dir_Home)
            Run_Dir = Run_Dir_Home
            
            datm_atm_in_file_name = "datm_atm_in"
            datm_streams_txt_file_name = "datm.streams.txt"
            presaero_stream_txt_file_name = "presaero.stream.txt"
            
            if fcomm_rank == 0:
                Write_datm_atm_in(datm_atm_in_file_name, datm_streams_txt_file_name, presaero_stream_txt_file_name, domain_file_path,domain_name, rdirc_name,align_year,first_year,last_year)
                
                Forcing_File_Path = Forcing_File_Path_Array[0]
                Write_datm_streams_txt(datm_streams_txt_file_name, Def_SpinUp, domain_file_path,domain_name, rdirc_name, Forcing_File_Path, start_ymd, stop_ymd)
                
                Write_presaero_stream_txt(presaero_stream_txt_file_name,aero_file_path,aero_file_name)
                
                lnd_in_file_name = "lnd_in"
                
                Write_lnd_in(Run_Dir, lnd_in_file_name,Model_Driver,dtime,rtm_nsteps,domain_file_lnd_path, domain_name, rdirc_name, fatmgrid_name,fatmlndfrc_name,fglcmask_name,finidat_name,\
                                flndtopo_name,fndepdat_name,fpftcon_name,frivinp_rtm_name,fsnowaging_name,fsnowoptics_name,fsurdat_name, popd_streams_name, light_streams_name,\
                             wrtdia,hist_nhtfrq,hist_mfilt,hist_crtinic,hist_dov2xy,hist_ndens,hist_type1d_pertape,hist_empty_htapes,hist_avgflag_pertape,hist_fincl1,hist_fexcl1,first_year)
                
                rof_in_file_name = "rof_in"
                Write_rof_in(rof_in_file_name,frivinp_rtm_path,frivinp_rtm_name)
                
                drv_flds_in_file_name = "drv_flds_in"
                Write_45_drv_flds_in(drv_flds_in_file_name, megan_factors_file_path, megan_factors_file_name)
                
                if Def_First_Run:
                    datm_in_file_name = "datm_in"
                    atm_modelio_file_name = "atm_modelio.nml"
                    cpl_modelio_file_name = "cpl_modelio.nml"
                    glc_modelio_file_name = "glc_modelio.nml"
                    ice_modelio_file_name = "ice_modelio.nml"
                    lnd_modelio_file_name = "lnd_modelio.nml"
                    ocn_modelio_file_name = "ocn_modelio.nml"
                    rof_modelio_file_name = "rof_modelio.nml"
                    wav_modelio_file_name = "wav_modelio.nml"
                    logfile_atm = "atm.log"
                    logfile_cpl = "cpl.log"
                    logfile_glc = "glc.log"
                    logfile_ice = "ice.log"
                    logfile_lnd = "lnd.log"
                    logfile_ocn = "ocn.log"
                    logfile_rof = "rof.log"
                    logfile_wav = "wav.log"
                    Run_Dir = Run_Dir_Home
                    Write_Config_Files(datm_in_file_name,datm_atm_in_file_name,atm_modelio_file_name,cpl_modelio_file_name,glc_modelio_file_name,ice_modelio_file_name,lnd_modelio_file_name,ocn_modelio_file_name,rof_modelio_file_name,
                                       wav_modelio_file_name, logfile_atm, logfile_cpl, logfile_glc, logfile_ice, logfile_lnd, logfile_ocn, logfile_rof, logfile_wav, Run_Dir)
                
        
        elif Ensemble_Number > 1:
            os.chdir(Run_Dir_Multi_Instance)
            
            if fcomm_rank == 0:
                for Ens_Index in range(Ensemble_Number):
                    Run_Dir = Run_Dir_Array[Ens_Index]
                    Ens_Index_String = str(Ens_Index+1).zfill(4)
                    
                    datm_atm_in_file_name = "datm_atm_in_"+Ens_Index_String
                    datm_streams_txt_file_name = "datm.streams.txt_"+Ens_Index_String
                    presaero_stream_txt_file_name = "presaero.stream.txt_"+Ens_Index_String
                    
                    Write_datm_atm_in(datm_atm_in_file_name, datm_streams_txt_file_name, presaero_stream_txt_file_name, domain_file_path,domain_name, rdirc_name,align_year,first_year,last_year)
                    
                    Forcing_File_Path = Forcing_File_Path_Array[Ens_Index]
                    Write_datm_streams_txt(datm_streams_txt_file_name, Def_SpinUp, domain_file_path,domain_name, rdirc_name, Forcing_File_Path, start_ymd, stop_ymd)
                    
                    Write_presaero_stream_txt(presaero_stream_txt_file_name,aero_file_path,aero_file_name)
                    
                    lnd_in_file_name = "lnd_in_"+Ens_Index_String
                        
                    Write_lnd_in(Run_Dir, lnd_in_file_name,Model_Driver,dtime,rtm_nsteps,domain_file_lnd_path, domain_name, rdirc_name, fatmgrid_name,fatmlndfrc_name,fglcmask_name,finidat_name,\
                                    flndtopo_name,fndepdat_name,fpftcon_name,frivinp_rtm_name,fsnowaging_name,fsnowoptics_name,fsurdat_name, popd_streams_name, light_streams_name,\
                                 wrtdia,hist_nhtfrq,hist_mfilt,hist_crtinic,hist_dov2xy,hist_ndens,hist_type1d_pertape,hist_empty_htapes,hist_avgflag_pertape,hist_fincl1,hist_fexcl1,first_year)
                    
                    rof_in_file_name = "rof_in_"+Ens_Index_String
                    Write_rof_in(rof_in_file_name,frivinp_rtm_path,frivinp_rtm_name)
                    
                    if Def_First_Run:
                        datm_in_file_name = "datm_in_"+Ens_Index_String
                        atm_modelio_file_name = "atm_modelio.nml_"+Ens_Index_String
                        cpl_modelio_file_name = "cpl_modelio.nml"
                        glc_modelio_file_name = "glc_modelio.nml_"+Ens_Index_String
                        ice_modelio_file_name = "ice_modelio.nml_"+Ens_Index_String
                        lnd_modelio_file_name = "lnd_modelio.nml_"+Ens_Index_String
                        ocn_modelio_file_name = "ocn_modelio.nml_"+Ens_Index_String
                        rof_modelio_file_name = "rof_modelio.nml_"+Ens_Index_String
                        wav_modelio_file_name = "wav_modelio.nml"
                        logfile_atm = "atm.log"
                        logfile_cpl = "cpl.log"
                        logfile_glc = "glc.log"
                        logfile_ice = "ice.log"
                        logfile_lnd = "lnd.log"
                        logfile_ocn = "ocn.log"
                        logfile_rof = "rof.log"
                        logfile_wav = "wav.log"
                        
                        Write_Config_Files(datm_in_file_name,datm_atm_in_file_name,atm_modelio_file_name,cpl_modelio_file_name,glc_modelio_file_name,ice_modelio_file_name,lnd_modelio_file_name,ocn_modelio_file_name,rof_modelio_file_name,
                                              wav_modelio_file_name, logfile_atm, logfile_cpl, logfile_glc, logfile_ice, logfile_lnd, logfile_ocn, logfile_rof, logfile_wav, Run_Dir)
                    
        
    else:
        
        history_file_name = Run_Dir + Region_Name + '.clm2' + '.h0.' + Stop_Year + '-' + Stop_Month + '-' + Stop_Day + '-' + stop_tod_string + '.nc'
        
        os.chdir(Run_Dir)
        
        datm_atm_in_file_name = "datm_atm_in"
        datm_streams_txt_file_name = "datm.streams.txt"
        presaero_stream_txt_file_name = "presaero.stream.txt"
        
        if fcomm_rank == 0:
            Write_datm_atm_in(datm_atm_in_file_name, datm_streams_txt_file_name, presaero_stream_txt_file_name, domain_file_path,domain_name, rdirc_name,align_year,first_year,last_year)
            
            Write_datm_streams_txt(datm_streams_txt_file_name, Def_SpinUp, domain_file_path,domain_name, rdirc_name, Forcing_File_Path, start_ymd, stop_ymd)
            
            Write_presaero_stream_txt(presaero_stream_txt_file_name,aero_file_path,aero_file_name)
            
            lnd_in_file_name = "lnd_in"
            
            Write_lnd_in(Run_Dir, lnd_in_file_name,Model_Driver,dtime,rtm_nsteps,domain_file_lnd_path, domain_name, rdirc_name, fatmgrid_name,fatmlndfrc_name,fglcmask_name,finidat_name,\
                            flndtopo_name,fndepdat_name,fpftcon_name,frivinp_rtm_name,fsnowaging_name,fsnowoptics_name,fsurdat_name, popd_streams_name, light_streams_name,\
                         wrtdia,hist_nhtfrq,hist_mfilt,hist_crtinic,hist_dov2xy,hist_ndens,hist_type1d_pertape,hist_empty_htapes,hist_avgflag_pertape,hist_fincl1,hist_fexcl1,first_year)
            
            if Def_First_Run:
                rof_in_file_name = "rof_in"
                Write_rof_in(rof_in_file_name,frivinp_rtm_path,frivinp_rtm_name)
            
                datm_in_file_name = "datm_in"
                atm_modelio_file_name = "atm_modelio.nml"
                cpl_modelio_file_name = "cpl_modelio.nml"
                glc_modelio_file_name = "glc_modelio.nml"
                ice_modelio_file_name = "ice_modelio.nml"
                lnd_modelio_file_name = "lnd_modelio.nml"
                ocn_modelio_file_name = "ocn_modelio.nml"
                rof_modelio_file_name = "rof_modelio.nml"
                wav_modelio_file_name = "wav_modelio.nml"
                logfile_atm = "atm.log"
                logfile_cpl = "cpl.log"
                logfile_glc = "glc.log"
                logfile_ice = "ice.log"
                logfile_lnd = "lnd.log"
                logfile_ocn = "ocn.log"
                logfile_rof = "rof.log"
                logfile_wav = "wav.log"
                
                Write_Config_Files(datm_in_file_name,datm_atm_in_file_name,atm_modelio_file_name,cpl_modelio_file_name,glc_modelio_file_name,ice_modelio_file_name,lnd_modelio_file_name,ocn_modelio_file_name,rof_modelio_file_name,
                                     wav_modelio_file_name, logfile_atm, logfile_cpl, logfile_glc, logfile_ice, logfile_lnd, logfile_ocn, logfile_rof, logfile_wav, Run_Dir)
            
            
        if Def_First_Run:
            Write_seq_maps(seq_maps_file_name, DAS_Data_Path, Row_Numbers_String, Col_Numbers_String, Region_Name)
    
    if fcomm_rank == 0:
        Write_drv_in(Def_PP, Model_Driver,Def_CESM_Multi_Instance,Ensemble_Number,num_processors,case_name,hostname,orb_iyear_ad,start_type,username,
                     atm_cpl_dt,lnd_cpl_dt,ocn_cpl_dt,ice_cpl_dt,glc_cpl_dt,rof_cpl_dt, wav_cpl_dt,
                     end_restart,restart_option,start_tod,start_ymd,stop_tod,stop_ymd,ntasks_CLM,rootpe_CLM,nthreads_CLM)
    
        if os.path.exists(history_file_name):
            os.remove(history_file_name)
    
    if Def_PP == 2:
        fcomm.barrier()
        fcomm.Barrier()
    
    CLM_Output = open("CLM_Output.txt","w")
    subprocess.call(Model_Path, stdout=CLM_Output, stderr=CLM_Output, shell=True)
    CLM_Output.close()
     
    os.chdir(DasPy_Path)
    
    return

    