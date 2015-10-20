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
import numpy, socket, sys, os, subprocess, multiprocessing
sys.path.append('Utilities')
import pyper

def DAS_Initialize(Model_Driver, Def_Region, mpi4py_rank=0):
    
    if mpi4py_rank == 0:
        print "socket.gethostname()",socket.gethostname()
    
    HOME_Path="/lustre/jhome7/jicg41/jicg4128"
    
    DasPy_Path = HOME_Path+"/DasPy_Release/"
    DAS_Data_Path = HOME_Path+"/DasPy_Release/DAS_Data/"
    DAS_Output_Path = HOME_Path+"/DasPy_Release/DAS_Data/"
    DAS_Depends_Path = HOME_Path+"/DAS_Depends/"
    geog_data_path = ''
    WRF_WPS_Path = ""
    WRF_WRF_Path = ""
    
    if Def_Region == 3:
        PicHeight=1000.0
        PicWidth=1000.0
        RegionName="Rur River Basin"
        Row_Numbers = 100    # The Number of Rows of Output Data
        Col_Numbers = 100   # The Number of Cols of Output Data
        # Rur Boundary
        Grid_Resolution_CEA = 1000.0
        Grid_Resolution_CEA_String = "1km"
        mksrf_edgew = 5.9
        mksrf_edgee = 6.733
        mksrf_edges = 50.38
        mksrf_edgen = 51.213
        Region_Name = "Rur"
        Run_Dir_Home = DAS_Output_Path + "SysModel/CLM/"+Region_Name+"/3D"
        Hydraulic_File_Name = DAS_Data_Path + "DataBase/"+Region_Name+"_Hydraulic.nc"
        Mask_File = DAS_Data_Path + "Analysis/Data/"+Region_Name+"/VEG_"+Region_Name+"_1km.dat"
        Observation_Path = DAS_Data_Path + "Observation/"+Region_Name+"/"
        Forcing_Folder = "Bilinear_1km_1hour_"+Region_Name
        
        Station_XY = numpy.array([[6.2, 50.95],[6.36,50.5]])
    
    if mpi4py_rank == 0:
        print "mksrf_edgew,mksrf_edgee",mksrf_edgew,mksrf_edgee
        print "mksrf_edges,mksrf_edgen",mksrf_edges,mksrf_edgen
        print "Center Point x", mksrf_edgew+(mksrf_edgee-mksrf_edgew)/2.0
        print "Center Point y", mksrf_edges+(mksrf_edgen-mksrf_edges)/2.0
        
    Grid_Resolution_GEO = numpy.zeros(2)    # 0 : x-direction, 1: y-direction
    Grid_Resolution_GEO[0] = (mksrf_edgee - mksrf_edgew) / Col_Numbers
    Grid_Resolution_GEO[1] = (mksrf_edgen - mksrf_edges) / Row_Numbers
    if mpi4py_rank == 0:
        print "Grid_Resolution_GEO",Grid_Resolution_GEO
        print"mksrf_edgew,mksrf_edgee,mksrf_edges,mksrf_edgen",mksrf_edgew,mksrf_edgee,mksrf_edges,mksrf_edgen
        print Station_XY[0][0],Station_XY[0][1],numpy.size(Station_XY)
    
    Forcing_File_Path_Home = DAS_Data_Path + 'ForcingData/'+Forcing_Folder+'/CLM'
        
    Station_XY_Index = numpy.zeros(numpy.shape(Station_XY),dtype=numpy.integer)
    
    for Staton_Index in range(numpy.size(Station_XY)/2):
        Station_XY_Index[Staton_Index][0] = int((Station_XY[Staton_Index][0] - mksrf_edgew) / Grid_Resolution_GEO[0])
        Station_XY_Index[Staton_Index][1] = int((mksrf_edgen - Station_XY[Staton_Index][1]) / Grid_Resolution_GEO[1])
    
    if mpi4py_rank == 0:
        print "Col_Index,Row_Index"
        print Station_XY_Index
    
    Row_Numbers_String = str(Row_Numbers)
    Col_Numbers_String = str(Col_Numbers)
        
    if mpi4py_rank == 0:
        print DAS_Depends_Path+"lib64/R/bin/R"
        #print pyper.R(RCMD=DAS_Depends_Path+"lib64/R/bin/R",use_numpy=True)
        try:
            r = pyper.R(RCMD=DAS_Depends_Path+"lib64/R/bin/R",use_numpy=True)
        except:
            r = pyper.R(use_numpy=True)
        #r('.libPaths("~/Rlibs")')    # Add Library Path to R
        # because the area function can not be called after below line
        print r("options(menu.graphics=FALSE)")
        print r("suppressPackageStartupMessages(library(sp))")
        print r('suppressPackageStartupMessages(library(geoR))')
        print r('suppressPackageStartupMessages(library(intamap))')
        print r('suppressPackageStartupMessages(library(matlab))')
        print r('suppressPackageStartupMessages(library(rgdal))')
        #print r('suppressPackageStartupMessages(library(maptools))')
        #r('suppressPackageStartupMessages(library(CompModSA))')
        #r('suppressPackageStartupMessages(library(sensitivity))')
        #print r('suppressPackageStartupMessages(library(ncdf))')
        #r('suppressPackageStartupMessages(library(lhs))')
        #r('suppressPackageStartupMessages(library(coda))')
        #print r('suppressPackageStartupMessages(library(maptools))')
        #r('suppressPackageStartupMessages(library(GSIF))')
        print r('suppressPackageStartupMessages(library(R.utils))')
        
        r.assign("DasPy_Path",DasPy_Path)
        r.assign('Def_Region', Def_Region)
        r.assign('PicHeight', PicHeight)
        r.assign('PicWidth', PicWidth)
        r.assign('Row_Numbers', Row_Numbers)
        r.assign('Col_Numbers', Col_Numbers)
        r.assign('Grid_Resolution_CEA', Grid_Resolution_CEA)
        r.assign('Grid_Resolution_GEO', Grid_Resolution_GEO)
        r.assign('mksrf_edgew', mksrf_edgew)
        r.assign('mksrf_edges', mksrf_edges)
        r.assign('mksrf_edgee', mksrf_edgee)
        r.assign('mksrf_edgen', mksrf_edgen)
        
        r('xy <- cbind(mksrf_edgew,mksrf_edges)')
        r('xy_meter <- project(xy,"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +datum=WGS84 +no_defs")')
        #print  r['xy'],r['xy_meter']
        mksrf_edgew_CEA = r['xy_meter'][0][0]
        mksrf_edges_CEA = r['xy_meter'][0][1]
        
        r('xy <- cbind(mksrf_edgee,mksrf_edgen)')
        r('xy_meter <- project(xy,"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +datum=WGS84 +no_defs")')
        #print  r['X_Coordiates'],r['Y_Coordiates'],r['xy'],r['xy_meter']
        mksrf_edgee_CEA = r['xy_meter'][0][0]
        mksrf_edgen_CEA = r['xy_meter'][0][1]
        
        #print mksrf_edgew_CEA,mksrf_edgee_CEA,mksrf_edges_CEA,mksrf_edgen_CEA
        
        #============================== Projection Parameters
        r('xy <- cbind(mksrf_edgew,mksrf_edges)')
        r('xy_meter <- project(xy,"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +datum=WGS84 +no_defs")')
        #print  r['xy_meter']
        xllcenter = r['xy_meter'][0][0] + Grid_Resolution_CEA / 2.0
        yllcenter = r['xy_meter'][0][1] + Grid_Resolution_CEA / 2.0
        
        Z_Resolution = 1.0
        
        Proj_String = "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +datum=WGS84 +no_defs"    
        r.assign('Proj_String', Proj_String)
        r.assign('mksrf_edgew', mksrf_edgew)
        r.assign('mksrf_edges', mksrf_edges)
        r.assign('Prop_Grid_Array_Observation', numpy.zeros((Row_Numbers, Col_Numbers)))
        r('grid <- GridTopology(c(mksrf_edgew,mksrf_edges), c(Grid_Resolution_GEO[1],Grid_Resolution_GEO[2]), c(Col_Numbers,Row_Numbers))')
        r('data <- data.frame(value  = as.vector(fliplr(Prop_Grid_Array_Observation)))')
        r('Prop_Grid <- SpatialGridDataFrame(grid, data,  proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))')
        r('observations_grid <- as(Prop_Grid,"SpatialPointsDataFrame")')
        r('PRJ = spTransform(observations_grid,CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +datum=WGS84 +no_defs"))')
        
        MODEL_CEA_X = numpy.asarray(r['matrix(coordinates(PRJ)[,1],nrow=Row_Numbers,byrow=TRUE)'])
        MODEL_CEA_Y = numpy.asarray(r['matrix(coordinates(PRJ)[,2],nrow=Row_Numbers,byrow=TRUE)'])
        #print MODEL_CEA_X
		
        MODEL_X_Left = numpy.min(MODEL_CEA_X)
        MODEL_X_Right = numpy.max(MODEL_CEA_X)
        MODEL_Y_Lower = numpy.min(MODEL_CEA_Y)
        MODEL_Y_Upper = numpy.max(MODEL_CEA_Y)
                
        r.assign('MODEL_X_Left', MODEL_X_Left)
        r.assign('MODEL_X_Right', MODEL_X_Right)
        r.assign('MODEL_Y_Lower', MODEL_Y_Lower)
        r.assign('MODEL_Y_Upper', MODEL_Y_Upper)
        print "MODEL_X_Left,MODEL_X_Right,MODEL_Y_Lower,MODEL_Y_Upper"
        print MODEL_X_Left,MODEL_X_Right,MODEL_Y_Lower,MODEL_Y_Upper
        
        if Def_Region == 1:
            r.assign('XMin',mksrf_edgew+Grid_Resolution_GEO[0])
        if Def_Region == 3:
            r.assign('XMin',mksrf_edgew+Grid_Resolution_GEO[0])
            r.assign('YMin',mksrf_edges+Grid_Resolution_GEO[1])
        else:
            r.assign('XMin',mksrf_edgew)
        r.assign('YMin', mksrf_edges)
        r.assign('CellSize_X', Grid_Resolution_GEO[0])
        r.assign('CellSize_Y', Grid_Resolution_GEO[1])
        r.assign('XDim', Col_Numbers)
        r.assign('YDim', Row_Numbers)
        r.assign('RegionName', RegionName)
        r.assign('Region_Name', Region_Name)
        
        octave = []       
    else:
        r = None
        mksrf_edgew_CEA=mksrf_edgee_CEA=mksrf_edges_CEA=mksrf_edgen_CEA=None
        xllcenter=yllcenter=MODEL_X_Left=MODEL_X_Right=MODEL_Y_Lower=MODEL_Y_Upper=MODEL_CEA_X=MODEL_CEA_Y=Z_Resolution=Proj_String=None
        octave = None
    
    UTC_Zone = int(round(((mksrf_edgew+(mksrf_edgee-mksrf_edgew)/2.0) - 7.5) / 15.0))    # The Difference Between the UTM Time and the Local Time
    
    if mpi4py_rank == 0:
        print "UTC_Zone",UTC_Zone
        print "################################ Working on",RegionName
    
    return PicHeight, PicWidth, RegionName, Row_Numbers, Col_Numbers, Grid_Resolution_CEA, Grid_Resolution_GEO, \
            mksrf_edgee, mksrf_edgew, mksrf_edges, mksrf_edgen, Region_Name, Run_Dir_Home, DAS_Output_Path, Hydraulic_File_Name, \
            Mask_File, Observation_Path, DAS_Data_Path, DasPy_Path, Forcing_File_Path_Home, DAS_Depends_Path, geog_data_path, \
            WRF_WPS_Path, WRF_WRF_Path, Station_XY, Station_XY_Index, r, octave, Row_Numbers_String, Col_Numbers_String, Grid_Resolution_CEA_String,\
            xllcenter, yllcenter, MODEL_X_Left, MODEL_X_Right, MODEL_Y_Lower, MODEL_Y_Upper, MODEL_CEA_X, MODEL_CEA_Y, Z_Resolution, Proj_String, UTC_Zone
            
            
