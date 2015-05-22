import os, sys, time, datetime, math, re, unittest, gc, subprocess, glob, string, shutil, signal, warnings, multiprocessing, socket, getpass
import numpy, scipy, scipy.io.netcdf, netCDF4
path = os.getcwd()
sys.path.append(path)
sys.path.append(path + '/Utilities')

def Call_CMEM(Def_Print,CMEM_Work_Path,rowNum_Out,colNum_Out,Latitude,Longitude,Tair_Mat,Tair_Time,Clay_Mat,Sand_Mat,ECOCVL_Mat,ECOCVH_Mat,ECOTVL_Mat,ECOTVH_Mat,ECOWAT_Mat,ECOLAIL_Mat,
              RSN_Mat,SD_Mat,STL_Mat,SWVL_Mat,TSKIN_Mat,Z_Mat,Frequency,View_Angle):
    
    if Def_Print:
        print "CMEM is Running!"
    os.chdir(CMEM_Work_Path)
    
    Input_File = open("input",'w')
    Input_File.write("&NAMOPT CIDIEL='Dobson',\n")
    Input_File.write("       CITEFF='Wigneron',\n")
    Input_File.write("       CISMR='Wilheit',\n")
    Input_File.write("       CIRGHR='Wtexture',\n")
    Input_File.write("       CIVEG='Wigneron',\n")
    Input_File.write("       CIATM='Pellarin',\n")       
    Input_File.write("       CITVEG='Tair',\n")
    Input_File.write("       CIDVEG='Ecoclimap',\n")
    Input_File.write("       CITDIEL='Teff',\n") 
    Input_File.write("       /\n")
    Input_File.write("&NAMRAD FGHZ="+str(Frequency)+",\n") 
    Input_File.write("        THETA="+str(View_Angle)+",\n")
    Input_File.write("       /\n")
    Input_File.write("&NAMLEV NLAY_SOIL_MW=100,\n")
    Input_File.write("        NLAY_SOIL_LS=7,\n")
    Input_File.write("        NOBS_ATM=,\n")
    Input_File.write("       /\n")
    Input_File.write("&NAMDEF LOFFCMEM=.True.,\n")
    Input_File.write("        LOFIELDEXP=.False.,\n")
    Input_File.write("        LGPRINT=.False.,\n")
    Input_File.write("        JPHISTLEV=1,\n")
    Input_File.write("        CFINOUT='netcdf',\n")
    Input_File.write("        LOMASK_OCEAN=.False.,\n")
    Input_File.write("        LOMASK_AUTO=.True.,\n")
    Input_File.write("/\n")
    Input_File.close()
    
    #=========================================================Tair
    if Def_Print >= 2:
        print 'Write NetCDF File:',"2T.nc"
    
    Tair_2m = netCDF4.Dataset("2T.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    Tair_2m.createDimension('LATITUDE', rowNum_Out)
    Tair_2m.createDimension('LEV', 1)
    Tair_2m.createDimension('LONGITUDE', colNum_Out)
    Tair_2m.createDimension('TIME', None)
    
    LATITUDE = Tair_2m.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.axis = "2D"
    LATITUDE.units = "Degrees North"
    LATITUDE.point_spacing = "even"
    
    LEV = Tair_2m.createVariable('LEV','f4',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = Tair_2m.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.axis = "2D"
    LONGITUDE.units = "Degrees East"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    TAIR = Tair_2m.createVariable('TAIR','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    TAIR.missing_value = -1.e+34
    TAIR.long_name = "TAIR[L=2:2921:8]"
    TAIR.history = "From tair_exp1_2006_1"
    
    TIME = Tair_2m.createVariable('TIME','f4',('TIME',))
    TIME.axis = "T"
    TIME.units = "DoY"
    
    Tair_2m.history = "netcdf4-python"
    Tair_2m.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    TAIR[0,0,:,:] = numpy.flipud(Tair_Mat)
    TIME[:] = Tair_Time
    
    Tair_2m.close()

    #=========================================================CLAY
    if Def_Print >= 2:
        print 'Write NetCDF File:',"clay.nc"
    
    clay = netCDF4.Dataset("clay.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    clay.createDimension('LATITUDE', rowNum_Out)
    clay.createDimension('LONGITUDE', colNum_Out)
    
    LATITUDE = clay.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LONGITUDE = clay.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    CLAY = clay.createVariable('CLAY','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    CLAY.missing_value = -1.e+34
    CLAY.long_name = "CLAY*100"
    CLAY.history = "From ecoclimap_tiles-2_axes"
    
    clay.history = "netcdf4-python"
    clay.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    CLAY[::] = numpy.flipud(Clay_Mat)
    
    clay.close()
    
    #=========================================================SAND
    if Def_Print >= 2:
        print 'Write NetCDF File:',"sand.nc"
    
    sand = netCDF4.Dataset("sand.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    sand.createDimension('LATITUDE', rowNum_Out)
    sand.createDimension('LONGITUDE', colNum_Out)
    
    LATITUDE = sand.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LONGITUDE = sand.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    SAND = sand.createVariable('SAND','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    SAND.missing_value = -1.e+34
    SAND.long_name = "SAND*100"
    SAND.history = "From ecoclimap_tiles-2_axes"
    
    sand.history = "netcdf4-python"
    sand.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    SAND[::] = numpy.flipud(Sand_Mat)
    
    sand.close()
    
    #=========================================================ECOCVL
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOCVL.nc"
    
    ECOCVL = netCDF4.Dataset("ECOCVL.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    ECOCVL.createDimension('LATITUDE', rowNum_Out)
    ECOCVL.createDimension('LONGITUDE', colNum_Out)
    
    LATITUDE = ECOCVL.createVariable('LATITUDE','f8',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LONGITUDE = ECOCVL.createVariable('LONGITUDE','f8',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    CVL = ECOCVL.createVariable('CVL','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    CVL.missing_value = -1.e+34
    CVL.long_name = "F_C3_CROPS[D=1]+F_C4_CROPS[D=1]+F_IRR_CROPS[D=1]+F_GRASSLAND[D=1]+F_TROP_GRASS[D=1]+F_PARK_MARSHES[D=1]"
    
    ECOCVL.history = "netcdf4-python"
    ECOCVL.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    CVL[::] = numpy.flipud(ECOCVL_Mat)
    
    ECOCVL.close()
    
    #=========================================================ECOCVH
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOCVH.nc"
    
    ECOCVH = netCDF4.Dataset("ECOCVH.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    ECOCVH.createDimension('LATITUDE', rowNum_Out)
    ECOCVH.createDimension('LONGITUDE', colNum_Out)
    
    LATITUDE = ECOCVH.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LONGITUDE = ECOCVH.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    CVH = ECOCVH.createVariable('CVH','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    CVH.missing_value = -1.e+34
    CVH.long_name = "F_BROADLEAF_TREE[D=1] +F_CONIF_TREE[D=1]+F_TROPICAL_TREE[D=1]"
    
    ECOCVH.history = "netcdf4-python"
    ECOCVH.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    CVH[::] = numpy.flipud(ECOCVH_Mat)
    
    ECOCVH.close()
    
    #=========================================================ECOTVL
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOTVL.nc"
    
    ECOTVL = netCDF4.Dataset("ECOTVL.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    
    ECOTVL.createDimension('LONGITUDE', colNum_Out)
    ECOTVL.createDimension('LATITUDE', rowNum_Out)
    
    LONGITUDE = ECOTVL.createVariable('LONGITUDE','f8',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    LATITUDE = ECOTVL.createVariable('LATITUDE','f8',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    TVL = ECOTVL.createVariable('TVL','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    TVL.missing_value = -1.e+34
    TVL.long_name = "IF(TVL_BRUT GE 0) THEN TVL_BRUT ELSE 0"
    
    ECOTVL.history = "netcdf4-python"
    ECOTVL.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    TVL[::] = numpy.flipud(ECOTVL_Mat)
    
    ECOTVL.close()
    
    #=========================================================ECOTVH
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOTVH.nc"
    
    ECOTVH = netCDF4.Dataset("ECOTVH.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    ECOTVH.createDimension('LONGITUDE', colNum_Out)
    ECOTVH.createDimension('LATITUDE', rowNum_Out)
    
    LONGITUDE = ECOTVH.createVariable('LONGITUDE','f8',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    LATITUDE = ECOTVH.createVariable('LATITUDE','f8',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    
    TVH = ECOTVH.createVariable('TVH','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    TVH.missing_value = -1.e+34
    TVH.long_name = "IF(TVH_BRUT GE 0) THEN TVH_BRUT ELSE 0"
    
    ECOTVH.history = "netcdf4-python"
    ECOTVH.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    TVH[::] = numpy.flipud(ECOTVH_Mat)
    
    ECOTVH.close()
    
    #=========================================================ECOWAT
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOWAT.nc"
    
    ECOWAT = netCDF4.Dataset("ECOWAT.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    ECOWAT.createDimension('LATITUDE', rowNum_Out)
    ECOWAT.createDimension('LONGITUDE', colNum_Out)
    
    LATITUDE = ECOWAT.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LONGITUDE = ECOWAT.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    LONGITUDE.axis = "X"
    
    WATER = ECOWAT.createVariable('WATER','f4',('LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    WATER.missing_value = -1.e+34
    WATER.long_name = "MIN(1,OCEAN+LAKES)"
    WATER.history = "From floodplains"
    
    ECOWAT.history = "netcdf4-python"
    ECOWAT.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LONGITUDE[:] = Longitude
    WATER[::] = numpy.flipud(ECOWAT_Mat)
    
    ECOWAT.close()
    
    #=========================================================ECOLAIL
    if Def_Print >= 2:
        print 'Write NetCDF File:',"ECOLAIL.nc"
    
    ECOLAIL = netCDF4.Dataset("ECOLAIL.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    ECOLAIL.createDimension('TIME', None)
    ECOLAIL.createDimension('LEV', 1)
    ECOLAIL.createDimension('LATITUDE', rowNum_Out)
    ECOLAIL.createDimension('LONGITUDE', colNum_Out)
    
    LAI = ECOLAIL.createVariable('LAI','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    LAI.missing_value = -1.e+34
    LAI.long_name = "LAI_DECADE[G=GRIDTPS]"
    LAI.history = "From ecoclimap_tiles-2_axes"
    
    LATITUDE = ECOLAIL.createVariable('LATITUDE','f8',('LATITUDE',))
    LATITUDE.units = "degrees north"
    LATITUDE.point_spacing = "even"
    LATITUDE.axis = "Y"
    
    LEV = ECOLAIL.createVariable('LEV','f8',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = ECOLAIL.createVariable('LONGITUDE','f8',('LONGITUDE',))
    LONGITUDE.axis = "X"
    LONGITUDE.units = "degrees east"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    TIME = ECOLAIL.createVariable('TIME','f8',('TIME',))
    TIME.units = "TIME"
    TIME.axis = "T"
    
    ECOLAIL.history = "netcdf4-python"
    ECOLAIL.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    LAI[0,0,:,:] = numpy.flipud(ECOLAIL_Mat)
    TIME[:] = 1
    
    ECOLAIL.close()
    
    #=========================================================RSN
    if Def_Print >= 2:
        print 'Write NetCDF File:',"RSN.nc"
    
    RSN = netCDF4.Dataset("RSN.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    RSN.createDimension('LATITUDE', rowNum_Out)
    RSN.createDimension('LEV', 1)
    RSN.createDimension('LONGITUDE', colNum_Out)
    RSN.createDimension('TIME', None)
    
    LATITUDE = RSN.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.axis = "2D"
    LATITUDE.units = "Degrees North"
    LATITUDE.point_spacing = "even"
    
    LEV = RSN.createVariable('LEV','f4',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = RSN.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.axis = "2D"
    LONGITUDE.units = "Degrees East"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    RSN_Varialbe = RSN.createVariable('RSN','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    RSN_Varialbe.missing_value = -1.e+34
    RSN_Varialbe.long_name = "SD"
    RSN_Varialbe.history = "From ALMIP"
    
    TIME = RSN.createVariable('TIME','f4',('TIME',))
    TIME.axis = "T"
    TIME.units = "DoY"
    
    RSN.history = "netcdf4-python"
    RSN.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    RSN_Varialbe[0,0,:,:] = numpy.flipud(RSN_Mat)
    TIME[:] = Tair_Time
    
    RSN.close()
    
    #=========================================================SD
    if Def_Print >= 2:
        print 'Write NetCDF File:',"SD.nc"
    
    SD = netCDF4.Dataset("SD.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    SD.createDimension('LATITUDE', rowNum_Out)
    SD.createDimension('LEV', 1)
    SD.createDimension('LONGITUDE', colNum_Out)
    SD.createDimension('TIME', None)
    
    LATITUDE = SD.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.axis = "2D"
    LATITUDE.units = "Degrees North"
    LATITUDE.point_spacing = "even"
    
    LEV = SD.createVariable('LEV','f4',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = SD.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.axis = "2D"
    LONGITUDE.units = "Degrees East"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    SD_Varialbe = SD.createVariable('SD','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',),fill_value=-1.e+34)
    SD_Varialbe.missing_value = -1.e+34
    SD_Varialbe.long_name = "SNOW Dens"
    SD_Varialbe.history = "From ALMIP"
    
    TIME = SD.createVariable('TIME','f4',('TIME',))
    TIME.axis = "T"
    TIME.units = "DoY"
    
    SD.history = "netcdf4-python"
    SD.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    SD_Varialbe[0,0,:,:] = numpy.flipud(SD_Mat)
    TIME[:] = Tair_Time
    
    SD.close()
    
    #=========================================================STL
    for i in range(7):
        
        Index = i + 1
        String = "STL"+str(Index)
        if Def_Print >= 2:
            print 'Write NetCDF File:',String+".nc"
        
        STL = netCDF4.Dataset(String+".nc", 'w', format='NETCDF3_CLASSIC')
        # Dim the dimensions of NetCDF
        STL.createDimension('LATITUDE', rowNum_Out)
        STL.createDimension('LEV', 1)
        STL.createDimension('LONGITUDE', colNum_Out)
        STL.createDimension('TIME', None)
    
        LATITUDE = STL.createVariable('LATITUDE','f8',('LATITUDE',))
        LATITUDE.units = "degrees0north"
        LATITUDE.point0spacing = "even"
        LATITUDE.axis = "Y"
        
        LEV = STL.createVariable('LEV','f8',('LEV',))
        LEV.axis = "1D"
        LEV.units = "-"
        LEV.long0name = "Tile index"
        LEV.point0spacing = "even"
        
        LONGITUDE = STL.createVariable('LONGITUDE','f8',('LONGITUDE',))
        LONGITUDE.units = "degrees0east"
        LONGITUDE.modulo = 360.
        LONGITUDE.point0spacing = "even"
        LONGITUDE.axis = "X"
        
        STL_Varialbe = STL.createVariable(String,'f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',))
        STL_Varialbe.missing0value = -1.e+034
        STL_Varialbe.units = "-"
        STL_Varialbe.long0name = " "
        STL_Varialbe.history = "From ALMIP"
        
        TIME = STL.createVariable('TIME','f8',('TIME',))
        TIME.units = "DoY"
        TIME.axis = "T"
        
        STL.history = "netcdf4-python"
        STL.Conventions = "CF-1.0"
        
        LATITUDE[:] = Latitude
        LEV[:] = 1
        LONGITUDE[:] = Longitude
        STL_Varialbe[0,0,:,:] = numpy.flipud(STL_Mat[i,::])
        TIME[:] = Tair_Time
        
        STL.close()

    
    #=========================================================SWVL
    for i in range(7):
        
        Index = i + 1
        String = "SWVL"+str(Index)
        if Def_Print >= 2:
            print 'Write NetCDF File:',String+".nc"
    
        SWVL = netCDF4.Dataset(String+".nc", 'w', format='NETCDF3_CLASSIC')
        # Dim the dimensions of NetCDF
        SWVL.createDimension('LATITUDE', rowNum_Out)
        SWVL.createDimension('LEV', 1)
        SWVL.createDimension('LONGITUDE', colNum_Out)
        SWVL.createDimension('TIME', None)
        
        LATITUDE = SWVL.createVariable('LATITUDE','f8',('LATITUDE',))
        LATITUDE.units = "degrees0north"
        LATITUDE.point0spacing = "even"
        LATITUDE.axis = "Y"
        
        LEV = SWVL.createVariable('LEV','f8',('LEV',))
        LEV.axis = "1D"
        LEV.units = "-"
        LEV.long_name = "Tile index"
        LEV.point_spacing = "even"
        
        LONGITUDE = SWVL.createVariable('LONGITUDE','f8',('LONGITUDE',))
        LONGITUDE.units = "degrees0east"
        LONGITUDE.modulo = 360.
        LONGITUDE.point0spacing = "even"
        LONGITUDE.axis = "X"
        
        SWVL_Varialbe = SWVL.createVariable(String,'f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',))
        SWVL_Varialbe.missing0value = -1.e+34
        SWVL_Varialbe.units = "-"
        SWVL_Varialbe.long0name = " "
        SWVL_Varialbe.history = "From ALMIP"
        
        TIME = SWVL.createVariable('TIME','f4',('TIME',))
        TIME.axis = "T"
        TIME.units = "DoY"
        
        SWVL.history = "netcdf4-python"
        SWVL.Conventions = "CF-1.0"
        
        LATITUDE[:] = Latitude
        LEV[:] = 1
        LONGITUDE[:] = Longitude
        SWVL_Varialbe[0,0,:,:] = numpy.flipud(SWVL_Mat[i,:,:])
        TIME[:] = Tair_Time
        
        SWVL.close()
    
    
    #=========================================================TSKIN
    if Def_Print >= 2:
        print 'Write NetCDF File:',"TSKIN.nc"
    
    TSKIN = netCDF4.Dataset("TSKIN.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    TSKIN.createDimension('LATITUDE', rowNum_Out)
    TSKIN.createDimension('LEV', 1)
    TSKIN.createDimension('LONGITUDE', colNum_Out)
    TSKIN.createDimension('TIME', None)
    
    LATITUDE = TSKIN.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.axis = "Y"
    LATITUDE.units = "Degrees North"
    LATITUDE.point_spacing = "even"
    
    LEV = TSKIN.createVariable('LEV','f4',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = TSKIN.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.axis = "X"
    LONGITUDE.units = "Degrees East"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    TSKIN_Varialbe = TSKIN.createVariable('TSKIN','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',))
    TSKIN_Varialbe.missing0value = -1.e+34
    TSKIN_Varialbe.units = "-"
    TSKIN_Varialbe.long0name = " "
    TSKIN_Varialbe.history = "From ALMIP"
    
    TIME = TSKIN.createVariable('TIME','f4',('TIME',))
    TIME.axis = "T"
    TIME.units = "DoY"
    
    TSKIN.history = "netcdf4-python"
    TSKIN.Conventions = "CF-1.0"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    TSKIN_Varialbe[0,0,:,:] = numpy.flipud(TSKIN_Mat)
    TIME[:] = Tair_Time
    
    TSKIN.close()
    
    #=========================================================Z
    if Def_Print >= 2:
        print 'Write NetCDF File:',"Z.nc"
    
    Z = netCDF4.Dataset("Z.nc", 'w', format='NETCDF3_CLASSIC')
    # Dim the dimensions of NetCDF
    Z.createDimension('LATITUDE', rowNum_Out)
    Z.createDimension('LEV', 1)
    Z.createDimension('LONGITUDE', colNum_Out)
    Z.createDimension('TIME', 1)
    
    LATITUDE = Z.createVariable('LATITUDE','f4',('LATITUDE',))
    LATITUDE.axis = "Y"
    LATITUDE.units = "Degrees North"
    LATITUDE.point_spacing = "even"
    
    LEV = Z.createVariable('LEV','f4',('LEV',))
    LEV.axis = "1D"
    LEV.units = "-"
    LEV.long_name = "Tile index"
    LEV.point_spacing = "even"
    
    LONGITUDE = Z.createVariable('LONGITUDE','f4',('LONGITUDE',))
    LONGITUDE.axis = "X"
    LONGITUDE.units = "Degrees East"
    LONGITUDE.modulo = 360.
    LONGITUDE.point_spacing = "even"
    
    Z_Varialbe = Z.createVariable('Z','f4',('TIME', 'LEV', 'LATITUDE', 'LONGITUDE',))
    Z_Varialbe.long_name = "Geopotential"
    Z_Varialbe.short_name = "Z"
    Z_Varialbe.units = "KM"
    Z_Varialbe.missing_value = 999.9
    
    TIME = Z.createVariable('TIME','f4',('TIME',))
    TIME.units = "DoY"
    TIME.long_name = "time counter"
    TIME.short_name = "time"
    TIME.note = "Il peut etre n importe quoi"
    
    Z.title = "Geopotential"
    Z.model = "ECMWF model"
    
    LATITUDE[:] = Latitude
    LEV[:] = 1
    LONGITUDE[:] = Longitude
    Z_Varialbe[0,0,:,:] = numpy.flipud(Z_Mat)
    TIME[:] = Tair_Time
    
    Z.close()
    
    OutFile_Name = glob.glob('out*.nc')
    for OutFile_Name_Index in OutFile_Name:
        if os.path.exists(OutFile_Name_Index):
            os.remove(OutFile_Name_Index)
        
    #================================== Run CMEM
    if os.name == 'nt':
        subprocess.call("cmem.exe &> /dev/null",shell=True)
    elif os.name == 'posix':
        subprocess.call("./cmem &> /dev/null",shell=True)
    time.sleep(1)
    #================================== Read the Results
    OutFile_Name = glob.glob('out*.nc')
    
    TB = netCDF4.Dataset(OutFile_Name[0], 'r')
    TBH = TB.variables['TBH'][0,0,:,:]
    TBV = TB.variables['TBV'][0,0,:,:]
    EFFECTIVE_TEMP = TB.variables['EFFECTIVE_TEMP'][0,0,:,:]
    TB.close()
    
    return TBH,TBV,EFFECTIVE_TEMP
    