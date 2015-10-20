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
import os, sys, time, datetime, random, math, gc, subprocess, glob, string, shutil, warnings
import numpy, scipy, scipy.io.netcdf, netCDF4

def Read_Soil_Texture(Def_Region, DAS_Data_Path, Resolution_Name, Region_Name, Row_Numbers, Col_Numbers, Corner_Row_Index, Corner_Col_Index):
        
    if True:
        Top_Global_Sand_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Top_Sand.nc"
        Top_Global_Clay_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Top_Clay.nc"
        Top_Global_Organic_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Top_OC.nc"
        Top_Global_Bulk_Density_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Top_BD.nc"
        Sub_Global_Sand_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Sub_Sand.nc"
        Sub_Global_Clay_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Sub_Clay.nc"
        Sub_Global_Organic_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Sub_OC.nc"
        Sub_Global_Bulk_Density_Data_File_Name = DAS_Data_Path+"DataBase/Global_Soil_Data/"+Resolution_Name+"/"+Region_Name+"/Sub_BD.nc"
            
        Top_Global_Sand_Data = netCDF4.Dataset(Top_Global_Sand_Data_File_Name, 'r')
        Top_Global_Clay_Data = netCDF4.Dataset(Top_Global_Clay_Data_File_Name, 'r')
        Top_Global_Organic_Data = netCDF4.Dataset(Top_Global_Organic_Data_File_Name, 'r')
        Top_Global_Bulk_Density_Data = netCDF4.Dataset(Top_Global_Bulk_Density_Data_File_Name, 'r')
        Sub_Global_Sand_Data = netCDF4.Dataset(Sub_Global_Sand_Data_File_Name, 'r')
        Sub_Global_Clay_Data = netCDF4.Dataset(Sub_Global_Clay_Data_File_Name, 'r')
        Sub_Global_Organic_Data = netCDF4.Dataset(Sub_Global_Organic_Data_File_Name, 'r')
        Sub_Global_Bulk_Density_Data = netCDF4.Dataset(Sub_Global_Bulk_Density_Data_File_Name, 'r')
        
        print "Read the Global Top Clay Data"
        Sand_Top = Top_Global_Sand_Data.variables['Top_sand'][::]
        print "Read the Global Top Clay Data"
        Clay_Top = Top_Global_Clay_Data.variables['Top_Clay'][::]
        print "Read the Global Top OC Data"
        Organic_Top = Top_Global_Organic_Data.variables['Top_OC'][::]
        print "Read the Global Top BD Data"
        Bulk_Density_Top = Top_Global_Bulk_Density_Data.variables['Top_BD'][::]
        print "Read the Global Sub Sand Data"
        Sand_Sub = Sub_Global_Sand_Data.variables['Sub_Sand'][::]
        print "Read the Global Sub Clay Data"
        Clay_Sub = Sub_Global_Clay_Data.variables['Sub_Clay'][::]
        print "Read the Global Sub OC Data"
        Organic_Sub = Sub_Global_Organic_Data.variables['Sub_OC'][::]
        print "Read the Global Sub BD Data"
        Bulk_Density_Sub = Sub_Global_Bulk_Density_Data.variables['Sub_BD'][::]
    
        #Sand_Top_Region = Sand_Top[Corner_Row_Index:(Corner_Row_Index+Row_Numbers),Corner_Col_Index:(Corner_Col_Index+Col_Numbers)]
        #w,h = plt.figaspect(float(Row_Numbers)/Col_Numbers)
        #fig1 = plt.figure(figsize=(w,h))
        #ax1 = fig1.add_subplot(1,1,1)
        #im1 = ax1.imshow(Sand_Top_Region, cmap=cm.jet, interpolation='bilinear')
        #plt.colorbar(im1)
        #plt.show()
        
        print "Assign the new Soil Data in the research region"
        
            
        Sand_Top_Region = Sand_Top
        Clay_Top_Region = Clay_Top
        Organic_Top_Region = Organic_Top / 0.58 # From Organic Carbon to Organic Matter
        Bulk_Density_Top_Region = Bulk_Density_Top
        Sand_Sub_Region = Sand_Sub
        Clay_Sub_Region = Clay_Sub
        Organic_Sub_Region = Organic_Sub / 0.58 # From Organic Carbon to Organic Matter
        Bulk_Density_Sub_Region = Bulk_Density_Sub
    
    #print "Sand_Top_Region",Sand_Top_Region.min(),Sand_Top_Region.max()
    #print "Clay_Top_Region",Clay_Top_Region.min(),Clay_Top_Region.max()
    #print "Organic_Top_Region",Organic_Top_Region.min(),Organic_Top_Region.max()
    #print "Bulk_Density_Top_Region",Bulk_Density_Top_Region.min(),Bulk_Density_Top_Region.max()
    #print "Sand_Sub_Region",Sand_Sub_Region.min(),Sand_Sub_Region.max()
    #print "Clay_Sub_Region",Clay_Sub_Region.min(),Clay_Sub_Region.max()
    #print "Organic_Sub_Region",Organic_Sub_Region.min(),Organic_Sub_Region.max()
    #print "Bulk_Density_Sub_Region",Bulk_Density_Sub_Region.min(),Bulk_Density_Sub_Region.max()
    
    if True:
        Top_Global_Sand_Data.close()
        Top_Global_Clay_Data.close()
        Top_Global_Organic_Data.close()
        Top_Global_Bulk_Density_Data.close()
        Sub_Global_Sand_Data.close()
        Sub_Global_Clay_Data.close()
        Sub_Global_Organic_Data.close()
        Sub_Global_Bulk_Density_Data.close()
    
    return Sand_Top_Region,Clay_Top_Region,Organic_Top_Region,Bulk_Density_Top_Region,Sand_Sub_Region,Clay_Sub_Region,Organic_Sub_Region,Bulk_Density_Sub_Region
