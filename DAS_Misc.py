'''
# -*- coding: utf-8 -*- 
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
import numpy, netCDF4, gc, string


def Plot_Parameters(Def_Print, fm, legend, plt, cm, colors, r, Def_Region, DasPy_Path, Region_Name, Row_Numbers, Col_Numbers, DAS_Data_Path, Row_Numbers_String, Col_Numbers_String, Dim_Soil_Par, Dim_Veg_Par, Start_Month, DateString_Plot, Assim_Flag, Mask_Index, PFT_Dominant_Index,
                    Parameter_Range_Veg, Parameter_Range_PFT, Parameter_Range_Hard, Ensemble_Number, Soil_Par_Sens_Array, Veg_Par_Sens_Array, PFT_Par_Sens_Array, Hard_Par_Sens_Array,  Station_XY, Station_XY_Index, NC_FileName_Assimilation_2_Parameter, NC_FileName_Optimized_Parameter, NC_FileName_Parameter_Space_Single):
    
    Veg_Par_Name = ['roota_par', 'rootb_par', 'displar', 'z0mr', 'smpsc', 'smpso', 'rholnir', 'rholvis', 'taulnir', 'taulvis', 'rhosnir', 'rhosvis', 'tausnir', 'tausvis', 'fnitr']
    
    NC_File_Out_Assimilation_2_Parameter = netCDF4.Dataset(NC_FileName_Assimilation_2_Parameter, 'r')
    
    NC_File_Parameter_Space_Single = netCDF4.Dataset(NC_FileName_Parameter_Space_Single,'r')
    Parameter_Soil_Space_Single = NC_File_Parameter_Space_Single.variables['Parameter_Soil_Space_Single'][:,:,:]
    Parameter_PFT_Space_Single = NC_File_Parameter_Space_Single.variables['Parameter_PFT_Space_Single'][:,:,:]
    GaussRF_Array = NC_File_Parameter_Space_Single.variables['GaussRF_Array'][:,:]
    NC_File_Parameter_Space_Single.close()
    
    NC_File_Out_Optimized_Parameter = netCDF4.Dataset(NC_FileName_Optimized_Parameter, 'r')
    Parameter_Soil_Optimized_Array = NC_File_Out_Optimized_Parameter.variables['Parameter_Soil_Optimized'][:,:,:,:]
    Parameter_PFT_Optimized_Array = NC_File_Out_Optimized_Parameter.variables['Parameter_PFT_Optimized'][:,:,:,:]
    Par_Steps_Soil = len(NC_File_Out_Optimized_Parameter.dimensions['time_soil'])
    Par_Steps_PFT = len(NC_File_Out_Optimized_Parameter.dimensions['time_pft'])
    NC_File_Out_Optimized_Parameter.close()
    
    if Assim_Flag[0]:
        Mask = Mask_Index[0, ::]
    if Assim_Flag[1] :
        Mask = Mask_Index[1, ::]
    
    mksurfdata_NC_FileName_In = DAS_Data_Path + "SysModel/CLM/tools/surfdata_"+Row_Numbers_String+"x"+Col_Numbers_String+"_"+Region_Name+".nc"
    mksurfdata_NC_File_In = netCDF4.Dataset(mksurfdata_NC_FileName_In, 'r')
    Sand_Ratio = numpy.zeros((len(mksurfdata_NC_File_In.dimensions['nlevsoi']),Row_Numbers,Col_Numbers))
    Clay_Ratio = numpy.zeros((len(mksurfdata_NC_File_In.dimensions['nlevsoi']),Row_Numbers,Col_Numbers))
    Organic_Ratio = numpy.zeros((8,Row_Numbers,Col_Numbers))
    
    Sand_Ratio[0,:,:] = 1.0
    Clay_Ratio[0,:,:] = 1.0
    Organic_Ratio[0,:,:] = 1.0
    
    for i in range(1,len(mksurfdata_NC_File_In.dimensions['nlevsoi']),1):
        Sand_Ratio[i,:,:] = numpy.flipud(mksurfdata_NC_File_In.variables["PCT_SAND"][i,:,:])/numpy.flipud(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:])
        Sand_Ratio[i,:,:][numpy.where(mksurfdata_NC_File_In.variables["PCT_SAND"][0,:,:]==0.0)] = 1.0
        Clay_Ratio[i,:,:] = numpy.flipud(mksurfdata_NC_File_In.variables["PCT_CLAY"][i,:,:])/numpy.flipud(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:])
        Clay_Ratio[i,:,:][numpy.where(mksurfdata_NC_File_In.variables["PCT_CLAY"][0,:,:]==0.0)] = 1.0
    for i in range(1,8,1):
        Organic_Ratio[i,:,:] = numpy.flipud(mksurfdata_NC_File_In.variables["ORGANIC"][i,:,:]) / numpy.flipud(mksurfdata_NC_File_In.variables["ORGANIC"][0,:,:])
        Organic_Ratio[i,:,:][numpy.where(mksurfdata_NC_File_In.variables["ORGANIC"][0,:,:]==0.0)] = 1.0
    
    mksurfdata_NC_File_In.close()
    
    
    w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)
    
    
    #GaussRF_Array = numpy.zeros_like(GaussRF_Array)
    
    Sand_10cm_True = (Parameter_Soil_Space_Single[0,::] + 0.5 * numpy.reshape(GaussRF_Array[:,0],(Row_Numbers,Col_Numbers))) * 2.0 * Sand_Ratio[2,:,:]
    Sand_50cm_True = (Parameter_Soil_Space_Single[0,::] + 0.5 * numpy.reshape(GaussRF_Array[:,0],(Row_Numbers,Col_Numbers))) * 2.0 * Sand_Ratio[5,:,:]
    Clay_10cm_True = (Parameter_Soil_Space_Single[1,::] - 2 * numpy.reshape(GaussRF_Array[:,1],(Row_Numbers,Col_Numbers))) / 2.0 * Clay_Ratio[2,:,:]
    Clay_50cm_True = (Parameter_Soil_Space_Single[1,::] - 2 * numpy.reshape(GaussRF_Array[:,1],(Row_Numbers,Col_Numbers))) / 2.0 * Clay_Ratio[5,:,:]
    
    Organic_10cm_True = (Parameter_Soil_Space_Single[2,::] - 2 * numpy.reshape(GaussRF_Array[:,2],(Row_Numbers,Col_Numbers))) / 2.0 * Organic_Ratio[2,:,:]
    Organic_50cm_True = (Parameter_Soil_Space_Single[2,::] - 2 * numpy.reshape(GaussRF_Array[:,2],(Row_Numbers,Col_Numbers))) / 2.0 * Organic_Ratio[5,:,:]
    
    Sand_10cm_True[numpy.where(Sand_10cm_True < 6)] = 6
    Sand_10cm_True[numpy.where(Sand_10cm_True > 90)] = 90
    Clay_10cm_True[numpy.where(Clay_10cm_True < 3)] = 3
    Clay_10cm_True[numpy.where(Clay_10cm_True > 80)] = 80
    Organic_10cm_True[numpy.where(Organic_10cm_True < 1)] = 1
    Organic_10cm_True[numpy.where(Organic_10cm_True > 130)] = 130
    
    Sand_50cm_True[numpy.where(Sand_50cm_True < 6)] = 6
    Sand_50cm_True[numpy.where(Sand_50cm_True > 90)] = 90
    Clay_50cm_True[numpy.where(Clay_50cm_True < 3)] = 3
    Clay_50cm_True[numpy.where(Clay_50cm_True > 80)] = 80
    Organic_50cm_True[numpy.where(Organic_50cm_True < 1)] = 1
    Organic_50cm_True[numpy.where(Organic_50cm_True > 130)] = 13
    
    LAI_True = Parameter_PFT_Space_Single[0,::]
    SAI_True = Parameter_PFT_Space_Single[1,::]
    
    
    Sand_10cm_Background = Parameter_Soil_Space_Single[0,::] * Sand_Ratio[2,:,:]
    Sand_50cm_Background = Parameter_Soil_Space_Single[0,::] * Sand_Ratio[5,:,:]
    Clay_10cm_Background = Parameter_Soil_Space_Single[1,::] * Clay_Ratio[2,:,:]
    Clay_50cm_Background = Parameter_Soil_Space_Single[1,::] * Clay_Ratio[5,:,:]
    
    Organic_10cm_Background = Parameter_Soil_Space_Single[2,::] * Organic_Ratio[2,:,:]
    Organic_50cm_Background = Parameter_Soil_Space_Single[2,::] * Organic_Ratio[5,:,:]
    
    LAI_Background = Parameter_PFT_Space_Single[0,::]
    SAI_Background = Parameter_PFT_Space_Single[1,::]
    
    Parameter_PFT_Space_Single_True = Parameter_PFT_Space_Single * 2.0
    
    if numpy.size(numpy.where(numpy.asarray(Soil_Par_Sens_Array) == True)) > 0:   
        
    
        
        Variable_Min = 10.0
        Variable_Max = 60.0
        ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
        color_boun_list = []
        color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
        for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
            color_bound[0] += color_bound[2]
            color_boun_list.append(color_bound[0])
        
        fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
        fig1.suptitle(DateString_Plot, fontsize=16)
        ax = fig1.add_subplot(2, 3, 1)
        Sand_10cm_True = numpy.ma.masked_where(Mask, Sand_10cm_True)
        im1 = ax.imshow(Sand_10cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        ax.set_title('Sand_10cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 2)
        Sand_10cm_Background = numpy.ma.masked_where(Mask, Sand_10cm_Background)
        im1 = ax.imshow(Sand_10cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Sand_10cm_True).flatten()-numpy.asarray(Sand_10cm_Background.flatten()))**2)/numpy.sum(((Sand_10cm_True).flatten()-numpy.mean((Sand_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Sand_10cm_Background",NSE1
            ax.set_title('Sand_10cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Sand_10cm_Background')
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 3)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,0,::],axis=0) * Sand_Ratio[2,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Sand_10cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Sand_10cm_True).flatten()-numpy.mean((Sand_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Sand_10cm_Optimized",NSE1
            ax.set_title('Sand_10cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Sand_10cm_Optimized')   
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 4)
        Sand_50cm_True = numpy.ma.masked_where(Mask, Sand_50cm_True)
        im2 = ax.imshow(Sand_50cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        ax.set_title('Sand_50cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 5)
        Sand_50cm_Background = numpy.ma.masked_where(Mask, Sand_50cm_Background)
        im2 = ax.imshow(Sand_50cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Sand_50cm_True).flatten()-numpy.asarray(Sand_50cm_Background.flatten()))**2)/numpy.sum(((Sand_50cm_True).flatten()-numpy.mean((Sand_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Sand_50cm_Background",NSE1
            ax.set_title('Sand_50cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Sand_50cm_Background')    
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 6)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,0,::],axis=0) * Sand_Ratio[5,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Sand_50cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Sand_50cm_True).flatten()-numpy.mean((Sand_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Sand_50cm_Optimized",NSE1
            ax.set_title('Sand_50cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Sand_50cm_Optimized')  
        plt.grid(True)
        
        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Sand_Optimized.png")
        plt.show()
        
        Variable_Min = 5.0
        Variable_Max = 50.0
        ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
        color_boun_list = []
        color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
        for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
            color_bound[0] += color_bound[2]
            color_boun_list.append(color_bound[0])
        
        fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
        fig1.suptitle(DateString_Plot, fontsize=16)
        ax = fig1.add_subplot(2, 3, 1)
        Clay_10cm_True = numpy.ma.masked_where(Mask, Clay_10cm_True)
        im1 = ax.imshow(Clay_10cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        ax.set_title('Clay_10cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 2)
        Clay_10cm_Background = numpy.ma.masked_where(Mask, Clay_10cm_Background)
        im1 = ax.imshow(Clay_10cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Clay_10cm_True).flatten()-numpy.asarray(Clay_10cm_Background.flatten()))**2)/numpy.sum(((Clay_10cm_True).flatten()-numpy.mean((Clay_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Clay_10cm_Background",NSE1
            ax.set_title('Clay_10cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Clay_10cm_Background')    
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 3)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,1,::],axis=0) * Clay_Ratio[2,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Clay_10cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Clay_10cm_True).flatten()-numpy.mean((Clay_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Clay_10cm_Optimized",NSE1
            ax.set_title('Clay_10cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Clay_10cm_Optimized')  
        plt.grid(True)
        
        
        ax = fig1.add_subplot(2, 3, 4)
        Clay_50cm_True = numpy.ma.masked_where(Mask, Clay_50cm_True)
        im1 = ax.imshow(Clay_50cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        ax.set_title('Clay_50cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 5)
        Clay_50cm_Background = numpy.ma.masked_where(Mask, Clay_50cm_Background)
        im1 = ax.imshow(Clay_50cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Clay_50cm_True).flatten()-numpy.asarray(Clay_50cm_Background.flatten()))**2)/numpy.sum(((Clay_50cm_True).flatten()-numpy.mean((Clay_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Clay_50cm_Background",NSE1
            ax.set_title('Clay_50cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Clay_50cm_Background')
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 6)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,1,::],axis=0) * Clay_Ratio[5,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Clay_50cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Clay_50cm_True).flatten()-numpy.mean((Clay_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Clay_50cm_Optimized",NSE1
            ax.set_title('Clay_50cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Clay_50cm_Optimized')    
        plt.grid(True)
        
        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Clay_Optimized.png")
        plt.show()
        
        Variable_Min = 5.0
        Variable_Max = 50.0
        ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
        color_boun_list = []
        color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
        for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
            color_bound[0] += color_bound[2]
            color_boun_list.append(color_bound[0])
        
        fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
        fig1.suptitle(DateString_Plot, fontsize=16)
        ax = fig1.add_subplot(2, 3, 1)
        Organic_10cm_True = numpy.ma.masked_where(Mask, Organic_10cm_True)
        im1 = ax.imshow(Organic_10cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        ax.set_title('Organic_10cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 2)
        Organic_10cm_Background = numpy.ma.masked_where(Mask, Organic_10cm_Background)
        im1 = ax.imshow(Organic_10cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Organic_10cm_True).flatten()-numpy.asarray(Organic_10cm_Background.flatten()))**2)/numpy.sum(((Organic_10cm_True).flatten()-numpy.mean((Organic_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Organic_10cm_Background",NSE1
            ax.set_title('Organic_10cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Organic_10cm_Background')  
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 3)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,2,::],axis=0) * Organic_Ratio[2,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Organic_10cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Organic_10cm_True).flatten()-numpy.mean((Organic_10cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Organic_10cm_Optimized",NSE1
            ax.set_title('Organic_10cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Organic_10cm_Optimized')  
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 4)
        Organic_50cm_True = numpy.ma.masked_where(Mask, Organic_50cm_True)
        im1 = ax.imshow(Organic_50cm_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        ax.set_title('Organic_50cm_True')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 5)
        Organic_50cm_Background = numpy.ma.masked_where(Mask, Organic_50cm_Background)
        im1 = ax.imshow(Organic_50cm_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im1, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Organic_50cm_True).flatten()-numpy.asarray(Organic_50cm_Background.flatten()))**2)/numpy.sum(((Organic_50cm_True).flatten()-numpy.mean((Organic_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Organic_50cm_Background",NSE1
            ax.set_title('Organic_50cm_Background-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Organic_50cm_Background')   
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 6)
        Mean_Soil_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_Soil_Space_Ensemble'][:,2,::],axis=0) * Organic_Ratio[5,:,:]
        Mean_Soil_Ensemble = numpy.ma.masked_where(Mask, Mean_Soil_Ensemble)
        im2 = ax.imshow(Mean_Soil_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        if Def_Print:
            NSE1 = 1.0 - numpy.sum(((Organic_50cm_True).flatten()-numpy.asarray(Mean_Soil_Ensemble.flatten()))**2)/numpy.sum(((Organic_50cm_True).flatten()-numpy.mean((Organic_50cm_True).flatten()))**2)
            print "Nash-Sutcliffe efficiencies of Organic_50cm_Optimized",NSE1
            ax.set_title('Organic_50cm_Optimized-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('Organic_50cm_Optimized')
        plt.grid(True)
        
        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Organic_Optimized.png")
        plt.show()
        
        
    if numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) > 0:
                
        if PFT_Par_Sens_Array[0][0] or PFT_Par_Sens_Array[1][0]:
            Variable_Min = -0.5
            Variable_Max = 3.0
            ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
            color_boun_list = []
            color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
            for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                color_bound[0] += color_bound[2]
                color_boun_list.append(color_bound[0])
            
            fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax = fig1.add_subplot(1, 3, 1)
            LAI_True = numpy.ma.masked_where(Mask, LAI_True)
            im1 = ax.imshow(LAI_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('LAI')
            plt.grid(True)
            ax = fig1.add_subplot(1, 3, 2)
            LAI_Background = numpy.ma.masked_where(Mask, LAI_Background)
            im1 = ax.imshow(LAI_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            if Def_Print:
                NSE1 = 1.0 - numpy.sum(((LAI_True).flatten()-numpy.asarray(LAI_Background.flatten()))**2)/numpy.sum(((LAI_True).flatten()-numpy.mean((LAI_True).flatten()))**2)
                print "Nash-Sutcliffe efficiencies of LAI_Background",NSE1
                ax.set_title('LAI_Background-NSE:'+str(NSE1)[0:4])
            else:
                ax.set_title('LAI_Background')
            plt.grid(True)
            ax = fig1.add_subplot(1, 3, 3)
            Mean_LAI_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,0,::],axis=0)
            Mean_LAI_Ensemble = numpy.ma.masked_where(Mask, Mean_LAI_Ensemble)
            im2 = ax.imshow(Mean_LAI_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')
            if Def_Print:
                NSE1 = 1.0 - numpy.sum(((LAI_True).flatten()-numpy.asarray(Mean_LAI_Ensemble.flatten()))**2)/numpy.sum(((LAI_True).flatten()-numpy.mean((LAI_True).flatten()))**2)
                print "Nash-Sutcliffe efficiencies of LAI_Optimized",NSE1
                ax.set_title('LAI_Optimized-NSE:'+str(NSE1)[0:4])
            else:
                ax.set_title('LAI_Optimized')
            plt.grid(True)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/LAI_Optimized.png")
            plt.show()        
        
        if PFT_Par_Sens_Array[0][1] or PFT_Par_Sens_Array[1][1]:
            Variable_Min = -0.5
            Variable_Max = 3.0
            ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
            color_boun_list = []
            color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
            for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
                color_bound[0] += color_bound[2]
                color_boun_list.append(color_bound[0])
            
            fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax = fig1.add_subplot(1, 3, 1)
            SAI_True = numpy.ma.masked_where(Mask, SAI_True)
            im1 = ax.imshow(SAI_True, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('SAI')
            plt.grid(True)
            ax = fig1.add_subplot(1, 3, 2)
            SAI_Background = numpy.ma.masked_where(Mask, SAI_Background)
            im1 = ax.imshow(SAI_Background, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            if Def_Print:
                NSE1 = 1.0 - numpy.sum(((SAI_True).flatten()-numpy.asarray(SAI_Background.flatten()))**2)/numpy.sum(((SAI_True).flatten()-numpy.mean((SAI_True).flatten()))**2)
                print "Nash-Sutcliffe efficiencies of SAI_Background",NSE1
                ax.set_title('SAI_Background-NSE:'+str(NSE1)[0:4])
            else:
                ax.set_title('SAI_Background')
            plt.grid(True)
            ax = fig1.add_subplot(1, 3, 3)
            Mean_SAI_Ensemble = numpy.mean(NC_File_Out_Assimilation_2_Parameter.variables['Parameter_PFT_Space_Ensemble'][:,1,::],axis=0)
            Mean_SAI_Ensemble = numpy.ma.masked_where(Mask, Mean_SAI_Ensemble)
            im2 = ax.imshow(Mean_SAI_Ensemble, cmap=cm.jet, norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')
            if Def_Print:
                NSE1 = 1.0 - numpy.sum(((SAI_True).flatten()-numpy.asarray(Mean_SAI_Ensemble.flatten()))**2)/numpy.sum(((SAI_True).flatten()-numpy.mean((SAI_True).flatten()))**2)
                print "Nash-Sutcliffe efficiencies of SAI_Optimized",NSE1
                ax.set_title('SAI_Optimized-NSE:'+str(NSE1)[0:4])
            else:
                ax.set_title('SAI_Optimized')
            plt.grid(True)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/SAI_Optimized.png")
            plt.show()
            
        
    print "*************************************** Plot the Parameter Ensembles"
    
    #for Station_Index in range(numpy.size(Station_XY)/2):
    for Station_Index in range(1):
        
        if numpy.size(numpy.where(numpy.asarray(Soil_Par_Sens_Array) == True)) > 0:
            x=numpy.arange(Par_Steps_Soil)
            
            fig1 = plt.figure()
            fig1.subplots_adjust(hspace=0.5)
            ax1 = fig1.add_subplot(2, 1, 1)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,110.0)
            
            Parameter_Soil_Optimized_Sub = numpy.zeros(Par_Steps_Soil)
            Parameter_Soil_Optimized_Sub_Ens = numpy.zeros((Par_Steps_Soil,Ensemble_Number))
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[0,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Sand_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Sand_10cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
            
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,0,Station_Index] * numpy.flipud(Sand_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Sand_10cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            
            ax1 = fig1.add_subplot(2, 1, 2)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,110.0)
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[0,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Sand_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Sand_50cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
            
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,0,Station_Index] * numpy.flipud(Sand_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Sand_50cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Sand_Ensembles_Station_"+str(Station_Index+1)+".png")
            plt.show()
            
            fig1 = plt.figure()
            fig1.subplots_adjust(hspace=0.5)
            ax1 = fig1.add_subplot(2, 1, 1)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,110.0)
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[1,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Clay_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Clay_10cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
           
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,1,Station_Index] * numpy.flipud(Clay_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                 #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Clay_10cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
                
            ax1 = fig1.add_subplot(2, 1, 2)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,110.0)
            
            Parameter_Soil_Optimized_Sub = numpy.zeros(Par_Steps_Soil)
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[1,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Clay_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Clay_50cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
           
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,1,Station_Index] * numpy.flipud(Clay_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                 #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Clay_50cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Clay_Ensembles_Station_"+str(Station_Index+1)+".png")
            plt.show()
            
            fig1 = plt.figure()
            fig1.subplots_adjust(hspace=0.5)
            ax1 = fig1.add_subplot(2, 1, 1)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,130.0)
            
            Parameter_Soil_Optimized_Sub = numpy.zeros(Par_Steps_Soil)
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[2,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Organic_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Organic_10cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
            
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,2,Station_Index] * numpy.flipud(Organic_Ratio[2,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Organic_10cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            
            ax1 = fig1.add_subplot(2, 1, 2)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_Soil-1)
            ax1.set_ylim(0.0,140.0)
            
            Parameter_Soil_Optimized_Sub = numpy.zeros(Par_Steps_Soil)
            
            ax1.plot(x,numpy.repeat(Parameter_Soil_Space_Single[2,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]] * numpy.flipud(Organic_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(Organic_50cm_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_Soil),color='r',linewidth=3)
            for Ens_Index in range(Ensemble_Number):
                for Parameter_Soil_Optimized_Sub_Index in range(Par_Steps_Soil):
                    Parameter_Soil_Optimized_Sub[Parameter_Soil_Optimized_Sub_Index] = Parameter_Soil_Optimized_Array[Parameter_Soil_Optimized_Sub_Index,Ens_Index,2,Station_Index] * numpy.flipud(Organic_Ratio[5,:,:])[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]]
                Parameter_Soil_Optimized_Sub_Ens[:,Ens_Index] = Parameter_Soil_Optimized_Sub
                ax1.plot(x,Parameter_Soil_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_Soil_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_Soil_Optimized_Sub
            #print x,Parameter_Soil_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('Organic_50cm',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
                
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Organic_Ensembles_Station_"+str(Station_Index+1)+".png")
            plt.show()
            
            
        
        if numpy.size(numpy.where(numpy.asarray(PFT_Par_Sens_Array) == True)) > 0:   
            
            x=numpy.arange(Par_Steps_PFT)
                
            fig1 = plt.figure()
            fig1.subplots_adjust(hspace=0.5)
            ax1 = fig1.add_subplot(1, 1, 1)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_PFT-1)
            ax1.set_ylim(-0.5,3.0)
            
            Parameter_PFT_Optimized_Sub = numpy.zeros(Par_Steps_PFT)
            Parameter_PFT_Optimized_Sub_Ens = numpy.zeros((Par_Steps_PFT,Ensemble_Number))
            
            ax1.plot(x,numpy.repeat(Parameter_PFT_Space_Single[0,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_PFT),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(LAI_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_PFT),color='r',linewidth=3)
            
            for Ens_Index in range(Ensemble_Number):
                for Parameter_PFT_Optimized_Sub_Index in range(Par_Steps_PFT):
                    Parameter_PFT_Optimized_Sub[Parameter_PFT_Optimized_Sub_Index] = Parameter_PFT_Optimized_Array[Parameter_PFT_Optimized_Sub_Index,Ens_Index,0,Station_Index]
                Parameter_PFT_Optimized_Sub_Ens[:,Ens_Index] = Parameter_PFT_Optimized_Sub
                ax1.plot(x,Parameter_PFT_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_PFT_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_PFT_Optimized_Sub
            #print x,Parameter_PFT_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('LAI',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/LAI_Ensembles_Station_"+str(Station_Index+1)+".png")
            plt.show()
            
            x=numpy.arange(Par_Steps_PFT)
                
            fig1 = plt.figure()
            fig1.subplots_adjust(hspace=0.5)
            ax1 = fig1.add_subplot(1, 1, 1)
            fig1.suptitle(DateString_Plot, fontsize=16)
            ax1.set_xlim(0.0,Par_Steps_PFT-1)
            ax1.set_ylim(-0.5,3.0)
            
            Parameter_PFT_Optimized_Sub = numpy.zeros(Par_Steps_PFT)
            Parameter_PFT_Optimized_Sub_Ens = numpy.zeros((Par_Steps_PFT,Ensemble_Number))
            
            ax1.plot(x,numpy.repeat(Parameter_PFT_Space_Single[1,Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_PFT),color='b',linewidth=3)
            ax1.plot(x,numpy.repeat(SAI_True[Station_XY_Index[Station_Index][1],Station_XY_Index[Station_Index][0]],Par_Steps_PFT),color='r',linewidth=3)
            
            for Ens_Index in range(Ensemble_Number):
                for Parameter_PFT_Optimized_Sub_Index in range(Par_Steps_PFT):
                    Parameter_PFT_Optimized_Sub[Parameter_PFT_Optimized_Sub_Index] = Parameter_PFT_Optimized_Array[Parameter_PFT_Optimized_Sub_Index,Ens_Index,1,Station_Index]
                Parameter_PFT_Optimized_Sub_Ens[:,Ens_Index] = Parameter_PFT_Optimized_Sub
                ax1.plot(x,Parameter_PFT_Optimized_Sub,color='k', linewidth=3 ,alpha=min([1.0,4.0/float(Ensemble_Number)]))
            ax1.plot(x,numpy.mean(Parameter_PFT_Optimized_Sub_Ens,axis=1),color='k',linewidth=3)
                #print Parameter_PFT_Optimized_Sub
            #print x,Parameter_PFT_Space_Single[12,Station_XY_Index[0][1],Station_XY_Index[0][0]]
            ax1.set_ylabel('SAI',fontsize=18)
            ax1.set_xlabel('Time Step',fontsize=18)
            prop = fm.FontProperties(size=12)
            legend(('Backgroud', 'Truth', "Ensembles"), 'upper right', shadow=True,prop=prop, ncol=3)
        
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fontsize=18)
            
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/SAI_Ensembles_Station_"+str(Station_Index+1)+".png")
            plt.show()
        
            
    NC_File_Out_Assimilation_2_Parameter.close()
    
    plt.close('all')
    gc.collect()

    
def Plot_States(octave, fm, legend, plt, cm, colors, Def_Region, Region_Name, Plot_Analysis, DasPy_Path, Soil_Layer_Num, Ensemble_Number, \
                Row_Numbers, Col_Numbers, Dim_Obs_Type, Observation_Matrix, NAvalue, Def_Print,
                Prop_Grid_Array_Sys_Index, Observation_Matrix_Index, SensorType_Sub, SensorVariable_Sub, SensorQuantity_Sub, SensorResolution_Sub, Variable_ID_Sub, QC_ID_Sub, Variable_List,\
                NA_Index_Prop_Grid_Array_Sys, Mask_Index, Analysis_Grid, Analysis_Variable_Name, DateString_Plot, Variable_Assimilation_Flag, 
                NC_FileName_Assimilation_2_Constant, NC_FileName_Assimilation_2_Diagnostic, NC_FileName_Assimilation_2_Initial, NC_FileName_Assimilation_2_Initial_Copy, ObsModel_Mat_Masked, Observation_File_Name):
#                            Mean_For_NA = numpy.mean(Prop_Grid_Array_Sys[:,:,0][numpy.where(Prop_Grid_Array_Sys[:,:,0] != NAvalue)])
#                            Prop_Grid_Array_Sys_Temp[numpy.where(Prop_Grid_Array_Sys[:,:,0] == NAvalue)] = Mean_For_NA
    
#                            Mean_For_NA = numpy.mean(Analysis_Grid_Temp[numpy.where(Prop_Grid_Array_Sys[:,:,0] != NAvalue)])
#                            Analysis_Grid_Temp[numpy.where(Prop_Grid_Array_Sys[:,:,0] == NAvalue)] = Mean_For_NA

    print "--------------------------------- Plot for Checking",SensorVariable_Sub,Variable_ID_Sub
    
    NC_File_Out_Assimilation_2_Diagnostic = netCDF4.Dataset(NC_FileName_Assimilation_2_Diagnostic, 'r')
    NC_File_Out_Assimilation_2_Initial = netCDF4.Dataset(NC_FileName_Assimilation_2_Initial, 'r')
    CLM_Soil_Moisture_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Moisture_Ensemble_Mat'][:,:,:,:]
    CLM_Soil_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Temperature_Ensemble_Mat'][:,:,:,:]
    #CLM_Soil_Ice_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Soil_Ice_Ensemble_Mat'][:,:,:,:]
    CLM_Vegetation_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Vegetation_Temperature_Ensemble_Mat'][:,:,:]
    CLM_Ground_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Ground_Temperature_Ensemble_Mat'][:,:,:]
    #CLM_2m_Air_Temperature_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_2m_Air_Temperature_Ensemble_Mat'][:,:,:]
    #CLM_Snow_Depth_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Depth_Ensemble_Mat'][:,:,:]
    #CLM_Snow_Water_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_Snow_Water_Ensemble_Mat'][:,:,:]
    #CLM_ROOTFR_Ensemble_Mat = NC_File_Out_Assimilation_2_Initial.variables['CLM_ROOTFR_Ensemble_Mat'][:,:,:,:]
    
    Analysis_Grid_Array = NC_File_Out_Assimilation_2_Diagnostic.variables['Analysis_Grid_Array'][:,Prop_Grid_Array_Sys_Index,:,:]
    Innovation_State = NC_File_Out_Assimilation_2_Diagnostic.variables['Innovation_State'][:,Prop_Grid_Array_Sys_Index,:,:]
    Increments_State = NC_File_Out_Assimilation_2_Diagnostic.variables['Increments_State'][:,Prop_Grid_Array_Sys_Index,:,:]
    
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
    Analysis_Grid_Array_Copy = NC_File_Out_Assimilation_2_Initial_Copy.variables['Prop_Grid_Array_Sys'][:,Prop_Grid_Array_Sys_Index,:,:]
    NC_File_Out_Assimilation_2_Initial_Copy.close()
    
    CLM_Soil_Moisture_Ensemble = numpy.zeros((Soil_Layer_Num, Row_Numbers * Col_Numbers, Ensemble_Number),dtype=numpy.float32)
    CLM_Soil_Temperature_Ensemble = numpy.zeros((Soil_Layer_Num, Row_Numbers * Col_Numbers, Ensemble_Number),dtype=numpy.float32)
    CLM_Ground_Temperature_Ensemble = numpy.zeros((Row_Numbers * Col_Numbers, Ensemble_Number),dtype=numpy.float32)
    CLM_Vegetation_Temperature_Ensemble = numpy.zeros((Row_Numbers * Col_Numbers, Ensemble_Number),dtype=numpy.float32)
    
    for Ens_Index in range(Ensemble_Number):
        for Soil_Layer_Index in range(Soil_Layer_Num):
            
            CLM_Soil_Moisture_Ensemble[Soil_Layer_Index,:,Ens_Index] = CLM_Soil_Moisture_Ensemble_Mat_Copy[Soil_Layer_Index,:,:,Ens_Index].flatten()
            CLM_Soil_Temperature_Ensemble[Soil_Layer_Index,:,Ens_Index] = CLM_Soil_Temperature_Ensemble_Mat_Copy[Soil_Layer_Index,:,:,Ens_Index].flatten()
        #CLM_Ground_Temperature_Ensemble[:,Ens_Index] = CLM_Ground_Temperature_Ensemble_Mat[Sub_Block_Left:Sub_Block_Right,Sub_Block_North:Sub_Block_South,Ens_Index].flatten()
        #CLM_Vegetation_Temperature_Ensemble[:,Ens_Index] = CLM_Vegetation_Temperature_Ensemble_Mat[Sub_Block_Left:Sub_Block_Right,Sub_Block_North:Sub_Block_South,Ens_Index].flatten()
             
    
    if SensorVariable_Sub != "Irrigation_Scheduling":
        
        Observation_Matrix_Copy = numpy.copy(Observation_Matrix[Observation_Matrix_Index,::])
        Observation_Matrix_Copy = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue, Observation_Matrix_Copy)
        
        
        Variable_Min_Innovation= 0.1
        Variable_Max_Innovation = 0.4
        Variable_Min_Increments = 0.1
        Variable_Max_Increments = 0.4
        Variable_Min_Bias = 0
        Variable_Max_Bias = 1
        
        Innovation_State_Sub = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue,numpy.mean(Innovation_State[:, ::],axis=0))
        Increments_State_Sub = numpy.ma.masked_where(Observation_Matrix_Copy == NAvalue,numpy.mean(Increments_State[:, ::],axis=0))
        Innovation_State_Sub = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Innovation_State_Sub)
        Increments_State_Sub = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Increments_State_Sub)
        
        print "Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)],SensorVariable_Sub"
        print Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)],SensorVariable_Sub
        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Soil_Moisture':
            Variable_Min = 0.1
            Variable_Max = 0.4
            Variable_Min_Innovation= -0.03
            Variable_Max_Innovation = 0.03
            Variable_Min_Increments = -0.03
            Variable_Max_Increments = 0.03
            
            if SensorQuantity_Sub == 'K':
                Variable_Min_Innovation= -10.0
                Variable_Max_Innovation = 10.0
                ObsModel_Mat_Masked = numpy.ma.masked_where(ObsModel_Mat_Masked<50,ObsModel_Mat_Masked)
                ObsModel_Mat_Masked = numpy.ma.masked_where(ObsModel_Mat_Masked>350,ObsModel_Mat_Masked)
                Variable_Obs_Min = 150.0
                Variable_Obs_Max = 250.0
            elif SensorQuantity_Sub == 'Neutron_Count':
                Variable_Min_Innovation= -20.0
                Variable_Max_Innovation = 20.0
                ObsModel_Mat_Masked = numpy.ma.masked_where(ObsModel_Mat_Masked==NAvalue,ObsModel_Mat_Masked)
                if Def_Region == -1:
                    Variable_Obs_Min = 500.0
                    Variable_Obs_Max = 700.0
                if Def_Region == -2:
                    Variable_Obs_Min = 1600.0
                    Variable_Obs_Max = 2200.0
            else:
                Variable_Obs_Min = 0.0
                Variable_Obs_Max = 0.6
            
            NA_Index_Prop_Grid_Array_Sys_Index = numpy.where(Mask_Index[Prop_Grid_Array_Sys_Index, ::])
            NA_Index_Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index,::][NA_Index_Prop_Grid_Array_Sys_Index] = NAvalue
            Prop_Grid_Array_Sys_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], NA_Index_Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index,::])
            Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid[Prop_Grid_Array_Sys_Index,::])
                             
            Layer_0_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Moisture_Ensemble[0,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_1_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Moisture_Ensemble[1,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_2_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Moisture_Ensemble[2,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_3_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Moisture_Ensemble[3,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_4_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Moisture_Ensemble[4,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            
            Layer_0_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[0,:,:,:],axis=2))
            Layer_1_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[1,:,:,:],axis=2))
            Layer_2_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[2,:,:,:],axis=2))
            Layer_3_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[3,:,:,:],axis=2))
            Layer_4_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Moisture_Ensemble_Mat[4,:,:,:],axis=2))
            
            
        if Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Surface_Temperature':
            Variable_Min = 280.0
            Variable_Max = 310.0
            Variable_Min_Innovation= -5.0
            Variable_Max_Innovation = 5.0
            Variable_Min_Increments = -5.0
            Variable_Max_Increments = 5.0
            Variable_Obs_Min = 280.0
            Variable_Obs_Max = 310.0
            NA_Index_Prop_Grid_Array_Sys_Index = numpy.where(Mask_Index[Prop_Grid_Array_Sys_Index, ::])
            NA_Index_Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index,::][NA_Index_Prop_Grid_Array_Sys_Index] = NAvalue
            Prop_Grid_Array_Sys_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], NA_Index_Prop_Grid_Array_Sys[Prop_Grid_Array_Sys_Index,::])
            Analysis_Grid_Temp = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::], Analysis_Grid[Prop_Grid_Array_Sys_Index,::])
            
            Layer_0_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Temperature_Ensemble[0,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_1_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Temperature_Ensemble[1,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_2_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Temperature_Ensemble[2,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_3_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Temperature_Ensemble[3,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            Layer_4_Model = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.reshape(numpy.mean(CLM_Soil_Temperature_Ensemble[4,:,:],axis=1),((Row_Numbers,Col_Numbers))))
            
            Layer_0_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[0,:,:,:],axis=2))
            Layer_1_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[1,:,:,:],axis=2))
            Layer_2_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[2,:,:,:],axis=2))
            Layer_3_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[3,:,:,:],axis=2))
            Layer_4_Analysis = numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],numpy.mean(CLM_Soil_Temperature_Ensemble_Mat[4,:,:,:],axis=2))
            
            
                
        #Variable_Min = min([Observation_Matrix.min(),Prop_Grid_Array_Sys_Temp.min(),Analysis_Grid_Temp.min()])
        #Variable_Max = max([Observation_Matrix.max(),Prop_Grid_Array_Sys_Temp.max(),Analysis_Grid_Temp.max()])
        
        Observation_Matrix_Copy = numpy.ma.masked_where(numpy.ma.getmask(Analysis_Grid_Temp), Observation_Matrix_Copy)
        
        print Observation_Matrix_Copy.min(), Prop_Grid_Array_Sys_Temp.min(), Analysis_Grid_Temp.min() 
        print Observation_Matrix_Copy.max(), Prop_Grid_Array_Sys_Temp.max(), Analysis_Grid_Temp.max()
        print "Variable_Obs_Min, Variable_Obs_Max, Variable_Min, Variable_Max"
        print Variable_Obs_Min, Variable_Obs_Max, Variable_Min, Variable_Max
        #ticks=range(min(numpy.min(Observation_Matrix_Copy),numpy.min(Prop_Grid_Array_Sys_Temp),numpy.min(Analysis_Grid_Temp)),max(numpy.max(Observation_Matrix_Copy),numpy.max(Prop_Grid_Array_Sys_Temp),numpy.max(Analysis_Grid_Temp)))
        
        #**********************************************************************
        ticks_Obs = numpy.arange(Variable_Obs_Min, Variable_Obs_Max, (Variable_Obs_Max - Variable_Obs_Min) / 5.0)
        color_boun_list_Obs = []
        color_bound_Obs = [Variable_Obs_Min, Variable_Obs_Max, (Variable_Obs_Max - Variable_Obs_Min) / 100.0]
        for i in range(int((color_bound_Obs[1] - color_bound_Obs[0]) / color_bound_Obs[2])):
            color_bound_Obs[0] += color_bound_Obs[2]
            color_boun_list_Obs.append(color_bound_Obs[0])
        #**********************************************************************
        ticks = numpy.arange(Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 5.0)
        color_boun_list = []
        color_bound = [Variable_Min, Variable_Max, (Variable_Max - Variable_Min) / 100.0]
        for i in range(int((color_bound[1] - color_bound[0]) / color_bound[2])):
            color_bound[0] += color_bound[2]
            color_boun_list.append(color_bound[0])
       
        #********************************************************************
        ticks_Innovation = numpy.arange(Variable_Min_Innovation, Variable_Max_Innovation, (Variable_Max_Innovation - Variable_Min_Innovation) / 5.0)
        color_boun_list_Innovation = []
        color_bound_Innovation = [Variable_Min_Innovation, Variable_Max_Innovation, (Variable_Max_Innovation - Variable_Min_Innovation) / 100.0]
        for i in range(int((color_bound_Innovation[1] - color_bound_Innovation[0]) / color_bound_Innovation[2])):
            color_bound_Innovation[0] += color_bound_Innovation[2]
            color_boun_list_Innovation.append(color_bound_Innovation[0])
        #**********************************************************************
        ticks_Increments = numpy.arange(Variable_Min_Increments, Variable_Max_Increments, (Variable_Max_Increments - Variable_Min_Increments) / 5.0)
        color_boun_list_Increments = []
        color_bound_Increments = [Variable_Min_Increments, Variable_Max_Increments, (Variable_Max_Increments - Variable_Min_Increments) / 100.0]
        for i in range(int((color_bound_Increments[1] - color_bound_Increments[0]) / color_bound_Increments[2])):
            color_bound_Increments[0] += color_bound_Increments[2]
            color_boun_list_Increments.append(color_bound_Increments[0])
        
        w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)
        fig1 = plt.figure(figsize=(w*3, h*2), dpi=80)
        fig1.text(.5,.95, Analysis_Variable_Name[Prop_Grid_Array_Sys_Index]+"-"+DateString_Plot, fontsize=16, ha='center')
        
        ax = fig1.add_subplot(2, 3, 1)
        im1 = ax.imshow(Observation_Matrix_Copy, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Obs, ncolors=300))
        plt.colorbar(im1, ticks=ticks_Obs, orientation='horizontal')
        ax.set_title('Observation')
        plt.grid(True)
        ax = fig1.add_subplot(2, 3, 2)
        im2 = ax.imshow(Prop_Grid_Array_Sys_Temp, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im2, ticks=ticks, orientation='horizontal')
        
        Prop_Grid_Array_Sys_Temp = numpy.ma.masked_where(numpy.where(Observation_Matrix_Copy==NAvalue),Prop_Grid_Array_Sys_Temp)
        if Def_Print:
            NSE1 = 1.0 - numpy.sum((Observation_Matrix_Copy.flatten()-numpy.asarray(Prop_Grid_Array_Sys_Temp.flatten()))**2)/numpy.sum((Observation_Matrix_Copy.flatten()-numpy.mean(Observation_Matrix_Copy.flatten()))**2)
            print "Nash-Sutcliffe efficiencies of CLM",NSE1
            ax.set_title('CLM-NSE:'+str(NSE1)[0:4])
        else:
            ax.set_title('CLM')  
        plt.grid(True)   
                              
        ax = fig1.add_subplot(2, 3, 3)
        im3 = ax.imshow(Analysis_Grid_Temp, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
        plt.colorbar(im3, ticks=ticks, orientation='horizontal')
        
        Analysis_Grid_Temp = numpy.ma.masked_where(numpy.where(Observation_Matrix_Copy==NAvalue),Analysis_Grid_Temp)
        if Def_Print:
            NSE2 = 1.0 - numpy.sum((Observation_Matrix_Copy.flatten()-numpy.asarray(Analysis_Grid_Temp.flatten()))**2)/numpy.sum((Observation_Matrix_Copy.flatten()-numpy.mean(Observation_Matrix_Copy.flatten()))**2)
            print "Nash-Sutcliffe efficiencies Analysis",NSE2
            ax.set_title('Analysis-NSE:'+str(NSE2)[0:4])
        else:
            ax.set_title('Analysis')
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 4)
        im4 = ax.imshow(ObsModel_Mat_Masked, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Obs, ncolors=300))
        plt.colorbar(im4, ticks=ticks_Obs, orientation='horizontal')
        ax.set_title('ObsModel_Mat_Masked')
        plt.grid(True)
        
        ax = fig1.add_subplot(2, 3, 5)
        im4 = ax.imshow(Innovation_State_Sub, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Innovation, ncolors=300))
        plt.colorbar(im4, ticks=ticks_Innovation, orientation='horizontal')
        ax.set_title('Innovation_State')
        plt.grid(True)
        
        
        ax = fig1.add_subplot(2, 3, 6)
        im4 = ax.imshow(Increments_State_Sub, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
        plt.colorbar(im4, ticks=ticks_Increments, orientation='horizontal')
        ax.set_title('Increments_State')
        plt.grid(True)
        
        plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Analysis_"+str(Observation_Matrix_Index)+".png") 
        
        if Plot_Analysis >= 2:
            for Ens_Index in range(Ensemble_Number):
                w, h = plt.figaspect(float(Row_Numbers) / Col_Numbers)
                fig1 = plt.figure(figsize=(w*2, h*2), dpi=80)
                fig1.text(.5,.95, Analysis_Variable_Name[Prop_Grid_Array_Sys_Index]+"-"+DateString_Plot, fontsize=16, ha='center')
                
                ax = fig1.add_subplot(2, 2, 1)
                im4 = ax.imshow(numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Analysis_Grid_Array_Copy[Ens_Index,Prop_Grid_Array_Sys_Index, ::]), cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                plt.colorbar(im4, ticks=ticks, orientation='horizontal')
                ax.set_title('Analysis_Grid_Array_Copy_Ens'+str(Ens_Index+1))
                plt.grid(True)
                
                ax = fig1.add_subplot(2, 2, 2)
                im4 = ax.imshow(numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Analysis_Grid_Array[Ens_Index,Prop_Grid_Array_Sys_Index, ::]), cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
                plt.colorbar(im4, ticks=ticks, orientation='horizontal')
                ax.set_title('Analysis_Grid_Array_Ens'+str(Ens_Index+1))
                plt.grid(True)
                
                Analysis_Grid_Array
                
                ax = fig1.add_subplot(2, 2, 3)
                im4 = ax.imshow(numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Innovation_State[Ens_Index,Prop_Grid_Array_Sys_Index, ::]), cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Innovation, ncolors=300))
                plt.colorbar(im4, ticks=ticks_Innovation, orientation='horizontal')
                ax.set_title('Innovation_State_Ens'+str(Ens_Index+1))
                plt.grid(True)
                
                
                ax = fig1.add_subplot(2, 2, 4)
                im4 = ax.imshow(numpy.ma.masked_where(Mask_Index[Prop_Grid_Array_Sys_Index, ::],Increments_State[Ens_Index,Prop_Grid_Array_Sys_Index, ::]), cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
                plt.colorbar(im4, ticks=ticks_Increments, orientation='horizontal')
                ax.set_title('Increments_State_Ens'+str(Ens_Index+1))
                plt.grid(True)
                
                plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Ens/" + "Analysis_Grid_Array_Innovation_Increments_Ens"+str(Ens_Index+1) + ".png") 
        
        if (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Surface_Temperature') or (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Soil_Moisture') or (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Latent_Heat') or (Variable_Assimilation_Flag[Variable_List.index(SensorVariable_Sub)] == 1 and SensorVariable_Sub == 'Latent_Heat_Daily'):
            # Plot the deep soil
            
            
            fig1 = plt.figure(figsize=(w*4, h*4), dpi=80)
            fig1.text(.5,.95, Analysis_Variable_Name[Prop_Grid_Array_Sys_Index]+"-"+DateString_Plot, fontsize=16, ha='center')
            
            ax = fig1.add_subplot(3, 5, 1)
            im1 = ax.imshow(Layer_0_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_0_Model')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 2)
            im2 = ax.imshow(Layer_1_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')                            
            ax.set_title('Layer_1_Model')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 3)
            im3 = ax.imshow(Layer_2_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im3, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_2_Model')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 4)
            im2 = ax.imshow(Layer_3_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')                            
            ax.set_title('Layer_3_Model')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 5)
            im3 = ax.imshow(Layer_4_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im3, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_4_Model')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 6)
            im1 = ax.imshow(Layer_0_Analysis, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_0_Analysis')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 7)
            im2 = ax.imshow(Layer_1_Analysis, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')                            
            ax.set_title('Layer_1_Analysis')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 8)
            im3 = ax.imshow(Layer_2_Analysis, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im3, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_2_Analysis')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 9)
            im1 = ax.imshow(Layer_3_Analysis, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im1, ticks=ticks, orientation='horizontal')
            ax.set_title('Layer_3_Analysis')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 10)
            im2 = ax.imshow(Layer_4_Analysis, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list, ncolors=300))
            plt.colorbar(im2, ticks=ticks, orientation='horizontal')                            
            ax.set_title('Layer_4_Analysis')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 11)
            im1 = ax.imshow(Layer_0_Analysis-Layer_0_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
            plt.colorbar(im1, ticks=ticks_Increments, orientation='horizontal')
            ax.set_title('Increment')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 12)
            im2 = ax.imshow(Layer_1_Analysis-Layer_1_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
            plt.colorbar(im2, ticks=ticks_Increments, orientation='horizontal')                            
            ax.set_title('Increment')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 13)
            im3 = ax.imshow(Layer_2_Analysis-Layer_2_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
            plt.colorbar(im3, ticks=ticks_Increments, orientation='horizontal')
            ax.set_title('Increment')
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 14)
            im2 = ax.imshow(Layer_3_Analysis-Layer_3_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
            plt.colorbar(im2, ticks=ticks_Increments, orientation='horizontal')                            
            ax.set_title('Increment')     
            plt.grid(True)
            ax = fig1.add_subplot(3, 5, 15)
            im3 = ax.imshow(Layer_4_Analysis-Layer_4_Model, cmap=cm.jet, interpolation='bilinear', norm=colors.BoundaryNorm(color_boun_list_Increments, ncolors=300))
            plt.colorbar(im3, ticks=ticks_Increments, orientation='horizontal')
            ax.set_title('Increment')
            plt.grid(True)
            
            plt.savefig(DasPy_Path+"Analysis/DAS_Temp/"+Region_Name+"/Deep_Layers_"+str(Observation_Matrix_Index)+".png")
                              
                                      
        plt.show()
    
    plt.close('all')
    gc.collect()


    