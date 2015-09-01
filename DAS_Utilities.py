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
import shutil

def copyLargeFile(src, dest, buffer_size=1024*1024*1024):
    try:
        with open(src, 'rb') as fsrc:
            with open(dest, 'wb') as fdest:
                shutil.copyfileobj(fsrc, fdest, buffer_size)
    except:
        shutil.copyfile(src, dest)