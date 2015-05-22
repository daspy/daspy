import os, sys, time, datetime, math, gc, subprocess, glob, string, shutil
import warnings, multiprocessing, socket, getpass, pickle, ctypes, platform
import numpy, scipy, netCDF4

def Clumping_Index(LAI, View_Zenith_Angle, Azimuth_Angle, Canopy_Height_Width_Ratio):
    
    '''
    
    Anderson, Martha C., J. M. Norman, William P. Kustas, Fuqin Li, John H. Prueger, John R. Mecikalski, 2005 
    Effects of Vegetation Clumping on Two Source Model Estimates of Surface Energy Fluxes from an Agricultural Landscape during SMACEX. J. Hydrometeor, 6, 892 909.
    doi: http://dx.doi.org/10.1175/JHM465.1
    Chen, J.M., 1996. Optically-based methods for measuring seasonal variation of leaf area index in boreal conifer stands. Agricultural and Forest Meteorology 80, 135-163.
'''
    
    p = 3.8 - 0.46 * Canopy_Height_Width_Ratio
    Omega_Zero = 0.492 * (1.0 + numpy.exp(-0.52 * (LAI - 0.45)))
    
    k = -1.0 * (0.3 + (1.7 * numpy.multiply(Omega_Zero, numpy.sin(Azimuth_Angle / 180.0 * numpy.pi) ** 0.1)) ** 14.0)
    Omega_Max = Omega_Zero + numpy.multiply((1.0 - Omega_Zero), (numpy.sin(Azimuth_Angle / 180.0 * numpy.pi) ** 0.05))
    Omega = numpy.multiply(Omega_Zero, Omega_Max) / (Omega_Zero + numpy.multiply((Omega_Max - Omega_Zero), numpy.exp(k*(View_Zenith_Angle / 180.0 * numpy.pi ** p))))
    
    return Omega


#print Clumping_Index(numpy.array([2.0,1.0]), numpy.array([45.0,45.0]), numpy.array([10.0,45.0]), numpy.array([1.0,1.0]))