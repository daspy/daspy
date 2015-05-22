import os, sys, time, datetime, random, math, gc, subprocess, glob, string, shutil, warnings, commands, socket
import numpy, scipy, scipy.io.netcdf, netCDF4

def COSMIC_Py(COSMIC,n,nlyr,soil_moisture,layerz,bd,lattwat,nthreads):
    
    '''
    !!!
    integer,             intent(in)  :: nlyr ! Total number of soil layers - ( default 3000 each 1mm for 3m)
    real(8),             intent(in)  :: soil_moisture(nsteps,nlevels) ! soil_moisture vector V/V
    real(8),             intent(in)  :: layerz(nlevels) ! original soil layer depths, COSMIC needs them to be positive and in centimeters
    
    !real(8) :: bd     = 0.0 ! Dry soil bulk density (g/m3)
    !real(8) :: vwclat = 0.0 ! Volumetric "lattice" water content (m3/m3)
    !real(8) :: N      = 0.0 ! High energy neutron flux (-)
    !real(8) :: alpha  = 0.0 ! Ratio of Fast Neutron Creation Factor (Soil to Water), alpha (-)
    !real(8) :: L1     = 0.0 ! High Energy Soil   Attenuation Length (g/cm2)
    !real(8) :: L2     = 0.0 ! High Energy Water  Attenuation Length (g/cm2)
    !real(8) :: L3     = 0.0 ! Fast Neutron Soil  Attenuation Length (g/cm2)
    !real(8) :: L4     = 0.0 ! Fast Neutron Water Attenuation Length (g/cm2)
    !!!
    '''
    
    soil_moisture = numpy.asarray(soil_moisture,dtype=numpy.float32)
    layerz = numpy.asarray(layerz,dtype=numpy.float32)
    n = numpy.asarray(n,dtype=numpy.float32)
    bd = numpy.asarray(bd,dtype=numpy.float32)
    lattwat = numpy.asarray(lattwat,dtype=numpy.float32)
    
    alpha = 0.401 - 0.101 * bd
    l1 = numpy.ones_like(bd,dtype=numpy.float32)*162.0
    l2 = numpy.ones_like(bd,dtype=numpy.float32)*129.1
    l3 = -31.65 + 99.29 * bd
    l4 = numpy.ones_like(bd,dtype=numpy.float32)*3.16
    
#    print "soil_moisture",numpy.max(soil_moisture),numpy.min(soil_moisture)
#    print "layerz",numpy.max(layerz),numpy.min(layerz)
#    print "bd",numpy.max(bd),numpy.min(bd)
#    print "lattwat",numpy.max(lattwat),numpy.min(lattwat)
#    print "alpha",numpy.max(alpha),numpy.min(alpha)
    
    #numpy.savez("COSMIC.npz",nlyr,soil_moisture,layerz,bd,lattwat,n,alpha,l1,l2,l3,l4,nthreads)
    val,soil_moisture_weighted,sensor_depth_weighted,istatus = COSMIC.cosmic(nlyr,soil_moisture,layerz,bd,lattwat,n,alpha,l1,l2,l3,l4,nthreads)
    
    return numpy.asarray(val,dtype=numpy.float32),numpy.asarray(soil_moisture_weighted,dtype=numpy.float32),numpy.asarray(sensor_depth_weighted,dtype=numpy.float32)
