#!/bin/bash -x
#MSUB -l nodes=1:ppn=16:turbomode
#MSUB -l walltime=00:10:00
#MSUB -M x.han@fz-juelich.de
#MSUB -m abe
#MSUB -v tpt=1

#total no of processors used for this simulation
# NSLOTS = nodes * ppn / tpt
NSLOTS=16

module unload parastation
module unload intel
module unload mkl
module load GCC/4.8.2
module load parastation/mpi2-gcc482-5.0.29-1

##just make sure ppserver and parallel python paths are correct, if previously not set by user
export LD_LIBRARY_PATH=/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib:/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib64:/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib64/R/lib:$LD_LIBRARY_PATH
export PATH=/lustre/jhome7/jicg41/jicg4128/DAS_Depends/bin:$PATH
export C_INCLUDE_PATH=/lustre/jhome7/jicg41/jicg4128/DAS_Depends/include

export PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib64/python2.7/site-packages:/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib/python2.7/site-packages:/lustre/jhome7/jicg41/jicg4128/DAS_Depends/lib/python2.7/site-packages/ga4py:$PYTHONPATH
export PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/ReBEL:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/Savitzky_Golay:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/PSRF/:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/SysModel/CLM:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/SysModel/PFCLM:PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/ReBEL:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/Savitzky_Golay:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/PSRF/DAS:PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/ReBEL:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/Savitzky_Golay:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm/PSRF/DAS:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ForcingData:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ForcingData/IPOLATES:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ForcingData/GMAO_Perturb:$PYTHONPATH
export PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ObsModel/LST:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ObsModel/ASCAT:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ObsModel/MODIS:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ObsModel/CMEM:$PYTHONPATH
export PYTHONPATH=/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Utilities/Soil:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Utilities/Algorithm/Geostatistics/Scripts:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/Algorithm:/lustre/jhome7/jicg41/jicg4128/DasPy_Release/ForcingData:$PYTHONPATH

export OMP_NUM_THREADS=8
cd $PBS_O_WORKDIR 
echo "Current Workdir: $PBS_O_WORKDIR" 

echo " Starting DAS on $NSLOTS processors "
python2.7 DAS.py $NSLOTS | tee log_hiwater.txt 

##end of of simulation

