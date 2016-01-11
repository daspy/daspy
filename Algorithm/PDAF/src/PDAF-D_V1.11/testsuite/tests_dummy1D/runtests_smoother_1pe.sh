#!/bin/sh
#$Id: runtests_1pe.sh 750 2009-08-04 13:42:11Z lnerger $

#############################################################
### Script for testing different configuration options    ###
### of PDAF. This version is to test all subtypes of the  ###
### different filters without parallelization.            ###
#############################################################

#set -vx

# General configuration
NENS=50                 # Ensemble size in EnKF/SEIK/LSEIK
NEOF=`expr $NENS - 1`    # Number of EOFs in SEEK
CONF="-dim_state 300 -screen 1 -dim_lag 10"    # General configuration (state dimension)
EXE="./pdaf_dummy_online"  # Name of executable

TEST_SEEK=1   # (1) to perform tests with the SEEK filter
TEST_SEIK=1   # (1) to perform tests with the SEIK filter
TEST_ENKF=1   # (1) to perform tests with the Ensemble Kalman filter
TEST_LSEIK=1  # (1) to perform tests with the local SEIK filter
TEST_ETKF=1   # (1) to perform tests with the ETKF
TEST_LETKF=1  # (1) to perform tests with the LETKF
TEST_ESTKF=1  # (1) to perform tests with the ESTKF
TEST_LESTKF=1 # (1) to perform tests with the LESTKF

# Perform tests
echo "====================  Testing smoother in PDAF  ===================="

echo "Machine: " `uname -a`
echo "Date: " `date`

if [ $TEST_ENKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run EnKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 2 -subtype 0 -filename output_smoother_enkf0.dat
fi
if [ $TEST_ETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 0 -filename output_smoother_etkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 1 -filename output_smoother_etkf1.dat
$EXE $CONF -dim_ens $NENS -filtertype 4 -subtype 2 -filename output_smoother_etkf2.dat
fi
if [ $TEST_LETKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LETKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 0 -filename output_smoother_letkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 1 -filename output_smoother_letkf1.dat
$EXE $CONF -dim_ens $NENS -filtertype 5 -subtype 2 -filename output_smoother_letkf2.dat
fi
if [ $TEST_ESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run ESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 0 -filename output_smoother_estkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 6 -subtype 2 -filename output_smoother_estkf2.dat
fi
if [ $TEST_LESTKF -eq 1 ]
then
echo "--------------------------------------------------------"
echo "Run LESTKF tests"
echo "--------------------------------------------------------"
$EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 0 -filename output_smoother_lestkf0.dat
$EXE $CONF -dim_ens $NENS -filtertype 7 -subtype 2 -filename output_smoother_lestkf2.dat
fi

echo "PDAF tests completed: " `date`

