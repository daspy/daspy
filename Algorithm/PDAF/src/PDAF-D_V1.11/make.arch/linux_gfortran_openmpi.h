######################################################
# Include file with machine-specific definitions     #
# for building PDAF.                                 #
#                                                    #
# Variant for Linux with gfortran and OpenMPI        #
# at AWI.                                            #
#                                                    #
# In the case of compilation without MPI, a dummy    #
# implementation of MPI, like provided in the        #
# directory nullmpi/ has to be linked when building  #
# an executable.                                     #
######################################################
# $Id: linux_gfortran_openmpi.h 1537 2014-12-20 08:18:03Z lnerger $


# Compiler, Linker, and Archiver
FC = mpif90
LD = mpif90
AR = ar
RANLIB = ranlib 

# C preprocessor
# (only required, if preprocessing is not performed via the compiler)
CPP = /usr/bin/cpp

# Definitions for CPP
# Define USE_PDAF to include PDAF
# Define BLOCKING_MPI_EXCHANGE to use blocking MPI commands to exchange data between model and PDAF
# Define PDAF_NO_UPDATE to deactivate the analysis step of the filter
# (if the compiler does not support get_command_argument()
# from Fortran 2003 you should define F77 here.)
CPP_DEFS = -DUSE_PDAF

# Optimization specs for compiler
#   (You should explicitly define double precision for floating point
#   variables in the compilation)  
OPT = -O3 -fdefault-real-8 -ffree-line-length-none

# Optimization specifications for Linker
OPT_LNK = $(OPT)

# Linking libraries (BLAS, LAPACK, if required: MPI)
LINK_LIBS =-L/usr/lib -llapack  -L$(HOME_Path)/DAS_Depends/lib -lblas   -lm

# Specifications for the archiver
AR_SPEC = 

# Specifications for ranlib
RAN_SPEC =

# Include path for MPI header file
MPI_INC =  #-Idummympi

# Object for nullMPI - if compiled without MPI library
OBJ_MPI = #nullmpi.o

# NetCDF (only required for Lorenz96)
NC_LIB   = -L$(HOME_Path)/DAS_Depends/lib
NC_INC   = -L$(HOME_Path)/DAS_Depends/include
