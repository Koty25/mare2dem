#
# Make include file for compiling MARE2DEM on MacOS systems with OpenMPI built around the Intel compiler suite
#

                                  
# MPI-Fortran compiler and compiler flags:  
    FC           = /opt/openmpi/bin/mpifort   # Open MPI built around the Intel Fortran compiler
    FFLAGS       = -cxxlib -O2 -fpp    
   
# MPI C compiler and compiler flags:    
    CC           = /opt/openmpi/bin/mpicc    # Open MPI built around the Intel C compiler
    CFLAGS       = -O2    
 
# Options for Triangle.c: 
    TRICOPTS = -O2  -fp-model precise -fp-model source # optimizations for Triangle compiler on Intel icc requires the -fp args
    
# SuperLU solver 
    SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran objects, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
    BLASDEF = -DUSE_VENDOR_BLAS  
    
# For both SuperLU and ScaLAPACK, use these commands to build the libraries:                                                                      
    ARCH = xiar   
    ARCHFLAGS = ruv
    RANLIB = ranlib
   
# Use the Intel Math Kernel library for the BLASLIB since its highly optimized:
    # sequential version so that MARE's threads don't overload cpu:
    BLASLIB =  -L$(MKLROOT)/lib -I$(MKLROOT)/include -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
    
#  Note that this also assumes that you have correctly installed the Intel compiler so that shell variable MKLROOT 
#  exists and points to your intel library. You can check this by typing "echo $MKLROOT" into a terminal shell and 
#  seeing if a directory for the Intel MKL appears.