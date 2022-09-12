#
# Make include file for compiling MARE2DEM on MacOS systems with OpenMPI built around the Intel compiler suite
#
#Scorep tracing manual instrumentation commands
ifdef TRACING
    PREP = scorep
    SFLAGS = --user --nocompiler
    TRACE_FLAG = -DTRACE
endif
                                  
# MPI-Fortran compiler and compiler flags:  
    FC_NO_INSTR = $(INTEL_PATH)/compilers_and_libraries/linux/mpi/intel64/bin/mpif90
    FC = $(PREP) $(SFLAGS) $(FC_NO_INSTR) $(TRACE_FLAG)
    FFLAGS  = -cxxlib -O2 -fpp -fPIC -fc=$(INTEL_PATH)/compilers_and_libraries/linux/bin/intel64/ifort    
   
# MPI C compiler and compiler flags:    
    CC_NO_INSTR = $(INTEL_PATH)/compilers_and_libraries/linux/mpi/intel64/bin/mpiicc
#   CC = $(PREP) $(CC_NO_INSTR)
    CC = $(CC_NO_INSTR)
    CFLAGS = -O2 -fPIC

# Options for Triangle.c: 
    TRICOPTS = -O2 -fPIC -fp-model precise -fp-model source # optimizations for Triangle compiler on Intel icc requires the -fp args
    
# SuperLU solver 
    SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran objects, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
    BLASDEF = -DUSE_VENDOR_BLAS  
    
# For both SuperLU and ScaLAPACK, use these commands to build the libraries:                                                                      
    ARCH = $(INTEL_PATH)/compilers_and_libraries/linux/bin/intel64/xiar
    ARCHFLAGS = ruv
    RANLIB = ranlib
   
# Use the Intel Math Kernel library for the BLASLIB since its highly optimized:
    # sequential version so that MARE's threads don't overload cpu:
    BLASLIB =  -L$(MKLROOT)/lib -I$(MKLROOT)/include -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
    
#  Note that this also assumes that you have correctly installed the Intel compiler so that shell variable MKLROOT 
#  exists and points to your intel library. You can check this by typing "echo $MKLROOT" into a terminal shell and 
#  seeing if a directory for the Intel MKL appears.