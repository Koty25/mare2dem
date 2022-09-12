#
# Make include file for compiling MARE2DEM on the Habanero Cluster at Columbia University
#
#  Habanero is a linux cluster with the intel parallel studio compiler suite.
#
# Notes: 
# 
#       (1) Make sure you have loaded the intel-parallel-studio/2019 module using "module load intel-parallel-studio"
#       in the terminal or even better add that line to your .bashrc file so it is always loaded and ready to go
#       whenever you log onto Habanero.
# 
#       (2) Make sure to specify the mpiicc compiler (note the 2 i's) since that's the intel compiler.  If you use
#       mpicc (one i) you will get the GNU C compiler and the compilation will likely issue an error about the C and Fortran 
#       compilers being incompatible. One way to check you are using the intel compiler is to type "mpiicc --version"
#       in the terminal and it should return the Intel compiler name, version number and copyright info. If you mistakenly
#       use mpicc, if will instead return "gcc (GCC) ...".  
#
#       (3) Make sure to  specify the mpiifort compiler which is the Intel Fortran compiler and DO NOT USE mpif90, 
#       which is often the gfortran compiler. You can verify this using  "mpiifort --version".  
#    
#       I don't understand why Intel includes the GNU compilers in the parallel studio in addition to the intel compilers, 
#       since the Intel compilers almost always produce way faster code.
# 
#       Items (2) and (3) have been a source of headache for the uninitiated. 
#  
#       You have been warned!  
#
# 

# MPI-Fortran compiler and compiler flags:  
   FC           = mpiifort  
   FFLAGS       = -cxxlib -O2 -fpp -fPIC  # -fPIC needed for kdtree on habanero   
   
 # MPI C compiler and compiler flags:    
   CC           = mpiicc   
   CFLAGS       = -O2 -fPIC   
 
# Options for Triangle.c: 
    TRICOPTS = -O2  -fp-model precise -fp-model source -fPIC # optimizations for Triangle compiler on Intel icc requires the -fp args

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
     BLASLIB = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
 
 

