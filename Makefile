# Makefile for MARE2DEM
#
# Usage:   make INCLUDE=<make.inc> 
#            
#   <make.inc> should be the name of an include (text) file that specifies 
#   the MPI Fortran and C compilers and compilers flags for your 
#   computing system.  See ./include/make_macos.inc for an example. 
#   I recommend you create your own make.inc files (with a specific name such as 
#   make_MyMegaCluster.inc) for each system where you build MARE2DEM. 
#
#   Example make include files:
#
#    make_macos.inc    - MacOS system with the intel compiler suite.
#    make_habanero.inc - Linux cluster with the intel parallel studio compilers.
# 
#   Note that MARE2DEM currently requires the Intel Fortran and C 
#   C compilers and relies on the Intel Math Kernel Library (MKL).
#   Either OpenMPI or MPICH MPI implementations work.
#                                            
# Cleaning the code: make clean_all
#
#   This will remove all compiled codes, libraries and the MARE2DEM executable, 
#   so that only the originally downloaded source files are present.
#   You will need to do a make clean_all before compiling MARE2DEM if you 
#   update the compilers, update the mpi library etc.  	
#

#---------------------------------------------------------------------
# Check for input INCLUDE file argument:
#----------------------------------------------------------------------

INC_FILE := $(shell echo $(INCLUDE) | tr A-Z a-z)
PREP := $(shell echo $(PREP) | tr A-Z a-z)

ifndef INCLUDE
    MARE2DEM: display_Error_NoInclude  display_UsageMsg 
else ifeq (,$(wildcard $(INC_FILE))) # ifneq (,$(wildcard $(INC_FILE)))   
    MARE2DEM: display_Error_IncludeNotFound display_UsageMsg
else 
    include $(INC_FILE)
    MARE2DEM: disp build_mare2dem build_mare2dem_lib;
endif


#-------------------------------------------------------------------------------
# Libraries required by MARE2DEM:
#-------------------------------------------------------------------------------
 
# 
# SuperLU library:
# 
SUPERLU_dir     = ./libraries/SuperLU_4.3
LIBSUPERLU      = $(SUPERLU_dir)/libsuperlu_4.3.a
SUPERLU_HEADER = -I$(SUPERLU_dir)/SRC

        
# 
# ScaLAPACK Library"
#
SCALAPACK_dir =./libraries/scalapack-2.0.2
LIBSCALAPACK  = $(SCALAPACK_dir)/libscalapack.a

#-------------------------------------------------------------------------------
# Recipes:
#-------------------------------------------------------------------------------
 
# 
# Error message if INCLUDE variable not defined on input:
#
display_Error_NoInclude:   
	@printf "\n\n\n !!!    Error making MARE2DEM     !!!   \n\n";
	@printf "    INCLUDE argument not given \n\n";

display_UsageMsg:	 
	@printf "    Usage:   make INCLUDE=<./include/make.inc> \n\n";
	@printf "    <make.inc> should be the name of an include (text) file that specifies \n";
	@printf "    the MPI Fortran and C compilers and compilers flags for your \n";
	@printf "    computing system.  See ./include/macos.inc for an example. I recommend \n";
	@printf "    you create your own make.inc file (with a specific name such as \n";
	@printf "    make_MyMegaCluster.inc) for each system where you build MARE2DEM. \n";
	@printf "    \n";
	@printf "    Example make include files are in ./include/:\n";
	@printf "    \n";
	@printf "     macos.inc    - MacOS system with the intel compiler suite.\n";
	@printf "     habanero.inc - Linux cluster with the intel parallel studio compilers.\n";
	@printf "    \n";
	@printf "    Note that MARE2DEM currently requires the Intel Fortran and C \n";
	@printf "    C compilers and relies on the Intel Math Kernel Library (MKL).\n";
	@printf "    Either OpenMPI or MPICH MPI implementations work.\n\n\n";

 
display_Error_IncludeNotFound: 
	@printf "\n\n\n !!!    Error making MARE2DEM     !!!   \n\n";
	@printf "    INCLUDE file not found:       %s \n\n" $(INC_FILE) ;
	@printf "    Did you specify the correct file path and name? \n\n";
	display_UsageMsg;
 

# 
#  Cleaning functions:
# 
clean: clean_base

clean_all: clean_base clean_superlu clean_scalapack 

clean_base:
	@printf "\n#\n# Cleaning \n#\n\n"; 
	$(RM) *.o
	$(RM) *.mod
	$(RM) MARE2DEM;
	@printf "\n# All clean \n"; 



clean_scalapack:
	@printf "#\n#\n# Cleaning ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); pwd; make clean;  	
		
clean_superlu:
	@printf "#\n#\n# Cleaning SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); pwd; make clean; 
		
 
 
 
# 
# SuperLU Library build:
# 
# Note that these are only executed if libsuperlu.a and libmetis.a can't be found:
# For example, this is not called if $(LIBSUPERLU) points to
# your cluster's own superlu library that is outside MARE2DEM/Source
# 
# 
$(LIBSUPERLU): 
	@printf "#\n#\n# Making SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); \
	make superlulib CC=$(CC) CFLAGS='$(CFLAGS)' \
	FORTRAN=$(FC_NO_INSTR) FFLAGS='$(FFLAGS)' CDEFS=$(SUPERLU_CDEFS) BLASDEF=$(BLASDEF) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);

#
# ScaLAPACK library build:
#
$(LIBSCALAPACK):
	@printf "#\n#\n# Making ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); \
	make lib CC=$(CC) CCFLAGS='$(CFLAGS)' \
	FC=$(FC_NO_INSTR) FCFLAGS='$(FFLAGS)'  CDEFS=$(SUPERLU_CDEFS) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);


# 
# MARE2DEM build:
# 

disp:   
	@printf "#\n#\n# Making MARE2DEM: \n#\n#\n"


mare2dem_core = em_constants.o kdtree2.o fem2d_utilities.o  binarytree.o  call_triangle.o sort.o \
			  c_fortran_zgssv.o  superlu_zsolver.o intelmkl_solver.o quad.o\
			  string_helpers.o triangle.o mt1d.o kx_io.o em2dkx.o dc2dkx.o mare2dem_scalapack.o occam.o\
		      c_fortran_triangle.o filtermodules.o   \
		      mare2dem_common.o  spline_kx_module.o mare2dem_worker.o mare2dem_io.o \
		      mare2dem_mpi.o em2d.o
		      
mare2dem_exe =  $(mare2dem_core) runmare2dem.o 

mare2dem_lib = $(mare2dem_core)  mare2dem_lib_interface.o

build_mare2dem:	$(LIBSUPERLU)  $(LIBSCALAPACK) $(mare2dem_exe) 
	        $(FC) $(FFLAGS) $(mare2dem_exe)  $(LIBSUPERLU) $(LIBSCALAPACK) \
	        $(BLASLIB) -o MARE2DEM
	
build_mare2dem_lib:	$(LIBSCALAPACK) $(LIBSUPERLU) $(mare2dem_lib) 
	        $(FC) $(FFLAGS) $(mare2dem_lib)  $(LIBSUPERLU) $(LIBSCALAPACK) \
	        $(BLASLIB) -shared -o mare2dem_lib.o
	        	        
      	        
 

# 
# Test build:
# 
Test_Files	= binaryTree.o call_triangle.o triangle.o  c_fortran_triangle.o  \
		      TestForSlivers.o  
		
TestForSlivers:	 $(Test_Files) 
	$(FC) $(FFLAGS) $(Test_Files) -o TestForSlivers
	
test_lib_files = $(mare2dem_lib)  test_library.o

test_library:	 $(LIBSUPERLU) $(LIBSCALAPACK)  $(test_lib_files) 
	$(FC) $(FFLAGS) $(test_lib_files) $(LIBSUPERLU) $(LIBSCALAPACK) $(BLASLIB) -o test_library
	
# 
# Triangle build:
# 
TRILIBDEFS = -DTRILIBRARY    

triangle.o:  triangle.c  triangle.h
	$(CC)  $(TRILIBDEFS) $(TRICOPTS) -c -o $(BIN)triangle.o triangle.c
		
# 

#-------------------------------------------------------------------------------
# Compiling recipes:
#-------------------------------------------------------------------------------
#
# General Fortran compile:
%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $^

%.o: %.F90 
	$(FC) $(FFLAGS) -c -o $@ $^

# General C compile:
%.o : %.c
	$(CC) $(CFLAGS) $(CDEFS) $(SUPERLU_HEADER) -c -o $@ $< 


