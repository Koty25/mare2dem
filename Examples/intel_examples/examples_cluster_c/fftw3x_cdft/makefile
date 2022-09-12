#===============================================================================
# Copyright 2006-2020 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

##  Content:
##      Intel(R) Math Kernel Library
##      C-style Cluster CDFT wrappers MPI FFTW3 examples creation and run.
##******************************************************************************

help:
	@echo "Usage: make libintel64"
	@echo "		[mpi=intelmpi|mpich2|openmpi] [mpidir=path]"
	@echo "		[compiler=intel|gnu] [workdir=path]"
	@echo
	@echo "mpi=intelmpi - using Intel(R) MPI Library, default."
	@echo "mpi=mpich2   - using MPICH2."
	@echo "mpi=openmpi  - using Open MPI."
	@echo
	@echo "mpidir=path - path to MPI installation directory. MPI scripts are taken"
	@echo "from mpidir/bin. If this directory is in PATH you can omit mpidir parameter."
	@echo "When using MPICH, make sure that MPICH was compiled by same compiler"
	@echo "as used to build examples."
	@echo
	@echo "compiler=intel - using Intel(R) C Compiler, default."
	@echo "compiler=gnu   - using GNU C compiler."
	@echo
	@echo "workdir=path - path to work directory, which is accessible from any node."
	@echo "If parameter is omitted, executable files and results will be located in"
	@echo "current directory."
	@echo
	@echo "Set environment variables LD_LIBRARY_PATH and etc properly before testing."

##*****************************************************************************
## examples of using:
##
## make libintel64
##              - compile and run wrappers examples in current dir,
##                static linking, using Intel(R) MPI Library,
##                Intel(R) C Compiler, for double precision data for
##                Intel(R) 64 based applications
##
## make libintel64 mpi=mpich2 mpidir=/opt/mpich compiler=gnu
##              - compile and run wrappers examples in current dir,
##                static linking, using MPICH from /opt/mpich,
##                GNU C compiler, for double precision data for
##                Intel(R) 64 based applications
##*****************************************************************************

include cdft_wrappers.lst

example = $(WRAPPER)

ifndef mpi
   mpi = intelmpi
endif

ifndef compiler
   compiler = intel
endif

ifndef workdir
   workdir = .
endif

RES = $(addsuffix .res,$(example))

ifndef MKLROOT
   MKLROOT = ../..
endif

MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

INSTALL_DIR="$(PWD)"

COPTS = -Wall -Werror

IFACE_LIB = $(MKL_PATH)/libmkl_intel_lp64.$(EXT)
THREADING_LIB = $(MKL_PATH)/libmkl_intel_thread.$(EXT)
CORE_LIB = $(MKL_PATH)/libmkl_core.$(EXT)
BLACS_PART = _lp64
BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
OMP_LIB = -L$(CMPLR_PATH) -liomp5
MKL_LIBS = -Wl,--start-group $(MKL_PATH)/libmkl_cdft_core.$(EXT) $(IFACE_LIB) $(THREADING_LIB) \
	$(CORE_LIB) $(BLACS_LIB) -Wl,--end-group $(OMP_LIB)

RES_DIR = _results/$(compiler)_$(_IA)_$(RES_EXT)

ifeq ($(mpi),mpich2)
   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   ifeq ($(compiler),intel)
      _CS = $(_BD)mpicc -cc=icc $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      _CS = $(_BD)mpicc -cc=gcc $(COPTS)
   endif
   _RS = $(_BD)mpiexec -n 2
   Bs = _intelmpi
   export OMP_NUM_THREADS=1
endif

ifeq ($(mpi),intelmpi)
   ifdef mpidir
      _BD = $(mpidir)/$(_IA)/bin/
   endif
   ifeq ($(compiler),intel)
      _CS = $(_BD)mpiicc -mt_mpi $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      _CS = $(_BD)mpicc -mt_mpi $(COPTS)
   endif
   _RS = $(_BD)mpiexec -n 2
   Bs = _intelmpi
endif

ifeq ($(mpi),openmpi)
   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   ifeq ($(compiler),intel)
      _CS = OMPI_CC=icc $(_BD)mpicc $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      _CS = OMPI_CC=gcc $(_BD)mpicc $(COPTS)
   endif
   _RS = $(_BD)mpiexec -np 2
   Bs = _openmpi
   DelF = rm $(RES_DIR)/*.o
endif

#-------------------------------------------------------------------------------

vpath %.c source

libintel64:
	$(MAKE) wrap_lib _IA=intel64
	$(MAKE) $(RES)   _IA=intel64 EXT=a RES_EXT=lib

wrap_lib:
	cd "$(MKLROOT)"/interfaces/fftw3x_cdft && \
		$(MAKE) lib$(_IA) mpi=$(mpi) compiler=$(compiler) mpidir="$(mpidir)" \
		INSTALL_DIR="$(INSTALL_DIR)/lib/$(_IA)" INSTALL_LIBNAME="libfftw3x_cdft"

.NOTPARALLEL: # as written, this can only be run with make -j1

$(RES): %.res: %.c
	mkdir -p $(RES_DIR)
	$(_CS) $(FL) $< -I $(MKLROOT)/include -I $(MKLROOT)/include/fftw \
		-c -o $(RES_DIR)/$*.o
	$(_CS) $(FL) $(RES_DIR)/$*.o -O0 \
		-L"$(INSTALL_DIR)/lib/$(_IA)" -lfftw3x_cdft \
		-L$(MKLROOT)/lib/$(_IA) $(MKL_LIBS) \
		-lpthread -lm -ldl -o $(RES_DIR)/$*.exe
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":$(CMPLR_PATH); \
		$(_RS) $(RES_DIR)/$*.exe > $(RES_DIR)/$@;
	-$(DelF)

.PHONY: clean
clean:
	-rm -rf _results;
#-------------------------------------------------------------------------------
