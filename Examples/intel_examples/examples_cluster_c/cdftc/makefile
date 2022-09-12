#===============================================================================
# Copyright 2005-2020 Intel Corporation.
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
##     Intel(R) Math Kernel Library C-style Cluster DFT
##     examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {sointel64|libintel64} [example=name]"
	@echo "[mpi=intelmpi|mpich2|openmpi|custom] [compiler=intel|gnu]"
	@echo "[interface=interface_name] [workdir=path] [mpidir=path]"
	@echo "[threading=threading_name]"
	@echo
	@echo "example=name - name of example. Please see names in cdftc.lst file."
	@echo "If parameter is omitted, all examples will be compiled."
	@echo
	@echo "mpi=intelmpi - using Intel(R) MPI Library, default"
	@echo "mpi=mpich2   - using MPICH2"
	@echo "mpi=openmpi  - using Open MPI"
	@echo "mpi=custom   - using custom MPI (MPICC, MPIRUN should be set)"
	@echo
	@echo "compiler=intel - using Intel(R) C Compiler, default"
	@echo "compiler=gnu   - using GNU C compiler"
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "threading_name - can be parallel or sequential. Default value is parallel."
	@echo
	@echo "workdir=path - path to work directory, which is accessible from any node."
	@echo "If parameter is omitted, executable files and results will be located in"
	@echo "current directory."
	@echo
	@echo "mpidir=path - path to MPI installation directory. MPI scripts are taken"
	@echo "from mpidir/bin (or mpidir/bin64 for Intel(R) MPI Library and Intel(R) 64 architecture). If this directory"
	@echo "is in PATH you can omit mpidir parameter."
	@echo "MPICH should be compiled by the compiler with same version as compiler used to build examples."
	@echo
	@echo "Set environment variables LD_LIBRARY_PATH and etc properly before testing."

##------------------------------------------------------------------------------
## examples of using:
##
## make libintel64 example=dm_complex_2d_double_ex1 workdir=/share
##                 - build and run DM_COMPLEX_2D_DOUBLE_EX1
##                   example for 64-bit applications in
##                   dir /share, static linking, using
##                   Intel(R) MPI Library, Intel(R) C Compiler
##
## make sointel64 mpi=mpich2 compiler=gnu mpidir=/opt/mpich
##                 - build and run all examples of CDFT for
##                   Intel(R) 64 processor family applications
##                   in current dir, dynamic linking, using
##                   MPICH from /opt/mpich, GNU C compiler
##------------------------------------------------------------------------------

include cdftc.lst

SHELL = /bin/sh
ENV = $(shell which env)

ifndef example
   example = $(CDFT)
endif

ifndef mpi
   mpi = intelmpi
endif
mpi.$(mpi) = $(mpi)

ifndef compiler
   compiler = intel
endif

ifndef workdir
   workdir = .
endif

ifndef interface
   interface=lp64
endif

ifndef threading
   threading=parallel
endif

RES = $(addsuffix .res,$(example))

ifndef MKLROOT
   MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

COPTS = -Wall -Werror

ifeq ($(compiler),intel)
   IFACE_COMP_PART=intel
   IFACE_THREADING_PART=intel
   OMP_LIB = -liomp5
endif

LOPTS = 

ifeq ($(compiler),gnu)
   IFACE_COMP_PART=intel
   IFACE_THREADING_PART=gnu
   OMP_LIB = -fopenmp
   ifeq ($(RES_EXT),so)
       LOPTS = -Wl,--no-as-needed
   endif
   COPTS += -m64
endif

ifeq ($(_IA),intel64)
   ifeq ($(interface),ilp64)
      COPTS += -DMKL_ILP64
      IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)_ilp64.$(EXT)
      BLACS_PART=_ilp64
   else
      COPTS +=
      IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)_lp64.$(EXT)
      BLACS_PART=_lp64
   endif
endif

ifeq ($(threading),sequential)
   THREADING_LIB=$(MKL_PATH)/libmkl_sequential.$(EXT)
else
   THREADING_LIB=$(MKL_PATH)/libmkl_$(IFACE_THREADING_PART)_thread.$(EXT)
endif

CORE_LIB=$(MKL_PATH)/libmkl_core.$(EXT)

ifeq ($(EXT),a)
   # Required to import weak symbols if static version of the pthread library is linked
   PTHREAD_LIB += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
   MKL_LIBS = -Wl,--start-group $(MKL_PATH)/libmkl_cdft_core.$(EXT) $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) $(BLACS_LIB) -Wl,--end-group
else
   PTHREAD_LIB += -lpthread
   MKL_LIBS=$(LOPTS) $(MKL_PATH)/libmkl_cdft_core.$(EXT) $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) $(BLACS_LIB)
endif

ifeq ($(threading),parallel)
   MKL_LIBS += $(OMP_LIB)
endif
MKL_LIBS += $(PTHREAD_LIB)

RES_DIR=_results/$(compiler)_$(mpi)_$(interface)_$(_IA)_$(RES_EXT)_$(threading)$Z

libintel64:
	$(MAKE) $(RES) _IA=intel64 Ibin=bin64 EXT=a RES_EXT=lib
sointel64:
	$(MAKE) $(RES) _IA=intel64 Ibin=bin64 EXT=so RES_EXT=so

.PHONY: clean
clean:
	-rm -rf _results;

ifeq ($(mpi),mpich2)
   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   ifeq ($(compiler),intel)
      CS = $(_BD)mpicc -cc=icc $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = $(_BD)mpicc -cc=gcc $(COPTS)
   endif
   MPITHR =#
   RS = $(_BD)mpiexec -n 2
   Bs = _intelmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
endif

ifeq ($(mpi),intelmpi)
   ifdef mpidir
      _BD = $(mpidir)/$(Ibin)/
   endif
   ifeq ($(compiler),intel)
      CS = $(_BD)mpiicc $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = $(_BD)mpicc -cc=gcc $(COPTS)
   endif
   ifeq ($(threading),parallel)
      MPITHR = -mt_mpi
#      MPITHR =#
   else
      MPITHR =#
   endif
   RS = $(_BD)mpiexec -n 2
   Bs = _intelmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
endif

ifeq ($(mpi),openmpi)
   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   ifeq ($(compiler),intel)
      CS = OMPI_CC=icc $(_BD)mpicc $(COPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = OMPI_CC=gcc $(_BD)mpicc $(COPTS)
   endif
   MPITHR =#
   RS = $(_BD)mpiexec -np 2
   Bs = _openmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
endif

ifeq ($(mpi),custom)
   MPICC  ?= mpicc
   MPIRUN ?= mpirun
   export MPICC
   wraplibdir = $(RES_DIR)
   wraplibname = libmkl_blacs_custom_$(interface)
   fullwraplibdir = $(realpath .)/$(wraplibdir)

$(wraplibdir)/$(wraplibname):
	mkdir -p $(wraplibdir)
	cd $(MKLROOT)/interfaces/mklmpi && \
	$(MAKE) $(RES_EXT)$(_IA) INSTALL_LIBNAME=$(wraplibname) INSTALL_DIR=$(fullwraplibdir)

   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   CS = $(_BD)$(MPICC) $(COPTS)
   MPITHR =#
   RS = $(_BD)$(MPIRUN) -np 2
   BLACS_LIB=$(wraplibdir)/$(wraplibname).$(EXT)
endif

COMMON = source/cdft_example_support.c

#-------------------------------------------------------------------------------

vpath %.c source

.NOTPARALLEL: # as written, this can only safely be run with 1 thread (make -j1) due to creation of COMMON obj files for each target

$(RES): %.res: %.c $(if $(mpi.custom), $(wraplibdir)/$(wraplibname))
	mkdir -p $(RES_DIR)
	$(CS) $(SPEC_OPT) $(COMMON) -I "$(MKLROOT)/include" -c -o $(RES_DIR)/cdft_example_support.o;
	$(CS) $(SPEC_OPT) $< -I "$(MKLROOT)/include" -c -o $(RES_DIR)/$*.o
	$(CS) $(SPEC_OPT) $(RES_DIR)/$*.o $(RES_DIR)/cdft_example_support.o -L $(MKL_PATH) $(MKL_LIBS) $(MPITHR) -ldl -lm -o $(RES_DIR)/$*.exe
	$(RS) $(ENV) LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(CMPLR_PATH) $(RES_DIR)/$*.exe data/$*.dat > $(RES_DIR)/$@
	-rm $(RES_DIR)/$*.o

#-------------------------------------------------------------------------------
