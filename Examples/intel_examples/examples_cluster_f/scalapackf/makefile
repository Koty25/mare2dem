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
##     Intel(R) Math Kernel Library SCALAPACK examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {sointel64|libintel64} [function=name]"
	@echo "[mpi=intelmpi|mpich2|openmpi|custom] [compiler=intel|gnu]"
	@echo "[interface=interface_name] [workdir=path] [mpidir=path]"
	@echo "[threading=threading_name]"
	@echo
	@echo "function=name - name of function. Please see names in scalapack.lst file."
	@echo "If parameter is omitted, all examples will be compiled."
	@echo
	@echo "mpi=intelmpi - using Intel(R) MPI Library, default"
	@echo "mpi=mpich2   - using MPICH2"
	@echo "mpi=openmpi  - using Open MPI"
	@echo "mpi=custom   - using custom MPI (MPICC, MPIF90, MPIRUN should be set)"
	@echo
	@echo "compiler=intel - using Intel(R) Fortran Compiler, default"
	@echo "compiler=gnu   - using GNU gfortran compiler"
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
## make sointel64 mpi=intelmpi function=pdgetrf workdir=/share
##                 - build and run PDGETRF
##                   example for 64-bit applications in
##                   dir /share, dynamic linking, using
##                   Intel(R) MPI Library, Intel(R) Fortran Compiler
##
## make libintel64 mpi=mpich2 compiler=intel mpidir=/opt/mpich
##                 - build and run all examples of SCALAPACK for
##                   Intel(R) 64 processor family applications
##                   in current dir, static linking, using
##                   MPICH from /opt/mpich,
##                   Intel(R) Fortran Compiler
##------------------------------------------------------------------------------

include scalapack.lst

SHELL = /bin/sh
ENV = $(shell which env)

ifndef function
   function = $(SCALAPACK)
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

ifndef NPROC
   NPROC=4
endif

RES = $(addsuffix x.res,$(function))

ifndef MKLROOT
   MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

FOPTS =

ifeq ($(compiler),intel)
   IFACE_COMP_PART=intel
   IFACE_THREADING_PART=intel
   OMP_LIB = -qopenmp
   FOPTS += -warn all -warn errors -nogen-interfaces
endif

ifeq ($(compiler),gnu)
   IFACE_COMP_PART=gf
   IFACE_THREADING_PART=gnu
   OMP_LIB = -fopenmp
   FOPTS += -Wall
   ifeq ($(RES_EXT),so)
       LOPTS = -Wl,--no-as-needed
   endif
   FOPTS += -m64
endif

IFACE_LIB=$(MKL_PATH)/libmkl_$(IFACE_COMP_PART)_$(interface).$(EXT)
SCALAPACK_LIB=$(MKL_PATH)/libmkl_scalapack_$(interface).$(EXT)

ifeq ($(interface),ilp64)
   ifeq ($(compiler),intel)
       FOPTS += -i8
   endif
   ifeq ($(compiler),gnu)
       FOPTS += -fdefault-integer-8
   endif
   BLACS_PART=_ilp64
else
   FOPTS +=
   BLACS_PART=_lp64
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
   MKL_LIBS = -Wl,--start-group $(SCALAPACK_LIB) $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) $(BLACS_LIB) -Wl,--end-group
else
   PTHREAD_LIB += -lpthread
   MKL_LIBS = $(LOPTS) $(SCALAPACK_LIB) $(IFACE_LIB) $(THREADING_LIB) $(CORE_LIB) $(BLACS_LIB)
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
      CS = $(_BD)mpif90 -fc=ifort $(FOPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = $(_BD)mpif90 -fc=gfortran $(FOPTS)
   endif
   MPITHR =#
   FOPTS += -DMPI_KIND_=4
   RS = $(_BD)mpiexec -n $(NPROC)
   Bs = _intelmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)_$(interface).$(EXT)
endif

ifeq ($(mpi),intelmpi)
   ifdef mpidir
      _BD = $(mpidir)/$(Ibin)/
   endif
   ifeq ($(compiler),intel)
      CS = $(_BD)mpiifort $(FOPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = $(_BD)mpif90 -f90=gfortran $(FOPTS)
   endif
   ifeq ($(threading),parallel)
      MPITHR = -mt_mpi
   else
      MPITHR =#
   endif
   ifeq ($(interface),ilp64)
      ifeq ($(compiler),intel)
         # -i8 enables Intel MPI ilp64 interface
         FOPTS += -DMPI_KIND_=8
      else
         FOPTS += -DMPI_KIND_=4
      endif
   else
      FOPTS += -DMPI_KIND_=4
   endif
   RS = $(_BD)mpiexec -n $(NPROC)
   Bs = _intelmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
endif

ifeq ($(mpi),openmpi)
   ifdef mpidir
      _BD = $(mpidir)/bin/
   endif
   ifeq ($(compiler),intel)
      CS = OMPI_FC=ifort $(_BD)mpif90 $(FOPTS)
   endif
   ifeq ($(compiler),gnu)
      CS = OMPI_FC=gfortran $(_BD)mpif90 $(FOPTS)
   endif
   MPITHR =#
   FOPTS += -DMPI_KIND_=4
   RS = $(_BD)mpiexec -np 2
   Bs = _openmpi
   BLACS_LIB=$(MKL_PATH)/libmkl_blacs$(Bs)$(BLACS_PART).$(EXT)
endif

ifeq ($(mpi),custom)
   MPICC  ?= mpicc
   MPIF90 ?= mpif90
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
   MPITHR =#
   FOPTS += -DMPI_KIND_=4
   CS = $(_BD)$(MPIF90) $(FOPTS)
   RS = $(_BD)$(MPIRUN) -np $(NPROC)
   BLACS_LIB=$(wraplibdir)/$(wraplibname).$(EXT)
endif


#-------------------------------------------------------------------------------

vpath %.f source

.NOTPARALLEL: # as written, this can only be run with make -j1, to support creation of .mod files first

$(RES): %.res: %.f $(if $(mpi.custom), $(wraplibdir)/$(wraplibname))
	mkdir -p $(RES_DIR)
	$(CS) source/common/pcmatgen.f -c -o $(RES_DIR)/pcmatgen.o
	$(CS) source/common/pdmatgen.f -c -o $(RES_DIR)/pdmatgen.o
	$(CS) source/common/psmatgen.f -c -o $(RES_DIR)/psmatgen.o
	$(CS) source/common/pzmatgen.f -c -o $(RES_DIR)/pzmatgen.o
	$(CS) source/common/pmatgeninc.f -c -o $(RES_DIR)/pmatgeninc.o
	$(CS) source/common/average_time.f -c -o $(RES_DIR)/average_time.o
	$(CS) $< -c -o $(RES_DIR)/$*.o
	$(CS) $(RES_DIR)/$*.o $(RES_DIR)/pcmatgen.o $(RES_DIR)/pdmatgen.o $(RES_DIR)/psmatgen.o $(RES_DIR)/pzmatgen.o $(RES_DIR)/pmatgeninc.o $(RES_DIR)/average_time.o -o $(RES_DIR)/$*.exe -L$(MKL_PATH) $(MKL_LIBS) $(MPITHR) -ldl -lm
	$(RS) $(ENV) LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(CMPLR_PATH) $(RES_DIR)/$*.exe in/$*.in $(RES_DIR)/$@
	-rm $(RES_DIR)/*.o

#-------------------------------------------------------------------------------
