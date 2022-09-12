#===============================================================================
# Copyright 2001-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library PBLAS C examples creation and run
##
##******************************************************************************
ROOTdir=$(shell pwd)
ifndef MKLROOT
ifndef LIBdir
MKLROOT = $(shell pwd)/../..
else
MKLROOT = $(LIBdir)/../..
endif
endif

include make.inc
include pblasc.lst

help:
	@echo "Usage:"
	@echo
	@echo "make [help]"
	@echo "   print this help."
	@echo
	@echo "make <targets> [mpi=mpi_name] [compiler=compiler_name]"
	@echo "      [interface=interface_name] [LIBdir=<libdir>] [MPIdir=<mpidir>] [example=example_list]"
	@echo "   makes examples. Make sure environment is set properly before testing."
	@echo "   If variables are not set, the defaults are:"
	@echo "      mpi=intelmpi - test with Intel(R) MPI"
	@echo "      compiler=intel - test with Intel(R) Compiler"
	@echo "      LIBdir=../../lib/intel64 - path to Intel(R) MKL binaries"
	@echo "      MPIdir not set - let the system locate MPI scripts like mpirun,"
	@echo "          set PATH, LD_LIBRARY_PATH properly before testing."
	@echo "   Variables could be set to:"
	@echo "      mpi=intelmpi - test with Intel(R) MPI"
	@echo "          Don't forget to start Intel(R) MPI daemon before testing."
	@echo "      mpi=openmpi - test with Open MPI"
	@echo "      mpi=mpich2 - test with MPICH2"
	@echo "      compiler=gnu - test with GNU C compiler, ilp64 libraries shouldn't be tested with gnu compiler."
	@echo "      compiler=intel - test with Intel(R) Compiler"
	@echo "      LIBdir=<MKLdir> - path to Intel(R) MKL binaries, if they're not in"
	@echo "          ../../lib/intel64."
	@echo "      MPIdir=<mpidir> - path to MPI. Specify MPIdir if it's not in"
	@echo "          the PATH. MPI scripts are taken from <mpidir>/bin (or bin64"
	@echo "          in case of Intel(R) 64)."
	@echo "          When testing with MPICH, make sure that MPICH was"
	@echo "          compiled by the same compiler as used to build the examples."
	@echo "      interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "      example_list - list of examples to be performed. You can specify the list by taking"
	@echo "          several example names into double quotes and separating it by spaces. If example_list is"
	@echo "          not specified - all examples will be performed by default. Please type 'make list' to"
	@echo "          see list of all possible example names. "
	@echo "   Targets:"
	@echo "      libintel64,sointel64 - to link with static/dynamic Intel(R) MKL"
	@echo "                     for Intel(R) 64 processor family"
	@echo "      cleanall - clean all objects"
	@echo "   Examples:"
	@echo "      make sointel64 mpi=mpich2 "
	@echo "      make libintel64 mpi=intelmpi example=pblas1_s_example "
	@echo "      make cleanall"


ifdef example
ifeq ("$(filter-out $(PBLAS),$(example))","")
   testset = $(foreach example_name,$(example),$(example_name).txt)
else
   $(error Invalid name of test. Possible values of 'example': $(PBLAS))
endif
else
   testset = $(foreach example_name,$(PBLAS),$(example_name).txt)
endif

list: ; @echo "Possible values of 'example': $(PBLAS)"

ifdef arch
libintel64: ; ($(MAKE) wrap_lib; cd $(SRCdir); $(MAKE) arch=intel64 dynamic=no  $(testset) MKLROOT=$(MKLROOT))
sointel64:  ; ($(MAKE) wrap_lib; cd $(SRCdir); $(MAKE) arch=intel64 dynamic=yes $(testset) MKLROOT=$(MKLROOT))
else
libintel64: ; ($(MAKE) arch=intel64 dynamic=no  libintel64)
sointel64:  ; ($(MAKE) arch=intel64 dynamic=yes sointel64)
endif

ifeq ($(mpi), custom)
wraplibdir = $(ROOTdir)/source/wrap_lib/$(arch)
wraplibname = libmkl_blacs_$(mpi)_$(interface)
export wraplibname wraplibdir

wrap_lib: $(wraplibdir)/$(wraplibname)

$(wraplibdir)/$(wraplibname):
	mkdir -p $(wraplibdir)
	cd $(MKLROOT)/interfaces/mklmpi && \
	$(MAKE) $(plat_s)$(arch) interface=$(interface) MPICC=$(MPICC) \
        INSTALL_LIBNAME=$(wraplibname) INSTALL_DIR=$(wraplibdir)
else
wrap_lib:
endif

cleanall: ; (cd $(SRCdir); $(MAKE) clean)
