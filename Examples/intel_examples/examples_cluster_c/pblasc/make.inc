# vim: syntax=make
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
############################################################################
#  Identifiers
#
#  arch
#  mpi
#  comp
#  opt
#  IDE
#  interface
############################################################################
# Intel architecture (must be defined)
#arch=intel64

# MPI
ifndef mpi
mpi = intelmpi
endif
ifneq ("$(mpi)","intelmpi")
ifneq ("$(mpi)","openmpi")
ifneq ("$(mpi)","mpich2")
ifneq ("$(mpi)","custom")
$(error ERROR: You set mpi to "$(mpi)". Allowed values are: intelmpi, openmpi, mpich2 and custom)
endif
endif
endif
endif

# Compiler
comp=intel
ifeq ($(compiler),gnu)
comp=gnu
endif

ifndef interface
interface = lp64
endif

ifeq ($(interface),ilp64)
ADD_IFACE = _ilp64
IFACE_PART = _ilp64
COPTS = -DMKL_ILP64
FOPTS = -i8
else
ADD_IFACE = _lp64
IFACE_PART = _lp64
COPTS =
FOPTS =
endif

ifeq ($(compiler),gnu)
	FOPTS+= -fno-f2c
endif

# Optimization level
opt=noopt
# executable identifier
ifeq ($(dynamic),yes)
plat_s=so
else
plat_s=lib
endif
ifeq ($(arch),intel64)
plat_i=$(plat_s)intel64
endif
IDE=_$(plat_i)_$(mpi)_$(comp)_$(opt)$(ADD_IFACE)
ifeq ($(compiler),gfortran)
IDE=_$(plat_i)_$(mpi)_$(compiler)_$(opt)$(ADD_IFACE)
endif
############################################################################
#  Paths
#
#  ROOTdir
#  LIBdir
#  SRCdir
#  MKLINC
############################################################################
#ROOTdir

ifeq ($(LIBdir),)
LIBpath=-L$(MKLROOT)/lib/$(arch)
MKLINC=$(MKLROOT)/include
else
LIBpath=-L$(LIBdir)
MKLINC=$(LIBdir)/../../include
endif

SRCdir=source
MKL_PATH=$(patsubst -L%,%,$(LIBpath))
############################################################################
#  The libraries to be linked to
#
#  LIBS
############################################################################
# BLACS

ifeq ($(dynamic),yes)
EXT_LIB=so
else
EXT_LIB=a
endif
ifeq ($(mpi),intelmpi)
BLACSLIB=$(MKL_PATH)/libmkl_blacs_intelmpi$(IFACE_PART).$(EXT_LIB)
else
ifeq ($(mpi),openmpi)
BLACSLIB=$(MKL_PATH)/libmkl_blacs_openmpi$(IFACE_PART).$(EXT_LIB)
#BLACSLIB=$(MKL_PATH)/libmkl_blacs_openmpi$(IFACE_PART).a
else
ifeq ($(mpi),mpich2)
BLACSLIB=$(MKL_PATH)/libmkl_blacs_intelmpi$(IFACE_PART).$(EXT_LIB)
#BLACSLIB=$(MKL_PATH)/libmkl_blacs_intelmpi$(IFACE_PART).a
else
ifeq ($(mpi), custom)
BLACSLIB = $(wraplibdir)/$(wraplibname).$(EXT_LIB)
else
BLACSLIB=$(MKL_PATH)/libmkl_blacs$(IFACE_PART).$(EXT_LIB)
#BLACSLIB=$(MKL_PATH)/libmkl_blacs$(IFACE_PART).a
endif
endif
endif
endif


FBLACSLIB=$(BLACSLIB)
# ScaLAPACK
SCALAPACKLIB= $(MKL_PATH)/libmkl_scalapack$(IFACE_PART).$(EXT_LIB)

# LAPACK and BLAS
LAPACKLIB=
ifeq ($(dynamic),yes)
MKL_LIBS=-L$(MKL_PATH) \
    -Wl,-rpath,$(MKL_PATH) \
    -Wl,--start-group \
        $(SCALAPACKLIB) $(FBLACSLIB) \
        $(MKL_PATH)/libmkl_intel$(ADD_IFACE).so \
        $(MKL_PATH)/libmkl_sequential.so \
        $(MKL_PATH)/libmkl_core.so \
    -Wl,--end-group -lpthread
else
MKL_LIBS=-L$(MKL_PATH) \
    -Wl,--start-group \
        $(SCALAPACKLIB) $(FBLACSLIB) \
        $(MKL_PATH)/libmkl_intel$(ADD_IFACE).a \
        $(MKL_PATH)/libmkl_sequential.a \
        $(MKL_PATH)/libmkl_core.a \
    -Wl,--end-group -lpthread
endif

# All libraries to link
LIBS=$(LIBpath) $(TOOLSLIB) $(MKL_LIBS) $(INSLIB) $(MPILIB)
ifeq ($(comp),gnu)
LIBS+= -lm -ldl 
endif

############################################################################
#  Tools
#
#  MAKE
#  CC
#  CCFLAGS
#  F77
#  F77FLAGS
#  MPIRUN
#  CCLOADER
#  F77LOADER
############################################################################
# Fortran, mpirun
ifndef MPIRUNOPTS
MPIRUNOPTS=-n 4
endif

MPICC  = mpicc
MPIF77 = mpif77

ifeq ($(MPIdir),)
	ifeq ("$(findstring intelmpi,$(mpi))","intelmpi")
		ifeq ($(comp),gnu)
			F77=mpif77 $(FOPTS)
			CC =mpicc $(COPTS)
		endif
		ifeq ($(comp),intel)
			F77=mpiifort $(FOPTS)
			CC =mpiicc $(COPTS)
		endif
		MPIRUN=mpiexec $(MPIRUNOPTS)
	else
		F77= $(MPIF77) $(FOPTS)
		CC = $(MPICC) $(COPTS)
		ifndef MPIRUN
		MPIRUN=mpirun -np 4
		endif
	endif
else
	ifeq ("$(findstring intelmpi,$(mpi))","intelmpi")
		MPIBinOld = $(MPIdir)/bin64
		MPIBinNew = $(MPIdir)/$(arch)/bin
		MPIBinDir = $(shell if [ -d $(MPIBinNew) ] ; then echo $(MPIBinNew) ; elif [ -d $(MPIBinOld) ] ; then echo $(MPIBinOld) ; fi )
		ifeq ($(comp),gnu)
			F77=$(MPIBinDir)/mpif77 $(FOPTS)
			CC =$(MPIBinDir)/mpicc $(COPTS)
		endif
		ifeq ($(comp),intel)
			F77=$(MPIBinDir)/mpiifort $(FOPTS)
			CC =$(MPIBinDir)/mpiicc $(COPTS)
		endif
		MPIRUN=$(MPIBinDir)/mpiexec $(MPIRUNOPTS)
	else
		F77=$(MPIdir)/bin/$(MPIF77) $(FOPTS)
		CC =$(MPIdir)/bin/$(MPICC) $(COPTS)
		ifndef MPIRUN
		MPIRUN=$(MPIdir)/bin/mpirun -np 4
		endif
	endif
endif

ifeq ($(mpi),mpich2)
ifeq ($(comp),intel)
F77=mpif90 $(FOPTS)
CC =mpicc $(COPTS)
else
ifeq ($(comp),gnu)
ifeq ($(interface),ilp64)
   $(error "ilp64 libraries should not be tested with gnu compiler")
endif
F77=mpif77 $(FOPTS)
endif
endif
endif
F77LOADER=$(F77)
CCLOADER=$(CC)
# Optimization
ifeq ($(opt),opt)
OPTS=-O2
else	# noopt
OPTS=-O0
endif
# Compiler flags
F77FLAGS=-w -u $(OPTS)
CCFLAGS =-w -u $(OPTS)
ifeq ($(comp),gnu)
F77FLAGS=$(OPTS) -ff90 -fno-globals -Wno-globals
endif

# Make options
MAKE := $(MAKE) MPIRUN="$(MPIRUN)" MPIRUNOPTS="$(MPIRUNOPTS)" ROOTdir=$(ROOTdir) MPIdir=$(MPIdir) arch=$(arch) mpi=$(mpi) comp=$(comp) opt=$(opt) dynamic=$(dynamic) interface=$(interface)

RES_DIR = ../_results/$(IDE)
