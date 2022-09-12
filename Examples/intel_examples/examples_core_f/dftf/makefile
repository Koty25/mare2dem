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

## Content:
##      Build and run Intel(R) Math Kernel Library Fortran-style DFT examples.
##*****************************************************************************

help:
	@echo "To run DFTF examples:"
	@echo "  make {libia32|soia32|libintel64|sointel64}"
	@echo "       [compiler=<name>] [interface=<name>] [threading=<name>]"
	@echo "       [parallel=<name>] [omp=<name>] [function=<name>]"
	@echo
	@echo "To get report on run examples:"
	@echo "  make report"
	@echo
	@echo "To clean results:"
	@echo "  make clean"
	@echo
	@echo "To get help just run make or:"
	@echo "  make help"
	@echo
	@echo "Main options:"
	@echo "  targets lib%   use static linkage"
	@echo "  targets so%    use dynamic linkage"
	@echo
	@echo "  compiler=<name> selects the compiler to build the examples:"
	@echo "      compiler=gnu   - GNU gfortran (GNU g77 is not supported)"
	@echo "      compiler=pgi   - PGI pgf95. Available only for intel64."
	@echo "      compiler=intel - Intel(R) Fortran compiler ifort (DEFAULT)"
	@echo
	@echo "  interface=<name> selects kind of MKL_INT type for %intel64 targets:"
	@echo "      interface=lp64  - 32-bit integers (DEFAULT)"
	@echo "      interface=ilp64 - 64-bit integers"
	@echo
	@echo "  threading=<name> selects threading of MKL:"
	@echo "      threading=parallel   - multithreaded version (DEFAULT)"
	@echo "      threading=sequential - sequential version"
	@echo
	@echo "  parallel=<name> selects Intel(R) MKL threading layer for threading=parallel:"
	@echo "      parallel=intel - libmkl_intel_thread"
	@echo "      parallel=gnu   - libmkl_gnu_thread"
	@echo "      parallel=pgi   - libmkl_pgi_thread"
	@echo "      Default value depends on the setting of compiler=<name>"
	@echo
	@echo "  omp=<name> selects OpenMP runtime library for threading=parallel:"
	@echo "      omp=iomp5 - Intel OpenMP runtime"
	@echo "      omp=gomp  - GNU OpenMP runtime"
	@echo "      omp=pgmp  - PGI OpenMP runtime"
	@echo "      Default value depends on the setting of parallel=<name>"
	@echo
	@echo "  function=<name> selects examples to execute"
	@echo "      Default value: all examples listed in file dftf.lst"
	@echo
	@echo "Additional options:"
	@echo "  RES_DIR=<path> defines where to place the results"
	@echo "      Default value: ./_results"
	@echo
	@echo "  MKLROOT=<path> defines alternative Intel(R) MKL root directory"
	@echo "      Default value: ../.."
	@echo
	@echo "  MKLRUNLIB_PATH=<path> defines alternative runtime library directory"
	@echo "      for dynamic linkage."
	@echo '      Default value: $$(MKLROOT)/lib/$$(IA), with $$(IA) defined by the target'
	@echo
	@echo "  LIB_PATH=<path> defines location of OpenMP runtime libraries"
	@echo "      needed by GNU or PGI threading layer. Set this value if"
	@echo "      respective libraries (libgomp, libpgmp) cannot be found in"
	@echo "      LIBRARY_PATH and LD_LIBRARY_PATH"
	@echo
	@echo "  TARGET_ARCH=<flags> defines additional compiler flags"
	@echo "      Refer to the compiler documentation about the architecture specific"
	@echo "      flags. For example, latest Intel compilers need TARGET_ARCH=-mAVX"
	@echo "      to generate Intel(R) AVX instructions."

##-----------------------------------------------------------------------------
## Usage examples:
##
## make libia32 function=basic_dp_complex_dft_1d
##     Build with Intel(R) Fortran compiler and run this example
##     as a 32-bit application with static linkage
##
## make soia32 compiler=gnu
##     Build with GNU Fortran compiler and run all examples as 32-bit applications
##     with dynamic linkage
##
## make libintel64 compiler=gnu
##     Build with GNU Fortran compiler and run all examples as statically linked
##     applications for Intel(R) 64 processor family
##
## make sointel64
##     Build with Intel(R) Fortran compiler and run all examples as
##     dynamically linked applications for Intel(R) 64 processor family
##-----------------------------------------------------------------------------

MY_MAKEFILE := $(MAKEFILE_LIST)

ifeq (,$(RUN)$(IA))
  target := $(firstword $(filter lib% so%,$(MAKECMDGOALS)))
else
  target = $(RUN)$(IA)
endif

MY_WAY := $(findstring target-specific,$(.FEATURES))

ifndef MKLROOT
  MKLROOT = ../..
endif

ifndef MKLRUNLIB_PATH
  MKLRUNLIB_PATH = $(MKLROOT)/lib/$(IA)
endif

include dftf.lst

ifndef RES_DIR
  RES_DIR = _results
endif

ifndef function
  function = $(DFT)
endif

RES = $(addsuffix .res,$(function))

ifeq (,$(filter gnu pgi,$(compiler)))
  override compiler = intel
  override parallel = intel
endif

ifeq ($(compiler)_$(_IA),pgi_ia32)
   $(info ***)
   $(info *** PGI compiler 32-bit support was dropped off starting from Intel(R) MKL 2017.)
   $(info *** Please try >make help)
   $(info ***)
   $(error )
endif

ifeq (,$(filter lp64 ilp64,$(interface)))
  override interface = lp64
endif
ifeq (,$(findstring 32,$(target)))
  iface = $(interface)_
else
  override interface =
endif

ifeq (,$(filter parallel sequential,$(threading)))
  override threading = parallel
endif

ifneq ($(threading),sequential)
_parallel_intel = intel
_parallel_gnu = gnu
_parallel_pgi = pgi

_omp_intel = iomp5
_omp_gnu = gomp
_omp_pgi = pgmp

parallel = $(_parallel_$(compiler))
omp = $(_omp_$(parallel))
endif

ifeq (,$(filter gomp pgmp,$(omp)))
  override omp = iomp5
endif

ifeq (,$(filter gnu pgi,$(parallel)))
  override parallel = intel
  override omp = iomp5
else
  ifeq ($(compiler)-$(parallel),gnu-gnu)
    ifneq ($(omp),gomp)
      override omp = iomp5
    endif
  else
    ifeq ($(compiler)-$(parallel),pgi-pgi)
      override omp = pgmp
      override MKL_LIBS_LINK = mixed
    else
      override omp = iomp5
    endif
  endif
endif
LOPTS =
ifeq ($(compiler),gnu)
  FC = gfortran
  FOPTS.ia32 = -m32 -fcray-pointer
  FOPTS.intel64 = -m64 -fcray-pointer
  FOPTS.DIAG = -Wall -Werror
  FOPTS = -I$(res_dir) # not needed for gfortran 4.6
  FOPTS += -J$(res_dir)
  FOPTS += -Isource # needed for gfortran 4.0 or older
  IFACE_COMP_PART = gf
  LOPTS = -Wl,--no-as-needed
else
  ifeq ($(compiler),pgi)
    FC = pgf95
    FOPTS.ia32 = -tp px
    FOPTS.intel64 = -tp x64
    FOPTS.DIAG = 
  else
    FC = ifort
    FOPTS.DIAG = -warn all -warn errors
  endif
  FOPTS = -module $(res_dir)
  IFACE_COMP_PART = intel
endif

FOPTS += -I$(MKLROOT)/include

ifeq ($(compiler)-$(parallel),gnu-gnu)
  IFACE_THREADING_PART = gnu
else
  ifeq ($(compiler)-$(parallel),pgi-pgi)
    IFACE_THREADING_PART = pgi
  else
    IFACE_THREADING_PART = intel
  endif
endif

IFACE_LIB.ia32 = mkl_$(IFACE_COMP_PART)
ifeq ($(interface),ilp64)
  ifeq ($(compiler),gnu)
    FOPTS.intel64 += -fdefault-integer-8
  else
    FOPTS.intel64 += -i8
  endif
  IFACE_LIB.intel64 = mkl_$(IFACE_COMP_PART)_ilp64
else
  IFACE_LIB.intel64 = mkl_$(IFACE_COMP_PART)_lp64
endif
IFACE_LIB = $(IFACE_LIB.$(IA))

FOPTS += $(FOPTS.$(IA)) $(FOPTS.DIAG)

ifeq ($(compiler),intel)
  FOPTS += -fpp
endif

ifeq ($(compiler),gnu)
  FOPTS += -cpp
endif

ifeq ($(compiler),pgi)
  FOPTS += -Mnokeepobj -Mpreprocess
endif

ifeq ($(threading),sequential)
  threadname = $(threading)
  threadlayer = $(threading)
  THREADING_LIB = mkl_sequential
  OMP_LIB =
else
  threadname = $(threading)_$(omp)
  threadlayer = $(parallel)
  THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
  ifeq ($(omp),iomp5)
    OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
  else
    ifneq ($(LIB_PATH),)
      OMP_LIB = -L$(LIB_PATH)
    endif
    OMP_LIB += -l$(omp)
  endif
  
  ifeq ($(compiler)-$(IFACE_THREADING_PART),pgi-pgi)
    FOPTS += -mp
  endif

  ifeq ($(compiler),intel)
    FOPTS += -qopenmp
  endif

  ifeq ($(compiler),gnu)
    FOPTS += -fopenmp
  endif
endif

CORE_LIB = mkl_core


MKL_PATH = $(MKLROOT)/lib/$(IA)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(IA)

res_dir = $(RES_DIR)/$(compiler)_$(iface)$(threadname)_$(target)

ifdef VERBOSE
  $(info )
  $(info MAKELEVEL=$(MAKELEVEL) target=$(target))
  $(info MKLROOT=$(MKLROOT) RES_DIR=$(RES_DIR))
  $(info compiler=$(compiler) interface=$(interface) threading=$(threading) \
         parallel=$(parallel) omp=$(omp))
endif

##-----------------------------------------------------------------------------
## Rules
.PHONY: libia32 soia32 libintel64 sointel64 \
	lib32 so32 libem64t soem64t \
	help run mkdir clean report prepare_run \
	FORCE

.SUFFIXES:
.SUFFIXES: .f90 .o .res

vpath %.f90 $(MKLROOT)/include
vpath %.f90 source

ifeq ($(RUN),so)
  ifneq ($(LIB_PATH),)
    RUNENV = LD_LIBRARY_PATH=$(LIB_PATH):$(MKLRUNLIB_PATH):$(CMPLR_PATH):$(LD_LIBRARY_PATH)
  else
    RUNENV = LD_LIBRARY_PATH=$(MKLRUNLIB_PATH):$(CMPLR_PATH):$(LD_LIBRARY_PATH)
  endif
  ifneq ($(MKL_LIBS_LINK),mixed)
    MKL_LIBS = $(LOPTS) -L$(MKL_PATH) -lmkl_rt
    RUNENV += MKL_INTERFACE_LAYER=$(interface) MKL_THREADING_LAYER=$(threadlayer)
  else
    MKL_LIBS = $(LOPTS) -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
  endif
else
  ifneq ($(LIB_PATH),)
    RUNENV = LD_LIBRARY_PATH=$(LIB_PATH):$(CMPLR_PATH):$(LD_LIBRARY_PATH)
  else
    RUNENV = LD_LIBRARY_PATH=$(CMPLR_PATH):$(LD_LIBRARY_PATH)
  endif
  MKL_LIBS = $(MKL_PATH)/lib$(IFACE_LIB).a
  MKL_LIBS += -Wl,--start-group
  MKL_LIBS += $(MKL_PATH)/lib$(THREADING_LIB).a $(MKL_PATH)/lib$(CORE_LIB).a
  MKL_LIBS += -Wl,--end-group
endif

libia32 lib32:
	$(MAKE) -f $(MY_MAKEFILE) prepare_run IA=ia32 RUN=lib
	$(MAKE) -f $(MY_MAKEFILE) run IA=ia32 RUN=lib

soia32 so32:
	$(MAKE) -f $(MY_MAKEFILE) prepare_run IA=ia32 RUN=so
	$(MAKE) -f $(MY_MAKEFILE) run IA=ia32 RUN=so

libintel64 libem64t:
	$(MAKE) -f $(MY_MAKEFILE) prepare_run IA=intel64 RUN=lib
	$(MAKE) -f $(MY_MAKEFILE) run IA=intel64 RUN=lib

sointel64 soem64t:
	$(MAKE) -f $(MY_MAKEFILE) prepare_run IA=intel64 RUN=so
	$(MAKE) -f $(MY_MAKEFILE) run IA=intel64 RUN=so

prepare_run: mkdir $(res_dir)/mkl_dfti.mod

mkdir: FORCE
	@echo
	@echo See results in $(res_dir)
	@echo
	mkdir -p $(res_dir)

$(res_dir)/mkl_dfti.mod: $(res_dir)/mkl_dfti.o
$(res_dir)/%.o: %.f90
	# compile $*.f90
	$(FC) $(FOPTS) $(FFLAGS) $(TARGET_ARCH) \
	  -c $< -o $@

run: $(RES)

%.res: %.f90 FORCE
	@echo
	# build $*.out
	$(FC) $(FOPTS) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH) \
	  $< \
	  $(res_dir)/mkl_dfti.o \
	  $(MKL_LIBS) \
	  $(OMP_LIB) -lpthread -lm -ldl -o $(res_dir)/$*.out
	-rm -f $(res_dir)/$@ # remove res-file
	@echo
	# run $*.out
	$(RUNENV) \
	  $(res_dir)/$*.out > $(res_dir)/$@

report: FORCE
	@echo
	@echo "--------------------------------------------------"
	@echo There are `ls source/*.f90 | wc -l` \
	  examples in source
	@echo
	@for D in `ls -d $(RES_DIR)/* 2>/dev/null`; do \
	  echo -n "In   $$D	Total run examples:" \
	    "`ls $$D/*.res 2>/dev/null | wc -l`"; \
	  echo -n "	and Passed: "; \
	  cat /dev/null $$D/*.res 2>/dev/null | grep -c PASSED | cat; \
	done
	@echo "--------------------------------------------------"
	@echo

clean: FORCE
	-rm -rf $(RES_DIR) # clean all results

FORCE: ;

##-----------------------------------------------------------------------------
