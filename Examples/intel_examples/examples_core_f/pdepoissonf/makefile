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
##      Intel(R) Math Kernel Library Poisson Solver examples creation
##      and run
##
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name][gnu_path=gnu_lib_path]"
	@echo
	@echo "name           - function name. Please see pdepoissonf.lst file"
	@echo
	@echo "compiler_name  - can be gnu, pgi, or intel. Default value is intel."
	@echo "                 Intel (R) Fortran Compiler as default."
	@echo "                 If compiler=gnu, then GNU gfortran compiler will be used."
	@echo "                 Note: GNU g77 compiler is not supported."
	@echo "                 If compiler=pgi, then PGI Fortran 95 compiler will be used. Available only for intel64."
	@echo
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo
	@echo "threading_name - can be omp, tbb or sequential. Default value is omp."
	@echo "parallel_name  - can be intel, pgi (only if compiler=pgi), or gnu (only if compiler=gnu)."
	@echo "                 Default value is intel."
	@echo "omp_name       - can be iomp5 if parallel=intel, or"
	@echo "                 can be iomp5 or gomp if parallel=gnu, or"
	@echo "                 can be pgmp if parallel=pgi"
	@echo "                 Default value is iomp5."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."

##------------------------------------------------------------------------------
## examples of usage:
##
## make libia32 function=d_Helmholtz_2D_f  - build by Intel(R) Fortran Compiler (as default) and
##                                           run d_Helmholtz_2D_f example for 32-bit applications,
##                                           static linking
##
## make libintel64                         - build by Intel(R) Fortran Compiler (as default) and run all
##                                           examples of Intel(R) MKL PDE POISSON for Intel(R) 64 processor family
##                                           applications, static linking
##
##------------------------------------------------------------------------------

include pdepoissonf.lst

function = $(PDE_POISSON_F)

ifeq (,$(filter gnu pgi,$(compiler)))
   override compiler=intel
   override parallel=intel
endif

ifeq ($(compiler)_$(_IA),pgi_ia32)
   $(info ***)
   $(info *** PGI compiler 32-bit support was dropped off starting from Intel(R) MKL 2017.)
   $(info *** Please try >make help)
   $(info ***)
   $(error )
endif

omp=iomp5

interface=lp64
threading=omp
gnu_path=/usr/lib

RES = $(addsuffix .res ,$(function))
MOD = $(addsuffix .mod ,$(modules))

ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"
TBB_PATH = $(shell if [ -d "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8" ]; then \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8"; \
	else \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.4"; \
	fi)

LOPTS = 
FOPTS = 

ifeq ($(compiler),gnu)
   FC = gfortran
   IFACE_COMP_PART=gf
   FOPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
   ifeq ($(RES_EXT),so)
       LOPTS = -Wl,--no-as-needed
   endif
else
   ifeq ($(compiler),pgi)
      FC=pgf95
   else
      override FC=ifort
   endif
   IFACE_COMP_PART=intel
endif

ifeq ($(parallel),pgi)
   override omp=pgmp
endif

ifeq ($(parallel),pgi)
   IFACE_THREADING_PART=pgi
else
ifeq ($(parallel),gnu)
   IFACE_THREADING_PART=gnu
else
   IFACE_THREADING_PART=intel
endif
endif


ifeq ($(interface),ilp64)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_ilp64
   FOPTS += $(if $(filter gnu, $(compiler)), -fdefault-integer-8, -i8)
else
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_lp64
endif

ifeq ($(_IA),ia32)
   ifeq ($(compiler),intel)
       SPEC_OPT=-mia32
#This option is required by Intel(R) 11.0 compiler to produce workable binaries for Pentium(R) III.
#If you don't need it, you can remove this option.
   endif
IFACE_LIB = mkl_$(IFACE_COMP_PART)
endif

ifeq ($(parallel),pgi)
   FOPTS += -mp
endif

ifeq ($(omp),gomp)
      LL_THREAD_LIB = -L$(gnu_path) -l$(omp)
else
   ifeq ($(omp),pgmp)
      LL_THREAD_LIB =
   else
      LL_THREAD_LIB = -L$(CMPLR_PATH) -l$(omp)
   endif
endif

ifeq ($(threading),sequential)
   THREADING_LIB = mkl_sequential
   LL_THREAD_LIB =
   GFORTRAN_LIB=
else
   ifeq ($(threading),tbb)
      THREADING_LIB = mkl_tbb_thread
      LL_THREAD_LIB.a  = -L$(TBB_PATH) -ltbb -lstdc++
      LL_THREAD_LIB.so =
      LL_THREAD_LIB   := $(LL_THREAD_LIB.$(EXT))
   else
      THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
      LL_THREAD_LIB = -L$(CMPLR_PATH) -l$(omp)
   endif
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group $(LL_THREAD_LIB)
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB) $(LL_THREAD_LIB)
endif

ifeq ($(_IA),ia32)
RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
else
RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
endif

libia32 lib32: clean_files
	$(MAKE) $(MOD) $(RES) _IA=ia32 EXT=a RES_EXT=lib
soia32 so32: clean_files
	$(MAKE) $(MOD) $(RES) _IA=ia32 EXT=so RES_EXT=so
libintel64 libem64t: clean_files
	$(MAKE) $(MOD) $(RES) _IA=intel64 EXT=a RES_EXT=lib
sointel64 soem64t: clean_files
	$(MAKE) $(MOD) $(RES) _IA=intel64 EXT=so RES_EXT=so
clean_files:
	rm -rf *.mod *.o;

.PHONY: clean
clean:
	-rm -rf _results;
#-------------------------------------------------------------------------------

vpath %.f90 source
vpath %.f90 $(MKLROOT)/include

.NOTPARALLEL: # as written, this can only be run with make -j1, to support creation of .mod files first

$(MOD): %.mod: %.f90
	mkdir -p ./$(RES_DIR)
	$(FC) $(SPEC_OPT) $(FOPTS) -w -c $<

$(RES): %.res: %.f90
	mkdir -p ./$(RES_DIR)
	$(FC) $(SPEC_OPT) $(FOPTS) -w $< $(LOPTS) $(MKL_LIBS) -lpthread -ldl -o $(RES_DIR)/$*.out
	export LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(CMPLR_PATH):$(TBB_PATH); $(RES_DIR)/$*.out >$(RES_DIR)/$@
#-------------------------------------------------------------------------------
