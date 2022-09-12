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
##      Intel(R) Math Kernel Library LAPACK examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
	@echo "name     - function name. Please see lapack.lst file."
	@echo "compiler_name     - can be gnu or intel. Default value is intel."
	@echo "                    Intel (R) Fortran Compiler as default."
	@echo "                    If compiler=gnu then GNU gfortran compiler will be used."
	@echo "                    To use g77 add FC=g77. FC=g77 cannot be used with ilp64 interface."
	@echo "interface_name    - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "threading_name    - can be omp, tbb or sequential. Default value is omp."
	@echo "parallel_name  - can be intel or gnu. Default value is intel."
	@echo "                 Relevant if threading=omp."
	@echo "omp_name       - can be iomp5 if parallel=intel or"
	@echo "                 can be iomp5 or gomp if parallel=gnu."
	@echo "                 Default value is iomp5."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=dgebrd
##              - build  by  Intel(R) Fortran Compiler  (as default)
##                and run DGEBRD example for 32-bit applications,
##                static linking
##
## make soia32 compiler=gnu FC=g77
##              - build by g77 GNU Fortran  compiler  and  run
##                all examples of Intel(R) MKL for 32-bit applications,
##                dynamic linking
##
## make libintel64 compiler=gnu
##              - build by gfortran GNU  Fortran  compiler  and  run  all
##                examples of Intel(R) MKL for Intel(R) 64 processor family applications,
##                static linking
##
## make sointel64
##              - build by Intel(R) Fortran Compiler (as default) and run all
##                examples of Intel(R) MKL for Intel(R) 64 processor family applications,
##                dynamic linking
##------------------------------------------------------------------------------

include lapack.lst

ifndef function
   function = $(LAPACK)
endif

ifneq ($(compiler),gnu)
   override compiler=intel
endif

ifneq ($(interface),ilp64)
   override interface=lp64
endif

ifneq ($(threading),sequential)
   ifneq ($(threading),tbb)
      override threading=omp
   endif
endif

ifneq ($(parallel),gnu)
   override parallel=intel
endif

ifeq ($(parallel),gnu)
   ifneq ($(omp),gomp)
      override omp=iomp5
   endif
else
   override omp=iomp5
endif

ifndef gnu_path
   gnu_path = /usr/lib
endif

ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = $(MKLROOT)/lib/$(_IA)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(_IA)
TBB_PATH = $(shell if [ -d "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8" ]; then \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8"; \
	else \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.4"; \
	fi)

FOPTS =
LOPTS =

ifeq ($(compiler),gnu)
   override FC = gfortran
   IFACE_COMP_PART=gf
   FOPTS += $(if $(filter ia32, $(_IA)), -m32, -m64) -Wall -Werror -finit-local-zero
   ifeq ($(RES_EXT),so)
       LOPTS = -Wl,--no-as-needed
   endif
else
   override FC=ifort
   IFACE_COMP_PART=intel
   FOPTS += -warn all -warn errors -nogen-interfaces
endif

ifeq ($(interface),ilp64)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_ilp64
   FOPTS += $(if $(filter gnu, $(compiler)), -fdefault-integer-8, -i8)
else
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_lp64
endif

GFORTRAN_LIB =

ifeq ($(threading),sequential)
   THREAD_LIB = mkl_sequential
   LL_THREAD_LIB =
else
   ifeq ($(threading),omp)
      THREAD_LIB = mkl_$(parallel)_thread
      ifeq ($(parallel),intel)
         LL_THREAD_LIB = -L$(CMPLR_PATH)
      else
         GFORTRAN_LIB = -L$(gnu_path) -lgfortran
         LL_THREAD_LIB = -L$(if $(filter gomp, $(omp)), $(gnu_path), $(CMPLR_PATH))
      endif
      LL_THREAD_LIB += -l$(omp)
   else
      THREAD_LIB = mkl_tbb_thread
      LL_THREAD_LIB.a  = -L$(TBB_PATH) -ltbb -lstdc++
      LL_THREAD_LIB.so =
      LL_THREAD_LIB   := $(LL_THREAD_LIB.$(EXT))
   endif
endif

ifeq ($(_IA),ia32)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)
endif

RES = $(addsuffix x.res ,$(function))

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREAD_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREAD_LIB) -l$(CORE_LIB)
endif

ifeq ($(_IA),ia32)
   ifeq ($(threading),omp)
      RES_DIR=_results/$(compiler)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
      RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
else
   ifeq ($(threading),omp)
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
endif

libia32 lib32:
	$(MAKE) $(RES) EXT=a _IA=ia32 PLATF=ia32 RES_EXT=lib
soia32 so32:
	$(MAKE) $(RES) EXT=so _IA=ia32 PLATF=ia32 RES_EXT=so
libintel64 libem64t:
	$(MAKE) $(RES) EXT=a _IA=intel64 PLATF=intel64 RES_EXT=lib
sointel64 soem64t:
	$(MAKE) $(RES) EXT=so _IA=intel64 PLATF=intel64 RES_EXT=so

.PHONY: clean
clean:
	-rm -rf _results;

#-------------------------------------------------------------------------------

vpath %.f source

$(RES): %.res: %.f 
	mkdir -p ./$(RES_DIR)
	$(FC) $(FOPTS) $< $(LOPTS) $(MKL_LIBS) $(LL_THREAD_LIB) $(GFORTRAN_LIB) -lpthread -ldl -o $(RES_DIR)/$*.out
	@if [ -f data/$*.d ]; then \
		export LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(gnu_path):$(CMPLR_PATH):$(TBB_PATH):$(LD_LIBRARY_PATH); $(RES_DIR)/$*.out <data/$*.d >$(RES_DIR)/$@ ; \
	else \
		export LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(gnu_path):$(CMPLR_PATH):$(TBB_PATH):$(LD_LIBRARY_PATH); $(RES_DIR)/$*.out >$(RES_DIR)/$@ ; \
	fi
	
#-------------------------------------------------------------------------------
