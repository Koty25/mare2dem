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
##      Intel(R) Math Kernel Library Trigonometric Transforms examples creation
##      and run
##
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name][gnu_path=gnu_lib_path]"
	@echo
	@echo "name           - function name. Please see pdettc.lst file"
	@echo
	@echo "compiler_name  - can be gnu, pgi (only for intel64), or intel. Default value is intel."
	@echo "                 Intel(R) C Compiler as default."
	@echo "parallel_name  - can be intel, pgi (only if compiler=pgi), or gnu (only if compiler=gnu). Default value is intel."
	@echo "omp_name       - can be iomp5 if parallel=intel, or"
	@echo "                 can be iomp5 or gomp if parallel=gnu, or"
	@echo "                 can be pgmp if parallel=pgi."
	@echo "                 Default value is iomp5."
	@echo
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo
	@echo "threading_name - can be parallel or sequential. Default value is parallel."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."

##------------------------------------------------------------------------------
## examples of usage:
##
## make libia32 function=d_trig_tforms_bvp - build by Intel(R) C Compiler (as default) and
##                                           run d_trig_tforms_bvp example for 32-bit applications,
##                                           static linking
##
## make libintel64 compiler=gnu            - build by GNU C compiler and run all
##                                           examples of Intel(R) MKL PDE TT for Intel(R) 64 processor family
##                                           applications, static linking
##
##------------------------------------------------------------------------------

include pdettc.lst

ifndef function
function = $(PDE_TT_C)
endif

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

ifndef interface
interface=lp64
endif

ifndef threading
threading=parallel
endif

ifndef $(gnu_path)
   gnu_path=/usr/lib
endif

ifeq (,$(filter gnu pgi,$(parallel)))
   override parallel=intel
   override omp=iomp5
else
   ifeq ($(parallel),gnu)
      ifneq ($(omp),gomp)
      override omp=iomp5
      endif
   else
      override omp=pgmp
   endif
endif

RES = $(addsuffix .res ,$(function))

ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

LOPTS = 
COPTS = 

WARN_FLAGS = -Wall -Werror

ifeq ($(compiler),gnu)
   override CC = gcc
   ifeq ($(RES_EXT),so)
       LOPTS = -Wl,--no-as-needed
   endif
   COPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
   IFACE_COMP_PART=intel
else
   ifeq ($(compiler),pgi)
      override CC=pgcc
      override WARN_FLAGS = 
   else
      override CC=icc
   endif
   IFACE_COMP_PART=intel
endif

ifeq ($(parallel),gnu)
   IFACE_THREADING_PART=gnu
   GFORTRAN_LIB=-L$(gnu_path) -lgfortran
else
      GFORTRAN_LIB=
   ifeq ($(parallel),pgi)
      IFACE_THREADING_PART=pgi
   else
      IFACE_THREADING_PART=intel
   endif
endif

ifeq ($(omp),gomp)
      OMP_LIB = -L$(gnu_path) -l$(omp)
else
   ifeq ($(omp),pgmp)
      OMP_LIB =
   else
      OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
   endif
endif

ifeq ($(interface),ilp64)
IFACE_LIB = mkl_$(IFACE_COMP_PART)_ilp64
COPTS += -DMKL_ILP64
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

ifeq ($(compiler)-$(threading),pgi-parallel)
COPTS += -mp -pgf90libs -Mnokeepobj
endif

ifeq ($(threading),sequential)
THREADING_LIB = mkl_sequential
OMP_LIB =
else
THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group $(OMP_LIB)
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB) $(OMP_LIB)
endif

ifeq ($(_IA),ia32)
RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
else
RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
endif

libia32 lib32:
	$(MAKE) $(RES) _IA=ia32 EXT=a RES_EXT=lib
soia32 so32:
	$(MAKE) $(RES) _IA=ia32 EXT=so RES_EXT=so
libintel64 libem64t:
	$(MAKE) $(RES) _IA=intel64 EXT=a RES_EXT=lib
sointel64 soem64t:
	$(MAKE) $(RES) _IA=intel64 EXT=so RES_EXT=so

.PHONY: clean
clean:
	-rm -rf _results;

#-------------------------------------------------------------------------------

vpath %.c source

$(RES): %.res: %.c
	mkdir -p ./$(RES_DIR)
	$(CC) $(SPEC_OPT) $(WARN_FLAGS) $(COPTS) -I$(MKLROOT)/include $< -L$(MKL_PATH) $(LOPTS) $(MKL_LIBS) -lpthread -lm -ldl -o $(RES_DIR)/$*.out
	export LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(CMPLR_PATH); $(RES_DIR)/$*.out >$(RES_DIR)/$@

#-------------------------------------------------------------------------------
