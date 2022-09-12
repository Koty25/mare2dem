#===============================================================================
# Copyright 2010-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library LAPACKE examples creation and run
##******************************************************************************
##
## Usage examples:
##
## make libia32 function="lapacke_dgels_row lapacke_zgesv_col"
##      Link 'lapacke_dgels_row' and 'lapacke_zgesv_col' examples against
##      parallel (omp) static Intel(R) MKL for applications running on IA-32.
##
## make sointel64 interface=ilp64 threading=sequential
##      Link all the examples against sequential dynamic Intel(R) MKL for applications
##      running on Intel(R) 64 using ILP64 interface.
##******************************************************************************

help:
	@echo $$'$(WRONG_OPTION)'
	@echo "Intel(R) MKL LAPACKE examples help."
	@echo ""
	@echo "Usage: make Target [Options...] [Function...]"
	@echo ""
	@echo "   Target:"
	@echo "       libia32    - statically linked example for processors that use the IA-32 architecture"
	@echo "       soia32     - dynamically linked example for processors that use the IA-32 architecture"
	@echo "       libintel64 - statically linked example for processors that use the Intel(R) 64 architecture"
	@echo "       sointel64  - dynamically linked example for processors that use the Intel(R) 64 architecture"
	@echo "       help       - print this help"
	@echo ""
	@echo "   Options:"
	@echo "       interface={lp64|ilp64}"
	@echo "           programming interface for Intel(R) 64."
	@echo "           Default: lp64."
	@echo "       threading={omp|tbb|sequential}"
	@echo "           Specifies whether to use Intel(R) MKL in the threaded or sequential mode."
	@echo "           Default: omp."
	@echo "       compiler={intel|gnu}"
	@echo "           Specifies the compiler used to compile the example."
	@echo "           Default: intel."
	@echo "       MKLROOT=<MKL_directory>"
	@echo "           Specifies the location of Intel(R) MKL libraries used to build this example."
	@echo "           Default: the Intel(R) MKL installation directory."
	@echo ""
	@echo "   Function:"
	@echo "       function=function_name1 function_name2 ..."
	@echo "           Specifies the list of required functions."
	@echo "           Please find the supported function list in the lapacke.lst file."
	@echo "           Default: all examples."
	@echo ""

include lapacke.lst

ifndef function
function = $(LAPACKE)
endif

TARGET = $(function)

# Defaults
compiler = intel
interface = lp64
threading = omp

ifneq ($(compiler),intel)
ifneq ($(compiler),gnu)
MSG2+= compiler=$(compiler)
endif
endif

ifneq ($(interface),ilp64)
ifneq ($(interface),lp64)
MSG2+= interface=$(interface)
endif
endif

ifneq ($(threading),parallel)
ifneq ($(threading),sequential)
ifneq ($(threading),omp)
ifneq ($(threading),tbb)
MSG2+= threading=$(threading)
endif
endif
endif
endif

ifeq ($(threading),parallel)
    override threading=omp
endif


ifneq ("$(MSG2)","")
WRONG_OPTION=\n\n*** COMMAND LINE ERROR: Wrong value of option(s):  $(MSG2)\n\n
TARGET=help
endif
###

ifdef _IA

ifeq ($(SD),static)
  EXT = a
  RES_EXT = lib
else
  EXT = so
  RES_EXT = so
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

LOPTS =

ifeq ($(compiler),gnu)
  override COMPILER = gcc
  IFACE_COMP_PART = intel
  IFACE_THREADING_PART = gnu
  ifeq ($(RES_EXT),so)
      LOPTS = -Wl,--no-as-needed
  endif
else
  override COMPILER = icc
  IFACE_COMP_PART = intel
  IFACE_THREADING_PART = intel
endif
OPTIONS = -w

ifeq ($(_IA),ia32)
  ifeq ($(interface),ilp64)
    $(warning  *** ILP64 interface is not available for IA-32 architecture)
    $(error Try >make help)
  endif
  ifeq ($(compiler),intel)
      # This option tells the compiler to generate optimized code
      # for Intel(R) Pentium or later processor.
      # If you don't need it, you can remove this option.
    OPTIONS += -mia32
  endif
  IFACE_SUFF =
  M32_64 = -m32 # This option tells compiler to generate code for IA-32 architecture.
else
  IFACE_SUFF = _$(interface)
  M32_64 = -m64 # This option tells compiler to generate code for Intel(R) 64 architecture.
endif

IFACE_LIB = mkl_$(IFACE_COMP_PART)$(IFACE_SUFF)

ifeq ($(interface),ilp64)
  OPTIONS += -DMKL_ILP64
endif

ifeq ($(threading),sequential)
   THREAD_LIB = mkl_sequential
   LL_THREAD_LIB =
else
   ifeq ($(threading),omp)
      THREAD_LIB = mkl_intel_thread
      LL_THREAD_LIB = -L$(CMPLR_PATH) -liomp5
   else
      THREAD_LIB = mkl_tbb_thread
      LL_THREAD_LIB.a  = -L$(TBB_PATH) -ltbb -lstdc++
      LL_THREAD_LIB.so =
      LL_THREAD_LIB   := $(LL_THREAD_LIB.$(EXT))
   endif
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREAD_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREAD_LIB) -l$(CORE_LIB)
endif

LABEL = $(compiler)$(IFACE_SUFF)_$(threading)_$(_IA)_$(RES_EXT)
RES_DIR = _results/$(LABEL)

ifneq ("$(out_to_screen)","")
OUT_TO_FILE =
else
OUT_TO_FILE = > $(RES_DIR)/$@.res
endif

endif # ifdef _IA

#-------------------------------------------------------------------------------

ia32: libia32 soia32
intel64: libintel64 sointel64

libia32 lib32:
	@$(MAKE) $(TARGET) --no-print-directory SD=static  _IA=ia32
soia32 so32:
	@$(MAKE) $(TARGET) --no-print-directory SD=dynamic _IA=ia32
libintel64 libem64t:
	@$(MAKE) $(TARGET) --no-print-directory SD=static  _IA=intel64
sointel64 soem64t:
	@$(MAKE) $(TARGET) --no-print-directory SD=dynamic _IA=intel64

.PHONY: clean
clean:
	-rm -rf _results;
#-------------------------------------------------------------------------------

$(LAPACKE):
	@mkdir -p ./$(RES_DIR)
	@echo $$'\n----- Compiling $(LABEL) ----- $@'
	$(COMPILER) $(M32_64) $(OPTIONS) -I"$(MKLROOT)/include" \
	./source/$@.c \
	$(LOPTS) \
	$(MKL_LIBS) \
	$(LL_THREAD_LIB) -lpthread -lm -ldl -o $(RES_DIR)/$@.out
	@echo ----- Execution $(LABEL) ----- $@
	export LD_LIBRARY_PATH="$(MKL_PATH)":$(TBB_PATH):$(CMPLR_PATH):$(LD_LIBRARY_PATH); \
	$(RES_DIR)/$@.out $(OUT_TO_FILE)

#-------------------------------------------------------------------------------
