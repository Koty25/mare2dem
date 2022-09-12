#===============================================================================
# Copyright 2004-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library Extended Eigensolvers create and run examples
##******************************************************************************

include solvers_eec.lst

help:
	@echo $$'$(WRONG_OPTION)'
	@echo "Intel(R) MKL Extended Eigensolvers example "
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
	@echo "           Programming interface for intel64."
	@echo "           Default: lp64."
	@echo "       threading={parallel|sequential}"
	@echo "           Specifies whether to use Intel(R) MKL in threaded or sequential mode."
	@echo "           Default: parallel."
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
	@echo "           The supported function list:"
	@echo "               $(SOLVERS_EE_SPARSE)"
	@echo "               $(SOLVERS_EE_BANDED)"
	@echo "               $(SOLVERS_EE_DENSE)"
	@echo "           Default: all examples."
	@echo ""
	@echo "Usage examples:"
	@echo ""
	@echo "   make libia32 function=\"dexample_sparse_c\""
	@echo "       Link 'dexample_sparse_c' examples against parallel static MKL-IA32."
	@echo ""
	@echo "   make sointel64 interface=ilp64 threading=sequential"
	@echo "       Link all the examples against sequential dynamic MKL-Intel64"
	@echo "       using ILP64 interface."
	@echo ""
#------------------------------------------------------------------------------


ifndef function
function = $(SOLVERS_EE)
endif

TARGET = $(function)

# Defaults
compiler=intel
interface=lp64
threading=parallel

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
MSG2+= threading=$(threading)
endif
endif

ifneq ("$(MSG2)","")
WRONG_OPTION=\n\n*** COMMAND LINE ERROR: Wrong value of option(s):  $(MSG2)\n\n
TARGET=help
endif
###

ifdef _IA

ifeq ($(SD),static)
  EXT=a
  RES_EXT=lib
  S_GRP=-Wl,--start-group
  E_GRP=-Wl,--end-group
else
  EXT=so
  RES_EXT=so
  S_GRP=
  E_GRP=
endif

ifndef MKLROOT
  MKLROOT = ../..
endif
MKL_PATH = $(MKLROOT)/lib/$(_IA)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(_IA)

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
    $(warning  *** ILP64 interface is not available for MKL-IA32)
    $(error Try >make help)
  endif
  ifeq ($(compiler),intel)
      # This option tells the compiler to generate optimized code
      # for Pentium or later processor.
      # If you don't need it, you can remove this option.
    OPTIONS += -mia32
  endif
  IFACE_SUFF=
  M32_64 = -m32 # This option tells compiler to generate code for IA-32 architecture.
else
  IFACE_SUFF=_$(interface)
  M32_64 = -m64 # This option tells compiler to generate code for Intel64 architecture.
endif


ifeq ($(EXT),so)
    MKL_LD_PATH=-L"$(MKL_PATH)"
    MKL_PREFIX=-l
    MKL_SUFFIX=
else
    MKL_LD_PATH=
    MKL_PREFIX="$(MKL_PATH)/lib
    MKL_SUFFIX=.a"
endif
IFACE_LIB=$(MKL_PREFIX)mkl_$(IFACE_COMP_PART)$(IFACE_SUFF)$(MKL_SUFFIX)

ifeq ($(interface),ilp64)
  OPTIONS += -DMKL_ILP64
endif

ifeq ($(threading),sequential)
  THREADING_LIB=$(MKL_PREFIX)mkl_sequential$(MKL_SUFFIX)
  OMP_LIB =
else
  THREADING_LIB=$(MKL_PREFIX)mkl_$(IFACE_THREADING_PART)_thread$(MKL_SUFFIX)
  OMP_LIB = -L"$(CMPLR_PATH)" -liomp5
endif

CORE_LIB=$(MKL_PREFIX)mkl_core$(MKL_SUFFIX)

LABEL = $(compiler)$(IFACE_SUFF)_$(threading)_$(_IA)_$(RES_EXT)
RES_DIR = _results/$(LABEL)

ifneq ("$(out_to_screen)","")
OUT_TO_FILE=
else
OUT_TO_FILE=> $(RES_DIR)/$@.res
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

$(SOLVERS_EE):
	@mkdir -p ./$(RES_DIR)
	@echo $$'\n----- Compiling $(LABEL) ----- $@'
	$(COMPILER) $(M32_64) $(OPTIONS) -I"$(MKLROOT)/include" \
	./source/$@.c $(S_GRP) \
	$(MKL_LD_PATH) $(LOPTS) $(IFACE_LIB) \
	$(THREADING_LIB) \
	$(CORE_LIB) \
	$(E_GRP) $(OMP_LIB) -lpthread -lm -ldl -o $(RES_DIR)/$@.out
	@echo ----- Execution $(LABEL) ----- $@
	export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
	$(RES_DIR)/$@.out $(OUT_TO_FILE)

#-------------------------------------------------------------------------------
