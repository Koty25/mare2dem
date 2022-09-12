#===============================================================================
# Copyright 2008-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library. versionquery example creation and run
##******************************************************************************

help:
	@echo $$'$(WRONG_OPTION)'
	@echo "Intel(R) MKL versionquery example help."
	@echo ""
	@echo "Usage: make Target [Options...]"
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
	@echo "           programming interface for intel64."
	@echo "           Default: lp64."
	@echo "       threading={parallel|sequential}"
	@echo "           Specifies whether to use Intel(R) MKL in the threaded or sequential mode."
	@echo "           Default: parallel."
	@echo "       compiler={intel|gnu}"
	@echo "           Specifies the compiler used to compile the example."
	@echo "           Default: intel."
	@echo "       MKLROOT=<MKL_directory>"
	@echo "           Specifies the location of Intel(R) MKL libraries used to build this example."
	@echo "           Default: the Intel(R) MKL installation directory."
	@echo ""
	@echo "Usage examples:"
	@echo ""
	@echo "   make libia32"
	@echo "       Link 'versionquery' example against parallel static MKL-IA32."
	@echo "   make sointel64 interface=ilp64 threading=sequential"
	@echo "       Link 'versionquery' example against sequential dynamic MKL-Intel64"
	@echo "       using ILP64 interface."
	@echo ""
#------------------------------------------------------------------------------

FUNCTION_LIST = getversion
TARGET = $(FUNCTION_LIST)

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
else
  EXT=so
  RES_EXT=so
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
OPTIONS =

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

IFACE_LIB = mkl_$(IFACE_COMP_PART)$(IFACE_SUFF)

ifeq ($(interface),ilp64)
  OPTIONS += -DMKL_ILP64
endif

ifeq ($(threading),sequential)
  THREADING_LIB = mkl_sequential
  OMP_LIB =
else
  THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
  OMP_LIB = -L"$(CMPLR_PATH)" -liomp5
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group  
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
endif

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

$(FUNCTION_LIST):
	@mkdir -p ./$(RES_DIR)
	@echo $$'\n----- Compiling $(LABEL) ----- $@'
	$(COMPILER) -Wall -Werror $(M32_64) $(OPTIONS) -I"$(MKLROOT)/include" \
	./source/$@.c \
	$(LOPTS) \
	$(MKL_LIBS) \
	$(OMP_LIB) -lpthread -lm -ldl -o $(RES_DIR)/$@.out
	@echo ----- Execution $(LABEL) ----- $@
	export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
	$(RES_DIR)/$@.out $(OUT_TO_FILE)

#-------------------------------------------------------------------------------
