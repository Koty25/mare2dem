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
##      Intel(R) Math Kernel Library. Service functions example creation and run
##******************************************************************************

help:
	@echo "Intel(R) MKL Service functions example help."
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
	@echo "       Link 'servicefuncs' example against parallel static MKL-IA32."
	@echo "   make sointel64 interface=ilp64 threading=sequential"
	@echo "       Link 'servicefuncs' example against sequential dynamic MKL-Intel64"
	@echo "       using ILP64 interface."
	@echo ""
#------------------------------------------------------------------------------

ifndef compiler
compiler=intel
else
ifneq ($(compiler),intel)
ifneq ($(compiler),gnu)
$(warning *** Wrong value compiler=$(compiler))
$(warning *** The correct values for compiler are intel or gnu)
$(error Try >make help)
endif
endif
endif

ifndef interface
interface=lp64
else
ifneq ($(interface),ilp64)
ifneq ($(interface),lp64)
$(warning *** Wrong value interface=$(interface))
$(warning *** The correct values for interface are lp64 or ilp64)
$(error Try >make help)
endif
endif
endif

ifndef threading
threading=parallel
else
ifneq ($(threading),parallel)
ifneq ($(threading),sequential)
$(warning *** Wrong value threading=$(threading))
$(warning *** The correct values for threading are parallel or sequential)
$(error Try >make help)
endif
endif
endif

RES = servicefuncs.res

ifdef _IA

EXT_S=a
EXT_D=so
EXT=$(EXT_$(SD))

RES_EXT_S=lib
RES_EXT_D=so

ifndef MKLROOT
  MKLROOT = ../..
endif
MKL_PATH = $(MKLROOT)/lib/$(_IA)
CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(_IA)

LOPTS = 

ifeq ($(compiler),gnu)
  override COMPILER = gfortran
  OPTIONS = -fcray-pointer -x f77-cpp-input
  OPTIONS.DIAG = -Wall -Werror -Wno-intrinsic-shadow
  FCOPTS_END = -x none
  IFACE_COMP_PART = gf
  IFACE_THREADING_PART = gnu
  ifeq ($(SD),D)
      LOPTS = -Wl,--no-as-needed
  endif
else
  override COMPILER = ifort
  OPTIONS = -fpp
  FCOPTS_END =
  IFACE_COMP_PART = intel
  IFACE_THREADING_PART = intel
endif

OPTIONS += $(OPTIONS.DIAG)

ifeq ($(_IA),ia32)
  ifeq ($(interface),ilp64)
    $(warning  *** ILP64 interface is not available for MKL-IA32)
    $(error Try >make help)
  endif
  OPTIONS += -m32 # This option tells compiler to generate code for IA-32 architecture.
  ifeq ($(compiler),intel)
      # This option tells the compiler to generate optimized code
      # for Pentium or later processor.
      # If you don't need it, you can remove this option.
    OPTIONS += -mia32
  endif
  OPTIONS += -D_IA32 # This definition is needed for the example itself on IA-32 architecture.
  IFACE_SUFF=
else
  IFACE_SUFF=_$(interface)
endif

RES_DIR=_results/$(compiler)$(IFACE_SUFF)_$(threading)_$(_IA)_$(RES_EXT_$(SD))
IFACE_LIB=mkl_$(IFACE_COMP_PART)$(IFACE_SUFF)

ifeq ($(interface),ilp64)
  ifeq ($(compiler),gnu)
    OPTIONS += -fdefault-integer-8
  else
    OPTIONS += -i8
  endif
endif

ifeq ($(threading),sequential)
  THREADING_LIB=mkl_sequential
  OMP_LIB =
else
  THREADING_LIB=mkl_$(IFACE_THREADING_PART)_thread
  OMP_LIB = -L"$(CMPLR_PATH)" -liomp5
endif

CORE_LIB=mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT_$(SD))$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT_$(SD))$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT_$(SD))$(CORE_LIB).$(EXT) -Wl,--end-group
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
endif

LABEL = $(compiler)$(IFACE_SUFF)_$(threading)_$(_IA)_$(RES_EXT_$(SD))
endif # ifdef _IA

#-------------------------------------------------------------------------------

ia32: libia32 soia32
intel64: libintel64 sointel64

libia32 lib32:
	@$(MAKE) $(RES) --no-print-directory SD=S _IA=ia32
soia32 so32:
	@$(MAKE) $(RES) --no-print-directory SD=D _IA=ia32
libintel64 libem64t:
	@$(MAKE) $(RES) --no-print-directory SD=S _IA=intel64
sointel64 soem64t:
	@$(MAKE) $(RES) --no-print-directory SD=D _IA=intel64

.PHONY: clean
clean:
	-rm -rf _results;
#-------------------------------------------------------------------------------

vpath %.f source

$(RES): %.res: %.f
	@mkdir -p ./$(RES_DIR)
	@echo ----- Compiling $(LABEL) -----
	$(COMPILER) $(OPTIONS) -I"$(MKLROOT)/include" \
	$< $(FCOPTS_END) \
	$(LOPTS) \
	$(MKL_LIBS) \
	$(OMP_LIB) \
	-lpthread -lm -ldl -o $(RES_DIR)/$*.out
	@echo ----- Execution $(LABEL) -----
	export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
	$(RES_DIR)/$*.out > $(RES_DIR)/$@

#-------------------------------------------------------------------------------
