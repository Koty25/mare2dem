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
##      Intel(R) Math Kernel Library. Single Dynamic Library Interface example.
##******************************************************************************

help:
	@echo $$'$(WRONG_OPTION)'
	@echo "Intel(R) MKL Single Dynamic Library Interface example help."
	@echo "    The example demonstrates the dynamical setting of Intel(R) MKL threading layer:"
	@echo "        1) by the specifying of environment variable MKL_THREADING_LAYER and"
	@echo "        2) by the functional call to MKL_Set_Threading_Layer()."
	@echo ""
	@echo "Usage: make Target [Options...]"
	@echo ""
	@echo "   Target:"
	@echo "       ia32    - dynamically linked example"
	@echo "                     for processors that use the IA-32 architecture"
	@echo "       intel64 - dynamically linked example"
	@echo "                     for processors that use the Intel(R) 64 architecture"
	@echo "       help    - print this help"
	@echo ""
	@echo "   Options:"
	@echo "       compiler={intel|gnu}"
	@echo "           Specifies the compiler used to compile the example."
	@echo "           Default: intel."
	@echo "       MKLROOT=<MKL_directory>"
	@echo "           Specifies the location of Intel(R) MKL libraries used to build this example."
	@echo "           Default: the Intel(R) MKL installation directory."
	@echo ""
	@echo "Usage examples:"
	@echo ""
	@echo "   make intel64"
	@echo "       Link 'Single Dynamic Library Interface' example"
	@echo "           for processors that use the Intel(R) 64 architecture"
	@echo "           using Intel compiler"
	@echo ""
	@echo "   make ia32 compiler=gnu"
	@echo "       Link 'Single Dynamic Library Interface' example"
	@echo "           for processors that use the IA-32 architecture"
	@echo "           using GNU compiler"
	@echo ""
#------------------------------------------------------------------------------

FUNCTION_LIST = sdli
TARGET = $(FUNCTION_LIST)

# Defaults
compiler=intel

ifneq ($(compiler),intel)
ifneq ($(compiler),gnu)
MSG2+= compiler=$(compiler)
endif
endif

ifneq ("$(MSG2)","")
WRONG_OPTION=\n\n*** COMMAND LINE ERROR: Wrong value of option(s):  $(MSG2)\n\n
TARGET=help
endif
###

ifdef _IA

ifndef MKLROOT
  MKLROOT = ../..
endif
MKL_PATH = $(MKLROOT)/lib/$(_IA)

LOPTS = 

ifeq ($(compiler),gnu)
  override COMPILER = gcc
  LOPTS = -Wl,--no-as-needed
else
  override COMPILER = icc
endif
OPTIONS =

ifeq ($(_IA),ia32)
  ifeq ($(compiler),intel)
      # This option tells the compiler to generate optimized code
      # for Pentium or later processor.
      # If you don't need it, you can remove this option.
    OPTIONS += -mia32
  endif
  M32_64 = -m32 # This option tells compiler to generate code for IA-32 architecture.
else
  M32_64 = -m64 # This option tells compiler to generate code for Intel64 architecture.
endif

# required libraries
IFACE_LIB = $(LOPTS) -L$(MKL_PATH) -lmkl_rt

CMPLR_PATH = $(MKLROOT)/../compiler/lib/$(_IA)
ifeq ($(threading),tbb)
TBB_PATH = $(shell if [ -d "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8" ]; then \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8"; \
	else \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.4"; \
	fi)
CMPLR_PATH += $(TBB_PATH) 
OMP_LIB = -L$(TBB_PATH) -lstdc++
else
OMP_LIB = -L$(CMPLR_PATH) -liomp5
endif

LABEL = $(compiler)_$(_IA)
RES_DIR = _results/$(LABEL)

OUT_TO_FILE=>> $(RES_DIR)/$@.res
PRINT_RESULT=@cat $(RES_DIR)/$@.res

endif # ifdef _IA

#-------------------------------------------------------------------------------

ia32:
	@$(MAKE) $(TARGET) --no-print-directory _IA=ia32
intel64:
	@$(MAKE) $(TARGET) --no-print-directory _IA=intel64


.PHONY: clean
clean:
	-rm -rf _results;
#-------------------------------------------------------------------------------

$(FUNCTION_LIST):
	@echo "----- Intel(R) MKL Single Dynamic Library Interface example -----"
	@mkdir -p ./$(RES_DIR)
	@echo $$'\n----- Compiling $(LABEL) ----- $@'
	@echo $$'      The single Intel(R) MKL library is linked: $(IFACE_LIB)\n'
	$(COMPILER) $(M32_64) $(OPTIONS) -I"$(MKLROOT)/include" \
	./source/$@.c \
	$(IFACE_LIB) \
	$(OMP_LIB) -lpthread -lm -ldl \
	-o $(RES_DIR)/$@.out

	@echo "Intel(R) MKL Single Dynamic Library Interface example" > $(RES_DIR)/$@.res
	@echo
	@echo "----- Run the example $@"
	@echo "Case 1: Run $@ with default behavior: no parameters, no Intel(R) MKL environment variables"
	@echo "Case 1: Run $@ with default behavior:" $(OUT_TO_FILE)
	@echo "        no parameters, no Intel(R) MKL environment variables" $(OUT_TO_FILE)
	@echo "        Intel(R) MKL Intel parallel threading layer will be used" $(OUT_TO_FILE)
	@export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
		$(RES_DIR)/$@.out $(OUT_TO_FILE)

	@echo "Case 2: Set environment variable MKL_THREADING_LAYER"
	@echo "Case 2: Set environment variable MKL_THREADING_LAYER" $(OUT_TO_FILE)
	@echo "        set MKL_THREADING_LAYER=SEQUENTIAL" $(OUT_TO_FILE)
	@echo "        Intel(R) MKL sequential threading layer will be used" $(OUT_TO_FILE)
	@export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
		export MKL_THREADING_LAYER=SEQUENTIAL; \
		$(RES_DIR)/$@.out $(OUT_TO_FILE)

	@echo "Case 3: Environment variable vs the function call"
	@echo "Case 3: Environment variable vs the function call" $(OUT_TO_FILE)
	@echo "        set MKL_THREADING_LAYER=SEQUENTIAL" $(OUT_TO_FILE)
	@echo "        Run $@ with the parameter 'par'"  $(OUT_TO_FILE)
	@echo "        The example will set PARALLEL threading Layer" $(OUT_TO_FILE)
	@echo "            by the call MKL_Set_Threading_Layer(MKL_THREADING_INTEL)" $(OUT_TO_FILE)
	@echo "        The function call has precedence over environment variable" $(OUT_TO_FILE)
	@echo "        Intel threading layer will be used" $(OUT_TO_FILE)
	@export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
		export MKL_THREADING_LAYER=SEQUENTIAL; \
		$(RES_DIR)/$@.out 'par' $(OUT_TO_FILE)

	@echo "Case 4: Select sequential Intel(R) MKL by the call MKL_Set_Threading_Layer(MKL_THREADING_SEQUENTIAL)"
	@echo "Case 4: Select sequential Intel(R) MKL by the call MKL_Set_Threading_Layer(MKL_THREADING_SEQUENTIAL)" $(OUT_TO_FILE)
	@echo "        Run $@ with the parameter 'seq'" $(OUT_TO_FILE)
	@echo "        The example will set SEQUENTIAL threading Layer" $(OUT_TO_FILE)
	@echo "            by the call MKL_Set_Threading_Layer(MKL_THREADING_SEQUENTIAL)" $(OUT_TO_FILE)
	@export LD_LIBRARY_PATH="$(MKL_PATH)":$(LD_LIBRARY_PATH):$(CMPLR_PATH); \
		$(RES_DIR)/$@.out 'seq' $(OUT_TO_FILE)

	@echo ----- The example '$@' completed
	@echo
	@echo ----- Results [ file $(RES_DIR)/$@.res ]:
	$(PRINT_RESULT)
	@echo ----- EOF -----
	@rm -f $(RES_DIR)/$@.o

#-------------------------------------------------------------------------------
