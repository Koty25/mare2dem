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
##      Intel(R) Math Kernel Library TRANS examples creation and run
##******************************************************************************
help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo
	@echo "name     - function name. Please see trans.lst file"
	@echo
	@echo "compiler_name     - can be gnu or intel. Default value is intel."
	@echo "           Intel(R) C Compiler as default"
	@echo
	@echo "interface_name    - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo
	@echo "threading_name    - can be parallel or sequential. Default value is parallel."
##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=simatcopy  - build  by  Intel(R) C Compiler  (as default)
##                                    and  run  simatcopy  example  for 32-bit
##                                    applications, static linking
##
## make soia32 compiler=gnu           - build  by  GNU  C  compiler  and  run  all
##                                    examples of Intel(R) MKL  for 32-bit  applications,
##                                    dynamic linking
##
## make libintel64 compiler=gnu       - build  by  GNU  C  compiler  and  run  all
##                                    examples  of Intel(R) MKL  for Intel(R) 64 processor
##                                    family applications, static linking
##
## make sointel64                     - build  by  Intel(R) C Compiler  (as default)
##                                    and run all examples of Intel(R) MKL for Intel(R) 64
##                                    processor  family  applications,   dynamic
##                                    linking
##------------------------------------------------------------------------------

include trans.lst

ifndef function
  function = $(TRANS)
endif

ifndef compiler
  compiler=intel
endif

ifndef interface
  interface=lp64
endif

ifndef threading
  threading=parallel
endif

RES = $(addsuffix .res ,$(function))

ifndef MKLROOT
  MKLROOT = ../..
endif

MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

LOPTS = 
COPTS = 

ifeq ($(compiler),gnu)
  override CC = gcc
  IFACE_COMP_PART=intel
  IFACE_THREADING_PART=gnu
ifeq ($(RES_EXT),so)
    LOPTS = -Wl,--no-as-needed
endif
  COPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
else
  override CC=icc
  IFACE_COMP_PART=intel
  IFACE_THREADING_PART=intel
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

ifeq ($(threading),sequential)
  THREADING_LIB = mkl_sequential
  OMP_LIB =
else
ifeq ($(compiler),gnu)
  COPTS += -fopenmp
  THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
  OMP_LIB = -lgomp
else
  COPTS += -qopenmp
  THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
  OMP_LIB = -L$(CMPLR_PATH) -liomp5
endif
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
  MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group $(OMP_LIB)
else
  MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB) $(OMP_LIB)
endif

ifeq ($(_IA),ia32)
  RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)
else
  RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)
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
	$(CC) -Wall -Werror -std=c99 $(SPEC_OPT) $(COPTS) -I"$(MKLROOT)/include" $< source/common_func.c $(LOPTS) $(MKL_LIBS) -lpthread -lm -ldl -o $(RES_DIR)/$*.out
	export LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(CMPLR_PATH); $(RES_DIR)/$*.out >$(RES_DIR)/$@
#-------------------------------------------------------------------------------
