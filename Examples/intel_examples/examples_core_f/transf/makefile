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
	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
	@echo
	@echo "name     - function name. Please see trans.lst file"
	@echo
	@echo "compiler_name     - can be gnu, pgi (only for intel64) or intel. Default value is intel."
	@echo "           Intel(R) Fortran Compiler as default"
	@echo
	@echo "interface_name    - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo
	@echo "threading_name    - can be parallel or sequential. Default value is parallel."
	@echo
	@echo "parallel_name     - can be intel, pgi (only if compiler=pgi)  or gnu (only if compiler=gnu). Default value is intel."
	@echo
	@echo "omp_name          - can be iomp5 if parallel=intel or"
	@echo "                    can be iomp5 or gomp if parallel=gnu or"
	@echo "                    can be pgmp if parallel=pgi."
	@echo "Default value is iomp5."
	@echo
	@echo "gnu_lib_path      - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                    such as libgomp and libgfortran, with gnu_path variable."
	@echo "                    Default value is /usr/lib."
##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=simatcopy  - build  by  Intel(R) Fortran Compiler  (as default)
##                                    and  run  simatcopy  example  for 32-bit
##                                    applications, static linking
##
## make soia32 compiler=gnu           - build  by  GNU Fortran 95 compiler  and  run  all
##                                    examples of Intel(R) MKL  for 32-bit  applications,
##                                    dynamic linking
##
## make libintel64 compiler=gnu       - build  by  GNU  Fortran 95 compiler  and  run  all
##                                    examples  of Intel(R) MKL  for Intel(R) 64 processor
##                                    family applications, static linking
##
## make sointel64                     - build  by  Intel(R) Fortran Compiler  (as default)
##                                    and run all examples of Intel(R) MKL for Intel(R) 64
##                                    processor  family  applications,   dynamic
##                                    linking
##------------------------------------------------------------------------------
include trans.lst

ifndef function
function = $(TRANS)
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

ifneq ($(interface),ilp64)
   override interface=lp64
endif

ifneq ($(threading),sequential)
   override threading=parallel
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

ifndef gnu_path
   gnu_path=/usr/lib
endif

RES = $(addsuffix .res ,$(function))

ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

FOPTS =
LOPTS = 

ifeq ($(compiler),gnu)
   override FC = gfortran
   FOPTS.DIAG = -Wall -Werror
   ifeq ($(RES_EXT),so)
     LOPTS = -Wl,--no-as-needed
   endif
   IFACE_COMP_PART=gf
else
   ifeq ($(compiler),pgi)
      FC=pgf95
      FOPTS.DIAG = 
   else
      override FC=ifort
      FOPTS.DIAG = -warn all -warn errors
   endif
   IFACE_COMP_PART=intel
endif

FOPTS += $(FOPTS.DIAG)

ifeq ($(interface),ilp64)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_ilp64
   FOPTS += $(if $(filter gnu, $(compiler)), -fdefault-integer-8, -i8)
else
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_lp64
endif

ifeq ($(_IA),ia32)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)
   ifeq ($(compiler),gnu)
       FOPTS += -m32
   endif
else
   ifeq ($(compiler),gnu)
       FOPTS += -m64
   endif
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

ifeq ($(_IA),ia32)
RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
MODPATH=-I $(F95ROOT)/include
else
RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
endif

ifeq ($(threading),sequential)
   THREADING_LIB = mkl_sequential
   OMP_LIB =
   GFORTRAN_LIB=
else
   THREADING_LIB = mkl_$(IFACE_THREADING_PART)_thread
   ifeq ($(parallel),pgi)
      FOPTS += -mp
   endif
   ifeq ($(omp),gomp)
      OMP_LIB = -L$(gnu_path) -l$(omp)
      CMPLR_PATH=
   else
      ifeq ($(omp),pgmp)
        OMP_LIB = 
      else
        OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
      endif
   endif
endif

ifndef F95ROOT
   F95ROOT = "$(MKLROOT)"
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group $(OMP_LIB)
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB) $(OMP_LIB)
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
vpath %.f90 source

$(RES): %.res: %.f90
	mkdir -p ./$(RES_DIR)
	$(FC) $(FOPTS) -I "$(MKLROOT)/include" $(MODPATH) $< $(LOPTS) $(MKL_LIBS) -lpthread -lm -ldl -o $(RES_DIR)/$*.out
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)":$(CMPLR_PATH); $(RES_DIR)/$*.out data/$*.d >$(RES_DIR)/$@

#-------------------------------------------------------------------------------
