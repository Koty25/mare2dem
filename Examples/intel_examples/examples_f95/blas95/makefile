#===============================================================================
# Copyright 2005-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library BLAS95 examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
	@echo "[F95ROOT=<f95_root_path>] [MKLROOT=<MKL_root_path>]"
	@echo "name           - function name. Please see blas95.lst file."
	@echo "compiler_name  - can be gnu, pgi or intel. Default value is intel."
	@echo "                 Intel (R) Fortran Compiler as default."
	@echo "                 If compiler=gnu then GNU gfortran compiler will be used."
	@echo "                 Note: GNU g77 compiler is not supported."
	@echo "                 If compiler=pgi then PGI Fortran 95 compiler will be used. Available only for intel64."
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "threading_name - can be parallel or sequential. Default value is parallel."
	@echo "parallel_name  - can be intel, pgi (only if compiler=pgi)  or gnu (only if compiler=gnu). Default value is intel."
	@echo "omp_name       - can be iomp5 if parallel=intel or"
	@echo "                 can be iomp5 or gomp if parallel=gnu or"
	@echo "                 can be pgmp if parallel=pgi."
	@echo "                 Default value is iomp5."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."
	@echo "f95_root_path  - Path to blas95 library and .mod files that are built on user side."
	@echo "                 Mandatory if compiler is not intel. To build library and .mod files,"
	@echo "                 go to the directory <MKLROOT>/interfaces/blas95 and run:"
	@echo "                 make <target> [options] FC=<f95compiler> INSTALL_DIR=<f95_root_path>"
	@echo "                 See description of the FC and INSTALL_DIR in the interface makefile."
	@echo "MKL_root_path  - Locate Intel(R) MKL libraries relative to <MKL_root_path>."
	@echo "                 Default value: ../.., unless defined in environment."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=dgemmx  - build by Intel(R) Fortran Compiler (as default) and
##                               run DGEMM  example  for  32-bit  applications,
##                               static linking
##
## make libintel64 compiler=gnu F95ROOT=f95_root_path
##                               - build by gfortran GNU  Fortran  compiler  and  run  all
##                               examples of Intel(R) MKL for Intel(R) 64 processor family
##                               applications, static linking
##
## make  sointel64               - build by Intel(R) Fortran Compiler (as default) and
##                               run all examples of Intel(R) MKL for Intel(R) 64 processor
##                               family applications, dynamic linking
##------------------------------------------------------------------------------

include blas95.lst

ifndef function
   function = $(SRC_F90)
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


ifndef $(gnu_path)
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
   override FC=gfortran
   FOPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
   FOPTS += -Werror
   ifeq ($(RES_EXT),so)
      LOPTS = -Wl,--no-as-needed
   endif
   IFACE_COMP_PART=gf
else
   ifeq ($(compiler),pgi)
      override FC=pgf95
   else
      override FC=ifort
      FOPTS += -warn errors -nogen-interfaces
   endif
   IFACE_COMP_PART=intel
endif

ifeq ($(interface),ilp64)
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_ilp64
   FOPTS += $(if $(filter gnu, $(compiler)), -fdefault-integer-8, -i8)
else
   IFACE_LIB = mkl_$(IFACE_COMP_PART)_lp64
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
   IFACE_LIB = mkl_$(IFACE_COMP_PART)
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
   else
      ifeq ($(omp),pgmp)
        OMP_LIB =
      else
        OMP_LIB = -L$(CMPLR_PATH) -l$(omp)
      endif
   endif
endif

CORE_LIB = mkl_core

ifeq ($(EXT),a)
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group $(OMP_LIB)
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB) $(OMP_LIB)
endif

ifeq ($(_IA),ia32)
   ifeq ($(threading),parallel)
      RES_DIR=_results/$(compiler)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
   RES_DIR=_results/$(compiler)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
else
   ifeq ($(threading),parallel)
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(parallel)_$(omp)_$(_IA)_$(RES_EXT)$Z
   else
      RES_DIR=_results/$(compiler)_$(interface)_$(threading)_$(_IA)_$(RES_EXT)$Z
   endif
endif

ifndef F95ROOT
   F95ROOT = $(MKLROOT)
   ifneq ($(compiler),intel)
     $(error F95ROOT is mandatory parameter for $(compiler))
   endif
endif

ifeq ($(_IA),ia32)
   MODPATH=$(F95ROOT)/include/$(_IA)
   LIB95=$(F95ROOT)/lib/$(_IA)/libmkl_blas95.a
else
   MODPATH=$(F95ROOT)/include/$(_IA)/$(interface)
   LIB95=$(F95ROOT)/lib/$(_IA)/libmkl_blas95_$(interface).a
endif

libia32 lib32:
	$(MAKE) $(RES) EXT=a _IA=ia32 RES_EXT=lib
soia32 so32:
	$(MAKE) $(RES) EXT=so _IA=ia32 RES_EXT=so
libintel64 libem64t:
	$(MAKE) $(RES) EXT=a _IA=intel64 RES_EXT=lib
sointel64 soem64t:
	$(MAKE) $(RES) EXT=so _IA=intel64 RES_EXT=so

.PHONY: clean
clean:
	-rm -rf _results;

#-------------------------------------------------------------------------------

vpath %.f90 source

$(RES): %.res: %.f90
	mkdir -p ./$(RES_DIR)
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)" 
	$(FC) $(SPEC_OPT) $(FOPTS) -I $(MODPATH) $< source/common_func.f \
	$(LIB95) $(LOPTS) $(MKL_LIBS) -lpthread -ldl $(GFORTRAN_LIB) \
	-o $(RES_DIR)/$*.out 
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)":$(CMPLR_PATH); $(RES_DIR)/$*.out <data/$*.d >$(RES_DIR)/$@

#-------------------------------------------------------------------------------
