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
##      Intel(R) Math Kernel Library BLAS examples creation and run
##******************************************************************************

help:
	@echo "Usage: make {libia32|soia32|libintel64|sointel64} [function=name]"
	@echo "[compiler=compiler_name] [interface=interface_name] [threading=threading_name]"
	@echo "[parallel=parallel_name] [omp=omp_name] [gnu_path=gnu_lib_path]"
	@echo "name           - function name. Please see blas.lst file."
	@echo "compiler_name  - can be gnu, pgi, or intel. Default value is intel."
	@echo "                 Intel (R) Fortran Compiler as default."
	@echo "                 If compiler=gnu then GNU gfortran compiler will be used."
	@echo "                 If compiler=pgi then PGI Fortran 95 compiler will be used. Available only for intel64."
	@echo "interface_name - can be lp64 or ilp64 for intel64. Default value is lp64."
	@echo "threading_name - can be parallel or sequential. Default value is parallel."
	@echo "parallel_name  - can be intel, pgi (only if compiler=pgi)  or gnu (only if compiler=gnu). Default value is intel."
	@echo "omp_name       - can be iomp5 if parallel=intel or"
	@echo "                 can be iomp5 or gomp if parallel=gnu or"
	@echo "                 can be pgmp if parallel=pgi."
	@echo "Default value is iomp5."
	@echo "gnu_lib_path   - If you are using gnu threading layer, specify path to GNU libraries,"
	@echo "                 such as libgomp and libgfortran, with gnu_path variable."
	@echo "                 Default value is /usr/lib."

##------------------------------------------------------------------------------
## examples of using:
##
## make libia32 function=dgemm     - build by Intel(R) Fortran Compiler (as default) and
##                                 run DGEMM  example  for  32-bit  applications,
##                                 static linking
##
## make soia32 compiler=gnu        - build by GNU gfortran compiler  and  run
##                                 all examples of Intel(R) MKL for 32-bit applications, dynamic
##                                 linking
##
## make libintel64 compiler=gnu    - build by gfortran GNU  Fortran  compiler  and  run  all
##                                 examples of Intel(R) MKL for Intel(R) 64 processor family
##                                 applications, static linking
##
## make sointel64                  - build by Intel(R) Fortran Compiler (as default) and
##                                 run all examples of Intel(R) MKL for Intel(R) 64 processor
##                                 family applications, dynamic linking
##------------------------------------------------------------------------------

include blas.lst

ifndef function
   function = $(BLAS) $(BLAS_90)
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

RES = $(addsuffix x.res ,$(function))
RES_90 = $(addsuffix x.res ,$(filter $(BLAS_90), $(function)))
RES_77 = $(addsuffix x.res ,$(filter-out $(BLAS_90), $(function)))

ifndef MKLROOT
MKLROOT = ../..
endif
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
CMPLR_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

FOPTS = -I$(MKLROOT)/include
LOPTS =
MODDIROPT = -module

ifeq ($(compiler),gnu)
    override FC = gfortran
    IFACE_COMP_PART=gf
    FOPTS += $(if $(filter ia32, $(_IA)), -m32, -m64)
    FOPTS += -fcray-pointer
    FOPTS += -Werror
    MODDIROPT = -J
    ifeq ($(RES_EXT),so)
        LOPTS = -Wl,--no-as-needed
    endif
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
      CMPLR_PATH=
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
    MKL_LIBS = -Wl,--start-group $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) -Wl,--end-group 
else
    MKL_LIBS = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
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

libia32 lib32 soia32 so32: RES := $(if $(filter intel, $(compiler)), $(RES_77), $(RES))
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

vpath %.f source
vpath %.f90 source

$(RES_77): %.res: %.f
	mkdir -p ./$(RES_DIR)/$*_modules
	$(FC) $(SPEC_OPT) $(FOPTS)  $< source/common_func.f \
	$(LOPTS) $(MKL_LIBS) \
	$(GFORTRAN_LIB) $(OMP_LIB) -lpthread  -ldl \
	-o $(RES_DIR)/$*.out $(MODDIROPT) $(RES_DIR)/$*_modules
	rm -rf ./$(RES_DIR)/$*_modules
	@echo ------------------
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)":$(CMPLR_PATH); \
	$(RES_DIR)/$*.out <data/$*.d >$(RES_DIR)/$@


$(RES_90): %.res: %.f90
	mkdir -p ./$(RES_DIR)/$*_modules
	$(FC) $(SPEC_OPT) $(FOPTS) $< source/common_func.f \
	$(LOPTS) $(MKL_LIBS) \
	$(GFORTRAN_LIB) $(OMP_LIB) -lpthread  -ldl \
	-o $(RES_DIR)/$*.out $(MODDIROPT) $(RES_DIR)/$*_modules
	rm -rf ./$(RES_DIR)/$*_modules
	@echo ------------------
	export LD_LIBRARY_PATH=$(MKL_PATH):"$(LD_LIBRARY_PATH)":"$(gnu_path)":$(CMPLR_PATH); \
	$(RES_DIR)/$*.out <data/$*.d >$(RES_DIR)/$@

#-------------------------------------------------------------------------------
