#===============================================================================
# Copyright 2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library Graph examples creation and run
##******************************************************************************

include graphc.lst

define usage

Usage:
------
  #> make Target                          \
          [compiler=intel|pgi|gnu]        \
          [interface=interface_name]      \
          [threading=omp|tbb|parallel|sequential]  \
          [MKLROOT=path]                  \
          [example=name]

  Targets:

    libia32    - statically linked example for processors that use the IA-32 architecture"
    soia32     - dynamically linked example for processors that use the IA-32 architecture"
    libintel64 - statically linked example for processors that use the Intel(R) 64 architecture"
    sointel64  - dynamically linked example for processors that use the Intel(R) 64 architecture"
    help       - print this help"

  Options:

    compiler=intel
        - Use Intel(R) C Compiler, default.
    compiler=pgi
        - Use PGI C compiler
    compiler=gnu
        - Use GNU C compiler

    interface=interface_name
                - Can be lp64 (default) or ilp64 for intel64.

    threading=omp
        - Use Intel(R) MKL in the threaded mode with OpenMP, default.
    threading=tbb
        - Use Intel(R) MKL in the threaded mode with TBB.
    threading=parallel
        - Use Intel(R) MKL in the default parallel mode (OpenMP currently).
    threading=sequential
        - Use Intel(R) MKL in the sequential mode.

    MKLROOT=<MKL_directory>
        - Specifies the location of Intel(R) MKL libraries used to build this example.
          Default: the Intel(R) MKL installation directory.

    example="example1 example2 ..."
        - Specifies the list of examples to run (source name without extension).
          Proper values are:
          $(GRAPH)
          All examples will be executed by default.

Usage examples:
--------------
  #> make libia32 function="graphc_mxv"
        -  Link 'graphc_mxv' example against parallel static MKL-IA32.

  #> make sointel64 interface=ilp64 threading=sequential
        - Link all the examples against sequential dynamic MKL-Intel64
          using ILP64 interface.

endef

help: export usage := $(usage)
help: ; @echo "$$usage"

# Main targets
libintel64: ; $(MAKE) ext=a  _ia=intel64 run
sointel64:  ; $(MAKE) ext=so _ia=intel64 run
rtintel64:  ; $(MAKE) ext=rt _ia=intel64 run
libia32:    ; $(MAKE) ext=a  _ia=ia32    run
soia32:     ; $(MAKE) ext=so _ia=ia32    run
rtia32:     ; $(MAKE) ext=rt _ia=ia32    run

.PHONY: clean
clean:
	-rm -rf _results;

# Check options and set default values
#   call args: $1 - option name, $2 - default value, $3 - other appropriate values for the option
check-opt = \
    $(if $(value $(strip $1)),,$(eval override $1 = $2)) \
    $(if $(filter $2 $3, $(value $(strip $1))),, \
        $(error [ERROR] Incorrect option value: '$(strip $1)=$(value $(strip $1))'! Should be some of '$(strip $2 $3)'))
$(call check-opt, compiler, intel, pgi gnu)
$(call check-opt, threading, omp, tbb parallel sequential)
$(call check-opt, example, $(GRAPH))
$(call check-opt, interface, lp64, ilp64)
$(call check-opt, MKLROOT, ../.., %)

ifeq ($(compiler)_$(_ia),pgi_ia32)
   $(info ***)
   $(info *** PGI compiler 32-bit support was dropped off starting from Intel(R) MKL 2017.)
   $(info *** Please try >make help)
   $(info ***)
   $(error )
endif

# Paths to Intel(R) MKL libraries
IMKL_PATH  = "$(MKLROOT)/lib/$(_ia)"
IOMP_PATH  = "$(MKLROOT)/../compiler/lib/$(_ia)"

# Path to TBB
TBB_PATH = $(shell if [ -d "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8" ]; then \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.8"; \
	else \
		echo "$(MKLROOT)/../tbb/lib/$(_IA)/gcc4.4"; \
	fi)

# Vendor specific compiler options
COPTS.gnu.ia32  = -m32

COPTS = $(COPTS.$(compiler).$(_ia)) $(COPTS.$(compiler))

omp_additional_libs = 

# Choosing correct OpenMP library (default is libiomp5)
ifeq ($(compiler),gnu)
  compiler_omp_lib = gomp
else
ifeq ($(compiler),pgi)
  compiler_omp_lib = pgmp
  omp_additional_libs = -pgf90libs
else
  compiler_omp_lib = iomp5
  COPTS += -qopenmp
endif
endif

# Macro to use in conditions
iface.$(interface) = $(interface)

# Linker options
LINK_OPTS = $(LINK_OPTS.$(compiler)) $(IMKL_LIBS.$(ext)) $(OMP_LIB) -lpthread -lm
LINK_OPTS.gnu = -Wl,--no-as-needed -ldl
LINK_OPTS.pgi = -Wl,--no-as-needed -ldl

# Intel(R) MKL libraries to use in linking (default is mkl_intel_thread)
ifeq ($(threading),sequential)
  THREADING_LIBNAME=mkl_sequential
  OMP_LIB = 
else
  ifeq ($(threading),omp)
    THREADING_LIBNAME=mkl_$(compiler)_thread
    OMP_LIB = -L"$(IOMP_PATH)" -l$(compiler_omp_lib) $(omp_additional_libs)
  else
    ifeq ($(threading),tbb)
      THREADING_LIBNAME=mkl_tbb_thread
      OMP_LIB.a  = -L"$(TBB_PATH)" -ltbb -lstdc++
      OMP_LIB.so =
      OMP_LIB.rt = -L"$(TBB_PATH)" -ltbb -lstdc++
      OMP_LIB := $(OMP_LIB.$(ext))
    else
      THREADING_LIBNAME=mkl_intel_thread
      OMP_LIB = -L"$(IOMP_PATH)" -liomp5
    endif
  endif
endif

IMKL_LIBS    = $(iface_lib.$(compiler))$(suffix_lib.$(_ia)) mkl_core $(THREADING_LIBNAME)

IMKL_LIBS.a  = -Wl,--start-group $(IMKL_LIBS:%=$(IMKL_PATH)/lib%.a) -Wl,--end-group
IMKL_LIBS.so = -L $(IMKL_PATH) $(IMKL_LIBS:%=-l%)
IMKL_LIBS.rt = -L $(IMKL_PATH) -lmkl_rt

iface_lib.$(compiler) = mkl_intel
iface_lib.gnu         = mkl_intel
suffix_lib.$(_ia)     = _$(interface)
suffix_lib.ia32       = 

# Compiler tools
CC.intel = icc
CC.pgi   = pgcc
CC.gnu   = gcc

ifeq ($(compiler),intel)
  warnflags = -Wall
else
ifeq ($(compiler),pgi)
  warnflags = 
else
  warnflags = -Wall
endif
endif

CC = $(CC.$(compiler)) $(COPTS) -I$(MKLROOT)/include $(warnflags) $(if $(iface.ilp64),-DMKL_ILP64) -I$(MKLROOT)/include

#-------------------------------------------------------------------------------

RES_DIR = _results/$(compiler)_$(interface)_$(_ia)_$(ext)
RES = $(example:%=$(RES_DIR)/%.res)

run   : $(RES)
$(RES): .FORCE
.FORCE:
.PHONY: run

$(RES_DIR)/%.res: $(RES_DIR)/%.exe ; $(shell which env) LD_LIBRARY_PATH=$(IMKL_PATH):$(IOMP_PATH):$(LD_LIBRARY_PATH) $< > $@

.SECONDARY:
$(RES_DIR)/%.exe: $(RES_DIR)/%.o ; $(CC) $^ -o $@ $(LINK_OPTS)

$(RES_DIR)/%.o: source/%.c | $(RES_DIR)/. ; $(CC) -c $< -o $@

%/. : ; mkdir -p $@
