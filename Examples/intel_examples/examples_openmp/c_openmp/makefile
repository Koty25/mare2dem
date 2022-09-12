#===============================================================================
# Copyright 2019-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library OpenMP offload C examples creation and run
##******************************************************************************

include all.lst

define usage

Usage:
------
  #> make Target                          \
          [interface=ilp64]               \
          [threading=sequential|omp]      \
          [MKLROOT=path]                  \
          [examples=name]                 \
          [domains=domain]

  Targets:

    libintel64 - statically linked example for processors that use the Intel(R) 64 architecture
    sointel64  - dynamically linked example for processors that use the Intel(R) 64 architecture
    # rtintel64  - single dynamic linking example for processors that use the Intel(R) 64 architecture
    clean      - remove all generated files and folders
    help       - print this help"

  Options:

    interface=ilp64
        - Use 64-bit integers and 64-bit pointers (default)
    interface=lp64
        - Use 32-bit integers and 64-bit pointers (not currently supported by lapack)

    threading=sequential
        - Use Intel(R) MKL in the sequential mode (default)
    threading=omp
        - Use Intel(R) MKL with the OpenMP threading mode with Intel OpenMP libraries.

    MKLROOT=<MKL_directory>
        - Specifies the location of Intel(R) MKL libraries used to build this example.
          Default location is the Intel(R) MKL installation directory (..\\..).

    examples="domain1/example1 domain2/example2 ..."
        - Specifies the list of domain/examples to run (source name without extension).
          Proper values are:
          $(list)
          All examples will be executed by default.

    domains="domain1 domain2 ..."
        - Specifies domain specific sets of examples to be run if "examples" parameter
          is not specified.
          Proper values are:
          $(domains_list)
          All domains examples will be executed by default.


Usage examples:
--------------
  #> make libintel64 examples="blas/sgemm"
        -  Statically link 'blas/sgemm'
           examples against Intel(R) MKL for Intel(R) 64-bit applications
           using the Intel (R) OpenMP offload C compiler, ILP64 interfaces and parallel threading.

  #> make sointel64 interface=lp64 threading=sequential
        - Dynamically link all the examples against Intel(R) MKL for
          Intel(R) 64-bit applications using the Intel (R) OpenMP offload C compiler, LP64
          interfaces and sequential threading.

  #> make libintel64 domains="BLAS"
        - Statically link all the BLAS examples against Intel(R) MKL for
          Intel(R) 64-bit applications using the Intel (R) OpenMP offload C compiler, ILP64 interfaces
          and sequential threading

  NOTE: If running with parallel GNU Make, and your GNU Make Version >= 4.0, we automatically
  add the "--output-sync=target" flag to your make command to have output more readable.

  #> make -j4 libintel64
        - Statically links all component examples against Intel(R) MKL for
		  Intel(R) 64 bit applications using the Intel (R) OpenMP offload C compiler, ILP64 interfaces
		  and sequential threading, using 3 parallel make job slots

  NOTE: the default value of environment variable OMP_TARGET_OFFLOAD is DEFAULT which if gpu is not found,
will execute without offloading. If you want to force it to offload or fail, export OMP_TARGET_OFFLOAD=MANDATORY
before calling this makefile.

endef

help: export usage := $(usage)
help: ; @echo "$$usage"

ifndef backend
  backend=opencl
endif

libintel64:
	$(MAKE) run _IA=intel64 EXT=a RES_EXT=lib
sointel64:
	$(MAKE) run _IA=intel64 EXT=so RES_EXT=so
rtintel64:
	$(MAKE) run _IA=intel64 EXT=rt RES_EXT=rt

.PHONY: clean
clean:
	-rm -rf _results;

# list of examples or domains to be run.  examples take precedence, otherwise domains are used
ifndef examples
  ifndef domains
    examples = $(list)
  else
    examples = $(foreach domain,$(domains),$($(domain)) )
  endif
endif

# lapack/getrf does not support lp64 interface
ifeq ($(interface), lp64)
  example_list := $(filter-out lapack/dgetrf, $(examples))
else
  example_list := $(examples)
endif

# list of targets to be run
RES = $(addsuffix .res, $(example_list))

# GNU Make specific flag with support starting in Make 4.0
# to synchronize output when run in parallel.
ifneq ($(findstring output-sync, $(.FEATURES)),)
GNUMAKEFLAGS += --output-sync=target
endif

# Check options and set default values
#   call args: $1 - option name, $2 - default value, $3 - other appropriate values for the option
check-opt = \
    $(if $(value $(strip $1)),,$(eval override $1 = $2)) \
    $(if $(filter $2 $3, $(value $(strip $1))),, \
        $(error [ERROR] Incorrect option value: '$(strip $1)=$(value $(strip $1))'! Should be some of '$(strip $2 $3)'))
$(call check-opt, compiler, intel)
$(call check-opt, threading, sequential, omp)
$(call check-opt, interface, ilp64, lp64)
$(call check-opt, MKLROOT, ../.., %)

# other compile and link time parameters
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
OMP_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"
OMP_PATH = "$(MKLROOT)/../compiler/lib/$(_IA)"

# set MKL_LIBS based on compiler, threading, interface
IFACE_LIB = mkl_intel_$(interface)
CORE_LIB = mkl_core

ifeq ($(threading),sequential)
  THREADING_LIB = mkl_sequential
  THREAD_SUPP_LIBS=
endif

ifeq ($(threading),omp)
  THREADING_LIB = mkl_$(compiler)_thread
  ifeq ($(compiler),intel)
    THREAD_SUPP_LIBS=-liomp5
  endif
endif

# determine MKL library linking
MKL_LIBS.a = $(MKL_PATH)/libmkl_sycl.a -Wl,--start-group \
             $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) \
             $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) \
             $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) \
             -Wl,--end-group
MKL_LIBS.so = -L$(MKL_PATH) -lmkl_sycl -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
MKL_LIBS.rt = -L$(MKL_PATH) -lmkl_rt
MKL_LIBS = $(MKL_LIBS.$(EXT))

SUPP_LIBS = -lpthread -ldl -lm -lstdc++

OMP_OPTS = -fiopenmp -fopenmp-targets=spir64 -fsycl

# determine compiler settings
COPTS.$(compiler) = $(compiler)
COPTS.$(interface) = $(interface)
COPTS = $(OMP_OPTS) $(if $(COPTS.ilp64),-DMKL_ILP64 -m64) -I$(MKLROOT)/include -Icommon

# determine linker setting
LOPTS.intel =
LOPTS.gnu   = -Wl,--no-as-needed
LOPTS       = $(LOPTS.$(compiler)) $(MKL_LIBS) -lsycl -lOpenCL $(THREAD_SUPP_LIBS) $(SUPP_LIBS)
CC.intel = icx
CC       = $(CC.$(compiler))

RES_DIR=_results/$(_IA)_$(interface)_$(threading)_$(RES_EXT)

#-------------------------------------------------------------------------------

run: $(RES)

%.res: %.c
	@echo ""
	@echo Building/Running $*
	mkdir -p ./$(RES_DIR)/$(@D)
	$(CC) $(COPTS) -c $< -o $(RES_DIR)/$*.o
	$(CC) -fsycl-device-code-split=per_kernel $(OMP_OPTS) $(RES_DIR)/$*.o $(LOPTS) -o $(RES_DIR)/$*.out
	LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH):$(OMP_PATH) $(RES_DIR)/$*.out 2>&1 | tee $(RES_DIR)/$@
