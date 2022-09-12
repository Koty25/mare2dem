#===============================================================================
# Copyright 2018-2020 Intel Corporation.
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
##      Intel(R) Math Kernel Library SYCL examples creation and run
##******************************************************************************

include all.lst

define usage

Usage:
------
  #> make Target                                       \
          [threading=sequential|tbb]                   \
          [MKLROOT=path]                               \
          [examples=name]                              \
          [domains=domain(s)]                          \
          [fail_on_missing_devices=yes|no]             \
          [sycl_devices=device(s)]

  Targets:

    libintel64   - statically linked examples for processors that use the Intel(R) 64 architecture
    sointel64    - dynamically linked examples for processors that use the Intel(R) 64 architecture
    # rtintel64  (currently not supported) - single dynamic linking example for processors that use the Intel(R) 64 architecture"
    clean        - remove all generated files and folders
    help         - print this help

  Options:

    compiler=intel
        - Use Intel(R) SYCL Compiler, default.

    threading=sequential
        - Use Intel(R) MKL in the sequential mode.
    threading=tbb
        - Use Intel(R) MKL in the tbb threaded mode, default.


    MKLROOT=<MKL_directory>
        - Specifies the location of Intel(R) MKL libraries used to build this example.
          Default: the Intel(R) MKL installation directory (../..).

    examples="domain1/example1 domain2/example2 ..."
        - Specifies the list of domain/examples to run (source name without extension).
          Proper values are:
            Host device:
              $(list_host)
            CPU device:
              $(list_cpu)
            GPU device:
              $(list_gpu)
            All devices:
              $(list_all)
          All examples will be executed by default.

    domains="domain1 domain2 ..."
        - Specifies domain specific sets of examples to be run if "examples" parameter
          is not specified.
          Proper values are:
            $(domains_list)
          All domain's examples will be executed by default.

    fail_on_missing_devices=yes
        - SYCL examples run on host, and when supported, cpu and gpu devices. By default,
          the examples just skip attempt at offloading to cpu or gpu device if they are
          not found, and returns success/failure only for those cases that were run.
          This option causes the example to fail if host and when supported: cpu, gpu
          devices either actually fail or were not found.
    fail_on_missing_devices=no
        - Default behaviour of skipping offload without failure when cpu or gpu device is
          not found. Example will still always run on host device if host is requested.

    sycl_devices="device1 device2 ..."
        - Specifies device list. SYCL examples might be run on host, and when supported,
          cpu and gpu devices. By default (the same as setting sycl_devices=all) examples
          will be run on every supported/founded device. To specify exact device(s) use
          device names as follows: host/cpu/gpu.
          For example: sycl_devices="host"
                       sycl_devices="host gpu"

Usage examples:
--------------
  #> make libintel64 examples="blas/gemm blas/gemv"
        -  Static link of 'blas/gemm' and 'blas/gemv' SYCL
           examples against Intel(R) MKL for Intel(R) 64-bit applications
           using the Intel(R) SYCL compiler, ILP64 interfaces and tbb threading.

  #> make libintel64 threading=sequential
        - Static link of all the examples against Intel(R) MKL for
          Intel(R) 64-bit applications using the Intel(R) SYCL compiler, ILP64
          interfaces and sequential threading.

  NOTE: If running with parallel GNU Make, and your GNU Make Version >= 4.0, we automatically
  add the "--output-sync=target" flag to your make command to have output more readable.

  #> make -j4 libintel64 domains=blas
        - Statically links all BLAS SYCL examples against Intel(R) MKL for
          Intel(R) 64 bit applications using the Intel(R) SYCL compiler, ILP64 interfaces
          and tbb threading, using 4 parallel make job slots.

  #> make -j4 libintel64 domains=lapack sycl_devices=gpu
        - Statically links all LAPACK SYCL examples that support SYCL GPU offload
          against Intel(R) MKL for Intel(R) 64 bit applications using the Intel(R) SYCL
          compiler, ILP64 interfaces and tbb threading, using 4 parallel make job
          slots.  Note that only offload to SYCL GPU device is enabled in example
          (skipping host and all other supported device types)

Note:
-----

Due to link times which are still longer than expected, the number of VM domain examples was reduced.
However, users may copy vml/vml_all.lst to vml/vml.lst in order to have all the examples available.


endef

help: export usage := $(usage)
help: ; @echo "$$usage"


mixedintel64:
	$(MAKE) run _IA=intel64 EXT=mix RES_EXT=mixed

libintel64:
	$(MAKE) run _IA=intel64 EXT=a RES_EXT=lib

sointel64:
	$(MAKE) run _IA=intel64 EXT=so RES_EXT=so

# rtintel64:
#	$(MAKE) run _IA=intel64 EXT=rt RES_EXT=rt

.PHONY: clean
clean:
	-rm -rf _results;

ifndef sycl_devices
  sycl_devices = all
endif

# list of examples or domains to be run.  examples take precedence, otherwise domains are used
ifndef examples
  ifndef domains
    examples = $(sort $(foreach sycl_device,$(sycl_devices),$(list_$(sycl_device))))
  else
    examples = $(sort $(foreach domain,$(domains),$(foreach sycl_device,$(sycl_devices),$($(domain)_$(sycl_device)))))
  endif
endif

# list of targets to be run
RES = $(addsuffix .res, $(examples))

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
$(call check-opt, compiler, intel, )
$(call check-opt, threading, tbb, sequential)
$(call check-opt, MKLROOT, ../.., %)
$(call check-opt, fail_on_missing_devices, no, yes)
$(call check-opt, sycl_devices, all, cpu gpu host)

# other compile and link time parameters
MKL_PATH = "$(MKLROOT)/lib/$(_IA)"
TBB_PATH = "$(TBBROOT)/lib/$(_IA)/gcc4.8"

# set MKL_LOPTS based on compiler, threading
IFACE_LIB = mkl_intel_ilp64
CORE_LIB  = mkl_core
SYCL_LIB  = mkl_sycl

ifeq ($(threading),sequential)
  THREADING_LIB = mkl_sequential
else
  ifeq ($(threading),tbb)
    THREADING_LIB = mkl_tbb_thread
  endif
endif

# threading library linking options
TOPT.sequential=
TOPT.tbb       = -L$(TBB_PATH) -ltbb

# determine MKL library linking
MKL_LOPTS.a = -Wl,-export-dynamic -Wl,--start-group \
             $(MKL_PATH)/$(RES_EXT)$(IFACE_LIB).$(EXT) \
             $(MKL_PATH)/$(RES_EXT)$(THREADING_LIB).$(EXT) \
             $(MKL_PATH)/$(RES_EXT)$(CORE_LIB).$(EXT) \
             -Wl,--end-group \
             $(TOPT.$(threading))

MKL_LOPTS.so  = -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
MKL_LOPTS.mix = -L$(MKL_PATH) -l$(IFACE_LIB) -l$(THREADING_LIB) -l$(CORE_LIB)
MKL_LOPTS.rt  = -L$(MKL_PATH) -lmkl_rt

MKL_LOPTS = $(MKL_LOPTS.$(EXT))

SUPP_LOPTS = -lsycl -lOpenCL -lpthread -lm -ldl

# determine compiler settings
COPTS.$(compiler) = $(compiler)
COPTS = -I$(MKLROOT)/include -Icommon -DMKL_ILP64

# pass flag to exmaple code to fail if no CPU or GPU is found instead of just skipping
ifneq ($(fail_on_missing_devices),no)
  COPTS += -DFAIL_ON_MISSING_DEVICES
endif

# pass flag to example code to only run for specified sycl_devices
COPTS += $(foreach sycl_device,$(sycl_devices),-DSYCL_DEVICES_$(sycl_device))

# determine linker setting
MKL_SYCL_LOPTS.mix  = $(MKL_PATH)/libmkl_sycl.a
MKL_SYCL_LOPTS.a    = $(MKL_PATH)/libmkl_sycl.a
MKL_SYCL_LOPTS.so   = -L$(MKL_PATH) -lmkl_sycl

# LOPTS.$(compiler)
LOPTS.intel = -fsycl-device-code-split=per_kernel
LOPTS       = $(LOPTS.$(compiler)) $(MKL_SYCL_LOPTS.$(EXT)) $(MKL_LOPTS) $(SUPP_LOPTS)

# set compile line
CXX.intel = dpcpp

RES_DIR=_results/$(_IA)_$(RES_EXT)_$(threading)/$(@D)

build_link_run_cmd = \
            mkdir -p ./$(RES_DIR) && \
            $(CXX.$(compiler)) $(COPTS) -c $< -o $(RES_DIR)/$(patsubst %.res,%.o,$(@F)) && \
            $(CXX.$(compiler)) $(RES_DIR)/$(patsubst %.res,%.o,$(@F)) $(LOPTS) -o $(RES_DIR)/$(patsubst %.res,%.out,$(@F)) && \
            LD_LIBRARY_PATH=$(MKL_PATH):$(LD_LIBRARY_PATH) $(RES_DIR)/$(patsubst %.res,%.out,$(@F)) > $(RES_DIR)/$(@F)

banner = @echo ""; \
         echo Building/Running $(@D)/$(patsubst %.res,%,$(@F)) for SYCL_DEVICES = $(sycl_devices)

#-------------------------------------------------------------------------------

run: $(RES)

%.res: %.cpp
	$(call banner)
	$(call build_link_run_cmd, "")
