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
##      Intel(R) Math Kernel Library SYCL device examples creation and run
##******************************************************************************

include all.lst

define usage

Usage:
------
  #> make                                              \
          [MKLROOT=path]                               \
          [examples=name]                              \
          [domains=domain(s)]                          \
          [fail_on_missing_devices=yes|no]             \
          [sycl_devices=device(s)]

  Options:

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
  #> make examples="rng/uniform"
        -  Run 'rng/uniform' SYCL examples using the Intel(R) DPC++ compiler.

  #> make -j4 domains=rng sycl_devices=gpu
        - Run all RNG SYCL examples using the Intel(R) DPC++ compiler.
          Note that only offload to SYCL GPU device is enabled in example
          (skipping host and all other supported device types)
endef

default:
	$(MAKE) run _IA=intel64

help: export usage := $(usage)
help: ; @echo "$$usage"

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
$(call check-opt, MKLROOT, ../.., %)
$(call check-opt, fail_on_missing_devices, no, yes)
$(call check-opt, sycl_devices, all, cpu gpu host)

# determine compiler settings
COPTS.$(compiler) = $(compiler)
COPTS = -fsycl-device-code-split=per_kernel -fno-sycl-early-optimizations -I$(MKLROOT)/include -Icommon -fsycl-unnamed-lambda

# pass flag to exmaple code to fail if no CPU or GPU is found instead of just skipping
ifneq ($(fail_on_missing_devices),no)
  COPTS += -DFAIL_ON_MISSING_DEVICES
endif

# pass flag to example code to only run for specified sycl_devices
COPTS += $(foreach sycl_device,$(sycl_devices),-DSYCL_DEVICES_$(sycl_device))

# set compile line
CXX.intel = dpcpp

RES_DIR=_results/$(_IA)/$(@D)

build_link_run_cmd = \
            mkdir -p ./$(RES_DIR) && \
            $(CXX.$(compiler)) $(COPTS) $< -o $(RES_DIR)/$(patsubst %.res,%.out,$(@F)) && \
            $(RES_DIR)/$(patsubst %.res,%.out,$(@F)) > $(RES_DIR)/$(@F)

banner = @echo ""; \
         echo Building/Running $(@D)/$(patsubst %.res,%,$(@F)) for SYCL_DEVICES = $(sycl_devices)

#-------------------------------------------------------------------------------

run: $(RES)

%.res: %.cpp
	$(call banner)
	$(call build_link_run_cmd, "")
