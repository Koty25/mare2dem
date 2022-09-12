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
##     Intel(R) Math Kernel Library Cluster PARDISO C examples makefile
##******************************************************************************

include cluster_sparse_solverc.lst

define usage

Usage:
------
  #> make {sointel64|libintel64}       \
          [mpi=intelmpi|mpich2|openmpi|custom] \
          [compiler=intel|gnu]         \
          [interface=interface_name]   \
          [ompthreads=omp_num_threads] \
          [mpidir=path]                \
          [example=name]

  Targets:

    sointel64  - Dynamically link example with Intel MKL
    libintel64 - Statically link example with Intel MKL

  Options:

    mpi=intelmpi - Use Intel(R) MPI Library 4.x, default
    mpi=mpich2   - Use MPICH2
    mpi=openmpi  - Use Open MPI
    mpi=custom   - Use custom MPI (set MPICC or make sure mpicc is in PATH,
                                   set MPIRUN or make sure mpirun is in PATH)

    compiler=intel - Use Intel(R) C Compiler, default.
    compiler=gnu   - Use GNU C compiler

    interface=interface_name
               - Can be lp64 (default) or ilp64 for intel64.

    ompthreads=omp_num_threads
               - Number of OpenMP threads, should be greater than 1. Default value is 2.

    mpidir=path
               - Path to MPI installation directory.
                 MPI scripts are taken from mpidir/bin (or mpidir/bin64
                 for Intel(R) MPI Library and Intel(R) 64 architecture).
                 If this directory is in PATH you can omit mpidir parameter.

    example="example1 example2 ..."
               - Specifies the list of examples to run (source name without extension).
                 Proper values are: $(CLUSTER_SPARSE_SOLVER)
                 All examples will be executed by default.

Usage example:
--------------
  #> make libintel64 mpi=mpich2 compiler=intel mpidir=/opt/mpich

  - Use Intel(R) C Compiler to compile all CLUSTER_SPARSE_SOLVER examples,
    statically link them against Intel MKL, and run using MPICH from /opt/mpich

endef

help: export usage := $(usage)
help: ; @echo "$$usage"

# Main targets
libintel64: ; $(MAKE) ext=a  run
sointel64:  ; $(MAKE) ext=so run

.PHONY: clean
clean:
	-rm -rf _results;

# Check options and set default values
#   call args: $1 - option name, $2 - default value, $3 - other appropriate values for the option
check-opt = \
    $(if $(value $(strip $1)),,$(eval override $1 = $2)) \
    $(if $(filter $2 $3, $(value $(strip $1))),, \
        $(error [ERROR] Incorrect option value: '$(strip $1)=$(value $(strip $1))'! Should be some of '$(strip $2 $3)'))
$(call check-opt, compiler, intel, gnu)
$(call check-opt, example, $(CLUSTER_SPARSE_SOLVER))
$(call check-opt, mpi, intelmpi, mpich2 openmpi custom)
$(call check-opt, interface, lp64, ilp64)
$(call check-opt, ompthreads, 2, %)
$(call check-opt, MKLROOT, ../.., %)

comma = ,

# Paths to Intel(R) MKL libraries
IMKL_PATH  = "$(MKLROOT)/lib/intel64"
IOMP_PATH  = "$(MKLROOT)/../compiler/lib/intel64"

# Macro to use in conditions
ext.$(ext)         = $(ext)
mpi.$(mpi)         = $(mpi)
iface.$(interface) = $(interface)

# Linker options
LINK_OPTS = $(LINK_OPTS.$(compiler)) $(IMKL_LIBS.$(ext)) -L $(IOMP_PATH) -liomp5 $(if $(mpi.intelmpi),-mt_mpi) $(PTHREAD.$(ext)) -ldl -lm
LINK_OPTS.gnu = -Wl,--no-as-needed
PTHREAD.a    = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
PTHREAD.so   = -lpthread

# Intel(R) MKL libraries to use in linking
IMKL_LIBS    = $(if $(mpi.custom),,$(blacs_lib.$(mpi))) $(iface_lib.$(compiler)) mkl_core mkl_intel_thread
IMKL_LIBS.a  = -Wl,--start-group $(BLACS_LIB.a) $(IMKL_LIBS:%=$(IMKL_PATH)/lib%.a) -Wl,--end-group
IMKL_LIBS.so = $(BLACS_LIB.so) -L $(IMKL_PATH) $(IMKL_LIBS:%=-l%)
BLACS.$(ext) = $(if $(mpi.custom),$(RES_DIR)/lib$(blacs_lib.$(mpi)).$(ext))
BLACS_LIB.a  = $(BLACS.a)
BLACS_LIB.so = $(if $(mpi.custom),-L $(RES_DIR) -Wl$(comma)-rpath$(comma)$(RES_DIR) -l$(blacs_lib.$(mpi)))
iface_lib.$(compiler) = mkl_intel_$(interface)
iface_lib.gnu         = mkl_gf_$(interface)
blacs_lib.$(mpi)   = mkl_blacs_intelmpi_$(interface)
blacs_lib.openmpi  = mkl_blacs_openmpi_$(interface)
blacs_lib.custom   = mkl_blacs_custom_$(interface)

# MPI tools
MPI_PATH = $(if $(mpidir),$(mpidir)/bin$(if $(mpi.intelmpi),64)/)
MPI_CPLR.$(mpi).intel   = $(MPI_PATH)mpicc -cc=icc
MPI_CPLR.$(mpi).gnu     = $(MPI_PATH)mpicc -cc=gcc
MPI_CPLR.openmpi.intel  = OMPI_CC=icc $(MPI_PATH)mpicc
MPI_CPLR.openmpi.gnu    = OMPI_CC=gcc $(MPI_PATH)mpicc
MPI_CPLR.intelmpi.intel = $(MPI_PATH)mpiicc
MPI_EXEC.$(mpi)         = $(MPI_PATH)mpiexec -np 2
MPI_EXEC.mpich2         = $(MPI_PATH)mpiexec -n 2
MPI_EXEC.intelmpi       = $(MPI_PATH)mpiexec -n 4

MPI_CPLR.custom.intel   = $(or $(MPICC),mpicc)
MPI_CPLR.custom.gnu     = $(or $(MPICC),mpicc)
MPI_EXEC.custom         = $(or $(MPIRUN),mpirun) -n 2

#-------------------------------------------------------------------------------

RES_DIR = _results/$(compiler)_$(mpi)_$(interface)_intel64_$(ext)
RES = $(example:%=$(RES_DIR)/%.res)

#-------------------------------------------------------------------------------

ifneq ($(mpi.custom),)
$(BLACS.$(ext)):
	mkdir -p $(@D) && \
	$(MAKE) -C ../../interfaces/mklmpi \
		$(if $(ext.a),lib,so)intel64 \
		INSTALL_LIBNAME=$(basename $(@F)) INSTALL_DIR=$(realpath .)/$(@D)
endif

#-------------------------------------------------------------------------------

run   : $(BLACS.$(ext)) $(RES)
$(RES): .FORCE
.FORCE:
.PHONY: run
.PRECIOUS: $(RES_DIR)/%.exe

$(RES_DIR)/%.res: $(RES_DIR)/%.exe
	$(MPI_EXEC.$(mpi)) $(shell which env) LD_LIBRARY_PATH=$(IMKL_PATH):$(IOMP_PATH):$(LD_LIBRARY_PATH) OMP_NUM_THREADS=$(ompthreads) $< > $@

$(RES_DIR)/%.exe: $(RES_DIR)/%.o $(MODULE)
	$(MPI_CPLR.$(mpi).$(compiler)) $^ -o $@ $(LINK_OPTS)

$(RES_DIR)/%.o: source/%.c | $(RES_DIR)/.
	$(MPI_CPLR.$(mpi).$(compiler)) -Wall $(if $(iface.ilp64),-DMKL_ILP64) -I$(MKLROOT)/include -c -o $@ $<

%/. : ; mkdir -p $@

