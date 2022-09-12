!===============================================================================
! Copyright 2011-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

mkl/examples/dftf directory contains examples of using Discrete Fourier
Transform Interface (DFTI) for Fortran to compute various FFT problems.

Each example is a self-contained Fortran program.

Examples with prefix 'basic_' demonstrate simple ways to compute complex and
real Fast Fourier transforms. The programs may be used as quick start
templates. Files containing '_dp_'/'_sp_' in their names use double and single
precision, respectively.

Examples with prefix 'config_' demonstrate usage of DFTI configuration
parameters.

For every FFT computed in the example an initialization and verification
function is provided. Initialization function shows how the input data is
indexed and what input will produce a unit peak in the result. Verification
function checks if the unit peak is produced by the computation.

Every example uses dynamically allocated arrays for the data.  Sizes of the
transforms and parameters for verification are selected randomly.

Note: In some examples that show in-place real-to-complex and complex-to-real
transforms the different type arrays are storage associated by use of Cray
pointers, which is not a standard Fortran feature.

Refer to "Fast Fourier Transforms" chapter in Intel(R) MKL Reference Manual for detailed
description of the functions used in these examples.

Your feedback on the examples is welcome at Intel(R) MKL Forum site:
http://software.intel.com/en-us/forums/intel-math-kernel-library

