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

mkl/examples/fftw3xf directory contains examples of using FFTW3 Fortran
interface to compute various FFT problems.

Each example is a self-contained Fortran program.

The examples are named by the name of the plan function they use.  Prefixes dp_
and sp_ indicate the floating point precision used in the example, double or
single precision, respectively.  For the FFT problems that
Intel(R) Math Kernel Library (Intel(R) MKL) does not support
examples are not provided.

For every FFT computed in the example an initialization and verification
function is provided. Initialization function shows how the input data is
indexed and what input will produce a unit peak in the result. Verification
function checks if the unit peak is produced by the computation.

Every example uses dynamically allocated arrays for the data.  Sizes of the
transforms and parameters for verification are selected randomly.

Note: In examples that show in-place real-to-complex and complex-to-real
transforms the different type arrays are storage associated by use of Cray
pointers, which is not a standard Fortran feature.

Refer to FFTW3 documentation for detailed description of the functions used in
these examples.  Refer to Intel MKL Reference Manual for limitations of the FFTW3
interface provided by Intel MKL.

Your feedback on the examples is welcome at Intel MKL Forum site:
http://software.intel.com/en-us/forums/intel-math-kernel-library


(set-fill-column 79)
