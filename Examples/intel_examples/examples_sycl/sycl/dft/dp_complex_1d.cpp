/*******************************************************************************
* Copyright 2019-2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*
*  Content:
*       This example demonstrates use of oneAPI Math Kernel Library (oneMKL)
*       API oneapi::mkl::dft to perform 1-D Double Precision Complex to Complex
*       Fast-Fourier Transform on a SYCL device (Host, CPU, GPU).
*
*       The supported floating point data types for data are:
*           double
*           std::complex<double>
*
*******************************************************************************/

#include <vector>
#include <iostream>
#include <CL/sycl.hpp>
#include "oneapi/mkl/dfti.hpp"

#include <stdexcept>
#include <cfloat>
#include <cstddef>
#include <limits>
#include <type_traits>
#include "mkl.h"

// local includes
#define NO_MATRIX_HELPERS
#include "common_for_examples.hpp"

typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::COMPLEX> descriptor_t;

constexpr int SUCCESS = 0;
constexpr int FAILURE = 1;
constexpr double TWOPI = 6.2831853071795864769;

// Compute (K*L)%M accurately
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

static void init(double *x, int N, int H)
{
    for (int n = 0; n < N; ++n) {
        double phase  = moda(n, H, N) / N;
        x[2*n+0] = cos(TWOPI * phase) / N;
        x[2*n+1] = sin(TWOPI * phase) / N;
    }
}

static int verify_fwd(double* x, int N, int H) {
    // Note: this simple error bound doesn't take into account error of
    //       input data
    double errthr = 5.0 * log((double) N) / log(2.0) * DBL_EPSILON;
    std::cout << "\t\tVerify the result, errthr = " << errthr << std::endl;

    double maxerr = 0.0;
    for (int n = 0; n < N; ++n) {
        double re_exp = ((n - H) % N == 0) ? 1.0 : 0.0;
        double im_exp = 0.0;

        double re_got = x[2*n+0];  // real component
        double im_got = x[2*n+1];  // imaginary component
        double err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            std::cout << "\t\tx[" << n << "]: "
                      << "expected (" << re_exp << "," << im_exp << "), "
                      << "got (" << re_got << "," << im_got << "), "
                      << "err " << err << std::endl;
            std::cout << "\t\tVerification FAILED" << std::endl;
            return FAILURE;
        }
    }
    std::cout << "\t\tVerified, maximum error was " << maxerr << std::endl;
    return SUCCESS;
}

static int verify_bwd(double* x, int N, int H) {
    // Note: this simple error bound doesn't take into account error of
    //       input data
    double errthr = 5.0 * log((double) N) / log(2.0) * DBL_EPSILON;
    std::cout << "\t\tVerify the result, errthr = " << errthr << std::endl;

    double maxerr = 0.0;
    for (int n = 0; n < N; ++n) {
        double phase  = moda(n, H, N) / N;
        double re_exp = cos(TWOPI * phase) / N;
        double im_exp = sin(TWOPI * phase) / N;

        int index = 2*n;
        double re_got = x[index+0];  // real component
        double im_got = x[index+1];  // imaginary component
        double err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            std::cout << "\t\tdata[" << n << "]: "
                      << "expected (" << re_exp << "," << im_exp << "), "
                      << "got (" << re_got << "," << im_got << "), "
                      << "err " << err << std::endl;
            std::cout << "\t\tVerification FAILED" << std::endl;
            return FAILURE;
        }
    }
    std::cout << "\t\tVerified, maximum error was " << maxerr << std::endl;
    return SUCCESS;
}

int run_dft_example(cl::sycl::device &dev) {
    //
    // Initialize data for DFT
    //
    int N = 16;
    int harmonic = 5;
    int result = FAILURE;

    //
    // Execute DFT
    //
    try {
        // Catch asynchronous exceptions
        auto exception_handler = [] (cl::sycl::exception_list exceptions) {
            for (std::exception_ptr const& e : exceptions) {
                try {
                    std::rethrow_exception(e);
                } catch(cl::sycl::exception const& e) {
                    std::cout << "Caught asynchronous SYCL exception:" << std::endl
                              << e.what() << std::endl;
                }
            }
        };

        // create execution queue with asynchronous error handling
        cl::sycl::queue queue(dev, exception_handler);

        // Setting up SYCL buffer and initialization
        double* x = (double*) mkl_calloc(N*2, sizeof(double), 64);
        init(x, N, harmonic);
        cl::sycl::buffer<double, 1> xBuffer(x, cl::sycl::range<1>(N*2));
        xBuffer.set_write_back(false);

        // Setting up USM and initialization
        double *x_usm = (double*) malloc_shared(N*2*sizeof(double), queue.get_device(), queue.get_context());
        init(x_usm, N, harmonic);

        descriptor_t desc(N);
        desc.set_value(oneapi::mkl::dft::config_param::BACKWARD_SCALE, (double)(1.0/N));
        desc.commit(queue);

        // Using SYCL buffer
        std::cout<<"\tUsing SYCL buffers"<<std::endl;
        oneapi::mkl::dft::compute_forward(desc, xBuffer);
        {
          auto xAcc = xBuffer.get_access<cl::sycl::access::mode::read>();
          result = verify_fwd(xAcc.get_pointer(), N, harmonic);
        }

        if (result != FAILURE) {
            oneapi::mkl::dft::compute_backward(desc, xBuffer);
            auto xAcc = xBuffer.get_access<cl::sycl::access::mode::read>();
            result = verify_bwd(xAcc.get_pointer(), N, harmonic);
        }

        if(result != SUCCESS) return result;
        result = FAILURE;

        // Using USM
        std::cout<<"\tUsing USM"<<std::endl;
        cl::sycl::event fwd, bwd;
        fwd = oneapi::mkl::dft::compute_forward(desc, x_usm);
        fwd.wait();
        result = verify_fwd(x_usm, N, harmonic);

        if (result != FAILURE) {
            bwd = oneapi::mkl::dft::compute_backward(desc, x_usm);
            bwd.wait();
            result = verify_bwd(x_usm, N, harmonic);
        }

        free(x_usm, queue.get_context());
        mkl_free(x);
    }
    catch(cl::sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during FFT" << std::endl;
        std::cout << "\t\t" << e.what() << std::endl;
        std::cout << "\t\tOpenCl status: " << e.get_cl_code() << std::endl;
    }
    catch(std::runtime_error const& e) {
        std::cout << "\t\truntime exception during FFT" << std::endl;
        std::cout << "\t\t" << e.what() << std::endl;
    }

    return result;
}

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# 1D FFT Complex-Complex Double-Precision Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   dft" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   double" << std::endl;
    std::cout << "#   std::complex<double>" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << std::endl;
}

//
// Main entry point for example.
//
// Dispatches to appropriate device types as set at build time with flag:
// -DSYCL_DEVICES_host -- only runs host implementation
// -DSYCL_DEVICES_cpu -- only runs SYCL CPU implementation
// -DSYCL_DEVICES_gpu -- only runs SYCL GPU implementation
// -DSYCL_DEVICES_all (default) -- runs on all: host, cpu and gpu devices
//
//  For each device selected and each supported data type, Basic_Dp_C2C_1D_FFTExample
//  runs is with all supported data types
//
int main() {
    print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    int returnCode = 0;
    for (auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {
        cl::sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

        if (my_dev_is_found) {
            if (!isDoubleSupported(my_dev)) {
                std::cout << "Double precision not supported on this device " << std::endl;
                std::cout << std::endl;
                continue;
            }
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

            std::cout << "\tRunning with double precision complex-to-complex 1-D FFT:" << std::endl;
            int status = run_dft_example(my_dev);
            if (status != SUCCESS) {
                std::cout << "\tTest Failed" << std::endl << std::endl;
                returnCode = status;
            } else {
                std::cout << "\tTest Passed" << std::endl << std::endl;
            }
        } else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[*it] << " devices found; Fail on missing devices is enabled." << std::endl;
            return 1;
#else
            std::cout << "No " << sycl_device_names[*it] << " devices found; skipping " << sycl_device_names[*it] << " tests." << std::endl << std::endl;
#endif
        }
    }

    return returnCode;
}
