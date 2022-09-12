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
*       API oneapi::mkl::dft to perform 3-D Double Precision Real to Complex
*       Fast-Fourier Transform on a SYCL device (Host, CPU).
*
*       The supported floating point data types for data are:
*           double
*           std::complex<double>
*
*******************************************************************************/

#define _USE_MATH_DEFINES
#include <cmath>
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

typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::DOUBLE, oneapi::mkl::dft::domain::REAL> descriptor_t;

constexpr int SUCCESS = 0;
constexpr int FAILURE = 1;
constexpr double TWOPI = 6.2831853071795864769;

// Compute (K*L)%M accurately
static double moda(int K, int L, int M)
{
    return (double)(((long long)K * L) % M);
}

// Initialize array data(N) to produce unit peaks at data(H) and data(N-H)
static void init_r(double *data, int N1, int H1)
{
    // Generalized strides for row-major addressing of data
    int S1 = 1;

    double factor = ((2 * (N1 - H1) % N1 == 0)) ? 1.0 : 2.0;

    for (int n1 = 0; n1 < N1; n1++) {
        double phase = moda(n1, H1, N1) / N1;
        int index = n1*S1;
        data[index] = factor * cos(TWOPI * phase) / N1;
    }
}

// Verify that data has unit peak at H
static int verify_c(double *data, int N1, int H1)
{
    // Note: this simple error bound doesn't take into account error of
    //       input data
    double errthr = 2.5 * log((double) N1) / log(2.0) * DBL_EPSILON;
    std::cout << "\t\tVerify the result, errthr = " << errthr << std::endl;

    // Generalized strides for row-major addressing of data
    int S1 = 1;

    double maxerr = 0.0;
    for (int n1 = 0; n1 < N1/2+1; n1++) {
        double re_exp = ((n1 - H1) % N1 == 0) ||
            ((-n1 - H1) % N1 == 0) ? 1.0 : 0.0;
        double im_exp = 0.0;

        int index = n1*S1;
        double re_got = data[index*2+0];
        double im_got = data[index*2+1];
        double err  = fabs(re_got - re_exp) + fabs(im_got - im_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            std::cout << "\t\tdata"
                      << "[" << n1 << "]: "
                      << " expected (" << re_exp << "," << im_exp << ")"
                      << " got (" << re_got << "," << im_got << ")"
                      << " err " << err << std::endl;
            std::cout << "\t\tVerification FAILED" << std::endl;
            return FAILURE;
        }
    }
    std::cout << "\t\tVerified, maximum error was " << maxerr << std::endl;
    return SUCCESS;
}

// Initialize array data(N) to produce unit peak at data(H)
static void init_c(double *data, int N1, int H1)
{
    // Generalized strides for row-major addressing of data
    int S1 = 1;

    for (int n1 = 0; n1 < N1/2+1; n1++) {
        double phase = moda(n1, H1, N1) / N1;
        int index = n1*S1;
        data[index*2+0] =  cos(TWOPI * phase) / N1;
        data[index*2+1] = -sin(TWOPI * phase) / N1;
    }
}

/* Verify that data has unit peak at H */
static int verify_r(double* data, int N1, int H1)
{
    // Generalized strides for row-major addressing of data
    int S1 = 1;

    // Note: this simple error bound doesn't take into account error of
    //       input data
    double errthr = 2.5 * log((double) N1) / log(2.0) * DBL_EPSILON;
    std::cout << "\t\tVerify the result, errthr = " << errthr << std::endl;

    double maxerr = 0.0;
    for (int n1 = 0; n1 < N1; n1++) {
        double re_exp = ((n1 - H1) % N1 == 0) ? 1.0 : 0.0;

        int index = n1*S1;
        double re_got = data[index];
        double err  = fabs(re_got - re_exp);
        if (err > maxerr) maxerr = err;
        if (!(err < errthr)) {
            std::cout << "\t\tdata"
                      << "[" << n1 << "]: "
                      << " expected (" << re_exp << ")"
                      << " got (" << re_got << ")"
                      << " err " << err << std::endl;
            std::cout << "\t\tVerification FAILED" << std::endl;
            return FAILURE;
        }
    }
    std::cout << "\t\tVerified, maximum error was " << maxerr << std::endl;
    return SUCCESS;
}

int run_dft_forward_example(cl::sycl::device &dev) {
    //
    // Initialize data for DFT
    //
    int N1 = 4;
    // Arbitrary harmonic used to verify FFT
    int H1 = -1;
    // Strides describe data layout in real and conjugate-even domain
    std::int64_t rs[2] = {0, 1};
    std::int64_t cs[2] = {0, 1};
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
        double* x = (double*) mkl_calloc(N1*2 + 2, sizeof(double), 64);
        init_r(x, N1, H1);
        cl::sycl::buffer<double, 1> xBuffer(x, cl::sycl::range<1>(N1*2 + 2));
        xBuffer.set_write_back(false);

        // Setting up USM and initialization
        double *x_usm = (double*) malloc_shared((N1*2 + 2)*sizeof(double), queue.get_device(), queue.get_context());
        init_r(x_usm, N1, H1);

        descriptor_t desc(N1);
        desc.set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, rs);
        desc.set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, cs);
        desc.commit(queue);

        // Using SYCL buffer
        std::cout<<"\tUsing SYCL buffers"<<std::endl;
        oneapi::mkl::dft::compute_forward(desc, xBuffer);

        auto xAcc = xBuffer.get_access<cl::sycl::access::mode::read>();
        result = verify_c(xAcc.get_pointer(), N1, H1);

        if(result != SUCCESS) return result;
        result = FAILURE;

        // Using USM
        std::cout<<"\tUsing USM"<<std::endl;
        auto fwd = oneapi::mkl::dft::compute_forward(desc, x_usm);
        fwd.wait();
        result = verify_c(x_usm, N1, H1);

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

int run_dft_backward_example(cl::sycl::device &dev) {
    //
    // Initialize data for DFT
    //
    int N1 = 4;
    // Arbitrary harmonic used to verify FFT
    int H1 = -1;
    // Strides describe data layout in real and conjugate-even domain
    std::int64_t rs[2] = {0, 1};
    std::int64_t cs[2] = {0, 1};
    int status = FAILURE;

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
        double* x = (double*) mkl_calloc(N1*2 + 2, sizeof(double), 64);
        init_c(x, N1, H1);
        cl::sycl::buffer<double, 1> xBuffer(x, cl::sycl::range<1>(N1*2 + 2));

        // Setting up USM and initialization
        double *x_usm = (double*) malloc_shared((N1*2 + 2)*sizeof(double), queue.get_device(), queue.get_context());
        init_c(x_usm, N1, H1);

        descriptor_t desc(N1);
        desc.set_value(oneapi::mkl::dft::config_param::INPUT_STRIDES, cs);
        desc.set_value(oneapi::mkl::dft::config_param::OUTPUT_STRIDES, rs);
        desc.commit(queue);


        // Using SYCL buffer
        std::cout<<"\tUsing SYCL buffers"<<std::endl;
        oneapi::mkl::dft::compute_backward(desc, xBuffer);

        auto xAcc = xBuffer.get_access<cl::sycl::access::mode::read>();
        status = verify_r(xAcc.get_pointer(), N1, H1);

        if(status != SUCCESS) return status;
        status = FAILURE;

        // Using USM
        std::cout<<"\tUsing USM"<<std::endl;
        auto bwd = oneapi::mkl::dft::compute_backward(desc, x_usm);
        bwd.wait();
        status = verify_r(x_usm, N1, H1);

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

    return status;
}


//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# 1D FFT Real-Complex Double-Precision Example: " << std::endl;
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
//  For each device selected and each supported data type, Basic_Dp_R2C_1D_FFTExample
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
            int status;
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

            std::cout << "\tRunning with double precision real-to-complex 1-D FFT:" << std::endl;
            status = run_dft_forward_example(my_dev);
            if (status != SUCCESS) {
                std::cout << "\tTest Forward Failed" << std::endl << std::endl;
                returnCode = status;
            } else {
                std::cout << "\tTest Forward Passed" << std::endl << std::endl;
            }
            status = run_dft_backward_example(my_dev);
            if (status != SUCCESS) {
                std::cout << "\tTest Backward Failed" << std::endl << std::endl;
                returnCode = status;
            } else {
                std::cout << "\tTest Backward Passed" << std::endl << std::endl;
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
