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
*       API oneapi::mkl::dft to perform Batch 1-D Single Precision Real to Complex
*       Fast-Fourier Transform on a SYCL device (Host, CPU, GPU).
*
*       The supported floating point data types for data are:
*           float
*           std::complex<float>
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

typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE, oneapi::mkl::dft::domain::REAL> descriptor_t;

constexpr int SUCCESS = 0;
constexpr int FAILURE = 1;

static void init(float *x, int N, int BATCH, int DIST)
{
    for (int k = 0; k < BATCH; k++) {
        for (int i = 0; i < N; i++) {
            x[i + k*DIST] = cosf(M_PI*(i%5)/19.0f)*sinf(M_PI*(i%7)/41.0f);
        }
    }
}

static int verify(cl::sycl::buffer<float, 1> xBuffer, float *in, int N, int BATCH, int DIST) {
    // Note: this simple error bound doesn't take into account error of
    //       input data
    float errthr = N * FLT_EPSILON;
    std::cout << "\t\tVerify the result, errthr = " << errthr << std::endl;

    auto x = xBuffer.template get_access<cl::sycl::access::mode::read>();
    float maxerr = 0.0f;
    for (int k = 0; k < BATCH; ++k) {
        for (int n = 0; n < N; ++n) {
            float re_exp = in[n+k*DIST];
            float re_got = x[n+k*DIST];

            float err = fabsf(re_got - re_exp);
            if (err > maxerr) maxerr = err;
            if (!(err < errthr)) {
                std::cout << "\t\tx[" << k << "][" << n << "]: "
                          << "expected (" << re_exp << "), "
                          << "got (" << re_got << "), "
                          << "err " << err << std::endl;
                std::cout << "\t\tVerification FAILED" << std::endl;
                return FAILURE;
            }
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
    int result  = FAILURE;
    // Note: forward  distance should be multiple by 8
    // Note: backward distance should be multiple by 4
    int PADDING = 4;
    int DIST    = ((N/2)+PADDING)*2;
    int BATCH   = 2;

    float* x = (float*) mkl_calloc(DIST*BATCH, sizeof(float), 64);
    float* expected = (float*) mkl_calloc(DIST*BATCH, sizeof(float), 64);
    init(x, N, BATCH, DIST);
    init(expected, N, BATCH, DIST);
    cl::sycl::buffer<float, 1> xBuffer(x, cl::sycl::range<1>(DIST*BATCH));
    xBuffer.set_write_back(false);

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

        descriptor_t desc(N);
        desc.set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, BATCH);
        desc.set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, DIST);
        desc.set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, DIST/2);
        desc.set_value(oneapi::mkl::dft::config_param::BACKWARD_SCALE,  (1.0/N));
        desc.commit(queue);

        oneapi::mkl::dft::compute_forward(desc, xBuffer);
        oneapi::mkl::dft::compute_backward(desc, xBuffer);
        result = verify(xBuffer, expected, N, BATCH, DIST);
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

    mkl_free(x);
    mkl_free(expected);
    return result;
}

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# 1D Batch FFT Real-Complex Single-Precision Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   dft" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   float" << std::endl;
    std::cout << "#   std::complex<float>" << std::endl;
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
//  For each device selected and each supported data type, Basic_Sp_R2C_1D_FFTExample
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
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

            std::cout << "\tRunning with single precision real-to-complex 1-D batch FFT:" << std::endl;
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
