/*******************************************************************************
* Copyright 2018-2020 Intel Corporation.
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
*       This example demonstrates use of oneapi::mkl::lapack::geqrf
*       to perform QR factorization on a SYCL device (HOST, CPU, GPU).
*
*
*       The supported floating point data types for matrix data are:
*           float
*           double
*           std::complex<float>
*           std::complex<double>
*
*
*******************************************************************************/

// stl includes
#include <iostream>
#include <complex>
#include <vector>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/lapack.hpp"
#include "mkl.h"

// local includes
#include "common_for_examples.hpp"


//
// Main example for QR consisting of initialization of
// general dense A matrix.
// Then the QR factorization
// A = Q * R
// is performed.
// Finally the results are post processed.
//
template <typename data_t>
void run_geqrf_example(const cl::sycl::device& dev) {
    // Matrix data sizes and leading dimension
    std::int64_t m   = 30;
    std::int64_t n   = 24;
    std::int64_t lda = 33;

    oneapi::mkl::side left_right;
    oneapi::mkl::transpose trans;

    // Variable holding status of calculations
    std::int64_t info = 0;

    // Asynchronous error handler
    auto error_handler = [&] (cl::sycl::exception_list exceptions) {
        for (auto const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(oneapi::mkl::lapack::exception const& e) {
                // Handle LAPACK related exceptions happened during asynchronous call
                info = e.info();
                std::cout << "Unexpected exception caught during asynchronous LAPACK operation:\n" << e.what() << "\ninfo: " << e.info() << std::endl;
            } catch(cl::sycl::exception const& e) {
                // Handle not LAPACK related exceptions happened during asynchronous call
                std::cout << "Unexpected exception caught during asynchronous operation:\n" << e.what() << std::endl;
                info = -1;
            }
        }
    };

    // Create execution queue for selected device
    cl::sycl::queue queue(dev, error_handler);
    cl::sycl::context context = queue.get_context();

    // Allocate matrices
    std::int64_t A_size = m * lda;
    std::int64_t tau_size = std::min(m,n);
    std::vector<data_t> A(A_size);
    std::vector<data_t> tau(tau_size);

    // Initialize matrix A
    rand_matrix(A, oneapi::mkl::transpose::nontrans, m, n, lda);

    try {
        // Get sizes of scratchpad for calculations
        std::int64_t geqrf_scratchpad_size = oneapi::mkl::lapack::geqrf_scratchpad_size<data_t>(queue, m, n, lda);

        // Define buffers for the data
        cl::sycl::buffer<data_t> A_buffer{A.data(), A_size};
        cl::sycl::buffer<data_t> tau_buffer{tau.data(), tau_size};
        cl::sycl::buffer<data_t> geqrf_scratchpad{geqrf_scratchpad_size};
    
        // Perform factorization
        oneapi::mkl::lapack::geqrf(queue, m, n, A_buffer, lda, tau_buffer, geqrf_scratchpad, geqrf_scratchpad_size);

        // Wait until the calculations are done
        queue.wait_and_throw();
    } catch(oneapi::mkl::lapack::exception const& e) {
        // Handle LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to LAPACK API:\nreason: " << e.what() << "\ninfo: " << e.info() << std::endl;
        info = e.info();
    } catch(cl::sycl::exception const& e) {
        // Handle not LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to SYCL API:\n" << e.what() << std::endl;
        info = -1;
    }

    std::cout << "geqrf " << ((info == 0) ? "ran OK" : "FAILED") << std::endl;

    return;
}


//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# QR Factorization Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Computes QA Factorization A = Q * R" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a general dense matrix." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   geqrf" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   float" << std::endl;
    std::cout << "#   double" << std::endl;
    std::cout << "#   std::complex<float>" << std::endl;
    std::cout << "#   std::complex<double>" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << std::endl;

}


//
// Main entry point for example.
//
// Dispatches to appropriate device types as set at build time with flag:
// -DSYCL_DEVICES_host -- only runs SYCL HOST device
// -DSYCL_DEVICES_cpu -- only runs SYCL CPU device
// -DSYCL_DEVICES_gpu -- only runs SYCL GPU device
// -DSYCL_DEVICES_all (default) -- runs on all: HOST, CPU and GPU devices
//
//  For each device selected and each data type supported, QR factorization example
//  runs with all supported data types
//
int main(int argc, char **argv) {

    print_example_banner();

    // find list of devices
    std::list<my_sycl_device_types> listOfDevices;
    set_list_of_devices(listOfDevices);

    for(auto &deviceType : listOfDevices) {
        cl::sycl::device myDev;
        bool myDevIsFound = false;
        get_sycl_device(myDev, myDevIsFound, deviceType);

        if(myDevIsFound) {
            std::cout << std::endl << "Running geqrf examples on " << sycl_device_names[deviceType] << "." << std::endl;

            std::cout << "\tRunning with single precision real data type:" << std::endl;
            run_geqrf_example<float>(myDev);

            std::cout << "\tRunning with single precision complex data type:" << std::endl;
            run_geqrf_example<std::complex<float>>(myDev);

            if (isDoubleSupported(myDev)) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                run_geqrf_example<double>(myDev);

                std::cout << "\tRunning with double precision complex data type:" << std::endl;
                run_geqrf_example<std::complex<double>>(myDev);
            } else {
                std::cout << "\tDouble precision not supported on this device " << std::endl;
                std::cout << std::endl;
            }
        }
        else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[deviceType] << " devices found; Fail on missing devices is enabled.\n";
            return 1;
#else
            std::cout << "No " << sycl_device_names[deviceType] << " devices found; skipping " << sycl_device_names[deviceType] << " tests.\n";
#endif
        }
    }
    return 0;
}
