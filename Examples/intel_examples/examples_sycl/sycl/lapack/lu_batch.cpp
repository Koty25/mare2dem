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
*       This example demonstrates use of oneapi::mkl::lapack::getrf_batch
*       to perform batched LU factorization on a SYCL device (HOST, CPU, GPU).
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
#include <cstdlib>
#include <complex>
#include <vector>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/lapack.hpp"
#include "mkl.h"

// local includes
#include "common_for_examples.hpp"

//
// Main example for batched LU consisting of initialization of
// multiple general dense A matrices.
// Then the LU factorization
// A = P * L * U
// is performed for each.
// Finally the results are post processed.
//
template <typename data_t>
void run_lu_batch_example(const cl::sycl::device& dev) {
    // Matrix data size and leading dimension
    std::int64_t m         = 29;
    std::int64_t n         = 29;
    std::int64_t lda       = 34;
    std::int64_t stride    = 1024;

    std::int64_t batch_size = 3;

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

    cl::sycl::queue main_queue(dev, error_handler);

    std::vector<data_t> A_initial(stride*batch_size);
    std::vector<std::int64_t> ipiv_initial(stride*batch_size);

    for (int i = 0; i < batch_size; i++)
        rand_matrix(A_initial, oneapi::mkl::transpose::nontrans, m, n, lda, i*stride);
    
    try {
        const std::int64_t getrf_batch_scratchpad_size = oneapi::mkl::lapack::getrf_batch_scratchpad_size<data_t>(main_queue, m, n, lda, stride, stride, batch_size);

        cl::sycl::buffer<data_t> A_buf{&A_initial[0], stride*batch_size};
        cl::sycl::buffer<std::int64_t> ipiv_buf{&ipiv_initial[0], stride*batch_size};
        cl::sycl::buffer<data_t> getrf_scratchpad{getrf_batch_scratchpad_size};

        oneapi::mkl::lapack::getrf_batch(main_queue, m, n, A_buf, lda, stride, ipiv_buf, stride, batch_size, getrf_scratchpad, getrf_batch_scratchpad_size);

    } catch(oneapi::mkl::lapack::exception const& e) {
        // Handle LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to LAPACK API:\nreason: " << e.what() << "\ninfo: " << e.info() << std::endl;
        info = e.info();
    } catch(cl::sycl::exception const& e) {
        // Handle not LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to SYCL API:\n" << e.what() << std::endl;
        info = -1;
    }
    
    std::cout << "getrf_batch " << ((info == 0) ? "ran OK" : "FAILED") << std::endl;
    
    return;
}


//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Batched LU Factorization Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Computes LU Factorization A = P * L * U" << std::endl;
    std::cout << "# for multiple general dense matrices A." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   getrf_batch" << std::endl;
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
//  For each device selected and each data type supported, lu_batch example
//  runs with all supported data types
//
int main(int argc, char **argv) {

    print_example_banner();

    // Find list of devices
    std::list<my_sycl_device_types> listOfDevices;
    set_list_of_devices(listOfDevices);

    for(auto &deviceType : listOfDevices) {
        cl::sycl::device myDev;
        bool myDevIsFound = false;
        get_sycl_device(myDev, myDevIsFound, deviceType);

        if(myDevIsFound) {
            std::cout << std::endl << "Running getrf_batch examples on " << sycl_device_names[deviceType] << "." << std::endl;

            std::cout << "\tRunning with single precision real data type:" << std::endl;
            run_lu_batch_example<float>(myDev);

            std::cout << "\tRunning with single precision complex data type:" << std::endl;
            run_lu_batch_example<std::complex<float>>(myDev);

            if (isDoubleSupported(myDev)) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                run_lu_batch_example<double>(myDev);

                std::cout << "\tRunning with double precision complex data type:" << std::endl;
                run_lu_batch_example<std::complex<double>>(myDev);
            } else {
                std::cout << "\tDouble precision not supported on this device " << std::endl;
                std::cout << std::endl;
            }

        } else {
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
