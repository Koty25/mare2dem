/*******************************************************************************
* Copyright 2020 Intel Corporation.
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
*       This example demonstrates use of oneapi::mkl::lapack::getrf and oneapi::mkl::lapack::getri
*       to perform LU factorization and compute inverse on a SYCL device (HOST, CPU, GPU).
*
*       The supported floating point data types for matrix data are:
*           float
*           double
*           std::complex<float>
*           std::complex<double>
*
*******************************************************************************/

// STL includes
#include <iostream>
#include <complex>
#include <vector>

// MKL/SYCL includes
#include "mkl.h"
#include "oneapi/mkl/lapack.hpp"

// local includes
#include "common_for_examples.hpp"

//
// Main example for LU consisting of initialization of
// a general dense A matrix.
// Then the LU factorization
// A = P * L * U
// is performed followed by computing the inverse
// inv(A) using the computed LU factorization.
// Finally the results are post processed.
//

template <typename data_t>
void run_getri_example(const cl::sycl::device& device)
{
    // Matrix sizes and leading dimensions
    std::int64_t n    = 23;
    std::int64_t lda  = 32;

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
    cl::sycl::queue queue(device, error_handler);
    cl::sycl::context context = queue.get_context();

    // Allocate matrices
    std::int64_t A_size = n * lda;
    std::int64_t ipiv_size = n; 
    std::vector<data_t> A(A_size);
    std::vector<int64_t> ipiv(ipiv_size);

    // Initialize matrix A
    rand_matrix(A, oneapi::mkl::transpose::nontrans, n, n, lda);

    try {
        // Get sizes of scratchpads for calculations
        std::int64_t getrf_scratchpad_size = oneapi::mkl::lapack::getrf_scratchpad_size<data_t>(queue, n, n, lda);
        std::int64_t getri_scratchpad_size = oneapi::mkl::lapack::getri_scratchpad_size<data_t>(queue, n, lda);

        // Define buffers for the data
        cl::sycl::buffer<data_t> A_buffer{A.data(), A_size};
        cl::sycl::buffer<int64_t> ipiv_buffer{ipiv.data(), ipiv_size};
        cl::sycl::buffer<data_t> getrf_scratchpad{getrf_scratchpad_size};
        cl::sycl::buffer<data_t> getri_scratchpad{getri_scratchpad_size};

        // Perform factorization
        oneapi::mkl::lapack::getrf(queue, n, n, A_buffer, lda, ipiv_buffer, getrf_scratchpad, getrf_scratchpad_size);
        oneapi::mkl::lapack::getri(queue, n, A_buffer, lda, ipiv_buffer, getri_scratchpad, getri_scratchpad_size);

        // Wait until calculations are done
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

    std::cout << "getri " << ((info == 0) ? "ran OK" : "FAILED") << std::endl;

    return;
}


//
// Description of example setup, apis used and supported floating point type precisions
//

void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# LU Factorization and Inverse Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Computes LU Factorization A = P * L * U" << std::endl;
    std::cout << "# and uses it to invert A: inv(A)" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a general dense matrix." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   getrf and getri" << std::endl;
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
//  For each device selected and each data type supported, matrix inversion example
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
          std::cout << std::endl << "Running getri examples on " << sycl_device_names[deviceType] << "." << std::endl;

          std::cout << "\tRunning with single precision real data type:" << std::endl;
          run_getri_example<float>(myDev);

          std::cout << "\tRunning with single precision complex data type:" << std::endl;
          run_getri_example<std::complex<float>>(myDev);

          if (isDoubleSupported(myDev)) {
              std::cout << "\tRunning with double precision real data type:" << std::endl;
              run_getri_example<double>(myDev);

              std::cout << "\tRunning with double precision complex data type:" << std::endl;
              run_getri_example<std::complex<double>>(myDev);
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
