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
*       This example demonstrates use of oneapi::mkl::lapack::potrf and oneapi::mkl::lapack::potrs
*       to perform Cholesky factorization and solve on a SYCL device (HOST, CPU, GPU).
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
// Example for Cholesky factorization consisting of initialization of
// symmetric positive definite matrix A and right-hand side B matrix.
// Then the Cholesky factorization A = L * L^T is performed by calling potrf,
// followed by solving L * L^T * X = B using potrs with computed factorization.
//
template <typename data_t>
void run_potrs_example(const cl::sycl::device& device)
{
    // Matrix sizes and leading dimensions
    std::int64_t n    = 32;
    std::int64_t lda  = 32;
    std::int64_t nrhs = 32;
    std::int64_t ldb  = 32;
    oneapi::mkl::uplo    uplo = oneapi::mkl::uplo::lower;

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

    // Allocators to use
    using allocator_t = cl::sycl::usm_allocator<data_t, cl::sycl::usm::alloc::shared>;
    allocator_t allocator(context, device);

    // Allocate matrices
    std::int64_t A_size = n * lda;
    std::int64_t B_size = nrhs * ldb;
    std::vector<data_t, allocator_t> A(A_size, allocator);
    std::vector<data_t, allocator_t> B(B_size, allocator);

    // Initialize positive definite matrix A
    for (int j=0; j<n; j++) {
        for (int i=0; i<n; i++) {
            // Initialize matrix value with 2*i+1 on diagonal and (i+j)/n for off-diagonal
            A[j*lda+i] = (i == j) ? 2*i+1 : static_cast<data_t>(i+j)/static_cast<decltype(abs(data_t()))>(n);
        }
    }

    // Initialize right hand side B of the system A * X = B
    for (int j=0; j<nrhs; j++) {
        for (int i=0; i<n; i++) {
            B[j*ldb+i] = j;
        }
    }

    try {
        // Allocate scratchpads for calculations
        std::int64_t potrf_scratchpad_size = oneapi::mkl::lapack::potrf_scratchpad_size<data_t>(queue, uplo, n, lda);
        std::vector<data_t, allocator_t> potrf_scratchpad(potrf_scratchpad_size, allocator);
        std::int64_t potrs_scratchpad_size = oneapi::mkl::lapack::potrs_scratchpad_size<data_t>(queue, uplo, n, nrhs, lda, ldb);
        std::vector<data_t, allocator_t> potrs_scratchpad(potrs_scratchpad_size, allocator);

        // Perform factorization
        auto potrf_done_event = oneapi::mkl::lapack::potrf(queue, uplo, n, A.data(), lda, potrf_scratchpad.data(), potrf_scratchpad_size);

        // Perform solve operation
        auto potrs_done_event = oneapi::mkl::lapack::potrs(queue, uplo, n, nrhs, A.data(), lda,  B.data(), ldb, potrs_scratchpad.data(), potrs_scratchpad_size, { potrf_done_event });

        // Wait until calculations are done
        potrs_done_event.wait_and_throw();
    } catch(oneapi::mkl::lapack::exception const& e) {
        // Handle LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to LAPACK API:\nreason: " << e.what() << "\ninfo: " << e.info() << std::endl;
        info = e.info();
    } catch(cl::sycl::exception const& e) {
        // Handle not LAPACK related exceptions happened during synchronous call
        std::cout << "Unexpected exception caught during synchronous call to SYCL API:\n" << e.what() << std::endl;
        info = -1;
    }

    std::cout << "potrs " << ((info == 0) ? "ran OK" : "FAILED") << std::endl;

    return;
}


//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Cholesky Factorization and Solve Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Computes Cholesky Factorization A = U^T * U" << std::endl;
    std::cout << "# and uses it to solve A * X = B" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A is a symmetric positive definite matrix and " << std::endl;
    std::cout << "# B contains the right hand sides." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   potrf and potrs" << std::endl;
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
//  For each device selected and each data type supported, Cholesky Example
//  runs with all supported data types
//
int main(int argc, char **argv) {

    print_example_banner();

    // Find list of devices
    std::list<my_sycl_device_types> listOfDevices;
    set_list_of_devices(listOfDevices);

    for (auto &deviceType: listOfDevices) {
        cl::sycl::device myDev;
        bool myDevIsFound = false;
        get_sycl_device(myDev, myDevIsFound, deviceType);

        if (myDevIsFound) {
          std::cout << std::endl << "Running potrs examples on " << sycl_device_names[deviceType] << "." << std::endl;

          std::cout << "Running with single precision real data type:" << std::endl;
          run_potrs_example<float>(myDev);

          std::cout << "Running with single precision complex data type:" << std::endl;
          run_potrs_example<std::complex<float>>(myDev);

          if (isDoubleSupported(myDev)) {
              std::cout << "Running with double precision real data type:" << std::endl;
              run_potrs_example<double>(myDev);

              std::cout << "Running with double precision complex data type:" << std::endl;
              run_potrs_example<std::complex<double>>(myDev);
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
