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
*       This example demonstrates use of DPCPP API oneapi::mkl::blas::gemm to perform General 
*       Matrix-Matrix Multiplication on a SYCL device (HOST, CPU, GPU).
*
*       C = alpha * op(A) * op(B) + beta * C
*
*       where op() is defined by one of oneapi::mkl::transpose::{nontrans,trans,conjtrans} 
*
*
*       The supported floating point data types for gemm matrix data are:
*           half
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
#include <limits>
#include <vector>
#include <algorithm>
#include <cstring>
#include <list>
#include <iterator>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/blas.hpp"
#include "mkl.h"

// local includes
#include "common_for_examples.hpp"


//
// Main example for Gemm consisting of 
// initialization of A, B and C matrices as well as 
// scalars alpha and beta.  Then the product
//
// C = alpha * op(A) * op(B) + beta * C
//
// is performed and finally the results are post processed.
//
template <typename fp>
void run_gemm_example(const cl::sycl::device &dev) {

    //
    // Initialize data for Gemm
    //
    // C = alpha * op(A) * op(B)  + beta * C 
    //

    oneapi::mkl::transpose transA = oneapi::mkl::transpose::trans;
    oneapi::mkl::transpose transB = oneapi::mkl::transpose::nontrans;

    // matrix data sizes
    int m = 79;
    int n = 83; 
    int k = 91;
    
    // leading dimensions of data
    int ldA = 103;
    int ldB = 105;
    int ldC = 106;

    // set scalar fp values     
    fp alpha = set_fp_value(fp(2.0), fp(-0.5)); 
    fp beta  = set_fp_value(fp(3.0), fp(-1.5));

    // prepare matrix data
    std::vector <fp, mkl_allocator<fp, 64>> A;
    std::vector <fp, mkl_allocator<fp, 64>> B;
    std::vector <fp, mkl_allocator<fp, 64>> C;

    rand_matrix(A, transA, m, k, ldA);
    rand_matrix(B, transB, k, n, ldB);
    rand_matrix(C, oneapi::mkl::transpose::nontrans, m, n, ldC);

    
    //
    // Execute Gemm
    //

    // Catch asynchronous exceptions
    auto exception_handler = [] (cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(cl::sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception during GEMM:\n"
                << e.what() << std::endl;
            }
        }
    };
    
    // create execution queue and buffers of matrix data
    cl::sycl::queue main_queue(dev, exception_handler);

    cl::sycl::buffer<fp, 1> A_buffer(A.data(), A.size());
    cl::sycl::buffer<fp, 1> B_buffer(B.data(), B.size());
    cl::sycl::buffer<fp, 1> C_buffer(C.data(), C.size());

    // add oneapi::mkl::blas::gemm to execution queue
    try {
        oneapi::mkl::blas::gemm(main_queue, transA, transB, m, n, k, alpha, A_buffer, ldA, B_buffer, ldB, beta, C_buffer, ldC);
    }
    catch(cl::sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during GEMM:\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
    }

    
    //
    // Post Processing
    //

    std::cout << "\n\t\tGEMM parameters:\n";
    std::cout << "\t\t\ttransA = " << ( transA == oneapi::mkl::transpose::nontrans ? "nontrans" : ( transA == oneapi::mkl::transpose::trans ? "trans" : "conjtrans")) 
              <<   ", transB = " << ( transB == oneapi::mkl::transpose::nontrans ? "nontrans" : ( transB == oneapi::mkl::transpose::trans ? "trans" : "conjtrans")) << std::endl;
    std::cout << "\t\t\tm = " << m << ", n = " << n << ", k = " << k << std::endl;
    std::cout << "\t\t\tlda = " << ldA << ", ldB = " << ldB << ", ldC = " << ldC << std::endl;
    std::cout << "\t\t\talpha = " << alpha << ", beta = " << beta << std::endl;

   
    std::cout << "\n\t\tOutputting 2x2 block of A,B,C matrices:" << std::endl;

    // output the top 2x2 block of A matrix
    auto A_accessor = A_buffer.template get_access<cl::sycl::access::mode::read>();
    print_2x2_matrix_values(A_accessor, ldA, "A");
    
    // output the top 2x2 block of B matrix
    auto B_accessor = B_buffer.template get_access<cl::sycl::access::mode::read>();
    print_2x2_matrix_values(B_accessor, ldB, "B");

    // output the top 2x2 block of C matrix
    auto C_accessor = C_buffer.template get_access<cl::sycl::access::mode::read>();
    print_2x2_matrix_values(C_accessor, ldC, "C");


}

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# General Matrix-Matrix Multiplication Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# C = alpha * A * B + beta * C" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where A, B and C are general dense matrices and alpha, beta are" << std::endl;
    std::cout << "# floating point type precision scalars." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   gemm" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   half" << std::endl;
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
//  For each device selected and each data type supported, Gemm Example 
//  runs with all supported data types
//
int main (int argc, char ** argv) {
   

    print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    for (auto dev_type : list_of_devices) {

        cl::sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, dev_type);

        if (my_dev_is_found) {
            std::cout << "Running tests on " << sycl_device_names[dev_type] << ".\n";
            
            std::cout << "\tRunning with half precision real data type:" << std::endl;
            run_gemm_example<half>(my_dev);

            std::cout << "\tRunning with single precision real data type:" << std::endl;
            run_gemm_example<float>(my_dev);

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                run_gemm_example<double>(my_dev);
            }

            std::cout << "\tRunning with single precision complex data type:" << std::endl;
            run_gemm_example<std::complex<float>>(my_dev);

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision complex data type:" << std::endl;
                run_gemm_example<std::complex<double>>(my_dev);
            }
        } else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[dev_type] << " devices found; Fail on missing devices is enabled.\n";
                return 1;
#else
            std::cout << "No " << sycl_device_names[dev_type] << " devices found; skipping " << sycl_device_names[dev_type] << " tests.\n";
#endif
        }


    }

    return 0;

}
