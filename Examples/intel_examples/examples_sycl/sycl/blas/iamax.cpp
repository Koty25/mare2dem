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
*       This example demonstrates use of DPCPP API oneapi::mkl::blas::iamax to find  
*       the index of the element with maximum absolute value in a vector
*       on a SYCL device (HOST, CPU, GPU).
*
*       result = argmax|x| (real)
*       result = argmax|Rx| + |Ix| (complex)
*
*       The supported floating point data types for iamax vector data are:
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
// Main example for Iamax consisting of 
// initialization x vector.  Then the computation 
//
// result = argmax|x| (real)
// result = argmax|Rx| + |Ix| (complex)
//
// is performed and finally the results are post processed.
//
template <typename fp>
void run_iamax_example(const cl::sycl::device &dev) {

    //
    // Initialize data for Iamax
    //
    // result = argmax|x| (real)
    // result = argmax|Rx| + |Ix| (complex)
    //

    // vector data size
    int n = 1357;
    
    // increment for x
    int incx = 2;

    // prepare vector data
    std::vector <fp, mkl_allocator<fp, 64>> x;

    // prepare result data
    int64_t result = -1;

    rand_vector(x, n, incx);
    
    //
    // Execute Iamax
    //
   
    // Catch asynchronous exceptions
    auto exception_handler = [] (cl::sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(cl::sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception during IAMAX:\n"
                << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
            }
        }
    };
    
    // create execution queue and buffers of matrix data
    cl::sycl::queue main_queue(dev, exception_handler);
 
    cl::sycl::buffer<fp, 1> x_buffer(x.data(), x.size());
    cl::sycl::buffer<int64_t, 1> result_buffer(&result, 1);

    // add oneapi::mkl::blas::iamax to execution queue
    try {
        oneapi::mkl::blas::iamax(main_queue, n, x_buffer, incx, result_buffer);
    }
    catch(cl::sycl::exception const& e) {
        std::cout << "\t\tCaught synchronous SYCL exception during IAMAX:\n"
                  << e.what() << std::endl;
    }

    
    //
    // Post Processing
    //

    std::cout << "\n\t\tIAMAX parameters:\n";
    std::cout << "\t\t\tn = " << n << std::endl;
    std::cout << "\t\t\tincx = " << incx << std::endl;
   
    std::cout << "\n\t\tOutputting 2x1 block of x vector and 1x1 block of result scalar:" << std::endl;

    // output the top 2x1 block of x vector
    auto x_accessor = x_buffer.template get_access<cl::sycl::access::mode::read>();
    print_2x1_vector_values(x_accessor, incx, "x");

    // output the 1x1 block of result scalar
    auto result_accessor = result_buffer.template get_access<cl::sycl::access::mode::read>();
    print_1x1_scalar_values(result_accessor, "result");


}

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Find Index of Maximum Absolute Argument Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# result = argmax|x| (real)" << std::endl;
    std::cout << "# result = argmax|Rx| + |Ix| (complex)" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# where x is general dense vector," << std::endl;
    std::cout << "# and result is the index of maximum absolute argument in vector x." << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   iamax" << std::endl;
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
//  For each device selected and each data type supported, Iamax Example 
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
            
            std::cout << "\tRunning with single precision real data type:" << std::endl;
            run_iamax_example<float>(my_dev);

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                run_iamax_example<double>(my_dev);
            }

            std::cout << "\tRunning with single precision complex data type:" << std::endl;
            run_iamax_example<std::complex<float>>(my_dev);

            if (my_dev.get_info<cl::sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision complex data type:" << std::endl;
                run_iamax_example<std::complex<double>>(my_dev);
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
