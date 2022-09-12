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
*            mkl::vm::exp example program text (SYCL interface)
*
*******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <limits>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstring>
#include <list>
#include <iterator>
#include <cstdint>
#include <random>
#include <CL/sycl.hpp>
#include "mkl_sycl_types.hpp"
#include "mkl_vml_sycl.hpp"
#include "mkl.h"
#include "common_for_examples.hpp"

#define VLEN 64*1024
#define BEG1 2.0
#define END1 +10.0

#include "vml_common.hpp"

#define MAX_STR 128

using namespace cl::sycl;
using std::vector;

//!
//! @brief Accuracy test
//!
template < typename A, typename R > bool vExpAccuracyLiteTest (const device & dev)
{

    // Catch asynchronous exceptions
    auto exception_handler = [] (cl::sycl::exception_list exceptions)
    {
        for (std::exception_ptr const& e : exceptions)
        {
            try
            {
                std::rethrow_exception(e);
            }
            catch(cl::sycl::exception const& e)
            {
                std::cout << "Caught asynchronous SYCL exception:\n" << e.what() << std::endl;
            }
        } // for (std::exception_ptr const& e : exceptions) 
    };

    // Create execution queue with asynchronous error handling
    cl::sycl::queue main_queue (dev, exception_handler);

    // Get device name
    std::string dev_name = main_queue.get_device ().get_info < cl::sycl::info::device::name > ();

    // *************************************************************
    // Variable declaraions
    // *************************************************************
    A  *varg = own_malloc < A > (main_queue, VLEN * sizeof (A));

    R  *vres11 = own_malloc < R > (main_queue, VLEN * sizeof (R));
    R  *vres12 = own_malloc < R > (main_queue, VLEN * sizeof (R));

    R  *vres21 = own_malloc < R > (main_queue, VLEN * sizeof (R));
    R  *vres22 = own_malloc < R > (main_queue, VLEN * sizeof (R));

    R  *vres = own_malloc < R > (main_queue, VLEN * sizeof (R));

    // Number of errors
    int errs = 0;
    // Number of printed errors
    int printed_errs = 0;

    std::random_device rdevice;
    std::mt19937 engine(rdevice());
    std::uniform_real_distribution<> distr(BEG1, END1);

    // *************************************************************
    // Vector input data initialization
    // *************************************************************
    for (size_t i = 0; i < VLEN; ++i)
    {
        
        varg[i] = distr(engine);
        
        vres11[i] = { 0 };
        vres12[i] = { 0 };
        vres21[i] = { 0 };
        vres22[i] = { 0 };
        vres[i]   = { 0 };
    }   // for (size_t i = 0; i < ACCURACY_LEN; ++i)

    mkl::vm::ln  (main_queue, VLEN, varg, vres11, nullptr, accuracy_mode[0]);
    mkl::vm::exp (main_queue, VLEN, varg, vres12, nullptr, accuracy_mode[0]);

    mkl::vm::exp (main_queue, VLEN, vres11, vres21, nullptr, accuracy_mode[0]);
    mkl::vm::ln (main_queue, VLEN,  vres12, vres22, nullptr, accuracy_mode[0]);

    mkl::vm::sub (main_queue, VLEN, vres21, vres22, vres, nullptr, accuracy_mode[0]);

    // Catch sycl exceptions
    try
    {
        main_queue.wait_and_throw ();
    }
    catch (cl::sycl::exception & e)
    {
        std::cerr << "SYCL exception during Accuracy Test\n" << e.what () << std::endl << "OpenCl status: " << e.get_cl_code () << std::endl;
        return false;
    }

    // *************************************************************
    // Compute ulp between computed and expected (reference)
    // values and check
    // *************************************************************
    for (size_t i = 0; i < VLEN; ++i)
    {
        // Check
        if (std::fabs(vres[i]) > 0.01 )
        {
            std::cout << "\t\t" << varg[i] << " -> " << vres[i] << std::endl;
            errs++;
        }
    }   // for (size_t i = 0; i < varg1.size (); ++i)

    std::cout << "\tResult: " << ((errs == 0) ? "PASS" : "FAIL") << std::endl;

    own_free < A > (main_queue, varg);

    own_free < R > (main_queue, vres11);
    own_free < R > (main_queue, vres12);

    own_free < R > (main_queue, vres21);
    own_free < R > (main_queue, vres22);

    own_free < R > (main_queue, vres);

    return (errs == 0);
}   // template <typename A, typename R> bool vFuncAccuracyText (const device &dev)

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner ()
{
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# General VM " << "exp" << " Function Example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   vm::" << "exp" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;

    std::cout << "#   float" << std::endl;

    std::cout << "#   double" << std::endl;

    std::cout << "#   std::complex<float>" << std::endl;

    std::cout << "#   std::complex<double>" << std::endl;

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
//  For each device selected and each data type supported, the example 
//  runs with all supported data types
//
int main (int argc, char **argv)
{
    // Number of errors occured
    int errs = 0;

    // Print standard banner for VM examples
    print_example_banner ();

    // List of available devices
    std::list < my_sycl_device_types > list_of_devices;
    set_list_of_devices (list_of_devices);

    // Loop by all available devices
  for (auto dev_type:list_of_devices)
    {
        cl::sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device (my_dev, my_dev_is_found, dev_type);


        // Run tests if the device is available
        if (my_dev_is_found)
        {
            std::cout << "Running tests on " << sycl_device_names[dev_type] << ".\n";

            std::cout << "\tRunning with single precision real data type:" << std::endl;
            errs += vExpAccuracyLiteTest < float, float >(my_dev);

            std::cout << "\tRunning with double precision real data type:" << std::endl;
            errs += vExpAccuracyLiteTest < double, double >(my_dev);
        }
        else
        {

            std::cout << "No " << sycl_device_names[dev_type] << " devices found; skipping " << sycl_device_names[dev_type] << " tests.\n";

        }
    }

    return 0;
}
