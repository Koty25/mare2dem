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
*       This example demonstrates usage of DPC++ API for MKL RNG service
*       functionality
*
*       Functions list:
*           oneapi::mkl::rng::skip_ahead
*           oneapi::mkl::rng::leapfrog
*
*
*******************************************************************************/

// stl includes
#include <iostream>
#include <vector>
#include <math.h>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl.hpp"

// local includes
#include "common_for_rng_examples.hpp"

// example parameters defines
#define N       100
#define S       10
#define NS      10
#define N_PRINT 100

template <typename RealType>
bool run_skip_ahead_example(sycl::queue& queue) {
    std::cout << "\n\tRun skip_ahead example" << std::endl;

    oneapi::mkl::rng::philox4x32x10 engine(queue);

    std::vector<oneapi::mkl::rng::philox4x32x10> engine_vec;

    for(int i = 0; i < S; i++) {
        // copy reference engine to engine_vec[i]
        engine_vec.push_back(oneapi::mkl::rng::philox4x32x10{engine});
        // skip ahead engine
        oneapi::mkl::rng::skip_ahead(engine_vec[i], i * NS);
    }

    // prepare array for random numbers
    std::vector<RealType, mkl_allocator<RealType, 64>> r_ref(N);
    std::vector<RealType, mkl_allocator<RealType, 64>> r(N);

    {
        sycl::buffer<RealType, 1> r_ref_buffer(r_ref.data(), r_ref.size());

        try {
            // fill r_ref with N random numbers
            oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine, N, r_ref_buffer);

            // fill r with random numbers by portions of NS
            for(int i = 0; i < S; i++) {
                sycl::buffer<RealType, 1> r_buffer(r.data() + i * NS, NS);
                oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine_vec[i], NS, r_buffer);
            }
            queue.wait_and_throw();
        }
        catch(sycl::exception const& e) {
            std::cout << "\t\tSYCL exception during generation\n"
                      << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
            return false;
        }
        sycl::buffer<RealType, 1> r_buffer(r.data(), r.size());
        std::cout << "\t\tOutput:" << std::endl;
        print_output(r_buffer, N_PRINT);
        std::cout << "\t\tReference output:" << std::endl;
        print_output(r_ref_buffer, N_PRINT);
    }

    // validation
    for(int i = 0; i < N; i++) {
        if(r[i] != r_ref[i]) {
            std::cout << "Fail at " << i << " element" << std::endl;
            return false;
        }
    }
    std::cout << "Success" << std::endl;
    return true;
}

template <typename RealType>
bool run_skip_ahead_ex_example(sycl::queue& queue) {
    std::cout << "\n\tRun skip_ahead extended example" << std::endl;

    oneapi::mkl::rng::mrg32k3a engine_1(queue);

    oneapi::mkl::rng::mrg32k3a engine_2(engine_1);

    // to skip 2^76 elements in the random engine with skip_ahead function should be called 2^14 times
    //    with nskip equal to 2^62
    std::uint64_t nskip = (std::uint64_t) pow(2,62);
    std::uint64_t skip_times = (std::uint64_t) pow(2,14);

    for(std::uint64_t i = 0; i < skip_times; i++) {
        oneapi::mkl::rng::skip_ahead(engine_1, nskip);
    }
    // skip 2^76 elements in the engine with advanced skip_ahead function should be called
    //    with nskip represented as
    //        nskip = 2^76 = 0 + 2^12 * 2^64
    //    in general case:
    //        nskip = params[0] + params[1] * 2^64 + params[2] * 2^128 + ...
    oneapi::mkl::rng::skip_ahead(engine_2, {0, (std::uint64_t) pow(2,12)});

    // prepare array for random numbers
    std::vector<RealType, mkl_allocator<RealType, 64>> r_ref(N);
    std::vector<RealType, mkl_allocator<RealType, 64>> r(N);

    {
        sycl::buffer<RealType, 1> r_ref_buffer(r_ref.data(), r_ref.size());
        sycl::buffer<RealType, 1> r_buffer(r.data(), r.size());

        try {
            oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine_1, N, r_ref_buffer);
            oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine_2, N, r_buffer);
            queue.wait_and_throw();
        }
        catch(sycl::exception const& e) {
            std::cout << "\t\tSYCL exception during generation\n"
                      << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
            return false;
        }
        std::cout << "\t\tOutput:" << std::endl;
        print_output(r_buffer, N_PRINT);
        std::cout << "\t\tReference output:" << std::endl;
        print_output(r_ref_buffer, N_PRINT);
    }

    // validation
    for(int i = 0; i < N; i++) {
        if(r[i] != r_ref[i]) {
            std::cout << "Fail at " << i << " element" << std::endl;
            return false;
        }
    }
    std::cout << "Success" << std::endl;
    return true;
}

template <typename RealType>
bool run_leapfrog_example(sycl::queue& queue) {
    std::cout << "\n\tRun leapfrog example:" << std::endl;

    oneapi::mkl::rng::mcg31m1 engine(queue);

    std::vector<oneapi::mkl::rng::mcg31m1> engine_vec;

    for(int i = 0; i < S; i++) {
        // copy reference engine to engine_vec[i]
        engine_vec.push_back(oneapi::mkl::rng::mcg31m1{engine});
        // skip ahead engine
        oneapi::mkl::rng::leapfrog(engine_vec[i], i, S);
    }

    // prepare array for random numbers
    std::vector<RealType, mkl_allocator<RealType, 64>> r_ref(N);
    std::vector<RealType, mkl_allocator<RealType, 64>> r(N);

    {
        sycl::buffer<RealType, 1> r_ref_buffer(r_ref.data(), r_ref.size());

        try {
            // fill r_ref with N random numbers
            oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine, N, r_ref_buffer);

            // fill r with random numbers by portions of NS
            for(int i = 0; i < S; i++) {
                sycl::buffer<RealType, 1> r_buffer(r.data() + i * NS, NS);
                oneapi::mkl::rng::generate(oneapi::mkl::rng::uniform<RealType>{}, engine_vec[i], NS, r_buffer);
            }
            queue.wait_and_throw();
        }
        catch(sycl::exception const& e) {
            std::cout << "\t\tSYCL exception during generation\n"
                      << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
            return false;
        }
        sycl::buffer<RealType, 1> r_buffer(r.data(), r.size());
        std::cout << "\t\tOutput:" << std::endl;
        print_output(r_buffer, N_PRINT);
        std::cout << "\t\tReference output:" << std::endl;
        print_output(r_ref_buffer, N_PRINT);
    } 

    // validation
    int j = 0;
    for(int i = 0; i < NS; i++) {
        for(int k =0; k < NS; k++) {
            if(r[j++] != r_ref[k * NS + i]) {
                std::cout << "Fail at " << i << " element" << std::endl;
                return false;
            }
        }
    }
    std::cout << "Success" << std::endl;
    return true;
}

void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Demonstrate service functionality usage for random number generators:" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   oneapi::mkl::rng::skip_ahead" << std::endl;
    std::cout << "#   oneapi::mkl::rng::leapfrog" << std::endl;
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

int main (int argc, char ** argv) {

    // catch asynchronous exceptions
    auto exception_handler = [] (sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception during generation:\n"
                << e.what() << std::endl;
            }
        }
    };

    print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    for(auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {

        sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

        if(my_dev_is_found) {
            std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";
            sycl::queue queue(my_dev, exception_handler);

            if(!run_skip_ahead_example<float>(queue) ||
                !run_skip_ahead_ex_example<float>(queue) ||
                !run_leapfrog_example<float>(queue)) {
                std::cout << "FAILED" << std::endl;
                return 1;
            }

        } else {
#ifdef FAIL_ON_MISSING_DEVICES
            std::cout << "No " << sycl_device_names[*it] << " devices found; Fail on missing devices is enabled.\n";
            std::cout << "FAILED" << std::endl;
            return 1;
#else
            std::cout << "No " << sycl_device_names[*it] << " devices found; skipping " << sycl_device_names[*it] << " tests.\n";
#endif
        }
    }
    std::cout << "PASSED" << std::endl;
    return 0;
}
