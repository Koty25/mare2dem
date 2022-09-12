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
*       This example demonstrates usage of DPC++ device API for RNG
*       continious distributions with oneapi::mkl::rng::philox4x32x10
*       random number generator
*
*       Continious distributions list:
*           oneapi::mkl::rng::gaussian
*           oneapi::mkl::rng::lognormal
*
*       The supported floating point data types for random numbers are:
*           float
*           double
*
*******************************************************************************/

// stl includes
#include <iostream>
#include <vector>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/rng/device.hpp"

// local includes
#include "common_for_rng_examples.hpp"

// example parameters defines
#define SEED    777
#define N       1000
#define N_PRINT 10

const int vec_size = 4;

template<typename RealType>
bool run_gaussian_example(sycl::queue& queue) {
    // prepare array for random numbers
    std::vector<RealType> r(N);

    RealType mean = static_cast<RealType>(2.0);
    RealType stddev = static_cast<RealType>(0.5);

    // submit a kernel to generate on device
    {
        sycl::buffer<RealType, 1> r_buf(r.data(), r.size());

        try {
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(sycl::range<1>(N / vec_size), [=](sycl::item<1> item) {
                    size_t id = item.get_id(0);
                    oneapi::mkl::rng::device::philox4x32x10<vec_size> engine(SEED, id * vec_size);
                    oneapi::mkl::rng::device::gaussian<RealType> distr(mean, stddev);

                    auto res = oneapi::mkl::rng::device::generate(distr, engine);

                    res.store(id, r_acc);
                });
            });
            queue.wait_and_throw();
        }
        catch(sycl::exception const& e) {
            std::cout << "\t\tSYCL exception\n"
                      << e.what() << std::endl;
            return 1;
        }

        std::cout << "\t\tOutput of generator:" << std::endl;

        std::cout << "\n\t\tOutput of generator with gaussian distribution:" << std::endl;
        print_output(r_buf, N_PRINT);
    } // buffer life-time ends

    // validation
    return statistics<RealType, oneapi::mkl::rng::device::gaussian<RealType>>{}
                .check(r, oneapi::mkl::rng::device::gaussian<RealType>{mean, stddev});
}

template<typename RealType>
bool run_lognormal_example(sycl::queue& queue) {
    // prepare array for random numbers
    std::vector<RealType> r(N);

    RealType mean = static_cast<RealType>(1.0);
    RealType stddev = static_cast<RealType>(0.5);

    // submit a kernel to generate on device
    {
        sycl::buffer<RealType, 1> r_buf(r.data(), r.size());

        try {
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(sycl::range<1>(N / vec_size), [=](sycl::item<1> item) {
                    size_t id = item.get_id(0);
                    oneapi::mkl::rng::device::philox4x32x10<vec_size> engine(SEED, id * vec_size);
                    oneapi::mkl::rng::device::lognormal<RealType> distr(mean, stddev);

                    auto res = oneapi::mkl::rng::device::generate(distr, engine);

                    res.store(id, r_acc);
                });
            });
            queue.wait_and_throw();
        }
        catch(sycl::exception const& e) {
            std::cout << "\t\tSYCL exception\n"
                      << e.what() << std::endl;
            return 1;
        }

        std::cout << "\t\tOutput of generator:" << std::endl;

        std::cout << "\n\t\tOutput of generator with lognormal distribution:" << std::endl;
        print_output(r_buf, N_PRINT);
    } // buffer life-time ends

    // validation
    return statistics<RealType, oneapi::mkl::rng::device::lognormal<RealType>>{}
            .check(r, oneapi::mkl::rng::device::lognormal<RealType>{mean, stddev});
}

void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Generate random numbers with continious rng distributions:" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "#   oneapi::mkl::rng::gaussian" << std::endl;
    std::cout << "#   oneapi::mkl::rng::lognormal" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Supported floating point type precisions:" << std::endl;
    std::cout << "#   float double" << std::endl;
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
            std::cout << "\tRunning with single precision real data type:" << std::endl;
            if(run_gaussian_example<float>(queue) ||
                run_lognormal_example<float>(queue)) {
                std::cout << "FAILED" << std::endl;
                return 1;
            }
            if(my_dev.get_info<sycl::info::device::double_fp_config>().size() != 0) {
                std::cout << "\tRunning with double precision real data type:" << std::endl;
                if(run_gaussian_example<double>(queue) ||
                run_lognormal_example<double>(queue)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
                }
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
