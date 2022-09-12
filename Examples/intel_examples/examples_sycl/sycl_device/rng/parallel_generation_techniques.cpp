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
*       This example demonstrates usage of different parallelizations techniques
*       for random number generators on a SYCL device (Host, CPU, GPU).
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
#define N       1024
#define N_PRINT 10

//
// examples show usage of rng device functionality to parallelize calculations
// on SYCL devices
//
template <typename Type>
int run_initial_skip_example(sycl::queue& queue) {
    std::cout << "\tRunning initial_skip example" << std::endl;
    // prepare array for random numbers
    std::vector<Type> r1(N);
    std::vector<Type> r2(N);

    // known amount of number to generate per each thread
    size_t num_per_thread = 32;

    // submit a kernel to generate on device
    {
        sycl::buffer<Type, 1> r1_buf(r1.data(), r1.size());
        sycl::buffer<Type, 1> r2_buf(r2.data(), r2.size());

        try {
            // the first kernel generates 32 scalar random numbers in a row
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r1_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(sycl::range<1>(N / num_per_thread), [=](sycl::item<1> item) {
                    // offset parameter = thread_id * num_per_thread
                    oneapi::mkl::rng::device::philox4x32x10<> engine(SEED, item.get_id(0) * num_per_thread);
                    oneapi::mkl::rng::device::uniform<Type> distr;

                    Type res;

                    for(int i = 0; i < num_per_thread; i++) {
                        res = oneapi::mkl::rng::device::generate(distr, engine);
                        r_acc[item.get_id(0) * num_per_thread + i] = res;
                    }
                });
            });
            // the second kernel generates 32 scalar random numbers, but state is skipped after each number
            // this aproach is less efficient, but it works if unknown amount of numbers required
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r2_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(sycl::range<1>(N / num_per_thread), [=](sycl::item<1> item) {
                    // offset parameter = thread_id
                    oneapi::mkl::rng::device::philox4x32x10<> engine(SEED, item.get_id(0));
                    oneapi::mkl::rng::device::uniform<Type> distr;

                    Type res;

                    for(int i = 0; i < num_per_thread; i++) {
                        res = oneapi::mkl::rng::device::generate(distr, engine);
                        // each engine generated a single number, need to skip on range - 1
                        oneapi::mkl::rng::device::skip_ahead(engine, N / num_per_thread - 1);
                        r_acc[item.get_id(0) + i * (N / num_per_thread)] = res;
                    }
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

        std::cout << "\n\t\tOutput of generator with uniform distribution:" << std::endl;
        print_output(r1_buf, N_PRINT);
    } // buffer life-time ends

    // compare results from two approaches
    int err = 0;
    for(int i = 0; i < N; i++) {
        if(r1[i] != r2[i]) {
            std::cout << "error in " << i << " element " << r1[i] << " " << r2[i] << std::endl;
            err++;
        }
    }
    return err;
}

template <typename Type>
int run_subsequence_example(sycl::queue& queue) {
    std::cout << "\tRunning subsequence example" << std::endl;
    // prepare array for random numbers
    std::vector<Type> r_dev(N);

    // known amount of threads
    size_t n_threads = 32;

    // submit a kernel to generate on device
    {
        sycl::buffer<Type, 1> r_buf(r_dev.data(), r_dev.size());

        try {
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(sycl::range<1>(n_threads), [=](sycl::item<1> item) {
                    size_t id = item.get_id(0);
                    // each engine skip sequence on 2^64 * id
                    // for this case each thread may produce up to 2^64 independent random numbers
                    oneapi::mkl::rng::device::philox4x32x10<> engine(SEED, {0, id});

                    oneapi::mkl::rng::device::gaussian<Type> distr;
                    Type res;

                    for(int i = 0; i < N / n_threads; i++) {
                        res = oneapi::mkl::rng::device::generate(distr, engine);
                        r_acc[i + id * (N / n_threads)] = res;
                    }
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

    return statistics<Type, oneapi::mkl::rng::device::gaussian<Type>>{}.check(r_dev, oneapi::mkl::rng::device::gaussian<Type>{});
}

//
// description of example setup, APIs used
//
void print_example_banner() {

    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Generate uniformly distributed random numbers with philox4x32x10\n# generator example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using APIs:" << std::endl;
    std::cout << "#   philox4x32x10 uniform gaussian skip_ahead" << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << std::endl;

}

//
// main entry point for example.
//
// Dispatches to appropriate device types as set at build time with flag:
// -DSYCL_DEVICES_host -- only runs host implementation
// -DSYCL_DEVICES_cpu -- only runs SYCL CPU implementation
// -DSYCL_DEVICES_gpu -- only runs SYCL GPU implementation
// -DSYCL_DEVICES_all (default) -- runs on all: host, cpu and gpu devices
//

int main (int argc, char ** argv) {

    print_example_banner();

    // handler to catch asynchronous exceptions
    auto exception_handler = [] (sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception:\n"
                << e.what() << std::endl;
            }
        }
    };

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    for (auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {

        sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

            if (my_dev_is_found) {
                std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";
    
                sycl::queue queue(my_dev, exception_handler);

                std::cout << "\n\tRunning with single precision real data type:" << std::endl;
                if(run_initial_skip_example<float>(queue) ||
                   run_subsequence_example<float>(queue)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
                }
                if(isDoubleSupported(my_dev)) {
                    std::cout << "\n\tRunning with double precision real data type:" << std::endl;
                    if(run_initial_skip_example<double>(queue) ||
                       run_subsequence_example<double>(queue)) {
                        std::cout << "FAILED" << std::endl;
                        return 1;
                    }
                }
                else {
                    std::cout << "Double precision is not supported for this device" << std::endl;
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
