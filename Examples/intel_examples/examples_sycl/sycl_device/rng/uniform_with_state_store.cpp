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
*       This example demonstrates usage of oneapi::mkl::rng::device::philox4x32x10
*       random number generator to produce random numbers on
*       a SYCL device (Host, CPU, GPU) with states, stored between kernels
*       via sycl::buffer or USM.
*
*******************************************************************************/

// stl includes
#include <iostream>
#include <vector>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl/rng/device.hpp"

// local includes
#include "common_for_examples.hpp"

// example parameters defines
#define SEED    777
#define N       1024
#define N_PRINT 10

//
// examples show usage of rng device functionality, which can be called from both
// host and device sides with state stored between kernels with buffer and USM api
//
template <typename Type>
int run_buffer_example(sycl::queue& queue) {
    std::cout << "\tRunning buffer example" << std::endl;
    // prepare array for random numbers
    std::vector<Type> r_dev(N);
    std::vector<Type> r_host(N);

    const int vec_size = 4;

    sycl::range<1> range(N / vec_size);

    // submit kernels to init engines and generate on device
    {
        sycl::buffer<Type, 1> r_buf(r_dev.data(), r_dev.size());
        sycl::buffer<oneapi::mkl::rng::device::philox4x32x10<vec_size>, 1> engine_buf(range);

        try {
            // kernel with initialization of rng engines
            queue.submit([&](sycl::handler& cgh) {
                // create accessor to sycl::buffer with engines to write initialized states
                auto engine_acc = engine_buf.template get_access<sycl::access::mode::write>(cgh);
                cgh.parallel_for(range, [=](sycl::item<1> item) {
                    size_t id = item.get_id(0);
                    oneapi::mkl::rng::device::philox4x32x10<vec_size> engine(SEED, id * vec_size);
                    engine_acc[id] = engine;
                });
            });
            // generate random numbers
            queue.submit([&](sycl::handler& cgh) {
                auto r_acc = r_buf.template get_access<sycl::access::mode::write>(cgh);
                // create accessor to sycl::buffer with engines to read initialized states
                auto engine_acc = engine_buf.template get_access<sycl::access::mode::read>(cgh);
                cgh.parallel_for(range, [=](sycl::item<1> item) {
                    size_t id = item.get_id(0);

                    auto engine = engine_acc[id];
                    oneapi::mkl::rng::device::uniform<Type> distr;
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

        auto r_acc = r_buf.template get_access<sycl::access::mode::read>();
        std::cout << "first "<< N_PRINT << " numbers of " << N << ": " << std::endl;
        for(int i = 0 ; i < N_PRINT; i++) {
            std::cout << r_acc[i] << " ";
        }
        std::cout << std::endl;
    } // buffer life-time ends

    // compare results with host-side generation
    oneapi::mkl::rng::device::philox4x32x10<> engine(SEED);
    oneapi::mkl::rng::device::uniform<Type> distr;

    int err = 0;
    for(int i = 0; i < N; i++) {
        r_host[i] = oneapi::mkl::rng::device::generate(distr, engine);
        if(r_host[i] != r_dev[i]) {
            std::cout << "error in " << i << " element " << r_host[i] << " " << r_dev[i] << std::endl;
            err++;
        }
    }
    return err;
}

template <typename Type>
int run_usm_example(sycl::queue& queue) {
    std::cout << "\tRunning USM example" << std::endl;
    // prepare array for random numbers
    std::vector<Type> r_host(N);
    Type* r_dev = (Type*)sycl::malloc_shared(sizeof(Type) * N, queue);

    const int vec_size = 4;

    sycl::range<1> range(N / vec_size);

    // create array with engines, allocated on the device
    oneapi::mkl::rng::device::philox4x32x10<vec_size>* engine_ptr = (oneapi::mkl::rng::device::philox4x32x10<vec_size>*)
            sycl::malloc_device(sizeof(oneapi::mkl::rng::device::philox4x32x10<vec_size>) * range.get(0), queue);

    // submit kernels to init engines and generate on device
    try {
        // kernel with initialization of rng engines
        auto event = queue.submit([&](sycl::handler& cgh) {
            cgh.parallel_for(range, [=](sycl::item<1> item) {
                size_t id = item.get_id(0);
                oneapi::mkl::rng::device::philox4x32x10<vec_size> engine(SEED, id * vec_size);
                engine_ptr[id] = engine;
            });
        });
        event.wait();
        // generate random numbers
        event = queue.submit([&](sycl::handler& cgh) {
            cgh.parallel_for(range, [=](sycl::item<1> item) {
                size_t id = item.get_id(0);
                oneapi::mkl::rng::device::uniform<Type> distr;
                auto res = oneapi::mkl::rng::device::generate(distr, engine_ptr[id]);

                for(int i = 0; i < vec_size; i++) {
                    r_dev[i + id * vec_size] = res[i];
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

    std::cout << "first "<< N_PRINT << " numbers of " << N << ": " << std::endl;
    for(int i = 0 ; i < N_PRINT; i++) {
        std::cout << r_dev[i] << " ";
    }
    std::cout << std::endl;

    // compare results with host-side generation
    oneapi::mkl::rng::device::philox4x32x10<> engine(SEED);
    oneapi::mkl::rng::device::uniform<Type> distr;

    int err = 0;
    for(int i = 0; i < N; i++) {
        r_host[i] = oneapi::mkl::rng::device::generate(distr, engine);
        if(r_host[i] != r_dev[i]) {
            std::cout << "error in " << i << " element " << r_host[i] << " " << r_dev[i] << std::endl;
            err++;
        }
    }
    sycl::free(engine_ptr, queue);
    return err;
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
    std::cout << "#   philox4x32x10 uniform" << std::endl;
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
                if(run_buffer_example<float>(queue) ||
                   run_usm_example<float>(queue)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
                }
                if(isDoubleSupported(my_dev)) {
                    std::cout << "\n\tRunning with double precision real data type:" << std::endl;
                    if(run_buffer_example<double>(queue) ||
                       run_usm_example<double>(queue)) {
                        std::cout << "FAILED" << std::endl;
                        return 1;
                    }
                }
                else {
                    std::cout << "Double precision is not supported for this device" << std::endl;
                }
                std::cout << "\n\tRunning with integer data type:" << std::endl;
                if(run_buffer_example<std::int32_t>(queue) ||
                   run_usm_example<std::int32_t>(queue)) {
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
