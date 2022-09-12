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
*       This example demonstrates use of DPC++ API for calculation of min/max
*       estimates
*
*
*******************************************************************************/

// stl includes
#include <iostream>
#include <vector>

// mkl/sycl includes
#include <CL/sycl.hpp>
#include "oneapi/mkl.hpp"

// local includes
#include "common_for_stats_examples.hpp"

template <typename RealType>
int run_min_example(sycl::queue& queue) {

    const size_t n_observations = 10;
    const size_t n_dims = 3;

    sycl::usm_allocator<RealType, sycl::usm::alloc::shared, 64> allocator(queue);
    std::vector<RealType, decltype(allocator)> x(n_observations * n_dims, allocator);
    random_vector(x, n_observations, n_dims);

    std::vector<RealType, decltype(allocator)> min_est(n_dims, allocator);

    for(int i = 0; i < n_dims; i++) {
        min_est[i] = x[i * n_observations];
    }
    try {
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x.data());
        oneapi::mkl::stats::min(queue, dataset, min_est.data());
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    // validation
    int num_errors = 0;
    for(int i = 0; i < n_dims; i++) {
        for(int j = 0; j < n_observations; j++) {
            if(x[j + i * n_observations] < min_est[i]) {
                num_errors++;
            }
        }
    }
    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Min for dimension " << i << ": " << min_est[i] << std::endl;
    }
    if(num_errors == 0) {
        std::cout << "All observations are within ranges for all dimensions" << std::endl;
    }
    else {
        std::cout << "There are " << num_errors << " observations beyond the ranges" << std::endl;
        return 1;
    }
    return 0;
}

template <typename RealType>
int run_max_example(sycl::queue& queue) {

    const size_t n_observations = 10;
    const size_t n_dims = 3;

    sycl::usm_allocator<RealType, sycl::usm::alloc::shared, 64> allocator(queue);
    std::vector<RealType, decltype(allocator)> x(n_observations * n_dims, allocator);
    random_vector(x, n_observations, n_dims);

    std::vector<RealType, decltype(allocator)> max_est(n_dims, allocator);
    for(int i = 0; i < n_dims; i++) {
        max_est[i] = x[i * n_observations];
    }
    try {
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x.data());
        oneapi::mkl::stats::max(queue, dataset, max_est.data());
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    // validation
    int num_errors = 0;
    for(int i = 0; i < n_dims; i++) {
        for(int j = 0; j < n_observations; j++) {
            if(x[j + i * n_observations] > max_est[i]) {
                num_errors++;
            }
        }
    }
    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Max for dimension " << i << ": " << max_est[i] << std::endl;
    }
    if(num_errors == 0) {
        std::cout << "All observations are within ranges for all dimensions" << std::endl;
    }
    else {
        std::cout << "There are " << num_errors << " observations beyond the ranges" << std::endl;
        return 1;
    }
    return 0;
}

template <typename RealType>
int run_min_max_example(sycl::queue& queue) {

    const size_t n_observations = 10;
    const size_t n_dims = 3;

    sycl::usm_allocator<RealType, sycl::usm::alloc::shared, 64> allocator(queue);
    std::vector<RealType, decltype(allocator)> x(n_observations * n_dims, allocator);
    random_vector(x, n_observations, n_dims);

    std::vector<RealType, decltype(allocator)> min_est(n_dims, allocator);
    std::vector<RealType, decltype(allocator)> max_est(n_dims, allocator);
    for(int i = 0; i < n_dims; i++) {
        min_est[i] = max_est[i] = x[i * n_observations];
    }
    try {
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x.data());
        // oneapi::mkl::stats::min() and oneapi::mkl::stats::max() functions may be called instead of oneapi::mkl::stats::min_max()
        oneapi::mkl::stats::min_max(queue, dataset, min_est.data(), max_est.data());
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    // validation
    int num_errors = 0;
    for(int i = 0; i < n_dims; i++) {
        for(int j = 0; j < n_observations; j++) {
            if(x[j + i * n_observations] < min_est[i] || x[j + i * n_observations] > max_est[i]) {
                num_errors++;
            }
        }
    }
    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Min for dimension " << i << ": " << min_est[i] << std::endl;
        std::cout << "Max for dimension " << i << ": " << max_est[i] << std::endl;
    }
    if(num_errors == 0) {
        std::cout << "All observations are within ranges for all dimensions" << std::endl;
    }
    else {
        std::cout << "There are " << num_errors << " observations beyond the ranges" << std::endl;
        return 1;
    }
    return 0;
}

//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Min max estimates example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "# oneapi::mkl::stats::min oneapi::mkl::stats::max oneapi::mkl::stats::min_max" << std::endl;
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
    auto exception_handler = [] (sycl::exception_list exceptions) {
        for (std::exception_ptr const& e : exceptions) {
            try {
                std::rethrow_exception(e);
            } catch(sycl::exception const& e) {
                std::cout << "Caught asynchronous SYCL exception during calculation:\n"
                << e.what() << std::endl;
            }
        }
    };

    print_example_banner();

    std::list<my_sycl_device_types> list_of_devices;
    set_list_of_devices(list_of_devices);

    for (auto it = list_of_devices.begin(); it != list_of_devices.end(); ++it) {

        sycl::device my_dev;
        bool my_dev_is_found = false;
        get_sycl_device(my_dev, my_dev_is_found, *it);

            if (my_dev_is_found) {
                sycl::queue queue(my_dev, exception_handler);

                std::cout << "Running tests on " << sycl_device_names[*it] << ".\n";

                std::cout << "\tRunning with single precision real data type:" << std::endl;
                if(run_min_example<float>(queue) ||
                    run_max_example<float>(queue) ||
                    run_min_max_example<float>(queue)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
                }
                if(isDoubleSupported(my_dev)) {
                    std::cout << "\tRunning with double precision real data type:" << std::endl;
                    if(run_min_example<double>(queue) ||
                        run_max_example<double>(queue) ||
                        run_min_max_example<double>(queue)) {
                        std::cout << "FAILED" << std::endl;
                        return 1;
                    }
                }
                else {
                    std::cout << "Double precision is not supported for this device" << std::endl;
                }
            }
            else {
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
