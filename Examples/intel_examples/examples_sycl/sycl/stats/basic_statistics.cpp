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
*       This example demonstrates use of DPC++ API for calculation of basic estimates:
*       raw and central sums
*       raw and central moments
*       mean, skewness, kurtosis, variation
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
int run_sum_example(sycl::queue& queue) {
    std::cout << "\nRunning raw and central sums example" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> sum(n_dims);
    std::vector<RealType> raw_sum_2(n_dims);
    std::vector<RealType> raw_sum_3(n_dims);
    std::vector<RealType> central_sum_2(n_dims);
    std::vector<RealType> central_sum_3(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> sum_buf(sum.data(), sum.size());
        sycl::buffer<RealType, 1> raw_sum_2_buf(raw_sum_2.data(), raw_sum_2.size());
        sycl::buffer<RealType, 1> raw_sum_3_buf(raw_sum_3.data(), raw_sum_3.size());
        sycl::buffer<RealType, 1> central_sum_2_buf(central_sum_2.data(), central_sum_2.size());
        sycl::buffer<RealType, 1> central_sum_3_buf(central_sum_3.data(), central_sum_3.size());
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::raw_sum(queue, dataset, sum_buf, raw_sum_2_buf, raw_sum_3_buf);
        oneapi::mkl::stats::central_sum(queue, dataset, central_sum_2_buf, central_sum_3_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Sum for dimension " << i << ": " << sum[i] << std::endl;
        std::cout << "raw sum 2 for dimension " << i << ": " << raw_sum_2[i] << std::endl;
        std::cout << "raw sum 3 for dimension " << i << ": " << raw_sum_3[i] << std::endl;
        std::cout << "central sum 2 for dimension " << i << ": " << central_sum_2[i] << std::endl;
        std::cout << "central sum 3 for dimension " << i << ": " << central_sum_3[i] << std::endl;
    }
    return 0;
}

template <typename RealType>
int run_central_sum_example_with_provided_mean(sycl::queue& queue) {
    std::cout << "\nRunning central sums example with provided mean" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> mean(n_dims);
    std::vector<RealType> central_sum_2(n_dims);
    std::vector<RealType> central_sum_3(n_dims);
    std::vector<RealType> central_sum_4(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> mean_buf(mean.data(), mean.size());
        sycl::buffer<RealType, 1> central_sum_2_buf(central_sum_2.data(), central_sum_2.size());
        sycl::buffer<RealType, 1> central_sum_3_buf(central_sum_3.data(), central_sum_3.size());
        sycl::buffer<RealType, 1> central_sum_4_buf(central_sum_4.data(), central_sum_4.size());
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::mean(queue, dataset, mean_buf);
        queue.wait_and_throw();

        oneapi::mkl::stats::central_sum(queue, mean_buf, dataset, central_sum_2_buf, central_sum_3_buf, central_sum_4_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Mean for dimension " << i << ": " << mean[i] << std::endl;
        std::cout << "central sum 2 for dimension " << i << ": " << central_sum_2[i] << std::endl;
        std::cout << "central sum 3 for dimension " << i << ": " << central_sum_3[i] << std::endl;
        std::cout << "central sum 4 for dimension " << i << ": " << central_sum_4[i] << std::endl;
    }
    return 0;
}

template <typename RealType>
int run_moment_example(sycl::queue& queue) {
    std::cout << "\nRunning raw and central moments example" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> mean(n_dims);
    std::vector<RealType> raw_moment_2(n_dims);
    std::vector<RealType> raw_moment_3(n_dims);
    std::vector<RealType> central_moment_2(n_dims);
    std::vector<RealType> central_moment_3(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> mean_buf(mean.data(), mean.size());
        sycl::buffer<RealType, 1> raw_moment_2_buf(raw_moment_2.data(), raw_moment_2.size());
        sycl::buffer<RealType, 1> raw_moment_3_buf(raw_moment_3.data(), raw_moment_3.size());
        sycl::buffer<RealType, 1> central_moment_2_buf(central_moment_2.data(), central_moment_2.size());
        sycl::buffer<RealType, 1> central_moment_3_buf(central_moment_3.data(), central_moment_3.size());
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::raw_moment(queue, dataset, mean_buf, raw_moment_2_buf, raw_moment_3_buf);
        oneapi::mkl::stats::central_moment(queue, dataset, central_moment_2_buf, central_moment_3_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Mean for dimension " << i << ": " << mean[i] << std::endl;
        std::cout << "raw moment 2 for dimension " << i << ": " << raw_moment_2[i] << std::endl;
        std::cout << "raw moment 3 for dimension " << i << ": " << raw_moment_3[i] << std::endl;
        std::cout << "central moment 2 for dimension " << i << ": " << central_moment_2[i] << std::endl;
        std::cout << "central moment 3 for dimension " << i << ": " << central_moment_3[i] << std::endl;
    }
    return 0;
}

template <typename RealType>
int run_central_moment_example_with_provided_mean(sycl::queue& queue) {
    std::cout << "\nRunning central moments example with provided mean" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> mean(n_dims);
    std::vector<RealType> central_moment_2(n_dims);
    std::vector<RealType> central_moment_3(n_dims);
    std::vector<RealType> central_moment_4(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> mean_buf(mean.data(), mean.size());
        sycl::buffer<RealType, 1> central_moment_2_buf(central_moment_2.data(), central_moment_2.size());
        sycl::buffer<RealType, 1> central_moment_3_buf(central_moment_3.data(), central_moment_3.size());
        sycl::buffer<RealType, 1> central_moment_4_buf(central_moment_4.data(), central_moment_4.size());
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::mean(queue, dataset, mean_buf);
        queue.wait_and_throw();

        oneapi::mkl::stats::central_moment(queue, mean_buf, dataset, central_moment_2_buf, central_moment_3_buf, central_moment_4_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }


    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Mean for dimension " << i << ": " << mean[i] << std::endl;
        std::cout << "central moment 2 for dimension " << i << ": " << central_moment_2[i] << std::endl;
        std::cout << "central moment 3 for dimension " << i << ": " << central_moment_3[i] << std::endl;
        std::cout << "central moment 4 for dimension " << i << ": " << central_moment_4[i] << std::endl;
    }
    return 0;
}

template <typename RealType>
int run_skewness_kurtosis_variation_example(sycl::queue& queue) {
    std::cout << "\nRunning skewness, kurtosis, variation example" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> variation(n_dims);
    std::vector<RealType> skewness(n_dims);
    std::vector<RealType> kurtosis(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> skewness_buf(skewness.data(), skewness.size());
        sycl::buffer<RealType, 1> kurtosis_buf(kurtosis.data(), kurtosis.size());
        sycl::buffer<RealType, 1> variation_buf(variation.data(), variation.size());

        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::skewness(queue, dataset, skewness_buf);
        oneapi::mkl::stats::kurtosis(queue, dataset, kurtosis_buf);
        oneapi::mkl::stats::variation(queue, dataset, variation_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Skewness for dimension " << i << ": " << skewness[i] << std::endl;
        std::cout << "Kurtosis for dimension " << i << ": " << kurtosis[i] << std::endl;
        std::cout << "Variation for dimension " << i << ": " << variation[i] << std::endl;
    }
    return 0;
}

template <typename RealType>
int run_skewness_kurtosis_variation_example_with_provided_mean(sycl::queue& queue) {
    std::cout << "\nRunning skewness, kurtosis, variation example with provided mean" << std::endl;

    const size_t n_observations = 100;
    const size_t n_dims = 3;

    std::vector<RealType> x;
    random_vector(x, n_observations, n_dims);

    std::vector<RealType> mean(n_dims);
    std::vector<RealType> variation(n_dims);
    std::vector<RealType> skewness(n_dims);
    std::vector<RealType> kurtosis(n_dims);

    try {
        sycl::buffer<RealType, 1> x_buf(x.data(), x.size());
        sycl::buffer<RealType, 1> mean_buf(mean.data(), mean.size());
        sycl::buffer<RealType, 1> skewness_buf(skewness.data(), skewness.size());
        sycl::buffer<RealType, 1> kurtosis_buf(kurtosis.data(), kurtosis.size());
        sycl::buffer<RealType, 1> variation_buf(variation.data(), variation.size());
        auto dataset = oneapi::mkl::stats::make_dataset<oneapi::mkl::stats::layout::row_major>(n_dims, n_observations, x_buf);
        oneapi::mkl::stats::mean(queue, dataset, mean_buf);
        queue.wait_and_throw();

        oneapi::mkl::stats::skewness(queue, mean_buf, dataset, skewness_buf);
        oneapi::mkl::stats::kurtosis(queue, mean_buf, dataset, kurtosis_buf);
        oneapi::mkl::stats::variation(queue, mean_buf, dataset, variation_buf);
        queue.wait_and_throw();
    }
    catch(sycl::exception const& e) {
        std::cout << "\t\tSYCL exception during calculation\n"
                  << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        return 1;
    }

    std::cout << "Number of dimensions: " << n_dims << std::endl;
    std::cout << "Number of observations: " << n_observations << std::endl;
    for(int i = 0; i < n_dims; i++) {
        std::cout << "Skewness for dimension " << i << ": " << skewness[i] << std::endl;
        std::cout << "Kurtosis for dimension " << i << ": " << kurtosis[i] << std::endl;
        std::cout << "Variation for dimension " << i << ": " << variation[i] << std::endl;
    }
    return 0;
}


//
// Description of example setup, apis used and supported floating point type precisions
//
void print_example_banner() {
    std::cout << "" << std::endl;
    std::cout << "########################################################################" << std::endl;
    std::cout << "# Basic statistics estimates example: " << std::endl;
    std::cout << "# " << std::endl;
    std::cout << "# Using apis:" << std::endl;
    std::cout << "# oneapi::mkl::stats::raw_sum oneapi::mkl::stats::central_sum" << std::endl;
    std::cout << "# oneapi::mkl::stats::raw_moment oneapi::mkl::stats::central_moment" << std::endl;
    std::cout << "# oneapi::mkl::stats::mean oneapi::mkl::stats::skewness oneapi::mkl::stats::kurtosis oneapi::mkl::stats::variance" << std::endl;
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
                if(run_sum_example<float>(queue) ||
                    run_central_sum_example_with_provided_mean<float>(queue) ||
                    run_moment_example<float>(queue) ||
                    run_central_moment_example_with_provided_mean<float>(queue) ||
                    run_skewness_kurtosis_variation_example<float>(queue) ||
                    run_skewness_kurtosis_variation_example_with_provided_mean<float>(queue)) {
                    std::cout << "FAILED" << std::endl;
                    return 1;
                }
                if(isDoubleSupported(my_dev)) {
                    std::cout << "\tRunning with double precision real data type:" << std::endl;
                    if(run_sum_example<double>(queue) ||
                       run_central_sum_example_with_provided_mean<double>(queue) ||
                       run_moment_example<double>(queue) ||
                       run_central_moment_example_with_provided_mean<double>(queue) ||
                       run_skewness_kurtosis_variation_example<double>(queue) ||
                       run_skewness_kurtosis_variation_example_with_provided_mean<double>(queue)) {
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
