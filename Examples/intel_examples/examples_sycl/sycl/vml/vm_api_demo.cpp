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
 *            demonstration of usage of VM APIs:
 *            for SYCL buffer,
 *            USM shared and device pointers,
 *            ordinary heap and stack pointers,
 *            error handler in replacement(fixup) mode
 *                          on
 *            generation of random normal variable N(0, 1)
 *            using inverse cumulative distribution function
 *
 *******************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>
#include <type_traits>

#include <CL/sycl.hpp>
#include "oneapi/mkl.hpp"

#include "common_for_examples.hpp"

namespace {

using std::int64_t;
using std::uint32_t;
using std::uint64_t;

void preamble(sycl::device & dev) {

    std::string dev_name       = dev.template get_info<sycl::info::device::name>();
    std::string driver_version = dev.template get_info<sycl::info::device::version>();

    std::cerr << std::endl
              << "running on:         " << std::endl
              << "       device name: " << dev_name << std::endl
              << "    driver version: " << driver_version << std::endl;
}

void async_sycl_error(sycl::exception_list el) {
    std::cerr << "async exceptions caught: " << std::endl;

    for (auto l = el.begin(); l != el.end(); ++l) {
        try {
            std::rethrow_exception(*l);
        } catch(const sycl::exception & e) {
            std::cerr << "SYCL exception occured with code " << e.get_cl_code() << " with " << e.what() << std::endl;
        }
    }
}

template <typename T>
struct UniformFiller {
    std::mt19937_64 rng_;
    UniformFiller(uint64_t seed): rng_ { seed } { }

    double gen64() {
        union {
            uint64_t w;
            double d;
        } arg = { .w = rng_() };

        arg.w &= UINT64_C(0x000F'FFFF'FFFF'FFFF);
        arg.w |= UINT64_C(0x3FF0'0000'0000'0000);
        arg.d -= 1.0;
        return arg.d;
    }

    float gen32() {
        union {
            uint32_t w32;
            float f;
        } arg = { .w32 = static_cast<uint32_t>(rng_() >> 32) };

        arg.w32 &= UINT32_C(0x007F'FFFF);
        arg.w32 |= UINT32_C(0x3F80'0000);
        arg.f -= 1.0f;
        return arg.f;
    }

    T operator()() {
        static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value, "UniformFiller is for floats and double's");
        if constexpr(std::is_same<T, float>::value) { return gen32(); }
        else if constexpr(std::is_same<T, double>::value) { return gen64(); }
    }

};

bool check_mean(double mean, double expected, double sigma, int64_t n) {
    double adiff = std::fabs(mean - expected);
    double err_estimate = sigma / std::sqrt(n * 1.0);

    return (adiff / err_estimate) < 3.0;
}

double check_stddev(double std_dev, double expected, double sigma, int64_t n) {
    double adiff = std::fabs(std_dev - expected);
    double err_estimate = sigma / std::sqrt(2 * (n - 1.0));
    return (adiff / err_estimate) < 3.0;
}

template<typename T>
void print_results(const char * method, int64_t n, T * y) {
    double s1 = std::accumulate(y, y + n, 0.0, [=](double s, double t) { return s + t; });
    double s2 = std::accumulate(y, y + n, 0.0, [=](double s, double t) { return s + t * t; });

    double mean     = s1 / n;
    double stddev   = std::sqrt((s2 - s1  * s1 / n) / (n - 1));

    std::string float_type_string { (sizeof(T) == 4) ? "float" : "double" };

    std::cout << std::setw(8) << float_type_string
              << std::setw(10) << method
              << "       n = " << std::setw(10) << n
              << " mean    = " << std::setw(16) << mean
              << std::setw(10) << ((check_mean(mean, 0.0, 1.0, n)) ? "( PASS )" : "( FAIL )")
              << " std.dev = " << std::setw(16) << stddev
              << std::setw(10) << ((check_stddev(stddev, 1.0, 1.0, n)) ? "( PASS )" : "( FAIL )")
              << std::endl;
}

template <typename T>
void run_usm(int64_t n, sycl::queue & queue) {
    namespace one = oneapi::mkl;

    // user can mix device and shared pointers
    T * a = sycl::malloc_shared<T>(n, queue);
    T * y = new T[n];
    T * dev_y = sycl::malloc_device<T>(n, queue);


    std::generate(a, a + n, UniformFiller<T>(88883)); // shared usm is accessible from host
    std::fill(y, y + n, std::nan(""));

    queue.memcpy(dev_y, y, n * sizeof(T)); // memcpy to device to clean device memory
    queue.wait_and_throw(); // memcpy is async too


    one::vm::cdfnorminv(queue, n, a, dev_y, { /* no dependent events */ }, one::vm::mode::ha, { one::vm::status::sing,  7.0f  });
    queue.wait_and_throw(); // USM call is asynchronous so wait is needed

    queue.memcpy(y, dev_y, n * sizeof(T)); // memcpy back to host
    queue.wait_and_throw(); // memcpy is async too

    print_results("on_usm", n, y);

    sycl::free(dev_y, queue);
    delete[] y;
    sycl::free(a, queue);

}

template <typename T>
void run_buffer(int64_t n, sycl::queue & queue) {
    namespace one = oneapi::mkl;

    T * a = new T[n];
    T * y = new T[n];

    std::generate(a, a + n, UniformFiller<T>(90001));

    {
        sycl::buffer<T, 1> buf_a { a , a + n }; // SYCL buffer which copies data from 'a', but does not copy back
        sycl::buffer<T, 1> buf_y { y , n };     // SYCL buffer which copy back to 'y' on destruction

        one::vm::cdfnorminv(queue, n, buf_a, buf_y, one::vm::mode::ha, { one::vm::status::sing,  7.0f  });

    } // buf_y destructed, data now in 'y'

    print_results("on_buffer", n, y);

    delete[] y;
    delete[] a;
}

template <typename T, size_t n>
void run_stack(sycl::queue & queue) {
    namespace one = oneapi::mkl;

    T a[n];
    T y[n];

    std::generate(a, a + n, UniformFiller<T>(80777));
    std::fill(y, y + n, std::nan(""));

    one::vm::cdfnorminv(queue, n, a, y, { /* no dependent events */ }, one::vm::mode::ha, { one::vm::status::sing,  7.0f  });
    // function returns with result ready when used with heap pointer as output

    print_results("on_stack", n, y);
}

template <typename T>
void run_heap(int64_t n, sycl::queue & queue) {
    namespace one = oneapi::mkl;

    T * a = new T[n];
    T * y = new T[n];

    std::generate(a, a + n, UniformFiller<T>(99001));
    std::fill(y, y + n, std::nan(""));

    one::vm::cdfnorminv(queue, n, a, y, { /* no dependent events */ }, one::vm::mode::ha, { one::vm::status::sing,  7.0f  });
    // function returns with result ready when used with heap pointer as output

    print_results("on_heap", n, y);

    delete[] y;
    delete[] a;
}

void run_on(sycl::device & dev) {
    double mean     = std::nan("");
    double std_dev  = std::nan("");

    constexpr int vector_stack_len   = 1024;
    int64_t       vector_heap_len    = 10'000'000;
    int64_t       vector_buffer_len  = 10'000'000;
    int64_t       vector_usm_len     = 10'000'000;

    preamble(dev);

    sycl::queue queue { dev, async_sycl_error };

    std::cout << std::fixed << std::setprecision(10);

    run_stack<float, vector_stack_len>(queue);
    run_heap<float>(vector_heap_len, queue);
    run_buffer<float>(vector_buffer_len, queue);
    run_usm<float>(vector_usm_len, queue);
    run_stack<double, vector_stack_len>(queue);
    run_heap<double>(vector_heap_len, queue);
    run_buffer<double>(vector_buffer_len, queue);
    run_usm<double>(vector_usm_len, queue);

    std::cout << std::endl << std::endl;
}


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
int main(int argc, char **argv) {
  // List of available devices
  std::list<my_sycl_device_types> list_of_devices;
  set_list_of_devices(list_of_devices);

  // Loop by all available devices
  for (auto dev_type : list_of_devices) {
    cl::sycl::device my_dev;
    bool my_dev_is_found = false;
    get_sycl_device(my_dev, my_dev_is_found, dev_type);

    // Run tests if the device is available
    if (my_dev_is_found) {
        std::cout << "Running tests on " << sycl_device_names[dev_type] << ".\n";
        run_on(my_dev);
    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
