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
 *            ordinary heap and stack pointers
 *            checking status and fixing results of numerical computation
 *            for SYCL-offloaded computation
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

namespace onemkl = oneapi::mkl;


constexpr onemkl::vm::status k_success   = onemkl::vm::status::success;
constexpr onemkl::vm::status k_errdom    = onemkl::vm::status::errdom;
constexpr onemkl::vm::status k_sing      = onemkl::vm::status::sing;
constexpr onemkl::vm::status k_underflow = onemkl::vm::status::underflow;
constexpr onemkl::vm::status k_overflow  = onemkl::vm::status::overflow;
constexpr onemkl::vm::status k_fix_all   = onemkl::vm::status::fix_all;




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

std::string decode_status(onemkl::vm::status st) {
    std::string r;

    if (has_any(st, onemkl::vm::status::errdom))     { r += "errdom ";    }
    if (has_any(st, onemkl::vm::status::sing))       { r += "sing ";      }
    if (has_any(st, onemkl::vm::status::overflow))   { r += "overflow ";  }
    if (has_any(st, onemkl::vm::status::underflow))  { r += "underflow "; }

    if (!st) { r = "success"; }

    return r;
}

void global_status(sycl::queue & queue) {
    std::vector<double>         a  { -2.0, -1.0, 0.0, 1.0 };
    std::vector<double>         y ( a.size() );

    onemkl::vm::set_status(queue, onemkl::vm::status::success);
    onemkl::vm::erfinv(queue, a.size(), a.data(), y.data(), {}, onemkl::vm::mode::global_status_report);
    auto st = onemkl::vm::get_status(queue);
    std::cout << decode_status(st) << std::endl;

    bool pass = has_only(st, onemkl::vm::status::errdom | onemkl::vm::status::sing);
    std::cout << "global_status: " << (pass ? "PASS" : "FAIL") << std::endl;
}

void local_status(sycl::queue & queue) {
    std::vector<double>         a  { -2.0, -1.0, 0.0, 1.0 };
    std::vector<double>         y ( a.size() );
    onemkl::vm::status          st = onemkl::vm::status::success;

    onemkl::vm::erfinv(queue, a.size(), a.data(), y.data(), {}, onemkl::vm::mode::la, onemkl::vm::error_handler<double> { &st });
    std::cout << decode_status(st) << std::endl;
    bool pass = !(st ^ (onemkl::vm::status::errdom | onemkl::vm::status::sing));
    std::cout << "local_status: " << (pass ? "PASS" : "FAIL") << std::endl;
}

void fix_numerical_results(sycl::queue & queue) {
    std::vector<double>         a  { -2.0, -1.0, 0.0, 1.0 };
    std::vector<double>         y0 ( a.size() );
    std::vector<double>         y1 ( a.size() );
    std::vector<double>         y2 ( a.size() );


    onemkl::vm::erfinv(queue, a.size(), a.data(), y1.data(), {}, onemkl::vm::mode::la);
    onemkl::vm::erfinv(queue, a.size(), a.data(), y1.data(), {}, onemkl::vm::mode::la, onemkl::vm::error_handler<double> { onemkl::vm::status::fix_all, 888.0, true });
    onemkl::vm::erfinv(queue, a.size(), a.data(), y2.data(), {}, onemkl::vm::mode::la, onemkl::vm::error_handler<double> { onemkl::vm::status::sing, 777.0, false });

    std::cout << std::setw(8) << "i" << " | "
        << std::setw(20) << "a[i]" << " | "
        << std::setw(20) << "y0[i]" << " | "
        << std::setw(20) << "y1[i]" << " | "
        << std::setw(20) << "y2[i]" << " |"
        << std::endl;


    for (int i = 0; i < a.size(); ++i) {
        std::cout << std::setw(8) << i << " | "
            << std::setw(20) << a[i] << " | "
            << std::setw(20) << y0[i] << " | "
            << std::setw(20) << y1[i] << " | "
            << std::setw(20) << y2[i] << " |"
            << std::endl;
    }

    bool pass = true;

    if (y1[0] != -888.0 || y1[1] != -888.0 || y1[3] != +888.0) { pass = false; }
    if (y2[1] != 777.0  || y2[3] != 777.0) { pass = false; }

    std::cout << "fix_numerical_results: " << (pass ? "PASS" : "FAIL") << std::endl;
}

void status_array(sycl::queue & queue) {
    std::mt19937_64 rng ( 90377 );

    std::vector<double>                      a  ( 1024 * 1024 );
    std::vector<double>                      y  ( a.size() );
    std::vector<onemkl::vm::status>          st ( a.size() );

    std::generate(a.begin(), a.end(), [&]() { return static_cast<double>(rng() >> 60) / 4.0 - 2.0; });

    onemkl::vm::erfinv(queue, a.size(), a.data(), y.data(), {}, onemkl::vm::mode::la, onemkl::vm::error_handler<double> { st.data(), static_cast<int64_t>(st.size()) });

    int errdom_found = 0;
    int sing_found = 0;

    int errdom = 0;
    int sing = 0;
    int underflow = 0;
    int overflow = 0;

    for (int i = 0; i < a.size(); ++i) {
        if (has_any(st[i], k_errdom))      { ++errdom; }
        if (has_any(st[i], k_sing))        { ++sing; }
        if (has_any(st[i], k_underflow))   { ++underflow; }
        if (has_any(st[i], k_overflow))    { ++overflow; }

        if (a[i] == -1.0 || a[i] == 1.0) { ++sing_found; }
        if (a[i] < -1.0 || a[i] > 1.0)   { ++errdom_found; }
    }

    std::cout << std::setw(12) << "status" << " | "
        << std::setw(20) << "number of elements" << " | "
        << std::setw(20) << "expected" << " | "
        << std::endl;

    std::cout << std::setw(12) << "errdom" << " | "
        << std::setw(20) << errdom << " | "
        << std::setw(20) << errdom_found << " | "
        << std::endl;

    std::cout << std::setw(12) << "sing" << " | "
        << std::setw(20) << sing << " | "
        << std::setw(20) << sing_found << " | "
        << std::endl;

    std::cout << std::setw(12) << "underflow" << " | "
        << std::setw(20) << underflow << " | "
        << std::setw(20) << 0 << " | "
        << std::endl;

    std::cout << std::setw(12) << "overflow" << " | "
        << std::setw(20) << overflow << " | "
        << std::setw(20) << 0 << " | "
        << std::endl;

    bool pass;
    pass = (errdom_found == errdom) && (sing_found && sing) && (0 == underflow) && (0 == overflow);
    std::cout << "status_array: " << (pass ? "PASS" : "FAIL") << std::endl;
}

void all_on(sycl::queue & queue) {
    std::mt19937_64 rng ( 90377 );

    std::vector<double>                      a  ( 1024 * 1024 );
    std::vector<double>                      y  ( a.size() );
    std::vector<onemkl::vm::status>          st ( a.size() );


    std::generate(a.begin(), a.end(), [&]() { return static_cast<double>(rng() >> 60) / 4.0 - 2.0; });
    onemkl::vm::erfinv(queue, a.size(), a.data(), y.data(), {}, onemkl::vm::mode::la, onemkl::vm::error_handler<double> { st.data(), static_cast<int64_t>(st.size()), k_fix_all, +888.0, false });

    int n_fixed = 0;
    int errdom_found = 0;
    int sing_found = 0;

    int errdom = 0;
    int sing = 0;
    int underflow = 0;
    int overflow = 0;

    for (int i = 0; i < a.size(); ++i) {
        if (has_any(st[i], k_errdom))      { ++errdom; }
        if (has_any(st[i], k_sing))        { ++sing; }
        if (has_any(st[i], k_underflow))   { ++underflow; }
        if (has_any(st[i], k_overflow))    { ++overflow; }

        if (a[i] == -1.0 || a[i] == 1.0) { ++sing_found; }
        if (a[i] < -1.0 || a[i] > 1.0)   { ++errdom_found; }

        n_fixed += (y[i] == +888.0);
    }


    std::cout << std::setw(12) << "status" << " | "
        << std::setw(20) << "number of elements" << " | "
        << std::setw(20) << "expected" << " | "
        << std::endl;

    std::cout << std::setw(12) << "errdom" << " | "
        << std::setw(20) << errdom << " | "
        << std::setw(20) << errdom_found << " | "
        << std::endl;

    std::cout << std::setw(12) << "sing" << " | "
        << std::setw(20) << sing << " | "
        << std::setw(20) << sing_found << " | "
        << std::endl;

    std::cout << std::setw(12) << "underflow" << " | "
        << std::setw(20) << underflow << " | "
        << std::setw(20) << 0 << " | "
        << std::endl;

    std::cout << std::setw(12) << "overflow" << " | "
        << std::setw(20) << overflow << " | "
        << std::setw(20) << 0 << " | "
        << std::endl;


    std::cout << "fixed " << n_fixed << " values " << std::endl;

    bool pass;
    pass = (errdom_found == errdom) && (sing_found == sing) && (0 == underflow) && (0 == overflow) && (n_fixed == (sing + errdom));
    std::cout << "all_on: " << (pass ? "PASS" : "FAIL") << std::endl;
}



void run_on(sycl::device & dev) {
    preamble(dev);

    sycl::queue queue { dev, async_sycl_error };

    std::cout << std::fixed << std::setprecision(10);
    global_status(queue);
    local_status(queue);
    fix_numerical_results(queue);
    status_array(queue);
    all_on(queue);
    std::cout << std::endl << std::endl;
}


} // namespace

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
