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
 *            oneapi::mkl::vm::erfinv example program text (SYCL interface)
 *
 *******************************************************************************/

#include <CL/sycl.hpp>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <vector>

#include "mkl.h"
#include "oneapi/mkl/types.hpp"
#include "oneapi/mkl/vm.hpp"

#include "common_for_examples.hpp"

#include "vml_common.hpp"

static constexpr int VLEN = 4;

static ulp_table_type ulp_table = {
    {MAX_HA_ULP_S, FLT_MAX}, {MAX_LA_ULP_S, FLT_MAX}, {MAX_EP_ULP_S, FLT_MAX},
    {MAX_HA_ULP_D, DBL_MAX}, {MAX_LA_ULP_D, DBL_MAX}, {MAX_EP_ULP_D, DBL_MAX},
    {MAX_HA_ULP_C, 2.0},     {MAX_LA_ULP_C, 5.0},     {MAX_EP_ULP_C, 5.0E3},
    {MAX_HA_ULP_Z, 2.0},     {MAX_LA_ULP_Z, 5.0},     {MAX_EP_ULP_Z, 7.0E7},
};

//!
//! @brief Accuracy test
//!
template <typename A, typename R>
bool vErfInvAccuracyLiteTest(const cl::sycl::device &dev) {
  static constexpr int ACCURACY_LEN = VLEN;

  // *************************************************************
  // Data table declaraion
  // *************************************************************
data_2_t data
{
    .i = 0
,

    .data_f32 = std::vector<data_2_f32_t>

    {

{ { UINT32_C(0x3F2C10F1) }, { UINT32_C(0x3F311CED) } }, //  0: vsErfInv ( 0.672133505     ) = ( 0.691847622     );
{ { UINT32_C(0xBE5BB8DC) }, { UINT32_C(0xBE4521D1) } }, //  1: vsErfInv ( -0.21457237     ) = ( -0.192511812    );
{ { UINT32_C(0x3F0EFF63) }, { UINT32_C(0x3F0B597E) } }, //  2: vsErfInv ( 0.558584392     ) = ( 0.544334292     );
{ { UINT32_C(0x3F16CF15) }, { UINT32_C(0x3F14DA73) } }, //  3: vsErfInv ( 0.589097321     ) = ( 0.581458271     );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x3FE5821DD706CAE3) }, { UINT64_C(0x3FE6239D393954E6) } }, //  0: vdErfInv ( 0.672133369420717108      ) = ( 0.69184743095925394       );
{ { UINT64_C(0xBFCB771B67504648) }, { UINT64_C(0xBFC8A43A11CD5897) } }, //  1: vdErfInv ( -0.214572358556823994     ) = ( -0.192511805241058126     );
{ { UINT64_C(0x3FE1DFEC4CAA263D) }, { UINT64_C(0x3FE16B2FAA0C0A6A) } }, //  2: vdErfInv ( 0.558584356055866871      ) = ( 0.544334251521218393      );
{ { UINT64_C(0x3FE2D9E26D4AF341) }, { UINT64_C(0x3FE29B4E1787A761) } }, //  3: vdErfInv ( 0.589097226583909728      ) = ( 0.581458135563689749      );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    { /* empty */ }

,
    .data_c64 = std::vector <data_2_c64_t>

    { /* empty */ }

};


  // *************************************************************
  // Variable declaraions
  // *************************************************************
  // Input arguments
  A arg1;
  std::vector<A> varg1;
  // Output results
  R ref1;
  std::vector<R> vref1;
  std::vector<R> vres1;

  // Number of errors
  int errs = 0;
  // Number of printed errors
  int printed_errs = 0;

  // *************************************************************
  // Vector input data initialization
  // *************************************************************
  for (int i = 0; i < ACCURACY_LEN; ++i) {
    // Getting values from reference data table

    data.get_values(arg1, ref1);

    // Pushing values into vectors
    varg1.push_back(arg1);
    vref1.push_back(ref1);

  } // for (int i = 0; i < ACCURACY_LEN; ++i)

  // Allocate memory for results
  vres1.assign(varg1.size(), 0);

  // Catch asynchronous exceptions
  auto exception_handler = [](cl::sycl::exception_list exceptions) {
    for (std::exception_ptr const &e : exceptions) {
      try {
        std::rethrow_exception(e);
      } catch (cl::sycl::exception const &e) {
        std::cout << "Caught asynchronous SYCL exception:\n"
                  << e.what() << std::endl;
      }
    } // for (std::exception_ptr const& e : exceptions)
  };

  // Create execution queue with asynchronous error handling
  cl::sycl::queue main_queue(dev, exception_handler);

  // Get device name
  std::string dev_name =
      main_queue.get_device().get_info<cl::sycl::info::device::name>();

  // *************************************************************
  // Loop by all 3 accuracy modes of VM: HA, LA, EP:
  // set computation mode, run VM and check results
  // *************************************************************
  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    // Clear result vectors
    std::fill(vres1.begin(), vres1.end(), 0);

    // Create sycl buffers
    cl::sycl::buffer<A, 1> in1(varg1.begin(), varg1.end());

    cl::sycl::buffer<R, 1> out1(vres1.begin(), vres1.end());

    // Run VM function
    oneapi::mkl::vm::erfinv(main_queue, varg1.size(), in1, out1, accuracy_mode[acc]);

    // Get results from sycl buffers
    auto host_vres1 = out1.template get_access<cl::sycl::access::mode::read>();

    for (int i = 0; i < vres1.size(); ++i) {
      vres1[i] = host_vres1[i];
    }

    // Catch sycl exceptions
    try {
      main_queue.wait_and_throw();
    } catch (cl::sycl::exception &e) {
      std::cerr << "SYCL exception during Accuracy Test\n"
                << e.what() << std::endl
                << "OpenCl status: " << e.get_cl_code() << std::endl;
      return false;
    }

    // *************************************************************
    // Compute ulp between computed and expected (reference)
    // values and check
    // *************************************************************
    for (int i = 0; i < ACCURACY_LEN; ++i) {
      // Compute ulp values for HA, LA and EP results
      R ulp = compute_ulp<R>(vres1[i], vref1[i]);

      std::cout << "\t\t" << accuracy_name[acc] << " vm::"
                << "erfinv"
                << "(" << varg1[i] << ") = " << vres1[i]

                << std::endl;

      // Check ulp
      if (check_ulp(ulp, false, accuracy_mode[acc], ulp_table) == false) {
        errs++;
      }
    } // for (int i = 0; i < varg1.size (); ++i)
  }   // for (int acc = 0; acc < 3; ++acc)

  std::cout << "\tResult: " << ((errs == 0) ? "PASS" : "FAIL") << std::endl;

  return (errs == 0);
} // template <typename A, typename R> bool vFuncAccuracyText (const device
  // &dev)
//
// Description of example setup, apis used and supported floating point type
// precisions
//
void print_example_banner() {
  std::cout << "" << std::endl;
  std::cout << "###############################################################"
               "#########"
            << std::endl;
  std::cout << "# General VM "
            << "erfinv"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "erfinv" << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Supported floating point type precisions:" << std::endl;

  std::cout << "#   float" << std::endl;

  std::cout << "#   double" << std::endl;

  std::cout << "###############################################################"
               "#########"
            << std::endl;
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
int main(int argc, char **argv) {
  // Number of errors occured
  int errs = 0;

  // Print standard banner for VM examples
  print_example_banner();

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

      std::cout << "\tRunning with single precision real data type:"
                << std::endl;
      errs += vErfInvAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vErfInvAccuracyLiteTest<double, double>(my_dev);
    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}