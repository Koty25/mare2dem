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
 *            oneapi::mkl::vm::powr example program text (SYCL interface)
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
bool vPowrAccuracyLiteTest(const cl::sycl::device &dev) {
  static constexpr int ACCURACY_LEN = VLEN;

  // *************************************************************
  // Data table declaraion
  // *************************************************************
data_3_t data
{
    .i = 0
,

    .data_f32 = std::vector<data_3_f32_t>

    {

{ { UINT32_C(0x41093E24) }, { UINT32_C(0x4093852F) }, { UINT32_C(0x469CE711) } }, //  0: vsPowr ( 8.57767105     , 4.61000776      ) = ( 20083.5332      );
{ { UINT32_C(0x41011D03) }, { UINT32_C(0x41034C40) }, { UINT32_C(0x4BD2F79E) } }, //  1: vsPowr ( 8.06958294     , 8.20611572      ) = ( 27651900        );
{ { UINT32_C(0x41136B29) }, { UINT32_C(0x4036ECDE) }, { UINT32_C(0x440EB888) } }, //  2: vsPowr ( 9.21366215     , 2.85820723      ) = ( 570.883301      );
{ { UINT32_C(0x4082ABE3) }, { UINT32_C(0x40FDFDE5) }, { UINT32_C(0x478A3D0D) } }, //  3: vsPowr ( 4.08348227     , 7.93724298      ) = ( 70778.1016      );
 }

,
    .data_f64 = std::vector<data_3_f64_t>

    {

{ { UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x401270A5F32DAE19) }, { UINT64_C(0x40D39CE2AABD156D) } }, //  0: vdPowr ( 8.57767068267691535      , 4.6100080486899282        ) = ( 20083.5416710576283       );
{ { UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40206988134D9FDD) }, { UINT64_C(0x417A5EF61F31C368) } }, //  1: vdPowr ( 8.06958309145159269      , 8.20611629793705255       ) = ( 27651937.9496492445       );
{ { UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x4006DD9BBAC0EE6B) }, { UINT64_C(0x4081D7109BAA7980) } }, //  2: vdPowr ( 9.21366148563738108      , 2.8582071867111174        ) = ( 570.883109409172903       );
{ { UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x401FBFBCBB737F7A) }, { UINT64_C(0x40F147A2E685447B) } }, //  3: vdPowr ( 4.0834825258625127       , 7.93724339382594657       ) = ( 70778.1812794375437       );
 }

,

    .data_c32 = std::vector <data_3_c32_t>

    { /* empty */ }

,
    .data_c64 = std::vector <data_3_c64_t>

    { /* empty */ }

};


  // *************************************************************
  // Variable declaraions
  // *************************************************************
  // Input arguments
  A arg1;
  std::vector<A> varg1;

  A arg2;

  std::vector<A> varg2;

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

    data.get_values(arg1, arg2, ref1);

    // Pushing values into vectors
    varg1.push_back(arg1);
    vref1.push_back(ref1);

    varg2.push_back(arg2);

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

    cl::sycl::buffer<A, 1> in2(varg2.begin(), varg2.end());

    cl::sycl::buffer<R, 1> out1(vres1.begin(), vres1.end());

    // Run VM function

    oneapi::mkl::vm::powr(main_queue, varg1.size(), in1, in2, out1, accuracy_mode[acc]);
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
                << "powr"
                << "(" << varg1[i]

                << "," << varg2[i] << ") = " << vres1[i] << std::endl;

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
            << "powr"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "powr" << std::endl;
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
      errs += vPowrAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vPowrAccuracyLiteTest<double, double>(my_dev);
    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
