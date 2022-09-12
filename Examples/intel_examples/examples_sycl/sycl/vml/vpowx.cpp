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
 *            oneapi::mkl::vm::powx example program text (SYCL interface)
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
    {MAX_HA_ULP_S, 4.5},     {MAX_LA_ULP_S, 5.0},     {MAX_EP_ULP_S, 5.0E3},
    {MAX_HA_ULP_D, 2.0},     {MAX_LA_ULP_D, 5.0},     {MAX_EP_ULP_D, 7.0E7},
    {MAX_HA_ULP_C, FLT_MAX}, {MAX_LA_ULP_C, FLT_MAX}, {MAX_EP_ULP_C, FLT_MAX},
    {MAX_HA_ULP_Z, DBL_MAX}, {MAX_LA_ULP_Z, DBL_MAX}, {MAX_EP_ULP_Z, DBL_MAX},
};

//!
//! @brief Accuracy test
//!
template <typename A, typename R>
bool vPowxAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x41093E24) }, { UINT32_C(0x4048F5C3) }, { UINT32_C(0x44552AC1) } }, //  0: vsPowx ( 8.57767105     , 3.1400001       ) = ( 852.66803       );
{ { UINT32_C(0x4093852F) }, { UINT32_C(0x4048F5C3) }, { UINT32_C(0x42F2B0D7) } }, //  1: vsPowx ( 4.61000776     , 3.1400001       ) = ( 121.34539       );
{ { UINT32_C(0x41011D03) }, { UINT32_C(0x4048F5C3) }, { UINT32_C(0x442FF9C3) } }, //  2: vsPowx ( 8.06958294     , 3.1400001       ) = ( 703.902527      );
{ { UINT32_C(0x41034C40) }, { UINT32_C(0x4048F5C3) }, { UINT32_C(0x44397EBA) } }, //  3: vsPowx ( 8.20611572     , 3.1400001       ) = ( 741.980103      );
 }

,
    .data_f64 = std::vector<data_3_f64_t>

    {

{ { UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0x408AA55778FCE709) } }, //  0: vdPowx ( 8.57767068267691535      , 3.14000000000000012       ) = ( 852.667711234856256       );
{ { UINT64_C(0x401270A5F32DAE19) }, { UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0x405E561AE4453C46) } }, //  1: vdPowx ( 4.6100080486899282       , 3.14000000000000012       ) = ( 121.345391337979066       );
{ { UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0x4085FF381E850103) } }, //  2: vdPowx ( 8.06958309145159269      , 3.14000000000000012       ) = ( 703.902401961415649       );
{ { UINT64_C(0x40206988134D9FDD) }, { UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0x40872FD74F3829E4) } }, //  3: vdPowx ( 8.20611629793705255      , 3.14000000000000012       ) = ( 741.980131567743683       );
 }

,

    .data_c32 = std::vector <data_3_c32_t>

    {

{ { UINT32_C(0x4093852F), UINT32_C(0x41093E24) }, { UINT32_C(0x4048F5C3), UINT32_C(0x4048F5C3) }, { UINT32_C(0xC19A8841), UINT32_C(0xC21A0136) } }, //  0: vcPowx ( 4.61000776      + i * 8.57767105     , 3.1400001       + i * 3.1400001       ) = ( -19.3165302     + i * -38.5011826     );
{ { UINT32_C(0x41034C40), UINT32_C(0x41011D03) }, { UINT32_C(0x4048F5C3), UINT32_C(0x4048F5C3) }, { UINT32_C(0xC310B7AE), UINT32_C(0xC2ED2BD0) } }, //  1: vcPowx ( 8.20611572      + i * 8.06958294     , 3.1400001       + i * 3.1400001       ) = ( -144.717499     + i * -118.585571     );
{ { UINT32_C(0x4036ECDE), UINT32_C(0x41136B29) }, { UINT32_C(0x4048F5C3), UINT32_C(0x4048F5C3) }, { UINT32_C(0x401FC609), UINT32_C(0xC1B5CAED) } }, //  2: vcPowx ( 2.85820723      + i * 9.21366215     , 3.1400001       + i * 3.1400001       ) = ( 2.49646211      + i * -22.7240849     );
{ { UINT32_C(0x40FDFDE5), UINT32_C(0x4082ABE3) }, { UINT32_C(0x4048F5C3), UINT32_C(0x4048F5C3) }, { UINT32_C(0xC2D4B709), UINT32_C(0x433D8524) } }, //  3: vcPowx ( 7.93724298      + i * 4.08348227     , 3.1400001       + i * 3.1400001       ) = ( -106.357491     + i * 189.520081      );
 }

,
    .data_c64 = std::vector <data_3_c64_t>

    {

{ { UINT64_C(0x401270A5F32DAE19), UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x40091EB851EB851F), UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0xC033510970136A6A), UINT64_C(0xC043402653EE7617) } }, //  0: vzPowx ( 4.6100080486899282        + i * 8.57767068267691535      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -19.3165502593423426      + i * -38.5011696733819733      );
{ { UINT64_C(0x40206988134D9FDD), UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40091EB851EB851F), UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0xC06216F614F1FE9F), UINT64_C(0xC05DA579773A3C38) } }, //  1: vzPowx ( 8.20611629793705255       + i * 8.06958309145159269      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -144.717539284368257      + i * -118.585538679952947      );
{ { UINT64_C(0x4006DD9BBAC0EE6B), UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x40091EB851EB851F), UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0x4003F8B937C10DF9), UINT64_C(0xC036B95D4F1E546A) } }, //  2: vzPowx ( 2.8582071867111174        + i * 9.21366148563738108      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( 2.49644702489763093       + i * -22.7240800332114432      );
{ { UINT64_C(0x401FBFBCBB737F7A), UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x40091EB851EB851F), UINT64_C(0x40091EB851EB851F) }, { UINT64_C(0xC05A96E0D33FD2E3), UINT64_C(0x4067B0A4892EA43F) } }, //  3: vzPowx ( 7.93724339382594657       + i * 4.0834825258625127       , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -106.357472240760714      + i * 189.520084944817398       );
 }

};


  // *************************************************************
  // Variable declaraions
  // *************************************************************
  // Input arguments
  A arg1;
  std::vector<A> varg1;

  A arg2;

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

    oneapi::mkl::vm::powx(main_queue, varg1.size(), in1, arg2, out1,
                  accuracy_mode[acc]);
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
                << "powx"
                << "(" << varg1[i]

                << "," << arg2 << ") = " << vres1[i] << std::endl;

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
            << "powx"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "powx" << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Supported floating point type precisions:" << std::endl;

  std::cout << "#   float" << std::endl;

  std::cout << "#   double" << std::endl;

  std::cout << "#   std::complex<float>" << std::endl;

  std::cout << "#   std::complex<double>" << std::endl;

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
      errs += vPowxAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vPowxAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vPowxAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vPowxAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
