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
 *            oneapi::mkl::vm::sqrt example program text (SYCL interface)
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
    {MAX_HA_ULP_S, 4.5},     {MAX_LA_ULP_S, 5.0},     {MAX_EP_ULP_S, FLT_MAX},
    {MAX_HA_ULP_D, 2.0},     {MAX_LA_ULP_D, 5.0},     {MAX_EP_ULP_D, 7.0E7},
    {MAX_HA_ULP_C, FLT_MAX}, {MAX_LA_ULP_C, FLT_MAX}, {MAX_EP_ULP_C, FLT_MAX},
    {MAX_HA_ULP_Z, DBL_MAX}, {MAX_LA_ULP_Z, DBL_MAX}, {MAX_EP_ULP_Z, DBL_MAX},
};

//!
//! @brief Accuracy test
//!
template <typename A, typename R>
bool vSqrtAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x41093E24) }, { UINT32_C(0x403B70E7) } }, //  0: vsSqrt ( 8.57767105      ) = ( 2.92876601      );
{ { UINT32_C(0x4093852F) }, { UINT32_C(0x400969F8) } }, //  1: vsSqrt ( 4.61000776      ) = ( 2.14709282      );
{ { UINT32_C(0x41011D03) }, { UINT32_C(0x4035CE0C) } }, //  2: vsSqrt ( 8.06958294      ) = ( 2.8407011       );
{ { UINT32_C(0x41034C40) }, { UINT32_C(0x40375621) } }, //  3: vsSqrt ( 8.20611572      ) = ( 2.86463189      );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x40076E1CE786FE27) } }, //  0: vdSqrt ( 8.57767068267691535       ) = ( 2.9287660682746437        );
{ { UINT64_C(0x401270A5F32DAE19) }, { UINT64_C(0x40012D3F0ED3A6F0) } }, //  1: vdSqrt ( 4.6100080486899282        ) = ( 2.14709292968188237       );
{ { UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x4006B9C187E1F54A) } }, //  2: vdSqrt ( 8.06958309145159269       ) = ( 2.84070116194076139       );
{ { UINT64_C(0x40206988134D9FDD) }, { UINT64_C(0x4006EAC429F8399F) } }, //  3: vdSqrt ( 8.20611629793705255       ) = ( 2.86463196553013644       );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0x4093852F), UINT32_C(0x41093E24) }, { UINT32_C(0x402B6B72), UINT32_C(0x3FCCF5B2) } }, //  0: vcSqrt ( 4.61000776      + i * 8.57767105      ) = ( 2.67843294      + i * 1.60124803      );
{ { UINT32_C(0x41034C40), UINT32_C(0x41011D03) }, { UINT32_C(0x4048F083), UINT32_C(0x3FA47E0B) } }, //  1: vcSqrt ( 8.20611572      + i * 8.06958294      ) = ( 3.13967967      + i * 1.28509653      );
{ { UINT32_C(0x4036ECDE), UINT32_C(0x41136B29) }, { UINT32_C(0x40200838), UINT32_C(0x3FEBD28B) } }, //  2: vcSqrt ( 2.85820723      + i * 9.21366215      ) = ( 2.50050163      + i * 1.84236276      );
{ { UINT32_C(0x40FDFDE5), UINT32_C(0x4082ABE3) }, { UINT32_C(0x4039D6BB), UINT32_C(0x3F34013F) } }, //  3: vcSqrt ( 7.93724298      + i * 4.08348227      ) = ( 2.90373111      + i * 0.703144014     );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0x401270A5F32DAE19), UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x40056D6E42511161), UINT64_C(0x3FF99EB631575346) } }, //  0: vzSqrt ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 2.67843295869731479       + i * 1.60124797128556073       );
{ { UINT64_C(0x40206988134D9FDD), UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40091E1072F5F2ED), UINT64_C(0x3FF48FC159BA032B) } }, //  1: vzSqrt ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 3.13967981160236898       + i * 1.28509650277573928       );
{ { UINT64_C(0x4006DD9BBAC0EE6B), UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x40040106EB21DA6F), UINT64_C(0x3FFD7A515924BC6F) } }, //  2: vzSqrt ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 2.50050147721349658       + i * 1.84236273595504563       );
{ { UINT64_C(0x401FBFBCBB737F7A), UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x40073AD76D49B864), UINT64_C(0x3FE68027E8B0FCB8) } }, //  3: vzSqrt ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 2.90373120671488216       + i * 0.703144030070595782      );
 }

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
    oneapi::mkl::vm::sqrt(main_queue, varg1.size(), in1, out1, accuracy_mode[acc]);

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
                << "sqrt"
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
            << "sqrt"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "sqrt" << std::endl;
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
      errs += vSqrtAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vSqrtAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vSqrtAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vSqrtAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
