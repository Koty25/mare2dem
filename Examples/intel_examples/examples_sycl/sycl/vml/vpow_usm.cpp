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
 *            oneapi::mkl::vm::pow example program text (SYCL interface)
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
bool vPowAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x41093E24) }, { UINT32_C(0x4093852F) }, { UINT32_C(0x469CE711) } }, //  0: vsPow ( 8.57767105     , 4.61000776      ) = ( 20083.5332      );
{ { UINT32_C(0x41011D03) }, { UINT32_C(0x41034C40) }, { UINT32_C(0x4BD2F79E) } }, //  1: vsPow ( 8.06958294     , 8.20611572      ) = ( 27651900        );
{ { UINT32_C(0x41136B29) }, { UINT32_C(0x4036ECDE) }, { UINT32_C(0x440EB888) } }, //  2: vsPow ( 9.21366215     , 2.85820723      ) = ( 570.883301      );
{ { UINT32_C(0x4082ABE3) }, { UINT32_C(0x40FDFDE5) }, { UINT32_C(0x478A3D0D) } }, //  3: vsPow ( 4.08348227     , 7.93724298      ) = ( 70778.1016      );
 }

,
    .data_f64 = std::vector<data_3_f64_t>

    {

{ { UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x401270A5F32DAE19) }, { UINT64_C(0x40D39CE2AABD156D) } }, //  0: vdPow ( 8.57767068267691535      , 4.6100080486899282        ) = ( 20083.5416710576283       );
{ { UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40206988134D9FDD) }, { UINT64_C(0x417A5EF61F31C368) } }, //  1: vdPow ( 8.06958309145159269      , 8.20611629793705255       ) = ( 27651937.9496492445       );
{ { UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x4006DD9BBAC0EE6B) }, { UINT64_C(0x4081D7109BAA7980) } }, //  2: vdPow ( 9.21366148563738108      , 2.8582071867111174        ) = ( 570.883109409172903       );
{ { UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x401FBFBCBB737F7A) }, { UINT64_C(0x40F147A2E685447B) } }, //  3: vdPow ( 4.0834825258625127       , 7.93724339382594657       ) = ( 70778.1812794375437       );
 }

,

    .data_c32 = std::vector <data_3_c32_t>

    {

{ { UINT32_C(0x4093852F), UINT32_C(0x41093E24) }, { UINT32_C(0x41034C40), UINT32_C(0x41011D03) }, { UINT32_C(0xC623D50C), UINT32_C(0x4693B11C) } }, //  0: vcPow ( 4.61000776      + i * 8.57767105     , 8.20611572      + i * 8.06958294      ) = ( -10485.2617     + i * 18904.5547      );
{ { UINT32_C(0x4036ECDE), UINT32_C(0x41136B29) }, { UINT32_C(0x40FDFDE5), UINT32_C(0x4082ABE3) }, { UINT32_C(0x489D1870), UINT32_C(0x482626A1) } }, //  1: vcPow ( 2.85820723      + i * 9.21366215     , 7.93724298      + i * 4.08348227      ) = ( 321731.5        + i * 170138.516      );
{ { UINT32_C(0x40C0F87C), UINT32_C(0x40649ED8) }, { UINT32_C(0x40D64D6C), UINT32_C(0x40AB29A5) }, { UINT32_C(0x4566C04A), UINT32_C(0x46CBEE93) } }, //  2: vcPow ( 6.03033257      + i * 3.57219505     , 6.69695091      + i * 5.34883356      ) = ( 3692.01807      + i * 26103.2871      );
{ { UINT32_C(0x40B56AA4), UINT32_C(0x408B1733) }, { UINT32_C(0x411410F1), UINT32_C(0x41193290) }, { UINT32_C(0x480FE542), UINT32_C(0xC714EC50) } }, //  3: vcPow ( 5.66926765      + i * 4.34658194     , 9.25413609      + i * 9.57484436      ) = ( 147349.031      + i * -38124.3125     );
 }

,
    .data_c64 = std::vector <data_3_c64_t>

    {

{ { UINT64_C(0x401270A5F32DAE19), UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x40206988134D9FDD), UINT64_C(0x402023A0651C4741) }, { UINT64_C(0xC0C47AA470F16D42), UINT64_C(0x40D27624C7C04344) } }, //  0: vzPow ( 4.6100080486899282        + i * 8.57767068267691535      , 8.20611629793705255       + i * 8.06958309145159269       ) = ( -10485.2846967490659      + i * 18904.5746918351069       );
{ { UINT64_C(0x4006DD9BBAC0EE6B), UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x401FBFBCBB737F7A), UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x4113A30DBF54F1C6), UINT64_C(0x4104C4D62BEAA20C) } }, //  1: vzPow ( 2.8582071867111174        + i * 9.21366148563738108      , 7.93724339382594657       + i * 4.0834825258625127        ) = ( 321731.436847474775       + i * 170138.771443620673       );
{ { UINT64_C(0x40181F0F82C50AEC), UINT64_C(0x400C93DAEEF5F483) }, { UINT64_C(0x401AC9AD9555935C), UINT64_C(0x40156534AD76CA6A) }, { UINT64_C(0x40ACD804D2E3A03F), UINT64_C(0x40D97DD363DA339D) } }, //  2: vzPow ( 6.03033260657933212       + i * 3.57219492614837142      , 6.69695123038112783       + i * 5.34883376157322665       ) = ( 3692.00942145663385       + i * 26103.3029695037876       );
{ { UINT64_C(0x4016AD549DF3C110), UINT64_C(0x401162E657685F66) }, { UINT64_C(0x4022821E20A96AA3), UINT64_C(0x40232652067FE63E) }, { UINT64_C(0x4101FCA94042F62A), UINT64_C(0xC0E29D899EA112CB) } }, //  3: vzPow ( 5.66926810074097887       + i * 4.34658180784740544      , 9.25413610523293606       + i * 9.57484455405494472       ) = ( 147349.156377719075       + i * -38124.3006139151621      );
 }

};

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
  // Variable declaraions
  // *************************************************************
  // Input arguments
  A arg1;

  A arg2;

  A *varg1 = own_malloc<A>(main_queue, VLEN * sizeof(A));

  A *varg2 = own_malloc<A>(main_queue, VLEN * sizeof(A));

  // Output results
  R ref1;
  R *vres1 = own_malloc<R>(main_queue, VLEN * sizeof(R));
  R *vref1 = own_malloc<R>(main_queue, VLEN * sizeof(R));

  // Number of errors
  int errs = 0;
  // Number of printed errors
  int printed_errs = 0;

  // *************************************************************
  // Vector input data initialization
  // *************************************************************
  for (int i = 0; i < VLEN; ++i) {
    // Getting values from reference data table

    data.get_values(arg1, arg2, ref1);

    // Pushing values into vectors
    varg1[i] = arg1;
    vres1[i] = R(0);
    vref1[i] = ref1;

    varg2[i] = arg2;

  } // for (int i = 0; i < ACCURACY_LEN; ++i)

  // *************************************************************
  // Loop by all 3 accuracy modes of VM: HA, LA, EP:
  // set computation mode, run VM and check results
  // *************************************************************
  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    cl::sycl::vector_class<cl::sycl::event> no_deps;
    // Run VM function

    oneapi::mkl::vm::pow(main_queue, VLEN, varg1, varg2, vres1, no_deps,
                 accuracy_mode[acc]);
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
    for (int i = 0; i < VLEN; ++i) {
      // Compute ulp values for HA, LA and EP results
      R ulp = compute_ulp<R>(vres1[i], vref1[i]);

      std::cout << "\t\t" << accuracy_name[acc] << " vm::"
                << "pow"
                << "(" << varg1[i]

                << "," << varg2[i] << ") = " << vres1[i] << std::endl;

      // Check ulp
      if (check_ulp(ulp, false, accuracy_mode[acc], ulp_table) == false) {
        errs++;
      }
    } // for (int i = 0; i < varg1.size (); ++i)
  }   // for (int acc = 0; acc < 3; ++acc)

  std::cout << "\tResult: " << ((errs == 0) ? "PASS" : "FAIL") << std::endl;

  own_free<A>(main_queue, varg1);

  own_free<A>(main_queue, varg2);

  own_free<R>(main_queue, vres1);
  own_free<R>(main_queue, vref1);

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
            << "pow"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "pow" << std::endl;
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
      errs += vPowAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vPowAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vPowAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vPowAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
