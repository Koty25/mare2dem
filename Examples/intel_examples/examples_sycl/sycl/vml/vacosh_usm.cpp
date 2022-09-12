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
 *            oneapi::mkl::vm::acosh example program text (SYCL interface)
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
bool vAcoshAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x41093E24) }, { UINT32_C(0x4035B072) } }, //  0: vsAcosh ( 8.57767105      ) = ( 2.83889437      );
{ { UINT32_C(0x4093852F) }, { UINT32_C(0x400D66CF) } }, //  1: vsAcosh ( 4.61000776      ) = ( 2.20939994      );
{ { UINT32_C(0x41011D03) }, { UINT32_C(0x4031C0B8) } }, //  2: vsAcosh ( 8.06958294      ) = ( 2.77738762      );
{ { UINT32_C(0x41034C40) }, { UINT32_C(0x4032D5B5) } }, //  3: vsAcosh ( 8.20611572      ) = ( 2.79429364      );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x4006B60E36BAC01A) } }, //  0: vdAcosh ( 8.57767068267691535       ) = ( 2.83889429814736349       );
{ { UINT64_C(0x401270A5F32DAE19) }, { UINT64_C(0x4001ACD9F73AF77A) } }, //  1: vdAcosh ( 4.6100080486899282        ) = ( 2.20940011166288475       );
{ { UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x40063816F3916B18) } }, //  2: vdAcosh ( 8.06958309145159269       ) = ( 2.77738752639323749       );
{ { UINT64_C(0x40206988134D9FDD) }, { UINT64_C(0x40065AB69C81EC96) } }, //  3: vdAcosh ( 8.20611629793705255       ) = ( 2.79429361602303583       );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0x4093852F), UINT32_C(0x41093E24) }, { UINT32_C(0x403E1EFD), UINT32_C(0x3F8A3801) } }, //  0: vcAcosh ( 4.61000776      + i * 8.57767105      ) = ( 2.97064137      + i * 1.0798341       );
{ { UINT32_C(0x41034C40), UINT32_C(0x41011D03) }, { UINT32_C(0x4048B869), UINT32_C(0x3F4765C9) } }, //  1: vcAcosh ( 8.20611572      + i * 8.06958294      ) = ( 3.1362555       + i * 0.778896868     );
{ { UINT32_C(0x4036ECDE), UINT32_C(0x41136B29) }, { UINT32_C(0x403D912A), UINT32_C(0x3FA2C0B4) } }, //  2: vcAcosh ( 2.85820723      + i * 9.21366215      ) = ( 2.96198511      + i * 1.27150583      );
{ { UINT32_C(0x40FDFDE5), UINT32_C(0x4082ABE3) }, { UINT32_C(0x403856E4), UINT32_C(0x3EF49846) } }, //  3: vcAcosh ( 7.93724298      + i * 4.08348227      ) = ( 2.88030338      + i * 0.477724254     );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0x401270A5F32DAE19), UINT64_C(0x402127C473A3E923) }, { UINT64_C(0x4007C3DFA89040E4), UINT64_C(0x3FF147001AC5D719) } }, //  0: vzAcosh ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 2.97064143839098627       + i * 1.07983408411150195       );
{ { UINT64_C(0x40206988134D9FDD), UINT64_C(0x402023A0651C4741) }, { UINT64_C(0x4009170D2ED1FCE7), UINT64_C(0x3FE8ECB91A991BB9) } }, //  1: vzAcosh ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 3.13625561312038625       + i * 0.778896858167050898      );
{ { UINT64_C(0x4006DD9BBAC0EE6B), UINT64_C(0x40226D6509CA7464) }, { UINT64_C(0x4007B22542341397), UINT64_C(0x3FF458166D8F9BF1) } }, //  2: vzAcosh ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 2.96198512765335975       + i * 1.27150576398139159       );
{ { UINT64_C(0x401FBFBCBB737F7A), UINT64_C(0x4010557C717977C6) }, { UINT64_C(0x40070ADC8F07EE90), UINT64_C(0x3FDE9308BD0B88AC) } }, //  3: vzAcosh ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 2.88030349486309234       + i * 0.477724251379309406      );
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

  A *varg1 = own_malloc<A>(main_queue, VLEN * sizeof(A));

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

    data.get_values(arg1, ref1);

    // Pushing values into vectors
    varg1[i] = arg1;
    vres1[i] = R(0);
    vref1[i] = ref1;

  } // for (int i = 0; i < ACCURACY_LEN; ++i)

  // *************************************************************
  // Loop by all 3 accuracy modes of VM: HA, LA, EP:
  // set computation mode, run VM and check results
  // *************************************************************
  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    cl::sycl::vector_class<cl::sycl::event> no_deps;
    // Run VM function
    oneapi::mkl::vm::acosh(main_queue, VLEN, varg1, vres1, no_deps, accuracy_mode[acc]);

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
                << "acosh"
                << "(" << varg1[i] << ") = " << vres1[i]

                << std::endl;

      // Check ulp
      if (check_ulp(ulp, false, accuracy_mode[acc], ulp_table) == false) {
        errs++;
      }
    } // for (int i = 0; i < varg1.size (); ++i)
  }   // for (int acc = 0; acc < 3; ++acc)

  std::cout << "\tResult: " << ((errs == 0) ? "PASS" : "FAIL") << std::endl;

  own_free<A>(main_queue, varg1);

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
            << "acosh"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "acosh" << std::endl;
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
      errs += vAcoshAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vAcoshAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vAcoshAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs +=
          vAcoshAccuracyLiteTest<std::complex<double>, std::complex<double>>(
              my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}