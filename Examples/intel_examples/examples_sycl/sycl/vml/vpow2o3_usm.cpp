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
 *            oneapi::mkl::vm::pow2o3 example program text (SYCL interface)
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
bool vPow2o3AccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x40D9B85C) }, { UINT32_C(0x4065CC38) } }, //  0: vsPow2o3 ( 6.80375481      ) = ( 3.59058952      );
{ { UINT32_C(0xC007309A) }, { UINT32_C(0x3FD2BA03) } }, //  1: vsPow2o3 ( -2.1123414      ) = ( 1.64630163      );
{ { UINT32_C(0x40B52EFA) }, { UINT32_C(0x404B4F67) } }, //  2: vsPow2o3 ( 5.66198444      ) = ( 3.17672133      );
{ { UINT32_C(0x40BF006A) }, { UINT32_C(0x40529713) } }, //  3: vsPow2o3 ( 5.96880054      ) = ( 3.29047084      );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x400CB986F5486496) } }, //  0: vdPow2o3 ( 6.80375434309419092       ) = ( 3.59058944346595421       );
{ { UINT64_C(0xC000E6134801CC26) }, { UINT64_C(0x3FFA574068B5519B) } }, //  1: vdPow2o3 ( -2.11234146361813924      ) = ( 1.64630165959997643       );
{ { UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x400969ECDC82CD04) } }, //  2: vdPow2o3 ( 5.66198447517211711       ) = ( 3.17672130846188061       );
{ { UINT64_C(0x4017E00D485FC01A) }, { UINT64_C(0x400A52E2644589BD) } }, //  3: vdPow2o3 ( 5.96880066952146571       ) = ( 3.29047087037289865       );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    { /* empty */ }

,
    .data_c64 = std::vector <data_2_c64_t>

    { /* empty */ }

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
    oneapi::mkl::vm::pow2o3(main_queue, VLEN, varg1, vres1, no_deps,
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
                << "pow2o3"
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
            << "pow2o3"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "pow2o3" << std::endl;
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
      errs += vPow2o3AccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vPow2o3AccuracyLiteTest<double, double>(my_dev);
    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
