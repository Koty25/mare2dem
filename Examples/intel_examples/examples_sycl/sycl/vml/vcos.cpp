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
 *            oneapi::mkl::vm::cos example program text (SYCL interface)
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
    {MAX_HA_ULP_S, 4.5}, {MAX_LA_ULP_S, 5.0}, {MAX_EP_ULP_S, 5.0E3},
    {MAX_HA_ULP_D, 2.0}, {MAX_LA_ULP_D, 5.0}, {MAX_EP_ULP_D, 7.0E7},
    {MAX_HA_ULP_C, 4.0}, {MAX_LA_ULP_C, 5.0}, {MAX_EP_ULP_C, 5.0E3},
    {MAX_HA_ULP_Z, 4.0}, {MAX_LA_ULP_Z, 5.0}, {MAX_EP_ULP_Z, 7.0E7},
};

//!
//! @brief Accuracy test
//!
template <typename A, typename R>
bool vCosAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x40D9B85C) }, { UINT32_C(0x3F5E16D8) } }, //  0: vsCos ( 6.80375481      ) = ( 0.867536068     );
{ { UINT32_C(0xC007309A) }, { UINT32_C(0xBF03F53A) } }, //  1: vsCos ( -2.1123414      ) = ( -0.51546061     );
{ { UINT32_C(0x40B52EFA) }, { UINT32_C(0x3F502C93) } }, //  2: vsCos ( 5.66198444      ) = ( 0.813180149     );
{ { UINT32_C(0x40BF006A) }, { UINT32_C(0x3F7373DF) } }, //  3: vsCos ( 5.96880054      ) = ( 0.950986803     );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x3FEBC2DB7AB89950) } }, //  0: vdCos ( 6.80375434309419092       ) = ( 0.867536296548488295      );
{ { UINT64_C(0xC000E6134801CC26) }, { UINT64_C(0xBFE07EA757C4010B) } }, //  1: vdCos ( -2.11234146361813924      ) = ( -0.51546065465666524      );
{ { UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x3FEA05925DBF776B) } }, //  2: vdCos ( 5.66198447517211711       ) = ( 0.813180144406698502      );
{ { UINT64_C(0x4017E00D485FC01A) }, { UINT64_C(0x3FEE6E7BF8882000) } }, //  3: vdCos ( 5.96880066952146571       ) = ( 0.950986848271895724      );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0xC007309A), UINT32_C(0x40D9B85C) }, { UINT32_C(0xC36845F3), UINT32_C(0x43C11152) } }, //  0: vcCos ( -2.1123414      + i * 6.80375481      ) = ( -232.273239     + i * 386.135315      );
{ { UINT32_C(0x40BF006A), UINT32_C(0x40B52EFA) }, { UINT32_C(0x4308CF67), UINT32_C(0x4231F100) } }, //  1: vcCos ( 5.96880054      + i * 5.66198444      ) = ( 136.810165      + i * 44.4853516      );
{ { UINT32_C(0xC0C1912F), UINT32_C(0x4103BA28) }, { UINT32_C(0x44E4C2CB), UINT32_C(0xC3DA5250) } }, //  2: vcCos ( -6.04897261     + i * 8.2329483       ) = ( 1830.08728      + i * -436.643066     );
{ { UINT32_C(0x40ABAABC), UINT32_C(0xC052EA36) }, { UINT32_C(0x41033D87), UINT32_C(0xC12B6150) } }, //  3: vcCos ( 5.3645916       + i * -3.2955451      ) = ( 8.20252132      + i * -10.7112579     );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0xC000E6134801CC26), UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0xC06D08BDB8FC7F5D), UINT64_C(0x407822296A7AAF2C) } }, //  0: vzCos ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -232.273159497412877      + i * 386.135111312137042       );
{ { UINT64_C(0x4017E00D485FC01A), UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x406119ECE50B1B29), UINT64_C(0x40463E1F6EEA30C8) } }, //  1: vzCos ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 136.810167810145941       + i * 44.4853342669971994       );
{ { UINT64_C(0xC0183225E080644C), UINT64_C(0x40207744D998EE8A) }, { UINT64_C(0x409C98572C8DB5F8), UINT64_C(0xC07B4A47E4EA8F4F) } }, //  2: vzCos ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 1830.08513089583539       + i * -436.642552295922485      );
{ { UINT64_C(0x4015755793FAEAB0), UINT64_C(0xC00A5D46A314BA8E) }, { UINT64_C(0x402067B107A38AEA), UINT64_C(0xC0256C29633824E0) } }, //  3: vzCos ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 8.20252250548715622       + i * -10.7112532621417245      );
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
    oneapi::mkl::vm::cos(main_queue, varg1.size(), in1, out1, accuracy_mode[acc]);

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
                << "cos"
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
            << "cos"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "cos" << std::endl;
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
      errs += vCosAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vCosAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vCosAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vCosAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
