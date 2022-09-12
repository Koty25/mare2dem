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
 *            oneapi::mkl::vm::tan example program text (SYCL interface)
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
bool vTanAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x40D9B85C) }, { UINT32_C(0x3F12C4FD) } }, //  0: vsTan ( 6.80375481      ) = ( 0.573318303     );
{ { UINT32_C(0xC007309A) }, { UINT32_C(0x3FD4CA42) } }, //  1: vsTan ( -2.1123414      ) = ( 1.66242242      );
{ { UINT32_C(0x40B52EFA) }, { UINT32_C(0xBF3739A7) } }, //  2: vsTan ( 5.66198444      ) = ( -0.715723455    );
{ { UINT32_C(0x40BF006A) }, { UINT32_C(0xBEA67C8E) } }, //  3: vsTan ( 5.96880054      ) = ( -0.325169027    );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x3FE2589E45D858D9) } }, //  0: vdTan ( 6.80375434309419092       ) = ( 0.573317657867643438      );
{ { UINT64_C(0xC000E6134801CC26) }, { UINT64_C(0x3FFA99480B7BBB50) } }, //  1: vdTan ( -2.11234146361813924      ) = ( 1.66242222295450759       );
{ { UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0xBFE6E734CC578392) } }, //  2: vdTan ( 5.66198447517211711       ) = ( -0.715723418336084771     );
{ { UINT64_C(0x4017E00D485FC01A) }, { UINT64_C(0xBFD4CF91224985AD) } }, //  3: vdTan ( 5.96880066952146571       ) = ( -0.325168879970159364     );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0xC007309A), UINT32_C(0x40D9B85C) }, { UINT32_C(0x3611FC00), UINT32_C(0x3F80000A) } }, //  0: vcTan ( -2.1123414      + i * 6.80375481      ) = ( 2.1753367e-06   + i * 1.00000119      );
{ { UINT32_C(0x40BF006A), UINT32_C(0x40B52EFA) }, { UINT32_C(0xB76E6472), UINT32_C(0x3F7FFEB8) } }, //  1: vcTan ( 5.96880054      + i * 5.66198444      ) = ( -1.42092922e-05 + i * 0.99998045      );
{ { UINT32_C(0xC0C1912F), UINT32_C(0x4103BA28) }, { UINT32_C(0x3388F265), UINT32_C(0x3F7FFFFE) } }, //  2: vcTan ( -6.04897261     + i * 8.2329483       ) = ( 6.37708482e-08  + i * 0.999999881     );
{ { UINT32_C(0x40ABAABC), UINT32_C(0xC052EA36) }, { UINT32_C(0xBB2DAE75), UINT32_C(0xBF801793) } }, //  3: vcTan ( 5.3645916       + i * -3.2955451      ) = ( -0.00265016896  + i * -1.00071943     );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0xC000E6134801CC26), UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x3EC23F813BAC1EE3), UINT64_C(0x3FF0000135BF0F45) } }, //  0: vzTan ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( 2.17533894664502934e-06   + i * 1.00000115389498601       );
{ { UINT64_C(0x4017E00D485FC01A), UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0xBEEDCC8D7ADCC53E), UINT64_C(0x3FEFFFD705FCDE2F) } }, //  1: vzTan ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( -1.42092866003163373e-05  + i * 0.999980460829595574      );
{ { UINT64_C(0xC0183225E080644C), UINT64_C(0x40207744D998EE8A) }, { UINT64_C(0x3E711E4F35CD0BC8), UINT64_C(0x3FEFFFFFBC562861) } }, //  2: vzTan ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 6.37709951070066822e-08   + i * 0.999999873967009845      );
{ { UINT64_C(0x4015755793FAEAB0), UINT64_C(0xC00A5D46A314BA8E) }, { UINT64_C(0xBF65B5CF78A52334), UINT64_C(0xBFF002F257CD6E64) } }, //  3: vzTan ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -0.00265017053348906052   + i * -1.00071939752424388      );
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
    oneapi::mkl::vm::tan(main_queue, VLEN, varg1, vres1, no_deps, accuracy_mode[acc]);

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
                << "tan"
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
            << "tan"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "tan" << std::endl;
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
      errs += vTanAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vTanAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vTanAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vTanAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
