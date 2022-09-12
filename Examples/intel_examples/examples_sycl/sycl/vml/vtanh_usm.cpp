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
 *            oneapi::mkl::vm::tanh example program text (SYCL interface)
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
bool vTanhAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x40D9B85C) }, { UINT32_C(0x3F7FFFD7) } }, //  0: vsTanh ( 6.80375481      ) = ( 0.999997556     );
{ { UINT32_C(0xC007309A) }, { UINT32_C(0xBF789E12) } }, //  1: vsTanh ( -2.1123414      ) = ( -0.971161962    );
{ { UINT32_C(0x40B52EFA) }, { UINT32_C(0x3F7FFE6B) } }, //  2: vsTanh ( 5.66198444      ) = ( 0.99997586      );
{ { UINT32_C(0x40BF006A) }, { UINT32_C(0x3F7FFF25) } }, //  3: vsTanh ( 5.96880054      ) = ( 0.999986947     );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x3FEFFFFAD5FE7BF2) } }, //  0: vdTanh ( 6.80375434309419092       ) = ( 0.999997537572083539      );
{ { UINT64_C(0xC000E6134801CC26) }, { UINT64_C(0xBFEF13C23C2086E3) } }, //  1: vdTanh ( -2.11234146361813924      ) = ( -0.971161954341564715     );
{ { UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x3FEFFFCD557AD24F) } }, //  2: vdTanh ( 5.66198447517211711       ) = ( 0.999975840523413484      );
{ { UINT64_C(0x4017E00D485FC01A) }, { UINT64_C(0x3FEFFFE491F87BD2) } }, //  3: vdTanh ( 5.96880066952146571       ) = ( 0.999986920451073624      );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0xC007309A), UINT32_C(0x40D9B85C) }, { UINT32_C(0xBF7C29D8), UINT32_C(0x3CCBCC3D) } }, //  0: vcTanh ( -2.1123414      + i * 6.80375481      ) = ( -0.985013485    + i * 0.0248776618    );
{ { UINT32_C(0x40BF006A), UINT32_C(0x40B52EFA) }, { UINT32_C(0x3F7FFFB9), UINT32_C(0xB74FB664) } }, //  1: vcTanh ( 5.96880054      + i * 5.66198444      ) = ( 0.999995768     + i * -1.23806276e-05 );
{ { UINT32_C(0xC0C1912F), UINT32_C(0x4103BA28) }, { UINT32_C(0xBF800044), UINT32_C(0xB7008037) } }, //  2: vcTanh ( -6.04897261     + i * 8.2329483       ) = ( -1.00000811     + i * -7.65924688e-06 );
{ { UINT32_C(0x40ABAABC), UINT32_C(0xC052EA36) }, { UINT32_C(0x3F7FFD44), UINT32_C(0xB75EA892) } }, //  3: vcTanh ( 5.3645916       + i * -3.2955451      ) = ( 0.999958277     + i * -1.32714795e-05 );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0xC000E6134801CC26), UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0xBFEF853AE5BB4742), UINT64_C(0x3F997986725FA047) } }, //  0: vzTanh ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -0.985013436026044298     + i * 0.0248776442821567988     );
{ { UINT64_C(0x4017E00D485FC01A), UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x3FEFFFF7272D333E), UINT64_C(0xBEE9F6CBFCEC114D) } }, //  1: vzTanh ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 0.999995781437611475      + i * -1.23806238696641625e-05  );
{ { UINT64_C(0xC0183225E080644C), UINT64_C(0x40207744D998EE8A) }, { UINT64_C(0xBFF000087C283048), UINT64_C(0xBEE0100456BD42E8) } }, //  2: vzTanh ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -1.00000809191534934      + i * -7.65922842274739419e-06  );
{ { UINT64_C(0x4015755793FAEAB0), UINT64_C(0xC00A5D46A314BA8E) }, { UINT64_C(0x3FEFFFA87AF27428), UINT64_C(0xBEEBD50EB5E0D319) } }, //  3: vzTanh ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 0.999958267336869433      + i * -1.32714537209671191e-05  );
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
    oneapi::mkl::vm::tanh(main_queue, VLEN, varg1, vres1, no_deps, accuracy_mode[acc]);

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
                << "tanh"
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
            << "tanh"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "tanh" << std::endl;
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
      errs += vTanhAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vTanhAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vTanhAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vTanhAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}