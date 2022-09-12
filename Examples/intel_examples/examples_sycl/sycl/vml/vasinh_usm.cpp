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
 *            oneapi::mkl::vm::asinh example program text (SYCL interface)
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
bool vAsinhAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x4106F102) }, { UINT32_C(0x40350CA0) } }, //  0: vsAsinh ( 8.4338398       ) = ( 2.82889557      );
{ { UINT32_C(0x40821418) }, { UINT32_C(0x40070FEC) } }, //  1: vsAsinh ( 4.06495285      ) = ( 2.11034679      );
{ { UINT32_C(0x40FBFADC) }, { UINT32_C(0x4030B06E) } }, //  2: vsAsinh ( 7.87437248      ) = ( 2.76076841      );
{ { UINT32_C(0x41006539) }, { UINT32_C(0x4031E3DE) } }, //  3: vsAsinh ( 8.02471256      ) = ( 2.77953291      );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x4020DE203A4CEF74) }, { UINT64_C(0x4006A19409F6875D) } }, //  0: vdAsinh ( 8.43383962811615362       ) = ( 2.82889564307781294       );
{ { UINT64_C(0x40104282F4C21EA0) }, { UINT64_C(0x4000E1FD71214ECF) } }, //  1: vdAsinh ( 4.06495268282711208       ) = ( 2.11034668333909492       );
{ { UINT64_C(0x401F7F5B79FEFEB7) }, { UINT64_C(0x4006160DBBDDDD1B) } }, //  2: vdAsinh ( 7.87437239283433765       ) = ( 2.7607683827478815        );
{ { UINT64_C(0x40200CA71821B2E8) }, { UINT64_C(0x40063C7BB36D8924) } }, //  3: vdAsinh ( 8.02471232806551882       ) = ( 2.77953281572367139       );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0x40821418), UINT32_C(0x4106F102) }, { UINT32_C(0x403B657B), UINT32_C(0x3F8F494C) } }, //  0: vcAsinh ( 4.06495285      + i * 8.4338398       ) = ( 2.92806888      + i * 1.11942434      );
{ { UINT32_C(0x41006539), UINT32_C(0x40FBFADC) }, { UINT32_C(0x40473A22), UINT32_C(0x3F462298) } }, //  1: vcAsinh ( 8.02471256      + i * 7.87437248      ) = ( 3.11292315      + i * 0.773965359     );
{ { UINT32_C(0x4008B448), UINT32_C(0x41122574) }, { UINT32_C(0x403B7891), UINT32_C(0x3FAB7EC8) } }, //  2: vcAsinh ( 2.13600349      + i * 9.13414383      ) = ( 2.92923379      + i * 1.33980656      );
{ { UINT32_C(0x40F7511A), UINT32_C(0x405F0D3D) }, { UINT32_C(0x40354EE7), UINT32_C(0x3ED793C5) } }, //  3: vcAsinh ( 7.72865009      + i * 3.485183        ) = ( 2.83294082      + i * 0.421049267     );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0x40104282F4C21EA0), UINT64_C(0x4020DE203A4CEF74) }, { UINT64_C(0x40076CAF627D52BE), UINT64_C(0x3FF1E929808BEB84) } }, //  0: vzAsinh ( 4.06495268282711208       + i * 8.43383962811615362       ) = ( 2.92806889481502619       + i * 1.11942434514523459       );
{ { UINT64_C(0x40200CA71821B2E8), UINT64_C(0x401F7F5B79FEFEB7) }, { UINT64_C(0x4008E7443A01B5E4), UINT64_C(0x3FE8C452F5B6F581) } }, //  1: vzAsinh ( 8.02471232806551882       + i * 7.87437239283433765       ) = ( 3.11292310064048827       + i * 0.773965339576236144      );
{ { UINT64_C(0x40011688F5E89379), UINT64_C(0x402244AE8957BC90) }, { UINT64_C(0x40076F1232EA82C3), UINT64_C(0x3FF56FD905B7041C) } }, //  2: vzAsinh ( 2.13600341907516311       + i * 9.13414410778048591       ) = ( 2.92923393037958268       + i * 1.33980657799134573       );
{ { UINT64_C(0x401EEA233BB5D447), UINT64_C(0x400BE1A7A0BAF683) }, { UINT64_C(0x4006A9DCD22D5BB3), UINT64_C(0x3FDAF278A674CEAD) } }, //  3: vzAsinh ( 7.72865002915666022       + i * 3.4851830060059128        ) = ( 2.83294071389124147       + i * 0.421049273066482155      );
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
    oneapi::mkl::vm::asinh(main_queue, VLEN, varg1, vres1, no_deps, accuracy_mode[acc]);

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
                << "asinh"
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
            << "asinh"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "asinh" << std::endl;
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
      errs += vAsinhAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vAsinhAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vAsinhAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs +=
          vAsinhAccuracyLiteTest<std::complex<double>, std::complex<double>>(
              my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
