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
 *            oneapi::mkl::vm::atanh example program text (SYCL interface)
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
bool vAtanhAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x3C82EB10) }, { UINT32_C(0x3C82EDEB) } }, //  0: vsAtanh ( 0.0159812272    ) = ( 0.0159825888    );
{ { UINT32_C(0x3D780F8E) }, { UINT32_C(0x3D785D5D) } }, //  1: vsAtanh ( 0.0605617091    ) = ( 0.0606359132    );
{ { UINT32_C(0x3CB1AF64) }, { UINT32_C(0x3CB1B687) } }, //  2: vsAtanh ( 0.0216900781    ) = ( 0.0216934811    );
{ { UINT32_C(0x3CA51E30) }, { UINT32_C(0x3CA523EA) } }, //  3: vsAtanh ( 0.0201559961    ) = ( 0.0201587267    );
 }

,
    .data_f64 = std::vector<data_2_f64_t>

    {

{ { UINT64_C(0x3F905D621353EDF8) }, { UINT64_C(0x3F905DBD64B240D8) } }, //  0: vdAtanh ( 0.0159812282845290532     ) = ( 0.0159825890264649606     );
{ { UINT64_C(0x3FAF01F1B0A46A4A) }, { UINT64_C(0x3FAF0BAB94755F2D) } }, //  1: vdAtanh ( 0.060561707318090699      ) = ( 0.0606359118198141825     );
{ { UINT64_C(0x3F9635EC782C6BD8) }, { UINT64_C(0x3F9636D0CCDE1024) } }, //  2: vdAtanh ( 0.0216900776241394089     ) = ( 0.0216934800187261884     );
{ { UINT64_C(0x3F94A3C609C2E128) }, { UINT64_C(0x3F94A47D42904F8C) } }, //  3: vdAtanh ( 0.020155996652392677      ) = ( 0.0201587268712298123     );
 }

,

    .data_c32 = std::vector <data_2_c32_t>

    {

{ { UINT32_C(0x3D780F8E), UINT32_C(0x3C82EB10) }, { UINT32_C(0x3D784D08), UINT32_C(0x3C836386) } }, //  0: vcAtanh ( 0.0605617091    + i * 0.0159812272    ) = ( 0.0606203377    + i * 0.0160386674    );
{ { UINT32_C(0x3CA51E30), UINT32_C(0x3CB1AF64) }, { UINT32_C(0x3CA51005), UINT32_C(0x3CB1BABB) } }, //  1: vcAtanh ( 0.0201559961    + i * 0.0216900781    ) = ( 0.0201492403    + i * 0.0216954853    );
{ { UINT32_C(0x3DA4576B), UINT32_C(0x3C10C1C8) }, { UINT32_C(0x3DA4AEBF), UINT32_C(0x3C11B0F3) } }, //  2: vcAtanh ( 0.0802448615    + i * 0.00883526355   ) = ( 0.0804114267    + i * 0.00889228564   );
{ { UINT32_C(0x3CBDDDC8), UINT32_C(0x3D88257A) }, { UINT32_C(0x3CBD1067), UINT32_C(0x3D8804D6) } }, //  3: vcAtanh ( 0.0231770426    + i * 0.0664777309    ) = ( 0.0230791103    + i * 0.0664154738    );
 }

,
    .data_c64 = std::vector <data_2_c64_t>

    {

{ { UINT64_C(0x3FAF01F1B0A46A4A), UINT64_C(0x3F905D621353EDF8) }, { UINT64_C(0x3FAF09A0E2B69FE1), UINT64_C(0x3F906C70CDA5C815) } }, //  0: vzAtanh ( 0.060561707318090699      + i * 0.0159812282845290532     ) = ( 0.060620334315274034      + i * 0.0160386682050060632     );
{ { UINT64_C(0x3F94A3C609C2E128), UINT64_C(0x3F9635EC782C6BD8) }, { UINT64_C(0x3F94A200AD66E706), UINT64_C(0x3F9637574ED5D99F) } }, //  1: vzAtanh ( 0.020155996652392677      + i * 0.0216900776241394089     ) = ( 0.020149241050353893      + i * 0.0216954843394546702     );
{ { UINT64_C(0x3FB48AED668F7C42), UINT64_C(0x3F821839168A96D8) }, { UINT64_C(0x3FB495D7D955D65A), UINT64_C(0x3F82361E797C0E88) } }, //  2: vzAtanh ( 0.0802448630706616151     + i * 0.00883526420632156639    ) = ( 0.0804114251712574613     + i * 0.00889228637925688903    );
{ { UINT64_C(0x3F97BBB8DC2F7770), UINT64_C(0x3FB104AF24553C92) }, { UINT64_C(0x3F97A20CB1D99626), UINT64_C(0x3FB1009AB2439241) } }, //  3: vzAtanh ( 0.0231770405188095885     + i * 0.0664777244285111035     ) = ( 0.023079107623195004      + i * 0.0664154706206057238     );
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
    oneapi::mkl::vm::atanh(main_queue, VLEN, varg1, vres1, no_deps, accuracy_mode[acc]);

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
                << "atanh"
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
            << "atanh"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "atanh" << std::endl;
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
      errs += vAtanhAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vAtanhAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vAtanhAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs +=
          vAtanhAccuracyLiteTest<std::complex<double>, std::complex<double>>(
              my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
