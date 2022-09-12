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
 *            oneapi::mkl::vm::sub example program text (SYCL interface)
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
    {MAX_HA_ULP_C, 2.0}, {MAX_LA_ULP_C, 5.0}, {MAX_EP_ULP_C, 5.0E3},
    {MAX_HA_ULP_Z, 2.0}, {MAX_LA_ULP_Z, 5.0}, {MAX_EP_ULP_Z, 7.0E7},
};

//!
//! @brief Accuracy test
//!
template <typename A, typename R>
bool vSubAccuracyLiteTest(const cl::sycl::device &dev) {
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

{ { UINT32_C(0x40D9B85C) }, { UINT32_C(0xC007309A) }, { UINT32_C(0x410EA854) } }, //  0: vsSub ( 6.80375481     , -2.1123414      ) = ( 8.91609573      );
{ { UINT32_C(0x40B52EFA) }, { UINT32_C(0x40BF006A) }, { UINT32_C(0xBE9D1700) } }, //  1: vsSub ( 5.66198444     , 5.96880054      ) = ( -0.306816101    );
{ { UINT32_C(0x4103BA28) }, { UINT32_C(0xC0C1912F) }, { UINT32_C(0x416482C0) } }, //  2: vsSub ( 8.2329483      , -6.04897261     ) = ( 14.2819214      );
{ { UINT32_C(0xC052EA36) }, { UINT32_C(0x40ABAABC) }, { UINT32_C(0xC10A8FEC) } }, //  3: vsSub ( -3.2955451     , 5.3645916       ) = ( -8.66013718     );
 }

,
    .data_f64 = std::vector<data_3_f64_t>

    {

{ { UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0xC000E6134801CC26) }, { UINT64_C(0x4021D50A8273AA16) } }, //  0: vdSub ( 6.80375434309419092      , -2.11234146361813924      ) = ( 8.91609580671233104       );
{ { UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0x4017E00D485FC01A) }, { UINT64_C(0xBFD3A2E0642745C0) } }, //  1: vdSub ( 5.66198447517211711      , 5.96880066952146571       ) = ( -0.306816194349348592     );
{ { UINT64_C(0x40207744D998EE8A) }, { UINT64_C(0xC0183225E080644C) }, { UINT64_C(0x402C9057C9D920B0) } }, //  2: vdSub ( 8.23294715873568705      , -6.04897261413232101      ) = ( 14.2819197728680081       );
{ { UINT64_C(0xC00A5D46A314BA8E) }, { UINT64_C(0x4015755793FAEAB0) }, { UINT64_C(0xC02151FD72C2A3FC) } }, //  3: vdSub ( -3.2955448857022196      , 5.36459189623808186       ) = ( -8.66013678194030234      );
 }

,

    .data_c32 = std::vector <data_3_c32_t>

    {

{ { UINT32_C(0xC007309A), UINT32_C(0x40D9B85C) }, { UINT32_C(0x40BF006A), UINT32_C(0x40B52EFA) }, { UINT32_C(0xC1014C5C), UINT32_C(0x3F922588) } }, //  0: vcSub ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( -8.08114243     + i * 1.14177036      );
{ { UINT32_C(0xC0C1912F), UINT32_C(0x4103BA28) }, { UINT32_C(0x40ABAABC), UINT32_C(0xC052EA36) }, { UINT32_C(0xC1369DF6), UINT32_C(0x413874B6) } }, //  1: vcSub ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -11.4135647     + i * 11.5284939      );
{ { UINT32_C(0x3F8A29C0), UINT32_C(0xC08E3964) }, { UINT32_C(0x4024F46C), UINT32_C(0xBEE77440) }, { UINT32_C(0xBFBFBF18), UINT32_C(0xC07F8440) } }, //  2: vcSub ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( -1.49801922     + i * -3.9924469      );
{ { UINT32_C(0x3E8939C0), UINT32_C(0xC02D136C) }, { UINT32_C(0x41052EB4), UINT32_C(0x4110B6A8) }, { UINT32_C(0xC100E4E6), UINT32_C(0xC13BFB83) } }, //  3: vcSub ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( -8.05588341     + i * -11.7489042     );
 }

,
    .data_c64 = std::vector <data_3_c64_t>

    {

{ { UINT64_C(0xC000E6134801CC26), UINT64_C(0x401B370B60E66E18) }, { UINT64_C(0x4017E00D485FC01A), UINT64_C(0x4016A5DF421D4BBE) }, { UINT64_C(0xC020298B76305316), UINT64_C(0x3FF244B07B248968) } }, //  0: vzSub ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( -8.08114213313960406      + i * 1.1417698679220738        );
{ { UINT64_C(0xC0183225E080644C), UINT64_C(0x40207744D998EE8A) }, { UINT64_C(0x4015755793FAEAB0), UINT64_C(0xC00A5D46A314BA8E) }, { UINT64_C(0xC026D3BEBA3DA77E), UINT64_C(0x40270E96825E1D2E) } }, //  1: vzSub ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -11.4135645103704029      + i * 11.5284920444379075       );
{ { UINT64_C(0x3FF1453801E28A70), UINT64_C(0xC011C72C86338E59) }, { UINT64_C(0x40049E8D96893D1C), UINT64_C(0xBFDCEE88B739DD20) }, { UINT64_C(0xBFF7F7E32B2FEFC8), UINT64_C(0xC00FF087F57FE10E) } }, //  2: vzSub ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( -1.49801937932987705      + i * -3.99244682117944816      );
{ { UINT64_C(0x3FD12735D3224E60), UINT64_C(0xC005A26D910B44DC) }, { UINT64_C(0x4020A5D666294BAC), UINT64_C(0x402216D5173C2DAA) }, { UINT64_C(0xC0201C9CB7903939), UINT64_C(0xC0277F707B7EFEE1) } }, //  3: vzSub ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( -8.05588315616170014      + i * -11.7489050476573897      );
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

    oneapi::mkl::vm::sub(main_queue, VLEN, varg1, varg2, vres1, no_deps,
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
                << "sub"
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
            << "sub"
            << " Function Example: " << std::endl;
  std::cout << "# " << std::endl;
  std::cout << "# Using apis:" << std::endl;
  std::cout << "#   vm::"
            << "sub" << std::endl;
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
      errs += vSubAccuracyLiteTest<float, float>(my_dev);

      std::cout << "\tRunning with double precision real data type:"
                << std::endl;
      errs += vSubAccuracyLiteTest<double, double>(my_dev);

      std::cout << "\tRunning with single precision complex data type:"
                << std::endl;
      errs += vSubAccuracyLiteTest<std::complex<float>, std::complex<float>>(
          my_dev);

      std::cout << "\tRunning with double precision complex data type:"
                << std::endl;
      errs += vSubAccuracyLiteTest<std::complex<double>, std::complex<double>>(
          my_dev);

    } else {

      std::cout << "No " << sycl_device_names[dev_type]
                << " devices found; skipping " << sycl_device_names[dev_type]
                << " tests.\n";
    }
  }

  return 0;
}
