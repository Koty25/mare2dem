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
 *            Log2 example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 2.0, 5.0, 5.0E3, }, // VM_COMPLEX8
    { 2.0, 5.0, 7.0E7, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_2_t data =
{

{

{ { 0x41093E24 }, { 0x40467000 } }, //  0: vsLog2 ( 8.57767105      ) = ( 3.10058594      );
{ { 0x4093852F }, { 0x400D1AF0 } }, //  1: vsLog2 ( 4.61000776      ) = ( 2.20476913      );
{ { 0x41011D03 }, { 0x4040CCB4 } }, //  2: vsLog2 ( 8.06958294      ) = ( 3.01249409      );
{ { 0x41034C40 }, { 0x40425949 } }, //  3: vsLog2 ( 8.20611572      ) = ( 3.03669953      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x4008CDFFFECB65A2 } }, //  0: vdLog2 ( 8.57767068267691535       ) = ( 3.1005859285184707        );
{ { 0x401270A5F32DAE19 }, { 0x4001A35E121D3680 } }, //  1: vdLog2 ( 4.6100080486899282        ) = ( 2.20476926948214214       );
{ { 0x402023A0651C4741 }, { 0x40081996870CC042 } }, //  2: vdLog2 ( 8.06958309145159269       ) = ( 3.01249413974440561       );
{ { 0x40206988134D9FDD }, { 0x40084B2928FBF51B } }, //  3: vdLog2 ( 8.20611629793705255       ) = ( 3.0366996004000959        );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vLog2AccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));
  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f32[i].v1.f;
      vref1[i] = data.data_f32[i].v2.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmsLog2(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Log2", acc);
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Double precision test
//!

int vLog2AccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));
  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f64[i].v1.f;
      vref1[i] = data.data_f64[i].v2.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmdLog2(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Log2", acc);
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  return errs;
}

int main(int argc, char **argv) {
  int errs = 0;
  int total_errs = 0;

  printf("Running %s functions:\n", "Log2");

  printf("\tRunning %s with single precision real data type:\n", "Log2");
  errs = vLog2AccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Log2",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Log2");
  errs = vLog2AccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Log2",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Log2",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}