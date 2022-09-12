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
 *            ExpInt1 example program text (OpenMP offload interface)
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

{ { 0x4086F102 }, { 0x3B3EA496 } }, //  0: vsExpInt1 ( 4.2169199       ) = ( 0.00290898001   );
{ { 0x40021418 }, { 0x3D3F82AE } }, //  1: vsExpInt1 ( 2.03247643      ) = ( 0.0467554852    );
{ { 0x407BFADC }, { 0x3B85A594 } }, //  2: vsExpInt1 ( 3.93718624      ) = ( 0.00407857634   );
{ { 0x40806539 }, { 0x3B740114 } }, //  3: vsExpInt1 ( 4.01235628      ) = ( 0.00372320879   );
}

,

{

{ { 0x4010DE203A4CEF74 }, { 0x3F67D492E5CCFD17 } }, //  0: vdExpInt1 ( 4.21691981405807681       ) = ( 0.00290898028325127121    );
{ { 0x40004282F4C21EA0 }, { 0x3FA7F055E5D8C495 } }, //  1: vdExpInt1 ( 2.03247634141355604       ) = ( 0.0467554896425862801     );
{ { 0x400F7F5B79FEFEB7 }, { 0x3F70B4B29A274FD2 } }, //  2: vdExpInt1 ( 3.93718619641716883       ) = ( 0.00407857672185581339    );
{ { 0x40100CA71821B2E8 }, { 0x3F6E8022BE9F9EE1 } }, //  3: vdExpInt1 ( 4.01235616403275941       ) = ( 0.00372320924815371145    );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vExpInt1AccuracyLiteTest_float() {
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
        vmsExpInt1(VLEN, (const float *)varg1, (float *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "ExpInt1", acc);
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

int vExpInt1AccuracyLiteTest_double() {
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
        vmdExpInt1(VLEN, (const double *)varg1, (double *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "ExpInt1", acc);
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

  printf("Running %s functions:\n", "ExpInt1");

  printf("\tRunning %s with single precision real data type:\n", "ExpInt1");
  errs = vExpInt1AccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "ExpInt1",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "ExpInt1");
  errs = vExpInt1AccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "ExpInt1",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "ExpInt1",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
