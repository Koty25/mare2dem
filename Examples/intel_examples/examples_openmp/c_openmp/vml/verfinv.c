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
 *            ErfInv example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { FLT_MAX, FLT_MAX, FLT_MAX, }, // float
    { DBL_MAX, DBL_MAX, DBL_MAX, }, // double
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

{ { 0x3F2C10F1 }, { 0x3F311CED } }, //  0: vsErfInv ( 0.672133505     ) = ( 0.691847622     );
{ { 0xBE5BB8DC }, { 0xBE4521D1 } }, //  1: vsErfInv ( -0.21457237     ) = ( -0.192511812    );
{ { 0x3F0EFF63 }, { 0x3F0B597E } }, //  2: vsErfInv ( 0.558584392     ) = ( 0.544334292     );
{ { 0x3F16CF15 }, { 0x3F14DA73 } }, //  3: vsErfInv ( 0.589097321     ) = ( 0.581458271     );
}

,

{

{ { 0x3FE5821DD706CAE3 }, { 0x3FE6239D393954E6 } }, //  0: vdErfInv ( 0.672133369420717108      ) = ( 0.69184743095925394       );
{ { 0xBFCB771B67504648 }, { 0xBFC8A43A11CD5897 } }, //  1: vdErfInv ( -0.214572358556823994     ) = ( -0.192511805241058126     );
{ { 0x3FE1DFEC4CAA263D }, { 0x3FE16B2FAA0C0A6A } }, //  2: vdErfInv ( 0.558584356055866871      ) = ( 0.544334251521218393      );
{ { 0x3FE2D9E26D4AF341 }, { 0x3FE29B4E1787A761 } }, //  3: vdErfInv ( 0.589097226583909728      ) = ( 0.581458135563689749      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vErfInvAccuracyLiteTest_float() {
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
        vmsErfInv(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "ErfInv", acc);
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

int vErfInvAccuracyLiteTest_double() {
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
        vmdErfInv(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "ErfInv", acc);
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

  printf("Running %s functions:\n", "ErfInv");

  printf("\tRunning %s with single precision real data type:\n", "ErfInv");
  errs = vErfInvAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "ErfInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "ErfInv");
  errs = vErfInvAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "ErfInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "ErfInv",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
