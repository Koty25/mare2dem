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
 *            Erf example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
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

{ { 0x3F2E2D16 }, { 0x3F29FF1C } }, //  0: vsErf ( 0.680375457     ) = ( 0.66404891      );
{ { 0xBE584DC4 }, { 0xBE707D97 } }, //  1: vsErf ( -0.211234152    ) = ( -0.234854087    );
{ { 0x3F10F262 }, { 0x3F13A33B } }, //  2: vsErf ( 0.566198468     ) = ( 0.576709449     );
{ { 0x3F18CD22 }, { 0x3F19F50C } }, //  3: vsErf ( 0.596880078     ) = ( 0.601395369     );
}

,

{

{ { 0x3FE5C5A2B3EB8B46 }, { 0x3FE53FE388571266 } }, //  0: vdErf ( 0.680375434309419047      ) = ( 0.664048925675683632      );
{ { 0xBFCB09B873361370 }, { 0xBFCE0FB2D80C1EAA } }, //  1: vdErf ( -0.211234146361813924     ) = ( -0.234854083530298852     );
{ { 0x3FE21E4C34E43C98 }, { 0x3FE2746759551790 } }, //  2: vdErf ( 0.566198447517211711      ) = ( 0.576709436871839287      );
{ { 0x3FE319A439E63348 }, { 0x3FE33EA175373E87 } }, //  3: vdErf ( 0.596880066952146571      ) = ( 0.601395348488907966      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vErfAccuracyLiteTest_float() {
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
        vmsErf(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Erf", acc);
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

int vErfAccuracyLiteTest_double() {
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
        vmdErf(VLEN, (const double *)varg1, (double *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Erf", acc);
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

  printf("Running %s functions:\n", "Erf");

  printf("\tRunning %s with single precision real data type:\n", "Erf");
  errs = vErfAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Erf",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Erf");
  errs = vErfAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Erf",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Erf",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
