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
 *            Acospi example program text (OpenMP offload interface)
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

{ { 0x3F2E2D16 }, { 0x3E860CD7 } }, //  0: vsAcospi ( 0.680375457     ) = ( 0.26181671      );
{ { 0xBE584DC4 }, { 0x3F1157F3 } }, //  1: vsAcospi ( -0.211234152    ) = ( 0.567748249     );
{ { 0x3F10F262 }, { 0x3E9DE862 } }, //  2: vsAcospi ( 0.566198468     ) = ( 0.308413565     );
{ { 0x3F18CD22 }, { 0x3E97C2A2 } }, //  3: vsAcospi ( 0.596880078     ) = ( 0.296406806     );
}

,

{

{ { 0x3FE5C5A2B3EB8B46 }, { 0x3FD0C19AFA2583F6 } }, //  0: vdAcospi ( 0.680375434309419047      ) = ( 0.261816734584555788      );
{ { 0xBFCB09B873361370 }, { 0x3FE22AFE6445C016 } }, //  1: vdAcospi ( -0.211234146361813924     ) = ( 0.567748256535199003      );
{ { 0x3FE21E4C34E43C98 }, { 0x3FD3BD0C52F35C07 } }, //  2: vdAcospi ( 0.566198447517211711      ) = ( 0.308413582807986975      );
{ { 0x3FE319A439E63348 }, { 0x3FD2F85444FAAB98 } }, //  3: vdAcospi ( 0.596880066952146571      ) = ( 0.296406810152512801      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vAcospiAccuracyLiteTest_float() {
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
        vmsAcospi(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Acospi", acc);
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

int vAcospiAccuracyLiteTest_double() {
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
        vmdAcospi(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Acospi", acc);
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

  printf("Running %s functions:\n", "Acospi");

  printf("\tRunning %s with single precision real data type:\n", "Acospi");
  errs = vAcospiAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Acospi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Acospi");
  errs = vAcospiAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Acospi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Acospi",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
