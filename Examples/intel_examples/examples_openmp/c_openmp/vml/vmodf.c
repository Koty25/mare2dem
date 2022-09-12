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
 *            Modf example program text (OpenMP offload interface)
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
data_3_t data =
{

{

{ { 0x40D9B85C }, { 0x40C00000 }, { 0x3F4DC2E0 } }, //  0: vsModf ( 6.80375481      ) = ( 6              , 0.803754807     );
{ { 0xC007309A }, { 0xC0000000 }, { 0xBDE61340 } }, //  1: vsModf ( -2.1123414      ) = ( -2             , -0.112341404    );
{ { 0x40B52EFA }, { 0x40A00000 }, { 0x3F2977D0 } }, //  2: vsModf ( 5.66198444      ) = ( 5              , 0.661984444     );
{ { 0x40BF006A }, { 0x40A00000 }, { 0x3F780350 } }, //  3: vsModf ( 5.96880054      ) = ( 5              , 0.968800545     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x4018000000000000 }, { 0x3FE9B85B073370C0 } }, //  0: vdModf ( 6.80375434309419092       ) = ( 6                        , 0.803754343094190915      );
{ { 0xC000E6134801CC26 }, { 0xC000000000000000 }, { 0xBFBCC269003984C0 } }, //  1: vdModf ( -2.11234146361813924      ) = ( -2                       , -0.112341463618139237     );
{ { 0x4016A5DF421D4BBE }, { 0x4014000000000000 }, { 0x3FE52EFA10EA5DF0 } }, //  2: vdModf ( 5.66198447517211711       ) = ( 5                        , 0.661984475172117115      );
{ { 0x4017E00D485FC01A }, { 0x4014000000000000 }, { 0x3FEF006A42FE00D0 } }, //  3: vdModf ( 5.96880066952146571       ) = ( 5                        , 0.968800669521465707      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vModfAccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));
  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  float * vres2 = (float *)malloc(VLEN * sizeof(float));
  float * vref2 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f32[i].v1.f;
      vref1[i] = data.data_f32[i].v2.f;
      vref2[i] = data.data_f32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) map(tofrom:vres2[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1, vres2) nowait
      {
        vmsModf(VLEN, (const float *)varg1, (float *)vres1, (float *)vres2,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG1_RES2, varg1[i], varg1[i], vres1[i],
                                 vres2[i], vref1[i], vref2[i], "Modf", acc);
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  free(vres2);
  free(vref2);

  return errs;
}
//!
//! @brief Double precision test
//!

int vModfAccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));
  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  double * vres2 = (double *)malloc(VLEN * sizeof(double));
  double * vref2 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f64[i].v1.f;
      vref1[i] = data.data_f64[i].v2.f;
      vref2[i] = data.data_f64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) map(tofrom:vres2[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1, vres2) nowait
      {
        vmdModf(VLEN, (const double *)varg1, (double *)vres1, (double *)vres2,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG1_RES2, varg1[i], varg1[i], vres1[i],
                                  vres2[i], vref1[i], vref2[i], "Modf", acc);
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  free(vres2);
  free(vref2);

  return errs;
}

int main(int argc, char **argv) {
  int errs = 0;
  int total_errs = 0;

  printf("Running %s functions:\n", "Modf");

  printf("\tRunning %s with single precision real data type:\n", "Modf");
  errs = vModfAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Modf",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Modf");
  errs = vModfAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Modf",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Modf",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}