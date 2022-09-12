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
 *            Pow2o3 example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x4065CC38 } }, //  0: vsPow2o3 ( 6.80375481      ) = ( 3.59058952      );
{ { 0xC007309A }, { 0x3FD2BA03 } }, //  1: vsPow2o3 ( -2.1123414      ) = ( 1.64630163      );
{ { 0x40B52EFA }, { 0x404B4F67 } }, //  2: vsPow2o3 ( 5.66198444      ) = ( 3.17672133      );
{ { 0x40BF006A }, { 0x40529713 } }, //  3: vsPow2o3 ( 5.96880054      ) = ( 3.29047084      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x400CB986F5486496 } }, //  0: vdPow2o3 ( 6.80375434309419092       ) = ( 3.59058944346595421       );
{ { 0xC000E6134801CC26 }, { 0x3FFA574068B5519B } }, //  1: vdPow2o3 ( -2.11234146361813924      ) = ( 1.64630165959997643       );
{ { 0x4016A5DF421D4BBE }, { 0x400969ECDC82CD04 } }, //  2: vdPow2o3 ( 5.66198447517211711       ) = ( 3.17672130846188061       );
{ { 0x4017E00D485FC01A }, { 0x400A52E2644589BD } }, //  3: vdPow2o3 ( 5.96880066952146571       ) = ( 3.29047087037289865       );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vPow2o3AccuracyLiteTest_float() {
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
        vmsPow2o3(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Pow2o3", acc);
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

int vPow2o3AccuracyLiteTest_double() {
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
        vmdPow2o3(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Pow2o3", acc);
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

  printf("Running %s functions:\n", "Pow2o3");

  printf("\tRunning %s with single precision real data type:\n", "Pow2o3");
  errs = vPow2o3AccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Pow2o3",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Pow2o3");
  errs = vPow2o3AccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Pow2o3",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Pow2o3",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
