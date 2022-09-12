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
 *            Sinpi example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3F140501 } }, //  0: vsSinpi ( 6.80375481      ) = ( 0.578201354     );
{ { 0xC007309A }, { 0xBEB0F8FB } }, //  1: vsSinpi ( -2.1123414      ) = ( -0.34564957     );
{ { 0x40B52EFA }, { 0xBF5F8FB0 } }, //  2: vsSinpi ( 5.66198444      ) = ( -0.873286247    );
{ { 0x40BF006A }, { 0xBDC86A5C } }, //  3: vsSinpi ( 5.96880054      ) = ( -0.0978591144   );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FE280A29F52D2F1 } }, //  0: vdSinpi ( 6.80375434309419092       ) = ( 0.578202544382831474      );
{ { 0xC000E6134801CC26 }, { 0xBFD61F2018124760 } }, //  1: vdSinpi ( -2.11234146361813924      ) = ( -0.345649741656520959     );
{ { 0x4016A5DF421D4BBE }, { 0xBFEBF1F5F0704FD0 } }, //  2: vdSinpi ( 5.66198447517211711       ) = ( -0.873286218268271242     );
{ { 0x4017E00D485FC01A }, { 0xBFB90D44ED33CC3C } }, //  3: vdSinpi ( 5.96880066952146571       ) = ( -0.0978587226016295886    );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vSinpiAccuracyLiteTest_float() {
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
        vmsSinpi(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Sinpi", acc);
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

int vSinpiAccuracyLiteTest_double() {
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
        vmdSinpi(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Sinpi", acc);
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

  printf("Running %s functions:\n", "Sinpi");

  printf("\tRunning %s with single precision real data type:\n", "Sinpi");
  errs = vSinpiAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Sinpi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Sinpi");
  errs = vSinpiAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Sinpi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Sinpi",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
