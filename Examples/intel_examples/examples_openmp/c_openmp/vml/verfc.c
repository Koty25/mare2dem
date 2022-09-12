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
 *            Erfc example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 256.0, 256.0, 256.0, }, // float
    { 1.0E10, 1.0E10, 1.0E10, }, // double
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

{ { 0x4042A1D0 }, { 0x378EC21C } }, //  0: vsErfc ( 3.04112625      ) = ( 1.70181083e-05  );
{ { 0x3EBB8B58 }, { 0x3F1ABCBB } }, //  1: vsErfc ( 0.366297483     ) = ( 0.604442298     );
{ { 0x402CB5CA }, { 0x390DFF0B } }, //  2: vsErfc ( 2.69859552      ) = ( 0.000135418188  );
{ { 0x403299DA }, { 0x38A643FA } }, //  3: vsErfc ( 2.79064035      ) = ( 7.92815845e-05  );
}

,

{

{ { 0x4008543A06F0A874 }, { 0x3EF1D8431328C9AA } }, //  0: vdErfc ( 3.04112630292825692       ) = ( 1.70181021507759823e-05   );
{ { 0x3FD7716B532EE2D8 }, { 0x3FE357973F1784C6 } }, //  1: vdErfc ( 0.366297560914558229      ) = ( 0.604442237116153747      );
{ { 0x400596B927AB2D72 }, { 0x3F21BFE29CC4BFA2 } }, //  2: vdErfc ( 2.69859534255163513       ) = ( 0.000135418331760115091   );
{ { 0x4006533B2B6CA676 }, { 0x3F14C8807F44637B } }, //  3: vdErfc ( 2.79064020085643971       ) = ( 7.92816570690645111e-05   );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vErfcAccuracyLiteTest_float() {
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
        vmsErfc(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Erfc", acc);
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

int vErfcAccuracyLiteTest_double() {
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
        vmdErfc(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Erfc", acc);
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

  printf("Running %s functions:\n", "Erfc");

  printf("\tRunning %s with single precision real data type:\n", "Erfc");
  errs = vErfcAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Erfc",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Erfc");
  errs = vErfcAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Erfc",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Erfc",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
