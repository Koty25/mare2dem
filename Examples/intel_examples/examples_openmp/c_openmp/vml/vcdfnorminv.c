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
 *            CdfNormInv example program text (OpenMP offload interface)
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

{ { 0x3C82EB10 }, { 0xC00945B6 } }, //  0: vsCdfNormInv ( 0.0159812272    ) = ( -2.14487982     );
{ { 0x3D780F8E }, { 0xBFC668DF } }, //  1: vsCdfNormInv ( 0.0605617091    ) = ( -1.55007541     );
{ { 0x3CB1AF64 }, { 0xC0014831 } }, //  2: vsCdfNormInv ( 0.0216900781    ) = ( -2.02003121     );
{ { 0x3CA51E30 }, { 0xC0033C02 } }, //  3: vsCdfNormInv ( 0.0201559961    ) = ( -2.05053759     );
}

,

{

{ { 0x3F905D621353EDF8 }, { 0xC00128B6C129E3D3 } }, //  0: vdCdfNormInv ( 0.0159812282845290532     ) = ( -2.14487982663238475      );
{ { 0x3FAF01F1B0A46A4A }, { 0xBFF8CD1BE28091B7 } }, //  1: vdCdfNormInv ( 0.060561707318090699      ) = ( -1.5500754211180785       );
{ { 0x3F9635EC782C6BD8 }, { 0xC0002906191BD5FA } }, //  2: vdCdfNormInv ( 0.0216900776241394089     ) = ( -2.02003116241644154      );
{ { 0x3F94A3C609C2E128 }, { 0xC0006780490CB477 } }, //  3: vdCdfNormInv ( 0.020155996652392677      ) = ( -2.05053765363714602      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vCdfNormInvAccuracyLiteTest_float() {
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
        vmsCdfNormInv(VLEN, (const float *)varg1, (float *)vres1,
                      accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i], vres1[i],
                             vref1[i], vref1[i], "CdfNormInv", acc);
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

int vCdfNormInvAccuracyLiteTest_double() {
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
        vmdCdfNormInv(VLEN, (const double *)varg1, (double *)vres1,
                      accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i], vres1[i],
                              vref1[i], vref1[i], "CdfNormInv", acc);
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

  printf("Running %s functions:\n", "CdfNormInv");

  printf("\tRunning %s with single precision real data type:\n", "CdfNormInv");
  errs = vCdfNormInvAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "CdfNormInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "CdfNormInv");
  errs = vCdfNormInvAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "CdfNormInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "CdfNormInv",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
