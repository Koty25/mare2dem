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
 *            TGamma example program text (OpenMP offload interface)
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

{ { 0x4083DF80 }, { 0x40DFFC58 } }, //  0: vsTGamma ( 4.12103271      ) = ( 6.99955368      );
{ { 0x3FD5A544 }, { 0x3F6734B2 } }, //  1: vsTGamma ( 1.66910601      ) = ( 0.903147817     );
{ { 0x4073A6A4 }, { 0x40977CDA } }, //  2: vsTGamma ( 3.80704594      ) = ( 4.73399067      );
{ { 0x40790D08 }, { 0x40A7CDD7 } }, //  3: vsTGamma ( 3.89142036      ) = ( 5.24387693      );
}

,

{

{ { 0x40107BEFEDD8F7E0 }, { 0x401BFF8A5E5FB880 } }, //  0: vdTGamma ( 4.12103244435090232       ) = ( 6.9995512720034867        );
{ { 0x3FFAB4A898656952 }, { 0x3FECE696513F6393 } }, //  1: vdTGamma ( 1.66910609750501182       ) = ( 0.903147848784202956      );
{ { 0x400E74D4645CE9A8 }, { 0x4012EF9AEA1EC930 } }, //  2: vdTGamma ( 3.80704573067233198       ) = ( 4.73398938954260018       );
{ { 0x400F21A0E7CE4342 }, { 0x4014F9BA8CE68C6E } }, //  3: vdTGamma ( 3.89142018411840329       ) = ( 5.24387569577366541       );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vTGammaAccuracyLiteTest_float() {
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
        vmsTGamma(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "TGamma", acc);
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

int vTGammaAccuracyLiteTest_double() {
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
        vmdTGamma(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "TGamma", acc);
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

  printf("Running %s functions:\n", "TGamma");

  printf("\tRunning %s with single precision real data type:\n", "TGamma");
  errs = vTGammaAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "TGamma",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "TGamma");
  errs = vTGammaAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "TGamma",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "TGamma",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
