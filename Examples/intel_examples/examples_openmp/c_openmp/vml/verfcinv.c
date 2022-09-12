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
 *            ErfcInv example program text (OpenMP offload interface)
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

{ { 0x3F57E803 }, { 0x3E0F0DBE } }, //  0: vsErfcInv ( 0.843383968     ) = ( 0.13970086      );
{ { 0x3ED02026 }, { 0x3F16428E } }, //  1: vsErfcInv ( 0.406495273     ) = ( 0.586953044     );
{ { 0x3F49957D }, { 0x3E433D7D } }, //  2: vsErfcInv ( 0.78743726      ) = ( 0.190664247     );
{ { 0x3F4D6EC1 }, { 0x3E3520C9 } }, //  3: vsErfcInv ( 0.80247122      ) = ( 0.176882878     );
}

,

{

{ { 0x3FEAFD005D47E586 }, { 0x3FC1E1B7D3A5B285 } }, //  0: vdErfcInv ( 0.843383962811615318      ) = ( 0.139700868934046124      );
{ { 0x3FDA0404BAD030FF }, { 0x3FE2C851C4CA7D32 } }, //  1: vdErfcInv ( 0.406495268282711153      ) = ( 0.586953052861565405      );
{ { 0x3FE932AF94CBFEF9 }, { 0x3FC867AFB9967916 } }, //  2: vdErfcInv ( 0.787437239283433787      ) = ( 0.190664258593593317      );
{ { 0x3FE9ADD8269C5173 }, { 0x3FC6A4190C656C3E } }, //  3: vdErfcInv ( 0.80247123280655186       ) = ( 0.176882868817161254      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vErfcInvAccuracyLiteTest_float() {
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
        vmsErfcInv(VLEN, (const float *)varg1, (float *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "ErfcInv", acc);
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

int vErfcInvAccuracyLiteTest_double() {
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
        vmdErfcInv(VLEN, (const double *)varg1, (double *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "ErfcInv", acc);
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

  printf("Running %s functions:\n", "ErfcInv");

  printf("\tRunning %s with single precision real data type:\n", "ErfcInv");
  errs = vErfcInvAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "ErfcInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "ErfcInv");
  errs = vErfcInvAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "ErfcInv",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "ErfcInv",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
