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
 *            CdfNorm example program text (OpenMP offload interface)
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

{ { 0x406150E8 }, { 0x3F7FF1E4 } }, //  0: vsCdfNorm ( 3.52056313      ) = ( 0.999784708     );
{ { 0x400BB8B6 }, { 0x3F7C48EA } }, //  1: vsCdfNorm ( 2.18314886      ) = ( 0.98548758      );
{ { 0x40565AE5 }, { 0x3F7FE574 } }, //  2: vsCdfNorm ( 3.34929776      ) = ( 0.999594927     );
{ { 0x40594CED }, { 0x3F7FE98A } }, //  3: vsCdfNorm ( 3.39532018      ) = ( 0.999657273     );
}

,

{

{ { 0x400C2A1D0378543A }, { 0x3FEFFE3C73412180 } }, //  0: vdCdfNorm ( 3.52056315146412846       ) = ( 0.999784684282573721      );
{ { 0x40017716B532EE2E }, { 0x3FEF891D3BBE3379 } }, //  1: vdCdfNorm ( 2.18314878045727934       ) = ( 0.985487572370046583      );
{ { 0x400ACB5C93D596B9 }, { 0x3FEFFCAE7A716C4F } }, //  2: vdCdfNorm ( 3.34929767127581757       ) = ( 0.999594916483497076      );
{ { 0x400B299D95B6533B }, { 0x3FEFFD3137C6318B } }, //  3: vdCdfNorm ( 3.39532010042821986       ) = ( 0.999657257970782864      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vCdfNormAccuracyLiteTest_float() {
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
        vmsCdfNorm(VLEN, (const float *)varg1, (float *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "CdfNorm", acc);
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

int vCdfNormAccuracyLiteTest_double() {
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
        vmdCdfNorm(VLEN, (const double *)varg1, (double *)vres1,
                   accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "CdfNorm", acc);
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

  printf("Running %s functions:\n", "CdfNorm");

  printf("\tRunning %s with single precision real data type:\n", "CdfNorm");
  errs = vCdfNormAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "CdfNorm",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "CdfNorm");
  errs = vCdfNormAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "CdfNorm",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "CdfNorm",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
