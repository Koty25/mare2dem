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
 *            LGamma example program text (OpenMP offload interface)
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

{ { 0x40866E17 }, { 0x40032FB5 } }, //  0: vsLGamma ( 4.2009387       ) = ( 2.04978681      );
{ { 0x3FFC67B3 }, { 0xBC3E5A32 } }, //  1: vsLGamma ( 1.97191465      ) = ( -0.0116181839   );
{ { 0x407A977D }, { 0x3FD7E3A1 } }, //  2: vsLGamma ( 3.91549611      ) = ( 1.68663418      );
{ { 0x407F8035 }, { 0x3FE4179D } }, //  3: vsLGamma ( 3.99220014      ) = ( 1.78197062      );
}

,

{

{ { 0x4010CDC2D8399B86 }, { 0x400065F67EE065B2 } }, //  0: vdLGamma ( 4.20093858577354773       ) = ( 2.04978655931764653       );
{ { 0x3FFF8CF65BFF19ED }, { 0xBF87CB4705A18EF5 } }, //  1: vdLGamma ( 1.97191463409546519       ) = ( -0.0116181896775695379    );
{ { 0x400F52EFA10EA5DF }, { 0x3FFAFC741D9972E9 } }, //  2: vdLGamma ( 3.91549611879302928       ) = ( 1.6866341739870967        );
{ { 0x400FF006A42FE00D }, { 0x3FFC82F39AFFB320 } }, //  3: vdLGamma ( 3.99220016738036643       ) = ( 1.78197060152451314       );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vLGammaAccuracyLiteTest_float() {
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
        vmsLGamma(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "LGamma", acc);
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

int vLGammaAccuracyLiteTest_double() {
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
        vmdLGamma(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "LGamma", acc);
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

  printf("Running %s functions:\n", "LGamma");

  printf("\tRunning %s with single precision real data type:\n", "LGamma");
  errs = vLGammaAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "LGamma",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "LGamma");
  errs = vLGammaAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "LGamma",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "LGamma",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
