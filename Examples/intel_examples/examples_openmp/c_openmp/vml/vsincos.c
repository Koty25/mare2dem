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
 *            SinCos example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3EFEA7D8 }, { 0x3F5E16D8 } }, //  0: vsSinCos ( 6.80375481      ) = ( 0.497374296    , 0.867536068     );
{ { 0xC007309A }, { 0xBF5B5EAB }, { 0xBF03F53A } }, //  1: vsSinCos ( -2.1123414      ) = ( -0.856913269   , -0.51546061     );
{ { 0x40B52EFA }, { 0xBF14FEBF }, { 0x3F502C93 } }, //  2: vsSinCos ( 5.66198444      ) = ( -0.582012117   , 0.813180149     );
{ { 0x40BF006A }, { 0xBE9E5396 }, { 0x3F7373DF } }, //  3: vsSinCos ( 5.96880054      ) = ( -0.30923146    , 0.950986803     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FDFD4F93E99B2E0 }, { 0x3FEBC2DB7AB89950 } }, //  0: vdSinCos ( 6.80375434309419092       ) = ( 0.497373877652348639     , 0.867536296548488295      );
{ { 0xC000E6134801CC26 }, { 0xBFEB6BD5549D70BC }, { 0xBFE07EA757C4010B } }, //  1: vdSinCos ( -2.11234146361813924      ) = ( -0.85691324735991925     , -0.51546065465666524      );
{ { 0x4016A5DF421D4BBE }, { 0xBFE29FD7C840E7D0 }, { 0x3FEA05925DBF776B } }, //  2: vdSinCos ( 5.66198447517211711       ) = ( -0.582012072677793313    , 0.813180144406698502      );
{ { 0x4017E00D485FC01A }, { 0xBFD3CA723281D19B }, { 0x3FEE6E7BF8882000 } }, //  3: vdSinCos ( 5.96880066952146571       ) = ( -0.309231328318924248    , 0.950986848271895724      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vSinCosAccuracyLiteTest_float() {
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
        vmsSinCos(VLEN, (const float *)varg1, (float *)vres1, (float *)vres2,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG1_RES2, varg1[i], varg1[i], vres1[i],
                                 vres2[i], vref1[i], vref2[i], "SinCos", acc);
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

int vSinCosAccuracyLiteTest_double() {
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
        vmdSinCos(VLEN, (const double *)varg1, (double *)vres1, (double *)vres2,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG1_RES2, varg1[i], varg1[i], vres1[i],
                                  vres2[i], vref1[i], vref2[i], "SinCos", acc);
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

  printf("Running %s functions:\n", "SinCos");

  printf("\tRunning %s with single precision real data type:\n", "SinCos");
  errs = vSinCosAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "SinCos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "SinCos");
  errs = vSinCosAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "SinCos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "SinCos",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
