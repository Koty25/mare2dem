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
 *            Cosd example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 5000.0, 5000.0, 5000.0, }, // float
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

{ { 0x442A1808 }, { 0x3F452E7B } }, //  0: vsCosd ( 680.375488      ) = ( 0.770240486     );
{ { 0xC3533BF0 }, { 0xBF5AE4EB } }, //  1: vsCosd ( -211.234131     ) = ( -0.855055511    );
{ { 0x440D8CB4 }, { 0xBF65B37D } }, //  2: vsCosd ( 566.198486      ) = ( -0.897270024    );
{ { 0x44153854 }, { 0xBF0BE061 } }, //  3: vsCosd ( 596.880127      ) = ( -0.5463925      );
}

,

{

{ { 0x40854300E3B40602 }, { 0x3FE8A5CE182B4917 } }, //  0: vdCosd ( 680.375434309419006       ) = ( 0.770239875035516941      );
{ { 0xC06A677E2082CEFC }, { 0xBFEB5C9D1A82A314 } }, //  1: vdCosd ( -211.234146361813941      ) = ( -0.85505538156312122      );
{ { 0x4081B1966BA6E32C }, { 0xBFECB67045A4F282 } }, //  2: vdCosd ( 566.198447517211662       ) = ( -0.897270332359383582     );
{ { 0x4082A70A608ACE14 }, { 0xBFE17C0DF1D13588 } }, //  3: vdCosd ( 596.880066952146535       ) = ( -0.546393368052734196     );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vCosdAccuracyLiteTest_float() {
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
        vmsCosd(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Cosd", acc);
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

int vCosdAccuracyLiteTest_double() {
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
        vmdCosd(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Cosd", acc);
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

  printf("Running %s functions:\n", "Cosd");

  printf("\tRunning %s with single precision real data type:\n", "Cosd");
  errs = vCosdAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Cosd",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Cosd");
  errs = vCosdAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Cosd",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Cosd",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
