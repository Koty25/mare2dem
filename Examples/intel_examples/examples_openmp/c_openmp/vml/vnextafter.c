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
 *            NextAfter example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0xC007309A }, { 0x40D9B85B } }, //  0: vsNextAfter ( 6.80375481     , -2.1123414      ) = ( 6.80375433      );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0x40B52EFB } }, //  1: vsNextAfter ( 5.66198444     , 5.96880054      ) = ( 5.66198492      );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0x4103BA27 } }, //  2: vsNextAfter ( 8.2329483      , -6.04897261     ) = ( 8.23294735      );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0xC052EA35 } }, //  3: vsNextAfter ( -3.2955451     , 5.3645916       ) = ( -3.29554486     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0x401B370B60E66E17 } }, //  0: vdNextAfter ( 6.80375434309419092      , -2.11234146361813924      ) = ( 6.80375434309419003       );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0x4016A5DF421D4BBF } }, //  1: vdNextAfter ( 5.66198447517211711      , 5.96880066952146571       ) = ( 5.661984475172118         );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0x40207744D998EE89 } }, //  2: vdNextAfter ( 8.23294715873568705      , -6.04897261413232101      ) = ( 8.23294715873568528       );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0xC00A5D46A314BA8D } }, //  3: vdNextAfter ( -3.2955448857022196      , 5.36459189623808186       ) = ( -3.29554488570221915      );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vNextAfterAccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));

  float * varg2 = (float *)malloc(VLEN * sizeof(float));
  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f32[i].v1.f;
      varg2[i] = data.data_f32[i].v2.f;
      vref1[i] = data.data_f32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmsNextAfter(VLEN, (const float *)varg1, (const float *)varg2,
                     (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i], vres1[i],
                             vref1[i], vref1[i], "NextAfter", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Double precision test
//!

int vNextAfterAccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));

  double * varg2 = (double *)malloc(VLEN * sizeof(double));
  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f64[i].v1.f;
      varg2[i] = data.data_f64[i].v2.f;
      vref1[i] = data.data_f64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmdNextAfter(VLEN, (const double *)varg1, (const double *)varg2,
                     (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i], vres1[i],
                              vref1[i], vref1[i], "NextAfter", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}

int main(int argc, char **argv) {
  int errs = 0;
  int total_errs = 0;

  printf("Running %s functions:\n", "NextAfter");

  printf("\tRunning %s with single precision real data type:\n", "NextAfter");
  errs = vNextAfterAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "NextAfter",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "NextAfter");
  errs = vNextAfterAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "NextAfter",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "NextAfter",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}