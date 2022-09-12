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
 *            Abs example program text (OpenMP offload interface)
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
data_2cf_t data =
{

{

{ { 0x40D9B85C }, { 0x40D9B85C } }, //  0: vsAbs ( 6.80375481      ) = ( 6.80375481      );
{ { 0xC007309A }, { 0x4007309A } }, //  1: vsAbs ( -2.1123414      ) = ( 2.1123414       );
{ { 0x40B52EFA }, { 0x40B52EFA } }, //  2: vsAbs ( 5.66198444      ) = ( 5.66198444      );
{ { 0x40BF006A }, { 0x40BF006A } }, //  3: vsAbs ( 5.96880054      ) = ( 5.96880054      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x401B370B60E66E18 } }, //  0: vdAbs ( 6.80375434309419092       ) = ( 6.80375434309419092       );
{ { 0xC000E6134801CC26 }, { 0x4000E6134801CC26 } }, //  1: vdAbs ( -2.11234146361813924      ) = ( 2.11234146361813924       );
{ { 0x4016A5DF421D4BBE }, { 0x4016A5DF421D4BBE } }, //  2: vdAbs ( 5.66198447517211711       ) = ( 5.66198447517211711       );
{ { 0x4017E00D485FC01A }, { 0x4017E00D485FC01A } }, //  3: vdAbs ( 5.96880066952146571       ) = ( 5.96880066952146571       );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40E3F8C8 } }, //  0: vcAbs ( -2.1123414      + i * 6.80375481      ) = ( 7.1241188       );
{ { 0x40BF006A, 0x40B52EFA }, { 0x4103A212 } }, //  1: vcAbs ( 5.96880054      + i * 5.66198444      ) = ( 8.22706795      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x412375B5 } }, //  2: vcAbs ( -6.04897261     + i * 8.2329483       ) = ( 10.2162371      );
{ { 0x40ABAABC, 0xC052EA36 }, { 0x40C978BB } }, //  3: vcAbs ( 5.3645916       + i * -3.2955451      ) = ( 6.29598761      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x401C7F18D50375CE } }, //  0: vzAbs ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( 7.12411816438311796       );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x402074424557113D } }, //  1: vzAbs ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 8.22706810653527931       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x40246EB68D8FA14C } }, //  2: vzAbs ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 10.2162365186528845       );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0x40192F176CA35836 } }, //  3: vzAbs ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 6.29598779437042388       );
}

};

//!
//! @brief Single precision test
//!

int vAbsAccuracyLiteTest_float() {
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
      { vmsAbs(VLEN, (const float *)varg1, vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1C_RES1R, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Abs", acc);
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

int vAbsAccuracyLiteTest_double() {
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
      { vmdAbs(VLEN, (const double *)varg1, vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1C_RES1R, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Abs", acc);
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex single precision test
//!

int vAbsAccuracyLiteTest_float_complex() {
  int errs = 0;
  VM_COMPLEX8 * varg1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_c32[i].v1.f;
      vref1[i] = data.data_c32[i].v2.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      { vmcAbs(VLEN, (const MKL_Complex8 *)varg1, vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      {
        VM_COMPLEX8 cres1 = vres1[i] + I * 0.0f, cref1 = vref1[i] + I * 0.0f;
        errs +=
            check_result_float_complex(ARG1C_RES1R, varg1[i], varg1[i], cres1,
                                       cres1, cref1, cref1, "Abs", acc);
      }
    }
  }

  free(varg1);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex double precision test
//!

int vAbsAccuracyLiteTest_double_complex() {
  int errs = 0;
  VM_COMPLEX16 * varg1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_c64[i].v1.f;
      vref1[i] = data.data_c64[i].v2.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      { vmzAbs(VLEN, (const MKL_Complex16 *)varg1, vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      {
        VM_COMPLEX16 cres1 = vres1[i] + I * 0.0, cref1 = vref1[i] + I * 0.0;
        errs +=
            check_result_double_complex(ARG1C_RES1R, varg1[i], varg1[i], cres1,
                                        cres1, cref1, cref1, "Abs", acc);
      }
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

  printf("Running %s functions:\n", "Abs");

  printf("\tRunning %s with single precision real data type:\n", "Abs");
  errs = vAbsAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Abs",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Abs");
  errs = vAbsAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Abs",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Abs");
  errs = vAbsAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Abs",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Abs");
  errs = vAbsAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Abs",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Abs",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
