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
 *            Cosh example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 4.0, 5.0, 5.0E3, }, // VM_COMPLEX8
    { 4.0, 5.0, 7.0E7, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_2_t data =
{

{

{ { 0x40D9B85C }, { 0x43E14E76 } }, //  0: vsCosh ( 6.80375481      ) = ( 450.612976      );
{ { 0xC007309A }, { 0x4086376C } }, //  1: vsCosh ( -2.1123414      ) = ( 4.19426537      );
{ { 0x40B52EFA }, { 0x430FDC7B } }, //  2: vsCosh ( 5.66198444      ) = ( 143.861252      );
{ { 0x40BF006A }, { 0x434384FB } }, //  3: vsCosh ( 5.96880054      ) = ( 195.519455      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x407C29CDF446D9F3 } }, //  0: vdCosh ( 6.80375434309419092       ) = ( 450.612781788600103       );
{ { 0xC000E6134801CC26 }, { 0x4010C6ED92DFB4E6 } }, //  1: vdCosh ( -2.11234146361813924      ) = ( 4.19426564684292735       );
{ { 0x4016A5DF421D4BBE }, { 0x4061FB8F749A7034 } }, //  2: vdCosh ( 5.66198447517211711       ) = ( 143.86126165546159        );
{ { 0x4017E00D485FC01A }, { 0x4068709F9995F450 } }, //  3: vdCosh ( 5.96880066952146571       ) = ( 195.51948241508444        );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x4068E013, 0xC001A955 } }, //  0: vcCosh ( -2.1123414      + i * 6.80375481      ) = ( 3.6386764       + i * -2.02596021     );
{ { 0x40BF006A, 0x40B52EFA }, { 0x431EFE17, 0xC2E3961F } }, //  1: vcCosh ( 5.96880054      + i * 5.66198444      ) = ( 158.992538      + i * -113.793205     );
{ { 0xC0C1912F, 0x4103BA28 }, { 0xC29CBEB1, 0xC344CEA2 } }, //  2: vcCosh ( -6.04897261     + i * 8.2329483       ) = ( -78.3724442     + i * -196.807159     );
{ { 0x40ABAABC, 0xC052EA36 }, { 0xC2D32E58, 0x41831431 } }, //  3: vcCosh ( 5.3645916       + i * -3.2955451      ) = ( -105.590515     + i * 16.3848591      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x400D1C030BF087F8, 0xC0003529C8009691 } }, //  0: vzCosh ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( 3.63867768600266217       + i * -2.02595859767718212      );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x4063DFC30F2B8DEF, 0xC05C72C3F7571452 } }, //  1: vzCosh ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 158.992560944621317       + i * -113.793210825956777      );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0xC05397D26E1FA502, 0xC06899D4FE6B7053 } }, //  2: vzCosh ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -78.3722186383274959      + i * -196.807250223008481      );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC05A65CB88661D7E, 0x40306284EF15E8D7 } }, //  3: vzCosh ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -105.590547656747702      + i * 16.3848409107675614       );
}

};

//!
//! @brief Single precision test
//!

int vCoshAccuracyLiteTest_float() {
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
        vmsCosh(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Cosh", acc);
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

int vCoshAccuracyLiteTest_double() {
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
        vmdCosh(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Cosh", acc);
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

int vCoshAccuracyLiteTest_float_complex() {
  int errs = 0;
  VM_COMPLEX8 * varg1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vres1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vref1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

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
      {
        vmcCosh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Cosh", acc);
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

int vCoshAccuracyLiteTest_double_complex() {
  int errs = 0;
  VM_COMPLEX16 * varg1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vres1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vref1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

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
      {
        vmzCosh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Cosh", acc);
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

  printf("Running %s functions:\n", "Cosh");

  printf("\tRunning %s with single precision real data type:\n", "Cosh");
  errs = vCoshAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Cosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Cosh");
  errs = vCoshAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Cosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Cosh");
  errs = vCoshAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Cosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Cosh");
  errs = vCoshAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Cosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Cosh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
