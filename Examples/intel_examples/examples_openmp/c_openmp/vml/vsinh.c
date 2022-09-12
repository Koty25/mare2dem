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
 *            Sinh example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x43E14E52 } }, //  0: vsSinh ( 6.80375481      ) = ( 450.611877      );
{ { 0xC007309A }, { 0xC0825890 } }, //  1: vsSinh ( -2.1123414      ) = ( -4.07331085     );
{ { 0x40B52EFA }, { 0x430FDB98 } }, //  2: vsSinh ( 5.66198444      ) = ( 143.857788      );
{ { 0x40BF006A }, { 0x43438454 } }, //  3: vsSinh ( 5.96880054      ) = ( 195.516907      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x407C29C968C677F1 } }, //  0: vdSinh ( 6.80375434309419092       ) = ( 450.611672187106763       );
{ { 0xC000E6134801CC26 }, { 0xC0104B1218DE4197 } }, //  1: vdSinh ( -2.11234146361813924      ) = ( -4.07331122261566403      );
{ { 0x4016A5DF421D4BBE }, { 0x4061FB72FBB708AE } }, //  2: vdSinh ( 5.66198447517211711       ) = ( 143.857786042678924       );
{ { 0x4017E00D485FC01A }, { 0x4068708AA6866883 } }, //  3: vdSinh ( 5.96880066952146571       ) = ( 195.516925108448135       );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xC06228DD, 0x400582FC } }, //  0: vcSinh ( -2.1123414      + i * 6.80375481      ) = ( -3.5337441      + i * 2.08611965      );
{ { 0x40BF006A, 0x40B52EFA }, { 0x431EFD8F, 0xC2E396E2 } }, //  1: vcSinh ( 5.96880054      + i * 5.66198444      ) = ( 158.990463      + i * -113.794693     );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x429CBE3E, 0x4344CF32 } }, //  2: vcSinh ( -6.04897261     + i * 8.2329483       ) = ( 78.3715668      + i * 196.809357      );
{ { 0x40ABAABC, 0xC052EA36 }, { 0xC2D32BFA, 0x418315A9 } }, //  3: vcSinh ( 5.3645916       + i * -3.2955451      ) = ( -105.585892     + i * 16.3855762      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xC00C451C45E4AF59, 0x4000B05EB8F0615E } }, //  0: vzSinh ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -3.533745332757388        + i * 2.08611816867430289       );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x4063DFB206091F94, 0xC05C72DC5A12846B } }, //  1: vzSinh ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 158.990481393641517       + i * -113.794699209292659      );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x405397C41F9D3258, 0x406899E6F517D13C } }, //  2: vzSinh ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 78.3713454280017459       + i * 196.809443041341979       );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC05A657FC5D4A250, 0x403062B3F5B1D862 } }, //  3: vzSinh ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -105.585923631334708      + i * 16.3855584677879804       );
}

};

//!
//! @brief Single precision test
//!

int vSinhAccuracyLiteTest_float() {
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
        vmsSinh(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Sinh", acc);
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

int vSinhAccuracyLiteTest_double() {
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
        vmdSinh(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Sinh", acc);
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

int vSinhAccuracyLiteTest_float_complex() {
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
        vmcSinh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Sinh", acc);
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

int vSinhAccuracyLiteTest_double_complex() {
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
        vmzSinh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Sinh", acc);
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

  printf("Running %s functions:\n", "Sinh");

  printf("\tRunning %s with single precision real data type:\n", "Sinh");
  errs = vSinhAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Sinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Sinh");
  errs = vSinhAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Sinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Sinh");
  errs = vSinhAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Sinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Sinh");
  errs = vSinhAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Sinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Sinh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
