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
 *            Acosh example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { FLT_MAX, FLT_MAX, FLT_MAX, }, // VM_COMPLEX8
    { DBL_MAX, DBL_MAX, DBL_MAX, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_2_t data =
{

{

{ { 0x41093E24 }, { 0x4035B072 } }, //  0: vsAcosh ( 8.57767105      ) = ( 2.83889437      );
{ { 0x4093852F }, { 0x400D66CF } }, //  1: vsAcosh ( 4.61000776      ) = ( 2.20939994      );
{ { 0x41011D03 }, { 0x4031C0B8 } }, //  2: vsAcosh ( 8.06958294      ) = ( 2.77738762      );
{ { 0x41034C40 }, { 0x4032D5B5 } }, //  3: vsAcosh ( 8.20611572      ) = ( 2.79429364      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x4006B60E36BAC01A } }, //  0: vdAcosh ( 8.57767068267691535       ) = ( 2.83889429814736349       );
{ { 0x401270A5F32DAE19 }, { 0x4001ACD9F73AF77A } }, //  1: vdAcosh ( 4.6100080486899282        ) = ( 2.20940011166288475       );
{ { 0x402023A0651C4741 }, { 0x40063816F3916B18 } }, //  2: vdAcosh ( 8.06958309145159269       ) = ( 2.77738752639323749       );
{ { 0x40206988134D9FDD }, { 0x40065AB69C81EC96 } }, //  3: vdAcosh ( 8.20611629793705255       ) = ( 2.79429361602303583       );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x403E1EFD, 0x3F8A3801 } }, //  0: vcAcosh ( 4.61000776      + i * 8.57767105      ) = ( 2.97064137      + i * 1.0798341       );
{ { 0x41034C40, 0x41011D03 }, { 0x4048B869, 0x3F4765C9 } }, //  1: vcAcosh ( 8.20611572      + i * 8.06958294      ) = ( 3.1362555       + i * 0.778896868     );
{ { 0x4036ECDE, 0x41136B29 }, { 0x403D912A, 0x3FA2C0B4 } }, //  2: vcAcosh ( 2.85820723      + i * 9.21366215      ) = ( 2.96198511      + i * 1.27150583      );
{ { 0x40FDFDE5, 0x4082ABE3 }, { 0x403856E4, 0x3EF49846 } }, //  3: vcAcosh ( 7.93724298      + i * 4.08348227      ) = ( 2.88030338      + i * 0.477724254     );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x4007C3DFA89040E4, 0x3FF147001AC5D719 } }, //  0: vzAcosh ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 2.97064143839098627       + i * 1.07983408411150195       );
{ { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0x4009170D2ED1FCE7, 0x3FE8ECB91A991BB9 } }, //  1: vzAcosh ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 3.13625561312038625       + i * 0.778896858167050898      );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x4007B22542341397, 0x3FF458166D8F9BF1 } }, //  2: vzAcosh ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 2.96198512765335975       + i * 1.27150576398139159       );
{ { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x40070ADC8F07EE90, 0x3FDE9308BD0B88AC } }, //  3: vzAcosh ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 2.88030349486309234       + i * 0.477724251379309406      );
}

};

//!
//! @brief Single precision test
//!

int vAcoshAccuracyLiteTest_float() {
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
        vmsAcosh(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Acosh", acc);
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

int vAcoshAccuracyLiteTest_double() {
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
        vmdAcosh(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Acosh", acc);
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

int vAcoshAccuracyLiteTest_float_complex() {
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
        vmcAcosh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float_complex(ARG1_RES1, varg1[i], varg1[i],
                                         vres1[i], vres1[i], vref1[i], vref1[i],
                                         "Acosh", acc);
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

int vAcoshAccuracyLiteTest_double_complex() {
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
        vmzAcosh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Acosh", acc);
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

  printf("Running %s functions:\n", "Acosh");

  printf("\tRunning %s with single precision real data type:\n", "Acosh");
  errs = vAcoshAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Acosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Acosh");
  errs = vAcoshAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Acosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Acosh");
  errs = vAcoshAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Acosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Acosh");
  errs = vAcoshAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Acosh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Acosh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
