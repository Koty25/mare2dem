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
 *            Powx example program text (OpenMP offload interface)
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
data_3_t data =
{

{

{ { 0x41093E24 }, { 0x4048F5C3 }, { 0x44552AC1 } }, //  0: vsPowx ( 8.57767105     , 3.1400001       ) = ( 852.66803       );
{ { 0x4093852F }, { 0x4048F5C3 }, { 0x42F2B0D7 } }, //  1: vsPowx ( 4.61000776     , 3.1400001       ) = ( 121.34539       );
{ { 0x41011D03 }, { 0x4048F5C3 }, { 0x442FF9C3 } }, //  2: vsPowx ( 8.06958294     , 3.1400001       ) = ( 703.902527      );
{ { 0x41034C40 }, { 0x4048F5C3 }, { 0x44397EBA } }, //  3: vsPowx ( 8.20611572     , 3.1400001       ) = ( 741.980103      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x40091EB851EB851F }, { 0x408AA55778FCE709 } }, //  0: vdPowx ( 8.57767068267691535      , 3.14000000000000012       ) = ( 852.667711234856256       );
{ { 0x401270A5F32DAE19 }, { 0x40091EB851EB851F }, { 0x405E561AE4453C46 } }, //  1: vdPowx ( 4.6100080486899282       , 3.14000000000000012       ) = ( 121.345391337979066       );
{ { 0x402023A0651C4741 }, { 0x40091EB851EB851F }, { 0x4085FF381E850103 } }, //  2: vdPowx ( 8.06958309145159269      , 3.14000000000000012       ) = ( 703.902401961415649       );
{ { 0x40206988134D9FDD }, { 0x40091EB851EB851F }, { 0x40872FD74F3829E4 } }, //  3: vdPowx ( 8.20611629793705255      , 3.14000000000000012       ) = ( 741.980131567743683       );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x4048F5C3, 0x4048F5C3 }, { 0xC19A8841, 0xC21A0136 } }, //  0: vcPowx ( 4.61000776      + i * 8.57767105     , 3.1400001       + i * 3.1400001       ) = ( -19.3165302     + i * -38.5011826     );
{ { 0x41034C40, 0x41011D03 }, { 0x4048F5C3, 0x4048F5C3 }, { 0xC310B7AE, 0xC2ED2BD0 } }, //  1: vcPowx ( 8.20611572      + i * 8.06958294     , 3.1400001       + i * 3.1400001       ) = ( -144.717499     + i * -118.585571     );
{ { 0x4036ECDE, 0x41136B29 }, { 0x4048F5C3, 0x4048F5C3 }, { 0x401FC609, 0xC1B5CAED } }, //  2: vcPowx ( 2.85820723      + i * 9.21366215     , 3.1400001       + i * 3.1400001       ) = ( 2.49646211      + i * -22.7240849     );
{ { 0x40FDFDE5, 0x4082ABE3 }, { 0x4048F5C3, 0x4048F5C3 }, { 0xC2D4B709, 0x433D8524 } }, //  3: vcPowx ( 7.93724298      + i * 4.08348227     , 3.1400001       + i * 3.1400001       ) = ( -106.357491     + i * 189.520081      );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x40091EB851EB851F, 0x40091EB851EB851F }, { 0xC033510970136A6A, 0xC043402653EE7617 } }, //  0: vzPowx ( 4.6100080486899282        + i * 8.57767068267691535      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -19.3165502593423426      + i * -38.5011696733819733      );
{ { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0x40091EB851EB851F, 0x40091EB851EB851F }, { 0xC06216F614F1FE9F, 0xC05DA579773A3C38 } }, //  1: vzPowx ( 8.20611629793705255       + i * 8.06958309145159269      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -144.717539284368257      + i * -118.585538679952947      );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x40091EB851EB851F, 0x40091EB851EB851F }, { 0x4003F8B937C10DF9, 0xC036B95D4F1E546A } }, //  2: vzPowx ( 2.8582071867111174        + i * 9.21366148563738108      , 3.14000000000000012       + i * 3.14000000000000012       ) = ( 2.49644702489763093       + i * -22.7240800332114432      );
{ { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x40091EB851EB851F, 0x40091EB851EB851F }, { 0xC05A96E0D33FD2E3, 0x4067B0A4892EA43F } }, //  3: vzPowx ( 7.93724339382594657       + i * 4.0834825258625127       , 3.14000000000000012       + i * 3.14000000000000012       ) = ( -106.357472240760714      + i * 189.520084944817398       );
}

};

//!
//! @brief Single precision test
//!

int vPowxAccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));

  float varg2;

  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f32[i].v1.f;
      varg2 = data.data_f32[i].v2.f;
      vref1[i] = data.data_f32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    float marg2;
    *((float *)&(marg2)) = varg2;
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmsPowx(VLEN, (const float *)varg1, marg2, (float *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2, vres1[i], vres1[i],
                                 vref1[i], vref1[i], "Powx", acc);
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

int vPowxAccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));

  double varg2;

  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f64[i].v1.f;
      varg2 = data.data_f64[i].v2.f;
      vref1[i] = data.data_f64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    double marg2;
    *((double *)&(marg2)) = varg2;
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmdPowx(VLEN, (const double *)varg1, marg2, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2, vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Powx", acc);
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

int vPowxAccuracyLiteTest_float_complex() {
  int errs = 0;
  VM_COMPLEX8 * varg1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  VM_COMPLEX8 varg2;

  VM_COMPLEX8 * vres1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vref1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c32[i].v1.f;
      varg2 = data.data_c32[i].v2.f;
      vref1[i] = data.data_c32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    MKL_Complex8 marg2;
    *((VM_COMPLEX8 *)&(marg2)) = varg2;
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmcPowx(VLEN, (const MKL_Complex8 *)varg1, marg2, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2, vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Powx", acc);
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

int vPowxAccuracyLiteTest_double_complex() {
  int errs = 0;
  VM_COMPLEX16 * varg1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  VM_COMPLEX16 varg2;

  VM_COMPLEX16 * vres1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vref1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c64[i].v1.f;
      varg2 = data.data_c64[i].v2.f;
      vref1[i] = data.data_c64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
    MKL_Complex16 marg2;
    *((VM_COMPLEX16 *)&(marg2)) = varg2;
#pragma omp target data map(to:varg1[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, vres1) nowait
      {
        vmzPowx(VLEN, (const MKL_Complex16 *)varg1, marg2,
                (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double_complex(ARG2_RES1, varg1[i], varg2, vres1[i],
                                          vres1[i], vref1[i], vref1[i], "Powx",
                                          acc);
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

  printf("Running %s functions:\n", "Powx");

  printf("\tRunning %s with single precision real data type:\n", "Powx");
  errs = vPowxAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Powx",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Powx");
  errs = vPowxAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Powx",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Powx");
  errs = vPowxAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Powx",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Powx");
  errs = vPowxAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Powx",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Powx",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
