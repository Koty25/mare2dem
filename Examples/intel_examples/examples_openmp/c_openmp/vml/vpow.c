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
 *            Pow example program text (OpenMP offload interface)
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

{ { 0x41093E24 }, { 0x4093852F }, { 0x469CE711 } }, //  0: vsPow ( 8.57767105     , 4.61000776      ) = ( 20083.5332      );
{ { 0x41011D03 }, { 0x41034C40 }, { 0x4BD2F79E } }, //  1: vsPow ( 8.06958294     , 8.20611572      ) = ( 27651900        );
{ { 0x41136B29 }, { 0x4036ECDE }, { 0x440EB888 } }, //  2: vsPow ( 9.21366215     , 2.85820723      ) = ( 570.883301      );
{ { 0x4082ABE3 }, { 0x40FDFDE5 }, { 0x478A3D0D } }, //  3: vsPow ( 4.08348227     , 7.93724298      ) = ( 70778.1016      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x401270A5F32DAE19 }, { 0x40D39CE2AABD156D } }, //  0: vdPow ( 8.57767068267691535      , 4.6100080486899282        ) = ( 20083.5416710576283       );
{ { 0x402023A0651C4741 }, { 0x40206988134D9FDD }, { 0x417A5EF61F31C368 } }, //  1: vdPow ( 8.06958309145159269      , 8.20611629793705255       ) = ( 27651937.9496492445       );
{ { 0x40226D6509CA7464 }, { 0x4006DD9BBAC0EE6B }, { 0x4081D7109BAA7980 } }, //  2: vdPow ( 9.21366148563738108      , 2.8582071867111174        ) = ( 570.883109409172903       );
{ { 0x4010557C717977C6 }, { 0x401FBFBCBB737F7A }, { 0x40F147A2E685447B } }, //  3: vdPow ( 4.0834825258625127       , 7.93724339382594657       ) = ( 70778.1812794375437       );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x41034C40, 0x41011D03 }, { 0xC623D50C, 0x4693B11C } }, //  0: vcPow ( 4.61000776      + i * 8.57767105     , 8.20611572      + i * 8.06958294      ) = ( -10485.2617     + i * 18904.5547      );
{ { 0x4036ECDE, 0x41136B29 }, { 0x40FDFDE5, 0x4082ABE3 }, { 0x489D1870, 0x482626A1 } }, //  1: vcPow ( 2.85820723      + i * 9.21366215     , 7.93724298      + i * 4.08348227      ) = ( 321731.5        + i * 170138.516      );
{ { 0x40C0F87C, 0x40649ED8 }, { 0x40D64D6C, 0x40AB29A5 }, { 0x4566C04A, 0x46CBEE93 } }, //  2: vcPow ( 6.03033257      + i * 3.57219505     , 6.69695091      + i * 5.34883356      ) = ( 3692.01807      + i * 26103.2871      );
{ { 0x40B56AA4, 0x408B1733 }, { 0x411410F1, 0x41193290 }, { 0x480FE542, 0xC714EC50 } }, //  3: vcPow ( 5.66926765      + i * 4.34658194     , 9.25413609      + i * 9.57484436      ) = ( 147349.031      + i * -38124.3125     );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0xC0C47AA470F16D42, 0x40D27624C7C04344 } }, //  0: vzPow ( 4.6100080486899282        + i * 8.57767068267691535      , 8.20611629793705255       + i * 8.06958309145159269       ) = ( -10485.2846967490659      + i * 18904.5746918351069       );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x4113A30DBF54F1C6, 0x4104C4D62BEAA20C } }, //  1: vzPow ( 2.8582071867111174        + i * 9.21366148563738108      , 7.93724339382594657       + i * 4.0834825258625127        ) = ( 321731.436847474775       + i * 170138.771443620673       );
{ { 0x40181F0F82C50AEC, 0x400C93DAEEF5F483 }, { 0x401AC9AD9555935C, 0x40156534AD76CA6A }, { 0x40ACD804D2E3A03F, 0x40D97DD363DA339D } }, //  2: vzPow ( 6.03033260657933212       + i * 3.57219492614837142      , 6.69695123038112783       + i * 5.34883376157322665       ) = ( 3692.00942145663385       + i * 26103.3029695037876       );
{ { 0x4016AD549DF3C110, 0x401162E657685F66 }, { 0x4022821E20A96AA3, 0x40232652067FE63E }, { 0x4101FCA94042F62A, 0xC0E29D899EA112CB } }, //  3: vzPow ( 5.66926810074097887       + i * 4.34658180784740544      , 9.25413610523293606       + i * 9.57484455405494472       ) = ( 147349.156377719075       + i * -38124.3006139151621      );
}

};

//!
//! @brief Single precision test
//!

int vPowAccuracyLiteTest_float() {
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
        vmsPow(VLEN, (const float *)varg1, (const float *)varg2, (float *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Pow", acc);
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

int vPowAccuracyLiteTest_double() {
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
        vmdPow(VLEN, (const double *)varg1, (const double *)varg2,
               (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Pow", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex single precision test
//!

int vPowAccuracyLiteTest_float_complex() {
  int errs = 0;
  VM_COMPLEX8 * varg1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  VM_COMPLEX8 * varg2 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vres1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vref1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c32[i].v1.f;
      varg2[i] = data.data_c32[i].v2.f;
      vref1[i] = data.data_c32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmcPow(VLEN, (const MKL_Complex8 *)varg1, (const MKL_Complex8 *)varg2,
               (MKL_Complex8 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Pow", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex double precision test
//!

int vPowAccuracyLiteTest_double_complex() {
  int errs = 0;
  VM_COMPLEX16 * varg1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  VM_COMPLEX16 * varg2 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vres1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vref1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c64[i].v1.f;
      varg2[i] = data.data_c64[i].v2.f;
      vref1[i] = data.data_c64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmzPow(VLEN, (const MKL_Complex16 *)varg1, (const MKL_Complex16 *)varg2,
               (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Pow", acc);
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

  printf("Running %s functions:\n", "Pow");

  printf("\tRunning %s with single precision real data type:\n", "Pow");
  errs = vPowAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Pow",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Pow");
  errs = vPowAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Pow",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Pow");
  errs = vPowAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Pow",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Pow");
  errs = vPowAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Pow",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Pow",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
