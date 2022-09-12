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
 *            Log10 example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 2.0, 5.0, 5.0E3, }, // VM_COMPLEX8
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

{ { 0x41093E24 }, { 0x3F6EF14C } }, //  0: vsLog10 ( 8.57767105      ) = ( 0.933369398     );
{ { 0x4093852F }, { 0x3F29E85A } }, //  1: vsLog10 ( 4.61000776      ) = ( 0.663701653     );
{ { 0x41011D03 }, { 0x3F682765 } }, //  2: vsLog10 ( 8.06958294      ) = ( 0.906851113     );
{ { 0x41034C40 }, { 0x3F6A04ED } }, //  3: vsLog10 ( 8.20611572      ) = ( 0.914137661     );
}

,

{

{ { 0x402127C473A3E923 }, { 0x3FEDDE297029A3ED } }, //  0: vdLog10 ( 8.57767068267691535       ) = ( 0.933369368617716355      );
{ { 0x401270A5F32DAE19 }, { 0x3FE53D0B503006B0 } }, //  1: vdLog10 ( 4.6100080486899282        ) = ( 0.663701683632288209      );
{ { 0x402023A0651C4741 }, { 0x3FED04EC97F0018B } }, //  2: vdLog10 ( 8.06958309145159269       ) = ( 0.906851097825027153      );
{ { 0x40206988134D9FDD }, { 0x3FED409DA344347E } }, //  3: vdLog10 ( 8.20611629793705255       ) = ( 0.914137667541254251      );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x3F7D0C5A, 0x3EEF9FB3 } }, //  0: vcLog10 ( 4.61000776      + i * 8.57767105      ) = ( 0.98846972      + i * 0.468015283     );
{ { 0x41034C40, 0x41011D03 }, { 0x3F87D028, 0x3EACC660 } }, //  1: vcLog10 ( 8.20611572      + i * 8.06958294      ) = ( 1.06103992      + i * 0.337450981     );
{ { 0x4036ECDE, 0x41136B29 }, { 0x3F7C0091, 0x3F0D3283 } }, //  2: vcLog10 ( 2.85820723      + i * 9.21366215      ) = ( 0.984383643     + i * 0.551551998     );
{ { 0x40FDFDE5, 0x4082ABE3 }, { 0x3F735E76, 0x3E534F92 } }, //  3: vcLog10 ( 7.93724298      + i * 4.08348227      ) = ( 0.95066011      + i * 0.206358224     );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x3FEFA18B2F6F5838, 0x3FDDF3F652C92C68 } }, //  0: vzLog10 ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 0.988469689031950871      + i * 0.468015271039524894      );
{ { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0x3FF0FA0504CCAB8B, 0x3FD598CBF8E47D48 } }, //  1: vzLog10 ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 1.06103994250108502       + i * 0.337450974520795643      );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x3FEF80121DF1D266, 0x3FE1A6505C169339 } }, //  2: vzLog10 ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 0.984383638845042652      + i * 0.551551990375265366      );
{ { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x3FEE6BCEC337BC1A, 0x3FCA69F241F4C86C } }, //  3: vzLog10 ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 0.950660115513417781      + i * 0.206358225064437462      );
}

};

//!
//! @brief Single precision test
//!

int vLog10AccuracyLiteTest_float() {
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
        vmsLog10(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Log10", acc);
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

int vLog10AccuracyLiteTest_double() {
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
        vmdLog10(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Log10", acc);
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

int vLog10AccuracyLiteTest_float_complex() {
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
        vmcLog10(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float_complex(ARG1_RES1, varg1[i], varg1[i],
                                         vres1[i], vres1[i], vref1[i], vref1[i],
                                         "Log10", acc);
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

int vLog10AccuracyLiteTest_double_complex() {
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
        vmzLog10(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Log10", acc);
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

  printf("Running %s functions:\n", "Log10");

  printf("\tRunning %s with single precision real data type:\n", "Log10");
  errs = vLog10AccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Log10",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Log10");
  errs = vLog10AccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Log10",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Log10");
  errs = vLog10AccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Log10",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Log10");
  errs = vLog10AccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Log10",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Log10",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
