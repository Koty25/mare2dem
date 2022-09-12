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
 *            Sqrt example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, FLT_MAX, }, // float
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

{ { 0x41093E24 }, { 0x403B70E7 } }, //  0: vsSqrt ( 8.57767105      ) = ( 2.92876601      );
{ { 0x4093852F }, { 0x400969F8 } }, //  1: vsSqrt ( 4.61000776      ) = ( 2.14709282      );
{ { 0x41011D03 }, { 0x4035CE0C } }, //  2: vsSqrt ( 8.06958294      ) = ( 2.8407011       );
{ { 0x41034C40 }, { 0x40375621 } }, //  3: vsSqrt ( 8.20611572      ) = ( 2.86463189      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x40076E1CE786FE27 } }, //  0: vdSqrt ( 8.57767068267691535       ) = ( 2.9287660682746437        );
{ { 0x401270A5F32DAE19 }, { 0x40012D3F0ED3A6F0 } }, //  1: vdSqrt ( 4.6100080486899282        ) = ( 2.14709292968188237       );
{ { 0x402023A0651C4741 }, { 0x4006B9C187E1F54A } }, //  2: vdSqrt ( 8.06958309145159269       ) = ( 2.84070116194076139       );
{ { 0x40206988134D9FDD }, { 0x4006EAC429F8399F } }, //  3: vdSqrt ( 8.20611629793705255       ) = ( 2.86463196553013644       );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x402B6B72, 0x3FCCF5B2 } }, //  0: vcSqrt ( 4.61000776      + i * 8.57767105      ) = ( 2.67843294      + i * 1.60124803      );
{ { 0x41034C40, 0x41011D03 }, { 0x4048F083, 0x3FA47E0B } }, //  1: vcSqrt ( 8.20611572      + i * 8.06958294      ) = ( 3.13967967      + i * 1.28509653      );
{ { 0x4036ECDE, 0x41136B29 }, { 0x40200838, 0x3FEBD28B } }, //  2: vcSqrt ( 2.85820723      + i * 9.21366215      ) = ( 2.50050163      + i * 1.84236276      );
{ { 0x40FDFDE5, 0x4082ABE3 }, { 0x4039D6BB, 0x3F34013F } }, //  3: vcSqrt ( 7.93724298      + i * 4.08348227      ) = ( 2.90373111      + i * 0.703144014     );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x40056D6E42511161, 0x3FF99EB631575346 } }, //  0: vzSqrt ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 2.67843295869731479       + i * 1.60124797128556073       );
{ { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0x40091E1072F5F2ED, 0x3FF48FC159BA032B } }, //  1: vzSqrt ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 3.13967981160236898       + i * 1.28509650277573928       );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x40040106EB21DA6F, 0x3FFD7A515924BC6F } }, //  2: vzSqrt ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 2.50050147721349658       + i * 1.84236273595504563       );
{ { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x40073AD76D49B864, 0x3FE68027E8B0FCB8 } }, //  3: vzSqrt ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 2.90373120671488216       + i * 0.703144030070595782      );
}

};

//!
//! @brief Single precision test
//!

int vSqrtAccuracyLiteTest_float() {
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
        vmsSqrt(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Sqrt", acc);
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

int vSqrtAccuracyLiteTest_double() {
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
        vmdSqrt(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Sqrt", acc);
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

int vSqrtAccuracyLiteTest_float_complex() {
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
        vmcSqrt(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Sqrt", acc);
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

int vSqrtAccuracyLiteTest_double_complex() {
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
        vmzSqrt(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Sqrt", acc);
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

  printf("Running %s functions:\n", "Sqrt");

  printf("\tRunning %s with single precision real data type:\n", "Sqrt");
  errs = vSqrtAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Sqrt",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Sqrt");
  errs = vSqrtAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Sqrt",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Sqrt");
  errs = vSqrtAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Sqrt",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Sqrt");
  errs = vSqrtAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Sqrt",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Sqrt",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
