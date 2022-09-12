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
 *            Ln example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 2.0, 5.0, 5.0E3, }, // VM_COMPLEX8
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

{ { 0x41093E24 }, { 0x40098BE1 } }, //  0: vsLn ( 8.57767105      ) = ( 2.14916253      );
{ { 0x4093852F }, { 0x3FC39D07 } }, //  1: vsLn ( 4.61000776      ) = ( 1.52822959      );
{ { 0x41011D03 }, { 0x4005A376 } }, //  2: vsLn ( 8.06958294      ) = ( 2.08810186      );
{ { 0x41034C40 }, { 0x4006B659 } }, //  3: vsLn ( 8.20611572      ) = ( 2.10487962      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x4001317C0DAF2E04 } }, //  0: vdLn ( 8.57767068267691535       ) = ( 2.14916239443641821       );
{ { 0x401270A5F32DAE19 }, { 0x3FF873A0E2559781 } }, //  1: vdLn ( 4.6100080486899282        ) = ( 1.52822960292675725       );
{ { 0x402023A0651C4741 }, { 0x4000B46EBA08EB66 } }, //  2: vdLn ( 8.06958309145159269       ) = ( 2.08810181941719275       );
{ { 0x40206988134D9FDD }, { 0x4000D6CB33EF951D } }, //  3: vdLn ( 8.20611629793705255       ) = ( 2.10487976622483908       );
}

,

{

{ { 0x4093852F, 0x41093E24 }, { 0x4011AA91, 0x3F89F046 } }, //  0: vcLn ( 4.61000776      + i * 8.57767105      ) = ( 2.27603555      + i * 1.07764506      );
{ { 0x41034C40, 0x41011D03 }, { 0x401C5C52, 0x3F46EA1A } }, //  1: vcLn ( 8.20611572      + i * 8.06958294      ) = ( 2.44313478      + i * 0.777009606     );
{ { 0x4036ECDE, 0x41136B29 }, { 0x4011106B, 0x3FA28F36 } }, //  2: vcLn ( 2.85820723      + i * 9.21366215      ) = ( 2.26662707      + i * 1.26999545      );
{ { 0x40FDFDE5, 0x4082ABE3 }, { 0x400C182E, 0x3EF347D4 } }, //  3: vcLn ( 7.93724298      + i * 4.08348227      ) = ( 2.18897581      + i * 0.47515738      );
}

,

{

{ { 0x401270A5F32DAE19, 0x402127C473A3E923 }, { 0x40023552232A5F81, 0x3FF13E08AB53D691 } }, //  0: vzLn ( 4.6100080486899282        + i * 8.57767068267691535       ) = ( 2.27603557084142993       + i * 1.07764498638917794       );
{ { 0x40206988134D9FDD, 0x402023A0651C4741 }, { 0x40038B8A3BF86018, 0x3FE8DD4333C089A0 } }, //  1: vzLn ( 8.20611629793705255       + i * 8.06958309145159269       ) = ( 2.44313475467425789       + i * 0.77700958354789762       );
{ { 0x4006DD9BBAC0EE6B, 0x40226D6509CA7464 }, { 0x4002220D62974699, 0x3FF451E6AFEA09CE } }, //  2: vzLn ( 2.8582071867111174        + i * 9.21366148563738108       ) = ( 2.26662709259182948       + i * 1.26999539104928116       );
{ { 0x401FBFBCBB737F7A, 0x4010557C717977C6 }, { 0x40018305BFEE26E9, 0x3FDE68FA78360A78 } }, //  3: vzLn ( 7.93724339382594657       + i * 4.0834825258625127        ) = ( 2.1889758104851933        + i * 0.47515737285008397       );
}

};

//!
//! @brief Single precision test
//!

int vLnAccuracyLiteTest_float() {
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
      { vmsLn(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Ln", acc);
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

int vLnAccuracyLiteTest_double() {
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
        vmdLn(VLEN, (const double *)varg1, (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Ln", acc);
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

int vLnAccuracyLiteTest_float_complex() {
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
        vmcLn(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
              accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Ln", acc);
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

int vLnAccuracyLiteTest_double_complex() {
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
        vmzLn(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
              accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Ln", acc);
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

  printf("Running %s functions:\n", "Ln");

  printf("\tRunning %s with single precision real data type:\n", "Ln");
  errs = vLnAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Ln",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Ln");
  errs = vLnAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Ln",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Ln");
  errs = vLnAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Ln",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Ln");
  errs = vLnAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Ln",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Ln",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
