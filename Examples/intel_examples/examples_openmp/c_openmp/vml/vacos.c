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
 *            Acos example program text (OpenMP offload interface)
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

{ { 0x3F2E2D16 }, { 0x3F5290C5 } }, //  0: vsAcos ( 0.680375457     ) = ( 0.822521508     );
{ { 0xBE584DC4 }, { 0x3FE44E1C } }, //  1: vsAcos ( -0.211234152    ) = ( 1.78363371      );
{ { 0x3F10F262 }, { 0x3F780A79 } }, //  2: vsAcos ( 0.566198468     ) = ( 0.9689098       );
{ { 0x3F18CD22 }, { 0x3F6E626E } }, //  3: vsAcos ( 0.596880078     ) = ( 0.931189418     );
}

,

{

{ { 0x3FE5C5A2B3EB8B46 }, { 0x3FEA5218ABED936B } }, //  0: vdAcos ( 0.680375434309419047      ) = ( 0.822521529957709219      );
{ { 0xBFCB09B873361370 }, { 0x3FFC89C38B81B281 } }, //  1: vdAcos ( -0.211234146361813924     ) = ( 1.78363375181939454       );
{ { 0x3FE21E4C34E43C98 }, { 0x3FEF014F38AD46CB } }, //  2: vdAcos ( 0.566198447517211711      ) = ( 0.968909846016879128      );
{ { 0x3FE319A439E63348 }, { 0x3FEDCC4DD5287A83 } }, //  3: vdAcos ( 0.596880066952146571      ) = ( 0.931189457249118724      );
}

,

{

{ { 0xBE584DC4, 0x3F2E2D16 }, { 0x3FDF6BAE, 0xBF252AF2 } }, //  0: vcAcos ( -0.211234152    + i * 0.680375457     ) = ( 1.74547362      + i * -0.645186543    );
{ { 0x3F18CD22, 0x3F10F262 }, { 0x3F8618E7, 0xBF1D4034 } }, //  1: vcAcos ( 0.596880078     + i * 0.566198468     ) = ( 1.04763496      + i * -0.614260912    );
{ { 0xBF1ADA8C, 0x3F52C372 }, { 0x400210E5, 0xBF52C236 } }, //  2: vcAcos ( -0.604897261    + i * 0.823294759     ) = ( 2.03228116      + i * -0.823275924    );
{ { 0x3F095564, 0xBEA8BB5E }, { 0x3F85D851, 0x3EBE9394 } }, //  3: vcAcos ( 0.536459208     + i * -0.329554498    ) = ( 1.04566395      + i * 0.372219682     );
}

,

{

{ { 0xBFCB09B873361370, 0x3FE5C5A2B3EB8B46 }, { 0x3FFBED75C72580BA, 0xBFE4A55E2B892514 } }, //  0: vzAcos ( -0.211234146361813924     + i * 0.680375434309419047      ) = ( 1.74547364989852705       + i * -0.645186505346972528     );
{ { 0x3FE319A439E63348, 0x3FE21E4C34E43C98 }, { 0x3FF0C31CDB3C5272, 0xBFE3A80667B42868 } }, //  1: vzAcos ( 0.596880066952146571      + i * 0.566198447517211711      ) = ( 1.04763494147223613       + i * -0.614260866686220375     );
{ { 0xBFE35B518066B6A3, 0x3FEA586E28F4B0DC }, { 0x4000421CA734509A, 0xBFEA5846A28C515F } }, //  2: vzAcos ( -0.604897261413232079     + i * 0.823294715873568617      ) = ( 2.03228121403124096       + i * -0.823275868870535166     );
{ { 0x3FE12AAC76625558, 0xBFD5176BB5AA2ED8 }, { 0x3FF0BB0A242AB476, 0x3FD7D272742FD332 } }, //  3: vzAcos ( 0.536459189623808008      + i * -0.32955448857022196      ) = ( 1.04566396835005326       + i * 0.372219670737922503      );
}

};

//!
//! @brief Single precision test
//!

int vAcosAccuracyLiteTest_float() {
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
        vmsAcos(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Acos", acc);
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

int vAcosAccuracyLiteTest_double() {
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
        vmdAcos(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Acos", acc);
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

int vAcosAccuracyLiteTest_float_complex() {
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
        vmcAcos(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Acos", acc);
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

int vAcosAccuracyLiteTest_double_complex() {
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
        vmzAcos(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Acos", acc);
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

  printf("Running %s functions:\n", "Acos");

  printf("\tRunning %s with single precision real data type:\n", "Acos");
  errs = vAcosAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Acos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Acos");
  errs = vAcosAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Acos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Acos");
  errs = vAcosAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Acos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Acos");
  errs = vAcosAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Acos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Acos",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
