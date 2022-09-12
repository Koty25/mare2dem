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
 *            Asin example program text (OpenMP offload interface)
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

{ { 0x3F2E2D16 }, { 0x3F3F8EF0 } }, //  0: vsAsin ( 0.680375457     ) = ( 0.748274803     );
{ { 0xBE584DC4 }, { 0xBE59F20E } }, //  1: vsAsin ( -0.211234152    ) = ( -0.212837428    );
{ { 0x3F10F262 }, { 0x3F1A153C } }, //  2: vsAsin ( 0.566198468     ) = ( 0.601886511     );
{ { 0x3F18CD22 }, { 0x3F23BD47 } }, //  3: vsAsin ( 0.596880078     ) = ( 0.639606893     );
}

,

{

{ { 0x3FE5C5A2B3EB8B46 }, { 0x3FE7F1DDFC9AC6C5 } }, //  0: vdAsin ( 0.680375434309419047      ) = ( 0.748274796837187339      );
{ { 0xBFCB09B873361370 }, { 0xBFCB3E41B9EC2B45 } }, //  1: vdAsin ( -0.211234146361813924     ) = ( -0.2128374250244979       );
{ { 0x3FE21E4C34E43C98 }, { 0x3FE342A76FDB1365 } }, //  2: vdAsin ( 0.566198447517211711      ) = ( 0.60188648077801743       );
{ { 0x3FE319A439E63348 }, { 0x3FE477A8D35FDFAD } }, //  3: vdAsin ( 0.596880066952146571      ) = ( 0.639606869545777834      );
}

,

{

{ { 0xBE584DC4, 0x3F2E2D16 }, { 0xBE32DE9D, 0x3F252AF2 } }, //  0: vcAsin ( -0.211234152    + i * 0.680375457     ) = ( -0.174677327    + i * 0.645186543     );
{ { 0x3F18CD22, 0x3F10F262 }, { 0x3F05EDE8, 0x3F1D4034 } }, //  1: vcAsin ( 0.596880078     + i * 0.566198468     ) = ( 0.523161411     + i * 0.614260912     );
{ { 0xBF1ADA8C, 0x3F52C372 }, { 0xBEEC47BF, 0x3F52C236 } }, //  2: vcAsin ( -0.604897261    + i * 0.823294759     ) = ( -0.461484879    + i * 0.823275924     );
{ { 0x3F095564, 0xBEA8BB5E }, { 0x3F066F13, 0xBEBE9394 } }, //  3: vcAsin ( 0.536459208     + i * -0.329554498    ) = ( 0.525132358     + i * -0.372219682    );
}

,

{

{ { 0xBFCB09B873361370, 0x3FE5C5A2B3EB8B46 }, { 0xBFC65BD3970A9D10, 0x3FE4A55E2B892514 } }, //  0: vzAsin ( -0.211234146361813924     + i * 0.680375434309419047      ) = ( -0.174677323103630489     + i * 0.645186505346972528      );
{ { 0x3FE319A439E63348, 0x3FE21E4C34E43C98 }, { 0x3FE0BDBCF20FB54C, 0x3FE3A80667B42868 } }, //  1: vzAsin ( 0.596880066952146571      + i * 0.566198447517211711      ) = ( 0.523161385322660433      + i * 0.614260866686220375      );
{ { 0xBFE35B518066B6A3, 0x3FEA586E28F4B0DC }, { 0xBFDD88F7E891D06C, 0x3FEA5846A28C515F } }, //  2: vzAsin ( -0.604897261413232079     + i * 0.823294715873568617      ) = ( -0.461484887236344177     + i * 0.823275868870535166      );
{ { 0x3FE12AAC76625558, 0xBFD5176BB5AA2ED8 }, { 0x3FE0CDE26032F145, 0xBFD7D272742FD332 } }, //  3: vzAsin ( 0.536459189623808008      + i * -0.32955448857022196      ) = ( 0.525132358444843406      + i * -0.372219670737922503     );
}

};

//!
//! @brief Single precision test
//!

int vAsinAccuracyLiteTest_float() {
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
        vmsAsin(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Asin", acc);
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

int vAsinAccuracyLiteTest_double() {
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
        vmdAsin(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Asin", acc);
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

int vAsinAccuracyLiteTest_float_complex() {
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
        vmcAsin(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Asin", acc);
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

int vAsinAccuracyLiteTest_double_complex() {
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
        vmzAsin(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Asin", acc);
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

  printf("Running %s functions:\n", "Asin");

  printf("\tRunning %s with single precision real data type:\n", "Asin");
  errs = vAsinAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Asin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Asin");
  errs = vAsinAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Asin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Asin");
  errs = vAsinAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Asin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Asin");
  errs = vAsinAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Asin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Asin",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
