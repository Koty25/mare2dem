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
 *            Asinh example program text (OpenMP offload interface)
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

{ { 0x4106F102 }, { 0x40350CA0 } }, //  0: vsAsinh ( 8.4338398       ) = ( 2.82889557      );
{ { 0x40821418 }, { 0x40070FEC } }, //  1: vsAsinh ( 4.06495285      ) = ( 2.11034679      );
{ { 0x40FBFADC }, { 0x4030B06E } }, //  2: vsAsinh ( 7.87437248      ) = ( 2.76076841      );
{ { 0x41006539 }, { 0x4031E3DE } }, //  3: vsAsinh ( 8.02471256      ) = ( 2.77953291      );
}

,

{

{ { 0x4020DE203A4CEF74 }, { 0x4006A19409F6875D } }, //  0: vdAsinh ( 8.43383962811615362       ) = ( 2.82889564307781294       );
{ { 0x40104282F4C21EA0 }, { 0x4000E1FD71214ECF } }, //  1: vdAsinh ( 4.06495268282711208       ) = ( 2.11034668333909492       );
{ { 0x401F7F5B79FEFEB7 }, { 0x4006160DBBDDDD1B } }, //  2: vdAsinh ( 7.87437239283433765       ) = ( 2.7607683827478815        );
{ { 0x40200CA71821B2E8 }, { 0x40063C7BB36D8924 } }, //  3: vdAsinh ( 8.02471232806551882       ) = ( 2.77953281572367139       );
}

,

{

{ { 0x40821418, 0x4106F102 }, { 0x403B657B, 0x3F8F494C } }, //  0: vcAsinh ( 4.06495285      + i * 8.4338398       ) = ( 2.92806888      + i * 1.11942434      );
{ { 0x41006539, 0x40FBFADC }, { 0x40473A22, 0x3F462298 } }, //  1: vcAsinh ( 8.02471256      + i * 7.87437248      ) = ( 3.11292315      + i * 0.773965359     );
{ { 0x4008B448, 0x41122574 }, { 0x403B7891, 0x3FAB7EC8 } }, //  2: vcAsinh ( 2.13600349      + i * 9.13414383      ) = ( 2.92923379      + i * 1.33980656      );
{ { 0x40F7511A, 0x405F0D3D }, { 0x40354EE7, 0x3ED793C5 } }, //  3: vcAsinh ( 7.72865009      + i * 3.485183        ) = ( 2.83294082      + i * 0.421049267     );
}

,

{

{ { 0x40104282F4C21EA0, 0x4020DE203A4CEF74 }, { 0x40076CAF627D52BE, 0x3FF1E929808BEB84 } }, //  0: vzAsinh ( 4.06495268282711208       + i * 8.43383962811615362       ) = ( 2.92806889481502619       + i * 1.11942434514523459       );
{ { 0x40200CA71821B2E8, 0x401F7F5B79FEFEB7 }, { 0x4008E7443A01B5E4, 0x3FE8C452F5B6F581 } }, //  1: vzAsinh ( 8.02471232806551882       + i * 7.87437239283433765       ) = ( 3.11292310064048827       + i * 0.773965339576236144      );
{ { 0x40011688F5E89379, 0x402244AE8957BC90 }, { 0x40076F1232EA82C3, 0x3FF56FD905B7041C } }, //  2: vzAsinh ( 2.13600341907516311       + i * 9.13414410778048591       ) = ( 2.92923393037958268       + i * 1.33980657799134573       );
{ { 0x401EEA233BB5D447, 0x400BE1A7A0BAF683 }, { 0x4006A9DCD22D5BB3, 0x3FDAF278A674CEAD } }, //  3: vzAsinh ( 7.72865002915666022       + i * 3.4851830060059128        ) = ( 2.83294071389124147       + i * 0.421049273066482155      );
}

};

//!
//! @brief Single precision test
//!

int vAsinhAccuracyLiteTest_float() {
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
        vmsAsinh(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Asinh", acc);
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

int vAsinhAccuracyLiteTest_double() {
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
        vmdAsinh(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Asinh", acc);
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

int vAsinhAccuracyLiteTest_float_complex() {
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
        vmcAsinh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float_complex(ARG1_RES1, varg1[i], varg1[i],
                                         vres1[i], vres1[i], vref1[i], vref1[i],
                                         "Asinh", acc);
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

int vAsinhAccuracyLiteTest_double_complex() {
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
        vmzAsinh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Asinh", acc);
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

  printf("Running %s functions:\n", "Asinh");

  printf("\tRunning %s with single precision real data type:\n", "Asinh");
  errs = vAsinhAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Asinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Asinh");
  errs = vAsinhAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Asinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Asinh");
  errs = vAsinhAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Asinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Asinh");
  errs = vAsinhAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Asinh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Asinh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
