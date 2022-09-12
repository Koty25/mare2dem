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
 *            Atanh example program text (OpenMP offload interface)
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

{ { 0x3C82EB10 }, { 0x3C82EDEB } }, //  0: vsAtanh ( 0.0159812272    ) = ( 0.0159825888    );
{ { 0x3D780F8E }, { 0x3D785D5D } }, //  1: vsAtanh ( 0.0605617091    ) = ( 0.0606359132    );
{ { 0x3CB1AF64 }, { 0x3CB1B687 } }, //  2: vsAtanh ( 0.0216900781    ) = ( 0.0216934811    );
{ { 0x3CA51E30 }, { 0x3CA523EA } }, //  3: vsAtanh ( 0.0201559961    ) = ( 0.0201587267    );
}

,

{

{ { 0x3F905D621353EDF8 }, { 0x3F905DBD64B240D8 } }, //  0: vdAtanh ( 0.0159812282845290532     ) = ( 0.0159825890264649606     );
{ { 0x3FAF01F1B0A46A4A }, { 0x3FAF0BAB94755F2D } }, //  1: vdAtanh ( 0.060561707318090699      ) = ( 0.0606359118198141825     );
{ { 0x3F9635EC782C6BD8 }, { 0x3F9636D0CCDE1024 } }, //  2: vdAtanh ( 0.0216900776241394089     ) = ( 0.0216934800187261884     );
{ { 0x3F94A3C609C2E128 }, { 0x3F94A47D42904F8C } }, //  3: vdAtanh ( 0.020155996652392677      ) = ( 0.0201587268712298123     );
}

,

{

{ { 0x3D780F8E, 0x3C82EB10 }, { 0x3D784D08, 0x3C836386 } }, //  0: vcAtanh ( 0.0605617091    + i * 0.0159812272    ) = ( 0.0606203377    + i * 0.0160386674    );
{ { 0x3CA51E30, 0x3CB1AF64 }, { 0x3CA51005, 0x3CB1BABB } }, //  1: vcAtanh ( 0.0201559961    + i * 0.0216900781    ) = ( 0.0201492403    + i * 0.0216954853    );
{ { 0x3DA4576B, 0x3C10C1C8 }, { 0x3DA4AEBF, 0x3C11B0F3 } }, //  2: vcAtanh ( 0.0802448615    + i * 0.00883526355   ) = ( 0.0804114267    + i * 0.00889228564   );
{ { 0x3CBDDDC8, 0x3D88257A }, { 0x3CBD1067, 0x3D8804D6 } }, //  3: vcAtanh ( 0.0231770426    + i * 0.0664777309    ) = ( 0.0230791103    + i * 0.0664154738    );
}

,

{

{ { 0x3FAF01F1B0A46A4A, 0x3F905D621353EDF8 }, { 0x3FAF09A0E2B69FE1, 0x3F906C70CDA5C815 } }, //  0: vzAtanh ( 0.060561707318090699      + i * 0.0159812282845290532     ) = ( 0.060620334315274034      + i * 0.0160386682050060632     );
{ { 0x3F94A3C609C2E128, 0x3F9635EC782C6BD8 }, { 0x3F94A200AD66E706, 0x3F9637574ED5D99F } }, //  1: vzAtanh ( 0.020155996652392677      + i * 0.0216900776241394089     ) = ( 0.020149241050353893      + i * 0.0216954843394546702     );
{ { 0x3FB48AED668F7C42, 0x3F821839168A96D8 }, { 0x3FB495D7D955D65A, 0x3F82361E797C0E88 } }, //  2: vzAtanh ( 0.0802448630706616151     + i * 0.00883526420632156639    ) = ( 0.0804114251712574613     + i * 0.00889228637925688903    );
{ { 0x3F97BBB8DC2F7770, 0x3FB104AF24553C92 }, { 0x3F97A20CB1D99626, 0x3FB1009AB2439241 } }, //  3: vzAtanh ( 0.0231770405188095885     + i * 0.0664777244285111035     ) = ( 0.023079107623195004      + i * 0.0664154706206057238     );
}

};

//!
//! @brief Single precision test
//!

int vAtanhAccuracyLiteTest_float() {
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
        vmsAtanh(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Atanh", acc);
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

int vAtanhAccuracyLiteTest_double() {
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
        vmdAtanh(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Atanh", acc);
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

int vAtanhAccuracyLiteTest_float_complex() {
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
        vmcAtanh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float_complex(ARG1_RES1, varg1[i], varg1[i],
                                         vres1[i], vres1[i], vref1[i], vref1[i],
                                         "Atanh", acc);
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

int vAtanhAccuracyLiteTest_double_complex() {
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
        vmzAtanh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Atanh", acc);
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

  printf("Running %s functions:\n", "Atanh");

  printf("\tRunning %s with single precision real data type:\n", "Atanh");
  errs = vAtanhAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Atanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Atanh");
  errs = vAtanhAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Atanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Atanh");
  errs = vAtanhAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Atanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Atanh");
  errs = vAtanhAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Atanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Atanh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
