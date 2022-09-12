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
 *            Tanh example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3F7FFFD7 } }, //  0: vsTanh ( 6.80375481      ) = ( 0.999997556     );
{ { 0xC007309A }, { 0xBF789E12 } }, //  1: vsTanh ( -2.1123414      ) = ( -0.971161962    );
{ { 0x40B52EFA }, { 0x3F7FFE6B } }, //  2: vsTanh ( 5.66198444      ) = ( 0.99997586      );
{ { 0x40BF006A }, { 0x3F7FFF25 } }, //  3: vsTanh ( 5.96880054      ) = ( 0.999986947     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FEFFFFAD5FE7BF2 } }, //  0: vdTanh ( 6.80375434309419092       ) = ( 0.999997537572083539      );
{ { 0xC000E6134801CC26 }, { 0xBFEF13C23C2086E3 } }, //  1: vdTanh ( -2.11234146361813924      ) = ( -0.971161954341564715     );
{ { 0x4016A5DF421D4BBE }, { 0x3FEFFFCD557AD24F } }, //  2: vdTanh ( 5.66198447517211711       ) = ( 0.999975840523413484      );
{ { 0x4017E00D485FC01A }, { 0x3FEFFFE491F87BD2 } }, //  3: vdTanh ( 5.96880066952146571       ) = ( 0.999986920451073624      );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xBF7C29D8, 0x3CCBCC3D } }, //  0: vcTanh ( -2.1123414      + i * 6.80375481      ) = ( -0.985013485    + i * 0.0248776618    );
{ { 0x40BF006A, 0x40B52EFA }, { 0x3F7FFFB9, 0xB74FB664 } }, //  1: vcTanh ( 5.96880054      + i * 5.66198444      ) = ( 0.999995768     + i * -1.23806276e-05 );
{ { 0xC0C1912F, 0x4103BA28 }, { 0xBF800044, 0xB7008037 } }, //  2: vcTanh ( -6.04897261     + i * 8.2329483       ) = ( -1.00000811     + i * -7.65924688e-06 );
{ { 0x40ABAABC, 0xC052EA36 }, { 0x3F7FFD44, 0xB75EA892 } }, //  3: vcTanh ( 5.3645916       + i * -3.2955451      ) = ( 0.999958277     + i * -1.32714795e-05 );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xBFEF853AE5BB4742, 0x3F997986725FA047 } }, //  0: vzTanh ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -0.985013436026044298     + i * 0.0248776442821567988     );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x3FEFFFF7272D333E, 0xBEE9F6CBFCEC114D } }, //  1: vzTanh ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 0.999995781437611475      + i * -1.23806238696641625e-05  );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0xBFF000087C283048, 0xBEE0100456BD42E8 } }, //  2: vzTanh ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -1.00000809191534934      + i * -7.65922842274739419e-06  );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0x3FEFFFA87AF27428, 0xBEEBD50EB5E0D319 } }, //  3: vzTanh ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 0.999958267336869433      + i * -1.32714537209671191e-05  );
}

};

//!
//! @brief Single precision test
//!

int vTanhAccuracyLiteTest_float() {
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
        vmsTanh(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Tanh", acc);
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

int vTanhAccuracyLiteTest_double() {
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
        vmdTanh(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Tanh", acc);
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

int vTanhAccuracyLiteTest_float_complex() {
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
        vmcTanh(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Tanh", acc);
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

int vTanhAccuracyLiteTest_double_complex() {
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
        vmzTanh(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Tanh", acc);
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

  printf("Running %s functions:\n", "Tanh");

  printf("\tRunning %s with single precision real data type:\n", "Tanh");
  errs = vTanhAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Tanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Tanh");
  errs = vTanhAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Tanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Tanh");
  errs = vTanhAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Tanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Tanh");
  errs = vTanhAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Tanh",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Tanh",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
