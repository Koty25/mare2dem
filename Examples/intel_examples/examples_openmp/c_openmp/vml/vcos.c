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
 *            Cos example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 4.0, 5.0, 5.0E3, }, // VM_COMPLEX8
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

{ { 0x40D9B85C }, { 0x3F5E16D8 } }, //  0: vsCos ( 6.80375481      ) = ( 0.867536068     );
{ { 0xC007309A }, { 0xBF03F53A } }, //  1: vsCos ( -2.1123414      ) = ( -0.51546061     );
{ { 0x40B52EFA }, { 0x3F502C93 } }, //  2: vsCos ( 5.66198444      ) = ( 0.813180149     );
{ { 0x40BF006A }, { 0x3F7373DF } }, //  3: vsCos ( 5.96880054      ) = ( 0.950986803     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FEBC2DB7AB89950 } }, //  0: vdCos ( 6.80375434309419092       ) = ( 0.867536296548488295      );
{ { 0xC000E6134801CC26 }, { 0xBFE07EA757C4010B } }, //  1: vdCos ( -2.11234146361813924      ) = ( -0.51546065465666524      );
{ { 0x4016A5DF421D4BBE }, { 0x3FEA05925DBF776B } }, //  2: vdCos ( 5.66198447517211711       ) = ( 0.813180144406698502      );
{ { 0x4017E00D485FC01A }, { 0x3FEE6E7BF8882000 } }, //  3: vdCos ( 5.96880066952146571       ) = ( 0.950986848271895724      );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xC36845F3, 0x43C11152 } }, //  0: vcCos ( -2.1123414      + i * 6.80375481      ) = ( -232.273239     + i * 386.135315      );
{ { 0x40BF006A, 0x40B52EFA }, { 0x4308CF67, 0x4231F100 } }, //  1: vcCos ( 5.96880054      + i * 5.66198444      ) = ( 136.810165      + i * 44.4853516      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x44E4C2CB, 0xC3DA5250 } }, //  2: vcCos ( -6.04897261     + i * 8.2329483       ) = ( 1830.08728      + i * -436.643066     );
{ { 0x40ABAABC, 0xC052EA36 }, { 0x41033D87, 0xC12B6150 } }, //  3: vcCos ( 5.3645916       + i * -3.2955451      ) = ( 8.20252132      + i * -10.7112579     );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xC06D08BDB8FC7F5D, 0x407822296A7AAF2C } }, //  0: vzCos ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -232.273159497412877      + i * 386.135111312137042       );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x406119ECE50B1B29, 0x40463E1F6EEA30C8 } }, //  1: vzCos ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 136.810167810145941       + i * 44.4853342669971994       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x409C98572C8DB5F8, 0xC07B4A47E4EA8F4F } }, //  2: vzCos ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 1830.08513089583539       + i * -436.642552295922485      );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0x402067B107A38AEA, 0xC0256C29633824E0 } }, //  3: vzCos ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 8.20252250548715622       + i * -10.7112532621417245      );
}

};

//!
//! @brief Single precision test
//!

int vCosAccuracyLiteTest_float() {
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
        vmsCos(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Cos", acc);
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

int vCosAccuracyLiteTest_double() {
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
        vmdCos(VLEN, (const double *)varg1, (double *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Cos", acc);
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

int vCosAccuracyLiteTest_float_complex() {
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
        vmcCos(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Cos", acc);
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

int vCosAccuracyLiteTest_double_complex() {
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
        vmzCos(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Cos", acc);
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

  printf("Running %s functions:\n", "Cos");

  printf("\tRunning %s with single precision real data type:\n", "Cos");
  errs = vCosAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Cos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Cos");
  errs = vCosAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Cos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Cos");
  errs = vCosAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Cos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Cos");
  errs = vCosAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Cos",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Cos",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
