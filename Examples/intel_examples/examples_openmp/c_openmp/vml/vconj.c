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
 *            Conj example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 2.0, 5.0, 5.0E3, }, // VM_COMPLEX8
    { 2.0, 5.0, 7.0E7, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_2_t data =
{

{ /* empty */ }

,

{ /* empty */}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xC007309A, 0xC0D9B85C } }, //  0: vcConj ( -2.1123414      + i * 6.80375481      ) = ( -2.1123414      + i * -6.80375481     );
{ { 0x40BF006A, 0x40B52EFA }, { 0x40BF006A, 0xC0B52EFA } }, //  1: vcConj ( 5.96880054      + i * 5.66198444      ) = ( 5.96880054      + i * -5.66198444     );
{ { 0xC0C1912F, 0x4103BA28 }, { 0xC0C1912F, 0xC103BA28 } }, //  2: vcConj ( -6.04897261     + i * 8.2329483       ) = ( -6.04897261     + i * -8.2329483      );
{ { 0x40ABAABC, 0xC052EA36 }, { 0x40ABAABC, 0x4052EA36 } }, //  3: vcConj ( 5.3645916       + i * -3.2955451      ) = ( 5.3645916       + i * 3.2955451       );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xC000E6134801CC26, 0xC01B370B60E66E18 } }, //  0: vzConj ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -2.11234146361813924      + i * -6.80375434309419092      );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A, 0xC016A5DF421D4BBE } }, //  1: vzConj ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 5.96880066952146571       + i * -5.66198447517211711      );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0xC0183225E080644C, 0xC0207744D998EE8A } }, //  2: vzConj ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -6.04897261413232101      + i * -8.23294715873568705      );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0, 0x400A5D46A314BA8E } }, //  3: vzConj ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 5.36459189623808186       + i * 3.2955448857022196        );
}

};

//!
//! @brief Complex single precision test
//!

int vConjAccuracyLiteTest_float_complex() {
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
        vmcConj(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Conj", acc);
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

int vConjAccuracyLiteTest_double_complex() {
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
        vmzConj(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Conj", acc);
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

  printf("Running %s functions:\n", "Conj");
  printf("\tRunning %s with single precision complex data type:\n", "Conj");
  errs = vConjAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Conj",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Conj");
  errs = vConjAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Conj",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Conj",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
