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
 *            CIS example program text (OpenMP offload interface)
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
data_2fc_t data =
{
{

{ { 0x40D9B85C }, { 0x3F5E16D8, 0x3EFEA7D8 } }, //  0: vcCIS ( 6.80375481      ) = ( 0.867536068     + i * 0.497374296     );
{ { 0xC007309A }, { 0xBF03F53A, 0xBF5B5EAB } }, //  1: vcCIS ( -2.1123414      ) = ( -0.51546061     + i * -0.856913269    );
{ { 0x40B52EFA }, { 0x3F502C93, 0xBF14FEBF } }, //  2: vcCIS ( 5.66198444      ) = ( 0.813180149     + i * -0.582012117    );
{ { 0x40BF006A }, { 0x3F7373DF, 0xBE9E5396 } }, //  3: vcCIS ( 5.96880054      ) = ( 0.950986803     + i * -0.30923146     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FEBC2DB7AB89950, 0x3FDFD4F93E99B2E0 } }, //  0: vzCIS ( 6.80375434309419092       ) = ( 0.867536296548488295      + i * 0.497373877652348639      );
{ { 0xC000E6134801CC26 }, { 0xBFE07EA757C4010B, 0xBFEB6BD5549D70BC } }, //  1: vzCIS ( -2.11234146361813924      ) = ( -0.51546065465666524      + i * -0.85691324735991925      );
{ { 0x4016A5DF421D4BBE }, { 0x3FEA05925DBF776B, 0xBFE29FD7C840E7D0 } }, //  2: vzCIS ( 5.66198447517211711       ) = ( 0.813180144406698502      + i * -0.582012072677793313     );
{ { 0x4017E00D485FC01A }, { 0x3FEE6E7BF8882000, 0xBFD3CA723281D19B } }, //  3: vzCIS ( 5.96880066952146571       ) = ( 0.950986848271895724      + i * -0.309231328318924248     );
}

};

//!
//! @brief Complex single precision test
//!

int vCISAccuracyLiteTest_float_complex() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));
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
      { vmcCIS(VLEN, varg1, (MKL_Complex8 *)vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      {
        VM_COMPLEX8 carg1 = varg1[i] + I * 0.0f;
        errs += check_result_float_complex(ARG1R_RES1C, carg1, carg1, vres1[i],
                                           vres1[i], vref1[i], vref1[i], "CIS",
                                           acc);
      }
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

int vCISAccuracyLiteTest_double_complex() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));
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
      { vmzCIS(VLEN, varg1, (MKL_Complex16 *)vres1, accuracy_mode[acc]); }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      {
        VM_COMPLEX16 carg1 = varg1[i] + I * 0.0;
        errs += check_result_double_complex(ARG1R_RES1C, carg1, carg1, vres1[i],
                                            vres1[i], vref1[i], vref1[i], "CIS",
                                            acc);
      }
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

  printf("Running %s functions:\n", "CIS");
  printf("\tRunning %s with single precision complex data type:\n", "CIS");
  errs = vCISAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "CIS",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "CIS");
  errs = vCISAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "CIS",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "CIS",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
