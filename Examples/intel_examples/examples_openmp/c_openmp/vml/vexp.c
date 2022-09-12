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
 *            Exp example program text (OpenMP offload interface)
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

{

{ { 0x40D9B85C }, { 0x44614E64 } }, //  0: vsExp ( 6.80375481      ) = ( 901.224854      );
{ { 0xC007309A }, { 0x3DF7B6F5 } }, //  1: vsExp ( -2.1123414      ) = ( 0.120954432     );
{ { 0x40B52EFA }, { 0x438FDC09 } }, //  2: vsExp ( 5.66198444      ) = ( 287.719025      );
{ { 0x40BF006A }, { 0x43C384A7 } }, //  3: vsExp ( 5.96880054      ) = ( 391.036346      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x408C29CBAE86A8F2 } }, //  0: vdExp ( 6.80375434309419092       ) = ( 901.224453975706865       );
{ { 0xC000E6134801CC26 }, { 0x3FBEF6DE805CD39C } }, //  1: vdExp ( -2.11234146361813924      ) = ( 0.120954424227262824      );
{ { 0x4016A5DF421D4BBE }, { 0x4071FB813828BC71 } }, //  2: vdExp ( 5.66198447517211711       ) = ( 287.719047698140514       );
{ { 0x4017E00D485FC01A }, { 0x40787095200E2E6A } }, //  3: vdExp ( 5.96880066952146571       ) = ( 391.036407523532603       );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x3DD6E6C3, 0x3D7669F0 } }, //  0: vcExp ( -2.1123414      + i * 6.80375481      ) = ( 0.104932331     + i * 0.0601596236    );
{ { 0x40BF006A, 0x40B52EFA }, { 0x439EFDD3, 0xC3639680 } }, //  1: vcExp ( 5.96880054      + i * 5.66198444      ) = ( 317.983002      + i * -227.587891     );
{ { 0xC0C1912F, 0x4103BA28 }, { 0xBA64E852, 0x3B0FB55F } }, //  2: vcExp ( -6.04897261     + i * 8.2329483       ) = ( -0.000873212819 + i * 0.0021928174    );
{ { 0x40ABAABC, 0xC052EA36 }, { 0xC3532D29, 0x420314ED } }, //  3: vcExp ( 5.3645916       + i * -3.2955451      ) = ( -211.176407     + i * 32.7704353      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x3FBADCD8C17B13EC, 0x3FAECD3C3BF2B352 } }, //  0: vzExp ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( 0.104932353245274335      + i * 0.0601595709971208953     );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x4073DFBA8A9A56C2, 0xC06C72D028B4CC5F } }, //  1: vzExp ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 317.983042338262862       + i * -227.58791003524945       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0xBF4C9D04E55315EA, 0x3F61F6AC60E93AF9 } }, //  2: vzExp ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -0.000873210325743554007  + i * 0.00219281833350437301    );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC06A65A5A71D5FE7, 0x4040629C7263E09C } }, //  3: vzExp ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -211.176471288082411      + i * 32.7703993785555383       );
}

};

//!
//! @brief Single precision test
//!

int vExpAccuracyLiteTest_float() {
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
        vmsExp(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Exp", acc);
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

int vExpAccuracyLiteTest_double() {
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
        vmdExp(VLEN, (const double *)varg1, (double *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Exp", acc);
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

int vExpAccuracyLiteTest_float_complex() {
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
        vmcExp(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Exp", acc);
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

int vExpAccuracyLiteTest_double_complex() {
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
        vmzExp(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Exp", acc);
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

  printf("Running %s functions:\n", "Exp");

  printf("\tRunning %s with single precision real data type:\n", "Exp");
  errs = vExpAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Exp",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Exp");
  errs = vExpAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Exp",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Exp");
  errs = vExpAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Exp",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Exp");
  errs = vExpAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Exp",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Exp",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
