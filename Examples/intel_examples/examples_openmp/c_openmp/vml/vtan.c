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
 *            Tan example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3F12C4FD } }, //  0: vsTan ( 6.80375481      ) = ( 0.573318303     );
{ { 0xC007309A }, { 0x3FD4CA42 } }, //  1: vsTan ( -2.1123414      ) = ( 1.66242242      );
{ { 0x40B52EFA }, { 0xBF3739A7 } }, //  2: vsTan ( 5.66198444      ) = ( -0.715723455    );
{ { 0x40BF006A }, { 0xBEA67C8E } }, //  3: vsTan ( 5.96880054      ) = ( -0.325169027    );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FE2589E45D858D9 } }, //  0: vdTan ( 6.80375434309419092       ) = ( 0.573317657867643438      );
{ { 0xC000E6134801CC26 }, { 0x3FFA99480B7BBB50 } }, //  1: vdTan ( -2.11234146361813924      ) = ( 1.66242222295450759       );
{ { 0x4016A5DF421D4BBE }, { 0xBFE6E734CC578392 } }, //  2: vdTan ( 5.66198447517211711       ) = ( -0.715723418336084771     );
{ { 0x4017E00D485FC01A }, { 0xBFD4CF91224985AD } }, //  3: vdTan ( 5.96880066952146571       ) = ( -0.325168879970159364     );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x3611FC00, 0x3F80000A } }, //  0: vcTan ( -2.1123414      + i * 6.80375481      ) = ( 2.1753367e-06   + i * 1.00000119      );
{ { 0x40BF006A, 0x40B52EFA }, { 0xB76E6472, 0x3F7FFEB8 } }, //  1: vcTan ( 5.96880054      + i * 5.66198444      ) = ( -1.42092922e-05 + i * 0.99998045      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x3388F265, 0x3F7FFFFE } }, //  2: vcTan ( -6.04897261     + i * 8.2329483       ) = ( 6.37708482e-08  + i * 0.999999881     );
{ { 0x40ABAABC, 0xC052EA36 }, { 0xBB2DAE75, 0xBF801793 } }, //  3: vcTan ( 5.3645916       + i * -3.2955451      ) = ( -0.00265016896  + i * -1.00071943     );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x3EC23F813BAC1EE3, 0x3FF0000135BF0F45 } }, //  0: vzTan ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( 2.17533894664502934e-06   + i * 1.00000115389498601       );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0xBEEDCC8D7ADCC53E, 0x3FEFFFD705FCDE2F } }, //  1: vzTan ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( -1.42092866003163373e-05  + i * 0.999980460829595574      );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x3E711E4F35CD0BC8, 0x3FEFFFFFBC562861 } }, //  2: vzTan ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 6.37709951070066822e-08   + i * 0.999999873967009845      );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xBF65B5CF78A52334, 0xBFF002F257CD6E64 } }, //  3: vzTan ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -0.00265017053348906052   + i * -1.00071939752424388      );
}

};

//!
//! @brief Single precision test
//!

int vTanAccuracyLiteTest_float() {
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
        vmsTan(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Tan", acc);
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

int vTanAccuracyLiteTest_double() {
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
        vmdTan(VLEN, (const double *)varg1, (double *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Tan", acc);
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

int vTanAccuracyLiteTest_float_complex() {
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
        vmcTan(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Tan", acc);
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

int vTanAccuracyLiteTest_double_complex() {
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
        vmzTan(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Tan", acc);
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

  printf("Running %s functions:\n", "Tan");

  printf("\tRunning %s with single precision real data type:\n", "Tan");
  errs = vTanAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Tan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Tan");
  errs = vTanAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Tan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Tan");
  errs = vTanAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Tan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Tan");
  errs = vTanAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Tan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Tan",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
