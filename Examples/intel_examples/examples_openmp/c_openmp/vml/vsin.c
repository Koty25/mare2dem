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
 *            Sin example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3EFEA7D8 } }, //  0: vsSin ( 6.80375481      ) = ( 0.497374296     );
{ { 0xC007309A }, { 0xBF5B5EAB } }, //  1: vsSin ( -2.1123414      ) = ( -0.856913269    );
{ { 0x40B52EFA }, { 0xBF14FEBF } }, //  2: vsSin ( 5.66198444      ) = ( -0.582012117    );
{ { 0x40BF006A }, { 0xBE9E5396 } }, //  3: vsSin ( 5.96880054      ) = ( -0.30923146     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FDFD4F93E99B2E0 } }, //  0: vdSin ( 6.80375434309419092       ) = ( 0.497373877652348639      );
{ { 0xC000E6134801CC26 }, { 0xBFEB6BD5549D70BC } }, //  1: vdSin ( -2.11234146361813924      ) = ( -0.85691324735991925      );
{ { 0x4016A5DF421D4BBE }, { 0xBFE29FD7C840E7D0 } }, //  2: vdSin ( 5.66198447517211711       ) = ( -0.582012072677793313     );
{ { 0x4017E00D485FC01A }, { 0xBFD3CA723281D19B } }, //  3: vdSin ( 5.96880066952146571       ) = ( -0.309231328318924248     );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xC3C11171, 0xC36845CE } }, //  0: vcSin ( -2.1123414      + i * 6.80375481      ) = ( -386.136261     + i * -232.272675     );
{ { 0x40BF006A, 0x40B52EFA }, { 0xC231F219, 0x4308CE8E } }, //  1: vcSin ( 5.96880054      + i * 5.66198444      ) = ( -44.4864235     + i * 136.806854      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x43DA5252, 0x44E4C2C8 } }, //  2: vcSin ( -6.04897261     + i * 8.2329483       ) = ( 436.643127      + i * 1830.08691      );
{ { 0x40ABAABC, 0xC052EA36 }, { 0xC12BD9EA, 0xC102E16D } }, //  3: vcSin ( 5.3645916       + i * -3.2955451      ) = ( -10.7407017     + i * -8.18003559     );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xC078222D4F7FC3B9, 0xC06D08B90982015C } }, //  0: vzSin ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -386.136062144356004      + i * -232.272587541500684      );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0xC0463E42A6B9A9E8, 0x406119D1D1680E9B } }, //  1: vzSin ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( -44.4864090353547112      + i * 136.806862548099929       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x407B4A4825964FE6, 0x409C9856E8CA5287 } }, //  2: vzSin ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( 436.642613970905927       + i * 1830.08487239960391       );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC0257B3CA2F0A873, 0xC0205C2DC648E14A } }, //  3: vzSin ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( -10.7406969946643809      + i * -8.18003673209809179      );
}

};

//!
//! @brief Single precision test
//!

int vSinAccuracyLiteTest_float() {
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
        vmsSin(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Sin", acc);
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

int vSinAccuracyLiteTest_double() {
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
        vmdSin(VLEN, (const double *)varg1, (double *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Sin", acc);
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

int vSinAccuracyLiteTest_float_complex() {
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
        vmcSin(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Sin", acc);
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

int vSinAccuracyLiteTest_double_complex() {
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
        vmzSin(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_double_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Sin", acc);
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

  printf("Running %s functions:\n", "Sin");

  printf("\tRunning %s with single precision real data type:\n", "Sin");
  errs = vSinAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Sin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Sin");
  errs = vSinAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Sin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Sin");
  errs = vSinAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Sin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Sin");
  errs = vSinAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Sin",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Sin",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
