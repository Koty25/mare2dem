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
 *            LinearFrac example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 1.0E6, 1.0E6, 1.0E6, }, // float
    { 1.0E6, 1.0E6, 1.0E6, }, // double
    { 2.0, 5.0, 5.0E3, }, // VM_COMPLEX8
    { 2.0, 5.0, 7.0E7, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_7_t data =
{

{

{ { 0x40D9B85C }, { 0xC007309A }, { 0x4048F5C3 }, { 0x40C8F5C3 }, { 0x4116B852 }, { 0x4148F5C3 }, { 0xC07117D3 } }, //  0: vsLinearFrac ( 6.80375481     , -2.1123414     , 3.1400001      , 6.28000021     , 9.42000008     , 12.5600004      ) = ( -3.76707911     );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0x4048F5C3 }, { 0x40C8F5C3 }, { 0x4116B852 }, { 0x4148F5C3 }, { 0x3EB313C1 } }, //  1: vsLinearFrac ( 5.66198444     , 5.96880054     , 3.1400001      , 6.28000021     , 9.42000008     , 12.5600004      ) = ( 0.349760085     );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0x4048F5C3 }, { 0x40C8F5C3 }, { 0x4116B852 }, { 0x4148F5C3 }, { 0xBF392C6C } }, //  2: vsLinearFrac ( 8.2329483      , -6.04897261    , 3.1400001      , 6.28000021     , 9.42000008     , 12.5600004      ) = ( -0.723334074    );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0x4048F5C3 }, { 0x40C8F5C3 }, { 0x4116B852 }, { 0x4148F5C3 }, { 0xBD840B70 } }, //  3: vsLinearFrac ( -3.2955451     , 5.3645916      , 3.1400001      , 6.28000021     , 9.42000008     , 12.5600004      ) = ( -0.0644749403   );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0x40091EB851EB851F }, { 0x40191EB851EB851F }, { 0x4022D70A3D70A3D7 }, { 0x40291EB851EB851F }, { 0xC00E22FA0DD4CC16 } }, //  0: vdLinearFrac ( 6.80375434309419092      , -2.11234146361813924     , 3.14000000000000012      , 6.28000000000000025      , 9.41999999999999993      , 12.5600000000000005       ) = ( -3.76707850270896483      );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0x40091EB851EB851F }, { 0x40191EB851EB851F }, { 0x4022D70A3D70A3D7 }, { 0x40291EB851EB851F }, { 0x3FD6627804818A52 } }, //  1: vdLinearFrac ( 5.66198447517211711      , 5.96880066952146571      , 3.14000000000000012      , 6.28000000000000025      , 9.41999999999999993      , 12.5600000000000005       ) = ( 0.349760059738547402      );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0x40091EB851EB851F }, { 0x40191EB851EB851F }, { 0x4022D70A3D70A3D7 }, { 0x40291EB851EB851F }, { 0xBFE7258D6AC8E3FB } }, //  2: vdLinearFrac ( 8.23294715873568705      , -6.04897261413232101     , 3.14000000000000012      , 6.28000000000000025      , 9.41999999999999993      , 12.5600000000000005       ) = ( -0.723334034503863577     );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0x40091EB851EB851F }, { 0x40191EB851EB851F }, { 0x4022D70A3D70A3D7 }, { 0x40291EB851EB851F }, { 0xBFB0816DEA262188 } }, //  3: vdLinearFrac ( -3.2955448857022196      , 5.36459189623808186      , 3.14000000000000012      , 6.28000000000000025      , 9.41999999999999993      , 12.5600000000000005       ) = ( -0.0644749352123935582    );
}

,
};

//!
//! @brief Single precision test
//!

int vLinearFracAccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));

  float * varg2 = (float *)malloc(VLEN * sizeof(float));

  float varg3;
  float varg4;
  float varg5;
  float varg6;

  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f32[i].v1.f;
      varg2[i] = data.data_f32[i].v2.f;
      varg3 = data.data_f32[i].v3.f;
      varg4 = data.data_f32[i].v4.f;
      varg5 = data.data_f32[i].v5.f;
      varg6 = data.data_f32[i].v6.f;
      vref1[i] = data.data_f32[i].v7.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmsLinearFrac(VLEN, (const float *)varg1, (const float *)varg2,
                      (const float)varg3, (const float)varg4,
                      (const float)varg5, (const float)varg6, (float *)vres1,
                      accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i], vres1[i],
                             vref1[i], vref1[i], "LinearFrac", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Double precision test
//!

int vLinearFracAccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));

  double * varg2 = (double *)malloc(VLEN * sizeof(double));

  double varg3;
  double varg4;
  double varg5;
  double varg6;

  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {
      varg1[i] = data.data_f64[i].v1.f;
      varg2[i] = data.data_f64[i].v2.f;
      varg3 = data.data_f64[i].v3.f;
      varg4 = data.data_f64[i].v4.f;
      varg5 = data.data_f64[i].v5.f;
      varg6 = data.data_f64[i].v6.f;
      vref1[i] = data.data_f64[i].v7.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {
#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmdLinearFrac(VLEN, (const double *)varg1, (const double *)varg2,
                      (const double)varg3, (const double)varg4,
                      (const double)varg5, (const double)varg6, (double *)vres1,
                      accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i], vres1[i],
                              vref1[i], vref1[i], "LinearFrac", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}

int main(int argc, char **argv) {
  int errs = 0;
  int total_errs = 0;

  printf("Running %s functions:\n", "LinearFrac");

  printf("\tRunning %s with single precision real data type:\n", "LinearFrac");
  errs = vLinearFracAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "LinearFrac",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "LinearFrac");
  errs = vLinearFracAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "LinearFrac",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "LinearFrac",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
