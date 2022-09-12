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
 *            Mul example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 32.0, 32.0, 32.0, }, // VM_COMPLEX8
    { 32.0, 32.0, 32.0, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_3_t data =
{

{

{ { 0x40D9B85C }, { 0xC007309A }, { 0xC165F31C } }, //  0: vsMul ( 6.80375481     , -2.1123414      ) = ( -14.3718529     );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0x42072E58 } }, //  1: vsMul ( 5.66198444     , 5.96880054      ) = ( 33.7952576      );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0xC247341A } }, //  2: vsMul ( 8.2329483      , -6.04897261     ) = ( -49.8008804     );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0xC18D6F1C } }, //  3: vsMul ( -3.2955451     , 5.3645916       ) = ( -17.6792526     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0xC02CBE63704FA37B } }, //  0: vdMul ( 6.80375434309419092      , -2.11234146361813924      ) = ( -14.3718524071898539      );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0x4040E5CAF8EF8918 } }, //  1: vdMul ( 5.66198447517211711      , 5.96880066952146571       ) = ( 33.7952567262274783       );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0xC048E682F866802F } }, //  2: vdMul ( 8.23294715873568705      , -6.04897261413232101      ) = ( -49.8008718967906745      );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0xC031ADE38CCD2028 } }, //  3: vdMul ( -3.2955448857022196      , 5.36459189623808186       ) = ( -17.6792533875269839      );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40BF006A, 0x40B52EFA }, { 0xC24C860A, 0x41E533A2 } }, //  0: vcMul ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( -51.1308975     + i * 28.6502113      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x40ABAABC, 0xC052EA36 }, { 0xC0AA2ED2, 0x428033BF } }, //  1: vcMul ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -5.31821537     + i * 64.1010666      );
{ { 0x3F8A29C0, 0xC08E3964 }, { 0x4024F46C, 0xBEE77440 }, { 0x3F45DBCD, 0xC13F17C4 } }, //  2: vcMul ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( 0.772885144     + i * -11.9433022     );
{ { 0x3E8939C0, 0xC02D136C }, { 0x41052EB4, 0x4110B6A8 }, { 0x41D585D7, 0xC1A0B0BB } }, //  3: vcMul ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( 26.6903515      + i * -20.0862942     );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0xC04990C13850811A, 0x403CA674173EC41B } }, //  0: vzMul ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( -51.1308966057860772      + i * 28.6502089050519366       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC01545DC22B5AAB8, 0x40500677CE4FD3E3 } }, //  1: vzMul ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -5.31822256311232167      + i * 64.1010623721663677       );
{ { 0x3FF1453801E28A70, 0xC011C72C86338E59 }, { 0x40049E8D96893D1C, 0xBFDCEE88B739DD20 }, { 0x3FE8BB786C331F6C, 0xC027E2F8AB9E2201 } }, //  2: vzMul ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( 0.772884570434127394      + i * -11.9433034544499623      );
{ { 0x3FD12735D3224E60, 0xC005A26D910B44DC }, { 0x4020A5D666294BAC, 0x402216D5173C2DAA }, { 0x403AB0BABC96CCD0, 0xC0341617A44203A2 } }, //  3: vzMul ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( 26.6903493755497152       + i * -20.0862982426803072      );
}

};

//!
//! @brief Single precision test
//!

int vMulAccuracyLiteTest_float() {
  int errs = 0;
  float * varg1 = (float *)malloc(VLEN * sizeof(float));

  float * varg2 = (float *)malloc(VLEN * sizeof(float));
  float * vres1 = (float *)malloc(VLEN * sizeof(float));
  float * vref1 = (float *)malloc(VLEN * sizeof(float));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f32[i].v1.f;
      varg2[i] = data.data_f32[i].v2.f;
      vref1[i] = data.data_f32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmsMul(VLEN, (const float *)varg1, (const float *)varg2, (float *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Mul", acc);
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

int vMulAccuracyLiteTest_double() {
  int errs = 0;
  double * varg1 = (double *)malloc(VLEN * sizeof(double));

  double * varg2 = (double *)malloc(VLEN * sizeof(double));
  double * vres1 = (double *)malloc(VLEN * sizeof(double));
  double * vref1 = (double *)malloc(VLEN * sizeof(double));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_f64[i].v1.f;
      varg2[i] = data.data_f64[i].v2.f;
      vref1[i] = data.data_f64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmdMul(VLEN, (const double *)varg1, (const double *)varg2,
               (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Mul", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex single precision test
//!

int vMulAccuracyLiteTest_float_complex() {
  int errs = 0;
  VM_COMPLEX8 * varg1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  VM_COMPLEX8 * varg2 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vres1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));
  VM_COMPLEX8 * vref1 = (VM_COMPLEX8 *)malloc(VLEN * sizeof(VM_COMPLEX8));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c32[i].v1.f;
      varg2[i] = data.data_c32[i].v2.f;
      vref1[i] = data.data_c32[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmcMul(VLEN, (const MKL_Complex8 *)varg1, (const MKL_Complex8 *)varg2,
               (MKL_Complex8 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Mul", acc);
    }
  }

  free(varg1);

  free(varg2);

  free(vres1);
  free(vref1);

  return errs;
}
//!
//! @brief Complex double precision test
//!

int vMulAccuracyLiteTest_double_complex() {
  int errs = 0;
  VM_COMPLEX16 * varg1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  VM_COMPLEX16 * varg2 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vres1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));
  VM_COMPLEX16 * vref1 = (VM_COMPLEX16 *)malloc(VLEN * sizeof(VM_COMPLEX16));

  {
    for (int i = 0; i < VLEN; ++i) {

      varg1[i] = data.data_c64[i].v1.f;
      varg2[i] = data.data_c64[i].v2.f;
      vref1[i] = data.data_c64[i].v3.f;
    }
  }

  for (int acc = 0; acc < ACCURACY_NUM; ++acc) {

#pragma omp target data map(to:varg1[0:VLEN]) map(to:varg2[0:VLEN]) map(tofrom:vres1[0:VLEN]) device(dnum)
    {
#pragma omp target variant dispatch device(dnum) use_device_ptr(varg1, varg2, vres1) nowait
      {
        vmzMul(VLEN, (const MKL_Complex16 *)varg1, (const MKL_Complex16 *)varg2,
               (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Mul", acc);
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

  printf("Running %s functions:\n", "Mul");

  printf("\tRunning %s with single precision real data type:\n", "Mul");
  errs = vMulAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Mul",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Mul");
  errs = vMulAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Mul",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Mul");
  errs = vMulAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Mul",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Mul");
  errs = vMulAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Mul",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Mul",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
