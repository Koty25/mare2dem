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
 *            Div example program text (OpenMP offload interface)
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
data_3_t data =
{

{

{ { 0x40D9B85C }, { 0xC007309A }, { 0xC04E241D } }, //  0: vsDiv ( 6.80375481     , -2.1123414      ) = ( -3.22095418     );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0x3F72D73C } }, //  1: vsDiv ( 5.66198444     , 5.96880054      ) = ( 0.948596716     );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0xBFAE36DB } }, //  2: vsDiv ( 8.2329483      , -6.04897261     ) = ( -1.36104906     );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0xBF1D43B3 } }, //  3: vsDiv ( -3.2955451     , 5.3645916       ) = ( -0.614314258    );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0xC009C483720BB98F } }, //  0: vdDiv ( 6.80375434309419092      , -2.11234146361813924      ) = ( -3.22095383737832419      );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0x3FEE5AE76A98B075 } }, //  1: vdDiv ( 5.66198447517211711      , 5.96880066952146571       ) = ( 0.948596676059891508      );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0xBFF5C6DB25F67A68 } }, //  2: vdDiv ( 8.23294715873568705      , -6.04897261413232101      ) = ( -1.36104883984776315      );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0xBFE3A876377717F0 } }, //  3: vdDiv ( -3.2955448857022196      , 5.36459189623808186       ) = ( -0.614314182596670477     );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40BF006A, 0x40B52EFA }, { 0x3EC407E7, 0x3F46D575 } }, //  0: vcDiv ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( 0.38287279      + i * 0.776694596     );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x40ABAABC, 0xC052EA36 }, { 0xBFC065C9, 0x3F1C7E64 } }, //  1: vcDiv ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -1.50310624     + i * 0.611303568     );
{ { 0x3F8A29C0, 0xC08E3964 }, { 0x4024F46C, 0xBEE77440 }, { 0x3F33205B, 0xBFCD03CA } }, //  2: vcDiv ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( 0.699712455     + i * -1.60167813     );
{ { 0x3E8939C0, 0xC02D136C }, { 0x41052EB4, 0x4110B6A8 }, { 0xBE16A63A, 0xBE28FD4F } }, //  3: vcDiv ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( -0.147118479    + i * -0.165028796    );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x3FD880FC9B51FFC7, 0x3FE8DAAE83F56877 } }, //  0: vzDiv ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( 0.382872726135243757      + i * 0.776694543582620578      );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xBFF80CB8F313B756, 0x3FE38FCC4B51C8A8 } }, //  1: vzDiv ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -1.50310606910666911      + i * 0.61130346976121519       );
{ { 0x3FF1453801E28A70, 0xC011C72C86338E59 }, { 0x40049E8D96893D1C, 0xBFDCEE88B739DD20 }, { 0x3FE6640B85CE8ECF, 0xBFF9A079189F1A3A } }, //  2: vzDiv ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( 0.699712525693451215      + i * -1.60167798631449765      );
{ { 0x3FD12735D3224E60, 0xC005A26D910B44DC }, { 0x4020A5D666294BAC, 0x402216D5173C2DAA }, { 0xBFC2D4C79C402529, 0xBFC51FA9A03BFCDC } }, //  3: vzDiv ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( -0.147118521970960786     + i * -0.1650287659067321       );
}

};

//!
//! @brief Single precision test
//!

int vDivAccuracyLiteTest_float() {
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
        vmsDiv(VLEN, (const float *)varg1, (const float *)varg2, (float *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Div", acc);
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

int vDivAccuracyLiteTest_double() {
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
        vmdDiv(VLEN, (const double *)varg1, (const double *)varg2,
               (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Div", acc);
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

int vDivAccuracyLiteTest_float_complex() {
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
        vmcDiv(VLEN, (const MKL_Complex8 *)varg1, (const MKL_Complex8 *)varg2,
               (MKL_Complex8 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Div", acc);
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

int vDivAccuracyLiteTest_double_complex() {
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
        vmzDiv(VLEN, (const MKL_Complex16 *)varg1, (const MKL_Complex16 *)varg2,
               (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Div", acc);
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

  printf("Running %s functions:\n", "Div");

  printf("\tRunning %s with single precision real data type:\n", "Div");
  errs = vDivAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Div",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Div");
  errs = vDivAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Div",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Div");
  errs = vDivAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Div",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Div");
  errs = vDivAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Div",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Div",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
