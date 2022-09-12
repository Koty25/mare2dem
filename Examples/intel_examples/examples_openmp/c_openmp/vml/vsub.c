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
 *            Sub example program text (OpenMP offload interface)
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
data_3_t data =
{

{

{ { 0x40D9B85C }, { 0xC007309A }, { 0x410EA854 } }, //  0: vsSub ( 6.80375481     , -2.1123414      ) = ( 8.91609573      );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0xBE9D1700 } }, //  1: vsSub ( 5.66198444     , 5.96880054      ) = ( -0.306816101    );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0x416482C0 } }, //  2: vsSub ( 8.2329483      , -6.04897261     ) = ( 14.2819214      );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0xC10A8FEC } }, //  3: vsSub ( -3.2955451     , 5.3645916       ) = ( -8.66013718     );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0x4021D50A8273AA16 } }, //  0: vdSub ( 6.80375434309419092      , -2.11234146361813924      ) = ( 8.91609580671233104       );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0xBFD3A2E0642745C0 } }, //  1: vdSub ( 5.66198447517211711      , 5.96880066952146571       ) = ( -0.306816194349348592     );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0x402C9057C9D920B0 } }, //  2: vdSub ( 8.23294715873568705      , -6.04897261413232101      ) = ( 14.2819197728680081       );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0xC02151FD72C2A3FC } }, //  3: vdSub ( -3.2955448857022196      , 5.36459189623808186       ) = ( -8.66013678194030234      );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40BF006A, 0x40B52EFA }, { 0xC1014C5C, 0x3F922588 } }, //  0: vcSub ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( -8.08114243     + i * 1.14177036      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x40ABAABC, 0xC052EA36 }, { 0xC1369DF6, 0x413874B6 } }, //  1: vcSub ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -11.4135647     + i * 11.5284939      );
{ { 0x3F8A29C0, 0xC08E3964 }, { 0x4024F46C, 0xBEE77440 }, { 0xBFBFBF18, 0xC07F8440 } }, //  2: vcSub ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( -1.49801922     + i * -3.9924469      );
{ { 0x3E8939C0, 0xC02D136C }, { 0x41052EB4, 0x4110B6A8 }, { 0xC100E4E6, 0xC13BFB83 } }, //  3: vcSub ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( -8.05588341     + i * -11.7489042     );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0xC020298B76305316, 0x3FF244B07B248968 } }, //  0: vzSub ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( -8.08114213313960406      + i * 1.1417698679220738        );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC026D3BEBA3DA77E, 0x40270E96825E1D2E } }, //  1: vzSub ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -11.4135645103704029      + i * 11.5284920444379075       );
{ { 0x3FF1453801E28A70, 0xC011C72C86338E59 }, { 0x40049E8D96893D1C, 0xBFDCEE88B739DD20 }, { 0xBFF7F7E32B2FEFC8, 0xC00FF087F57FE10E } }, //  2: vzSub ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( -1.49801937932987705      + i * -3.99244682117944816      );
{ { 0x3FD12735D3224E60, 0xC005A26D910B44DC }, { 0x4020A5D666294BAC, 0x402216D5173C2DAA }, { 0xC0201C9CB7903939, 0xC0277F707B7EFEE1 } }, //  3: vzSub ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( -8.05588315616170014      + i * -11.7489050476573897      );
}

};

//!
//! @brief Single precision test
//!

int vSubAccuracyLiteTest_float() {
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
        vmsSub(VLEN, (const float *)varg1, (const float *)varg2, (float *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Sub", acc);
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

int vSubAccuracyLiteTest_double() {
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
        vmdSub(VLEN, (const double *)varg1, (const double *)varg2,
               (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Sub", acc);
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

int vSubAccuracyLiteTest_float_complex() {
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
        vmcSub(VLEN, (const MKL_Complex8 *)varg1, (const MKL_Complex8 *)varg2,
               (MKL_Complex8 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Sub", acc);
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

int vSubAccuracyLiteTest_double_complex() {
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
        vmzSub(VLEN, (const MKL_Complex16 *)varg1, (const MKL_Complex16 *)varg2,
               (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Sub", acc);
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

  printf("Running %s functions:\n", "Sub");

  printf("\tRunning %s with single precision real data type:\n", "Sub");
  errs = vSubAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Sub",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Sub");
  errs = vSubAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Sub",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Sub");
  errs = vSubAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Sub",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Sub");
  errs = vSubAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Sub",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Sub",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
