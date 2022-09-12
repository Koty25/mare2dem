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
 *            Add example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0xC007309A }, { 0x4096200F } }, //  0: vsAdd ( 6.80375481     , -2.1123414      ) = ( 4.6914134       );
{ { 0x40B52EFA }, { 0x40BF006A }, { 0x413A17B2 } }, //  1: vsAdd ( 5.66198444     , 5.96880054      ) = ( 11.630785       );
{ { 0x4103BA28 }, { 0xC0C1912F }, { 0x400BC642 } }, //  2: vsAdd ( 8.2329483      , -6.04897261     ) = ( 2.1839757       );
{ { 0xC052EA36 }, { 0x40ABAABC }, { 0x40046B42 } }, //  3: vsAdd ( -3.2955451     , 5.3645916       ) = ( 2.0690465       );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xC000E6134801CC26 }, { 0x4012C401BCE58805 } }, //  0: vdAdd ( 6.80375434309419092      , -2.11234146361813924      ) = ( 4.69141287947605168       );
{ { 0x4016A5DF421D4BBE }, { 0x4017E00D485FC01A }, { 0x402742F6453E85EC } }, //  1: vdAdd ( 5.66198447517211711      , 5.96880066952146571       ) = ( 11.6307851446935828       );
{ { 0x40207744D998EE8A }, { 0xC0183225E080644C }, { 0x400178C7A562F190 } }, //  2: vdAdd ( 8.23294715873568705      , -6.04897261413232101      ) = ( 2.18397454460336604       );
{ { 0xC00A5D46A314BA8E }, { 0x4015755793FAEAB0 }, { 0x40008D6884E11AD2 } }, //  3: vdAdd ( -3.2955448857022196      , 5.36459189623808186       ) = ( 2.06904701053586226       );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40BF006A, 0x40B52EFA }, { 0x4076D03A, 0x414773AB } }, //  0: vcAdd ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( 3.85645914      + i * 12.4657393      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x40ABAABC, 0xC052EA36 }, { 0xBF2F3398, 0x409DFF35 } }, //  1: vcAdd ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -0.684381008    + i * 4.9374032       );
{ { 0x3F8A29C0, 0xC08E3964 }, { 0x4024F46C, 0xBEE77440 }, { 0x406A094C, 0xC09CB0A8 } }, //  2: vcAdd ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( 3.65681744      + i * -4.89656448     );
{ { 0x3E8939C0, 0xC02D136C }, { 0x41052EB4, 0x4110B6A8 }, { 0x41097882, 0x40CAE39A } }, //  3: vcAdd ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( 8.59192085      + i * 6.34028339      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x400EDA0748BDB40E, 0x4028EE755181DCEB } }, //  0: vzAdd ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( 3.85645920590332647       + i * 12.465738818266308        );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xBFE5E672642BCCE0, 0x4013BFE661A77FCD } }, //  1: vzAdd ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -0.684380717894239154     + i * 4.93740227303346746       );
{ { 0x3FF1453801E28A70, 0xC011C72C86338E59 }, { 0x40049E8D96893D1C, 0xBFDCEE88B739DD20 }, { 0x400D4129977A8254, 0xC013961511A72C2B } }, //  2: vzAdd ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( 3.65681761114709936       + i * -4.89656474669304043      );
{ { 0x3FD12735D3224E60, 0xC005A26D910B44DC }, { 0x4020A5D666294BAC, 0x402216D5173C2DAA }, { 0x40212F1014C25E1F, 0x40195C7365F2B8E6 } }, //  3: vzAdd ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( 8.59191956398632151       + i * 6.34028395933112243       );
}

};

//!
//! @brief Single precision test
//!

int vAddAccuracyLiteTest_float() {
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
        vmsAdd(VLEN, (const float *)varg1, (const float *)varg2, (float *)vres1,
               accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Add", acc);
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

int vAddAccuracyLiteTest_double() {
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
        vmdAdd(VLEN, (const double *)varg1, (const double *)varg2,
               (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Add", acc);
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

int vAddAccuracyLiteTest_float_complex() {
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
        vmcAdd(VLEN, (const MKL_Complex8 *)varg1, (const MKL_Complex8 *)varg2,
               (MKL_Complex8 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_float_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Add", acc);
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

int vAddAccuracyLiteTest_double_complex() {
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
        vmzAdd(VLEN, (const MKL_Complex16 *)varg1, (const MKL_Complex16 *)varg2,
               (MKL_Complex16 *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs +=
          check_result_double_complex(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                      vres1[i], vref1[i], vref1[i], "Add", acc);
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

  printf("Running %s functions:\n", "Add");

  printf("\tRunning %s with single precision real data type:\n", "Add");
  errs = vAddAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Add",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Add");
  errs = vAddAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Add",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Add");
  errs = vAddAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Add",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Add");
  errs = vAddAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Add",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Add",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
