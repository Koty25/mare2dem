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
 *            MulByConj example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { 4.5, 5.0, 5.0E3, }, // float
    { 2.0, 5.0, 7.0E7, }, // double
    { 16.0, 16.0, 16.0, }, // VM_COMPLEX8
    { 16.0, 16.0, 16.0, }, // VM_COMPLEX16
};

// device number
int dnum;

// *************************************************************
// Data table declaraion
// *************************************************************
data_3_t data =
{

{ /* empty */ }

,

{ /* empty */}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0x40BF006A, 0x40B52EFA }, { 0x41CF511F, 0x425247FD } }, //  0: vcMulByConj ( -2.1123414      + i * 6.80375481     , 5.96880054      + i * 5.66198444      ) = ( 25.9146099      + i * 52.5703011      );
{ { 0xC0C1912F, 0x4103BA28 }, { 0x40ABAABC, 0xC052EA36 }, { 0xC26E544C, 0x41C1DA9C } }, //  1: vcMulByConj ( -6.04897261     + i * 8.2329483      , 5.3645916       + i * -3.2955451      ) = ( -59.5823212     + i * 24.2317429      );
{ { 0x3F8A29C0, 0xC08E3964 }, { 0x4024F46C, 0xBEE77440 }, { 0x409951D8, 0xC12F7A77 } }, //  2: vcMulByConj ( 1.07939911      + i * -4.44450569    , 2.57741833      + i * -0.452058792    ) = ( 4.79124069      + i * -10.9673986     );
{ { 0x3E8939C0, 0xC02D136C }, { 0x41052EB4, 0x4110B6A8 }, { 0xC1B1D3D2, 0xC1C779EE } }, //  3: vcMulByConj ( 0.268018723     + i * -2.70431042    , 8.32390213      + i * 9.04459381      ) = ( -22.2284279     + i * -24.934536      );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x4039EA23A3CE3157, 0x404A48FF86CC45BE } }, //  0: vzMulByConj ( -2.11234146361813924      + i * 6.80375434309419092      , 5.96880066952146571       + i * 5.66198447517211711       ) = ( 25.9146063211822728       + i * 52.5702980515884377       );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0xC04DCA8957C19C5A, 0x40383B535E41D16B } }, //  1: vzMulByConj ( -6.04897261413232101      + i * 8.23294715873568705      , 5.36459189623808186       + i * -3.2955448857022196       ) = ( -59.5823163695683462      + i * 24.2317408476532528       );
{ { 0x3FF1453801E28A70, 0xC011C72C86338E59 }, { 0x40049E8D96893D1C, 0xBFDCEE88B739DD20 }, { 0x40132A3B52620689, 0xC025EF4EF83576E2 } }, //  2: vzMulByConj ( 1.07939911590861115       + i * -4.44450578393624429     , 2.57741849523848821       + i * -0.452058962756796134     ) = ( 4.79124191973972646       + i * -10.9673993649734633      );
{ { 0x3FD12735D3224E60, 0xC005A26D910B44DC }, { 0x4020A5D666294BAC, 0x402216D5173C2DAA }, { 0xC0363A7ABA8C5B9C, 0xC038EF3D5EDD987E } }, //  3: vzMulByConj ( 0.268018203912310682      + i * -2.70431054416313366     , 8.32390136007401082       + i * 9.04459450349425609       ) = ( -22.2284351914091616      + i * -24.9345301905636845      );
}

};

//!
//! @brief Complex single precision test
//!

int vMulByConjAccuracyLiteTest_float_complex() {
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
        vmcMulByConj(VLEN, (const MKL_Complex8 *)varg1,
                     (const MKL_Complex8 *)varg2, (MKL_Complex8 *)vres1,
                     accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float_complex(ARG2_RES1, varg1[i], varg2[i],
                                         vres1[i], vres1[i], vref1[i], vref1[i],
                                         "MulByConj", acc);
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

int vMulByConjAccuracyLiteTest_double_complex() {
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
        vmzMulByConj(VLEN, (const MKL_Complex16 *)varg1,
                     (const MKL_Complex16 *)varg2, (MKL_Complex16 *)vres1,
                     accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double_complex(ARG2_RES1, varg1[i], varg2[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "MulByConj", acc);
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

  printf("Running %s functions:\n", "MulByConj");
  printf("\tRunning %s with single precision complex data type:\n",
         "MulByConj");
  errs = vMulByConjAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "MulByConj",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n",
         "MulByConj");
  errs = vMulByConjAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "MulByConj",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "MulByConj",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
