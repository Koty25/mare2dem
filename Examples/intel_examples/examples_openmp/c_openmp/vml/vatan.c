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
 *            Atan example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0x3FB661ED } }, //  0: vsAtan ( 6.80375481      ) = ( 1.42486346      );
{ { 0xC007309A }, { 0xBF907785 } }, //  1: vsAtan ( -2.1123414      ) = ( -1.12864745     );
{ { 0x40B52EFA }, { 0x3FB2AF8F } }, //  2: vsAtan ( 5.66198444      ) = ( 1.39598262      );
{ { 0x40BF006A }, { 0x3FB3D07E } }, //  3: vsAtan ( 5.96880054      ) = ( 1.40480018      );
}

,

{

{ { 0x401B370B60E66E18 }, { 0x3FF6CC3DAB5A6D16 } }, //  0: vdAtan ( 6.80375434309419092       ) = ( 1.42486349997381501       );
{ { 0xC000E6134801CC26 }, { 0xBFF20EF0A5D69141 } }, //  1: vdAtan ( -2.11234146361813924      ) = ( -1.12864746838120333      );
{ { 0x4016A5DF421D4BBE }, { 0x3FF655F1DB631CA2 } }, //  2: vdAtan ( 5.66198447517211711       ) = ( 1.39598260591609646       );
{ { 0x4017E00D485FC01A }, { 0x3FF67A0FB53FB5DE } }, //  3: vdAtan ( 5.96880066952146571       ) = ( 1.40480013656939873       );
}

,

{

{ { 0xC007309A, 0x40D9B85C }, { 0xBFC3A3F8, 0x3E09DBCA } }, //  0: vcAtan ( -2.1123414      + i * 6.80375481      ) = ( -1.52844143     + i * 0.134627491     );
{ { 0x40BF006A, 0x40B52EFA }, { 0x3FBDB99A, 0x3DAA6186 } }, //  1: vcAtan ( 5.96880054      + i * 5.66198444      ) = ( 1.48222661      + i * 0.0831938237    );
{ { 0xC0C1912F, 0x4103BA28 }, { 0xBFC19B0F, 0x3DA15661 } }, //  2: vcAtan ( -6.04897261     + i * 8.2329483       ) = ( -1.51254451     + i * 0.0787780359    );
{ { 0x40ABAABC, 0xC052EA36 }, { 0x3FB7BA3C, 0xBDA78E76 } }, //  3: vcAtan ( 5.3645916       + i * -3.2955451      ) = ( 1.43537092      + i * -0.0818146914   );
}

,

{

{ { 0xC000E6134801CC26, 0x401B370B60E66E18 }, { 0xBFF8747EFBF5F2BF, 0x3FC13B795AD816CE } }, //  0: vzAtan ( -2.11234146361813924      + i * 6.80375434309419092       ) = ( -1.52844141409074985      + i * 0.13462750373599025       );
{ { 0x4017E00D485FC01A, 0x4016A5DF421D4BBE }, { 0x3FF7B7334500B661, 0x3FB54C30C3922F24 } }, //  1: vzAtan ( 5.96880066952146571       + i * 5.66198447517211711       ) = ( 1.48222662882053435       + i * 0.0831938245266284349     );
{ { 0xC0183225E080644C, 0x40207744D998EE8A }, { 0xBFF83361D7717724, 0x3FB42ACC31C9F775 } }, //  2: vzAtan ( -6.04897261413232101      + i * 8.23294715873568705       ) = ( -1.51254448087224436      + i * 0.0787780400805482978     );
{ { 0x4015755793FAEAB0, 0xC00A5D46A314BA8E }, { 0x3FF6F747899AD360, 0xBFB4F1CEA21834C6 } }, //  3: vzAtan ( 5.36459189623808186       + i * -3.2955448857022196       ) = ( 1.43537095786924596       + i * -0.0818146844614658642    );
}

};

//!
//! @brief Single precision test
//!

int vAtanAccuracyLiteTest_float() {
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
        vmsAtan(VLEN, (const float *)varg1, (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Atan", acc);
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

int vAtanAccuracyLiteTest_double() {
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
        vmdAtan(VLEN, (const double *)varg1, (double *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Atan", acc);
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

int vAtanAccuracyLiteTest_float_complex() {
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
        vmcAtan(VLEN, (const MKL_Complex8 *)varg1, (MKL_Complex8 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs +=
          check_result_float_complex(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                     vres1[i], vref1[i], vref1[i], "Atan", acc);
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

int vAtanAccuracyLiteTest_double_complex() {
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
        vmzAtan(VLEN, (const MKL_Complex16 *)varg1, (MKL_Complex16 *)vres1,
                accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double_complex(ARG1_RES1, varg1[i], varg1[i],
                                          vres1[i], vres1[i], vref1[i],
                                          vref1[i], "Atan", acc);
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

  printf("Running %s functions:\n", "Atan");

  printf("\tRunning %s with single precision real data type:\n", "Atan");
  errs = vAtanAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Atan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Atan");
  errs = vAtanAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Atan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with single precision complex data type:\n", "Atan");
  errs = vAtanAccuracyLiteTest_float_complex();
  printf("\t%s single precision complex result: %s\n\n", "Atan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision complex data type:\n", "Atan");
  errs = vAtanAccuracyLiteTest_double_complex();
  printf("\t%s double precision complex result: %s\n", "Atan",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("%s function result: %s\n\n", "Atan",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
