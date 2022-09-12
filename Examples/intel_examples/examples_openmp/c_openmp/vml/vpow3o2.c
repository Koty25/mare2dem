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
 *            Pow3o2 example program text (OpenMP offload interface)
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

{ { 0x4106AF8C }, { 0x41C362B2 } }, //  0: vsPow3o2 ( 8.41785812      ) = ( 24.4231911      );
{ { 0x408023F8 }, { 0x410035F8 } }, //  1: vsPow3o2 ( 4.00439072      ) = ( 8.01317596      );
{ { 0x40FB492C }, { 0x41B00AD4 } }, //  2: vsPow3o2 ( 7.85268211      ) = ( 22.0052872      );
{ { 0x410012AA }, { 0x41B52C8C } }, //  3: vsPow3o2 ( 8.00455666      ) = ( 22.6467514      );
}

,

{

{ { 0x4020D5F18943457D }, { 0x40386C564FFCEAF7 } }, //  0: vdPow3o2 ( 8.41785839983162454       ) = ( 24.4231920235133337       );
{ { 0x4010047F1160D5CB }, { 0x402006BF13594277 } }, //  1: vdPow3o2 ( 4.00439097550902101       ) = ( 8.01317654099078247       );
{ { 0x401F69258D86D24B }, { 0x4036015A817DC769 } }, //  2: vdPow3o2 ( 7.85268231521019811       ) = ( 22.005287259299994        );
{ { 0x40200255351CD177 }, { 0x4036A59172F8673B } }, //  3: vdPow3o2 ( 8.00455633141312539       ) = ( 22.6467506271794541       );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vPow3o2AccuracyLiteTest_float() {
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
        vmsPow3o2(VLEN, (const float *)varg1, (float *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Pow3o2", acc);
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

int vPow3o2AccuracyLiteTest_double() {
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
        vmdPow3o2(VLEN, (const double *)varg1, (double *)vres1,
                  accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Pow3o2", acc);
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

  printf("Running %s functions:\n", "Pow3o2");

  printf("\tRunning %s with single precision real data type:\n", "Pow3o2");
  errs = vPow3o2AccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Pow3o2",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Pow3o2");
  errs = vPow3o2AccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Pow3o2",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Pow3o2",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
