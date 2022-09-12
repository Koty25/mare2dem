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
 *            Powr example program text (OpenMP offload interface)
 *
 *******************************************************************************/

#define VLEN 4

#include "_vml_common.h"

max_ulp_table_t max_ulp_table =
{ //  HA   LA   EP
    { FLT_MAX, FLT_MAX, FLT_MAX, }, // float
    { DBL_MAX, DBL_MAX, DBL_MAX, }, // double
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

{ { 0x41093E24 }, { 0x4093852F }, { 0x469CE711 } }, //  0: vsPowr ( 8.57767105     , 4.61000776      ) = ( 20083.5332      );
{ { 0x41011D03 }, { 0x41034C40 }, { 0x4BD2F79E } }, //  1: vsPowr ( 8.06958294     , 8.20611572      ) = ( 27651900        );
{ { 0x41136B29 }, { 0x4036ECDE }, { 0x440EB888 } }, //  2: vsPowr ( 9.21366215     , 2.85820723      ) = ( 570.883301      );
{ { 0x4082ABE3 }, { 0x40FDFDE5 }, { 0x478A3D0D } }, //  3: vsPowr ( 4.08348227     , 7.93724298      ) = ( 70778.1016      );
}

,

{

{ { 0x402127C473A3E923 }, { 0x401270A5F32DAE19 }, { 0x40D39CE2AABD156D } }, //  0: vdPowr ( 8.57767068267691535      , 4.6100080486899282        ) = ( 20083.5416710576283       );
{ { 0x402023A0651C4741 }, { 0x40206988134D9FDD }, { 0x417A5EF61F31C368 } }, //  1: vdPowr ( 8.06958309145159269      , 8.20611629793705255       ) = ( 27651937.9496492445       );
{ { 0x40226D6509CA7464 }, { 0x4006DD9BBAC0EE6B }, { 0x4081D7109BAA7980 } }, //  2: vdPowr ( 9.21366148563738108      , 2.8582071867111174        ) = ( 570.883109409172903       );
{ { 0x4010557C717977C6 }, { 0x401FBFBCBB737F7A }, { 0x40F147A2E685447B } }, //  3: vdPowr ( 4.0834825258625127       , 7.93724339382594657       ) = ( 70778.1812794375437       );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vPowrAccuracyLiteTest_float() {
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
        vmsPowr(VLEN, (const float *)varg1, (const float *)varg2,
                (float *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_float(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Powr", acc);
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

int vPowrAccuracyLiteTest_double() {
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
        vmdPowr(VLEN, (const double *)varg1, (const double *)varg2,
                (double *)vres1, accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {

      errs += check_result_double(ARG2_RES1, varg1[i], varg2[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Powr", acc);
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

  printf("Running %s functions:\n", "Powr");

  printf("\tRunning %s with single precision real data type:\n", "Powr");
  errs = vPowrAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Powr",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Powr");
  errs = vPowrAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Powr",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Powr",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}
