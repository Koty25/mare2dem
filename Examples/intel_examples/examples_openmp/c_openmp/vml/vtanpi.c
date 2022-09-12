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
 *            Tanpi example program text (OpenMP offload interface)
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

{ { 0x40D9B85C }, { 0xBF356B88 } }, //  0: vsTanpi ( 6.80375481      ) = ( -0.708672047    );
{ { 0xC007309A }, { 0xBEBC98D5 } }, //  1: vsTanpi ( -2.1123414      ) = ( -0.368353516    );
{ { 0x40B52EFA }, { 0xBFE56E6A } }, //  2: vsTanpi ( 5.66198444      ) = ( -1.79243207     );
{ { 0x40BF006A }, { 0xBDC961CD } }, //  3: vsTanpi ( 5.96880054      ) = ( -0.0983310714   );
}

,

{

{ { 0x401B370B60E66E18 }, { 0xBFE6AD759D62B54F } }, //  0: vdTanpi ( 6.80375434309419092       ) = ( -0.708674247163495141     );
{ { 0xC000E6134801CC26 }, { 0xBFD7931B8E7CDD8B } }, //  1: vdTanpi ( -2.11234146361813924      ) = ( -0.368353737972449025     );
{ { 0x4016A5DF421D4BBE }, { 0xBFFCADCCD759EC08 } }, //  2: vdTanpi ( 5.66198447517211711       ) = ( -1.79243167993149832      );
{ { 0x4017E00D485FC01A }, { 0xBFB92C330946EDEC } }, //  3: vdTanpi ( 5.96880066952146571       ) = ( -0.0983306787156192352    );
}

,
{ /* empty */ }

,

{ /* empty */ }

};

//!
//! @brief Single precision test
//!

int vTanpiAccuracyLiteTest_float() {
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
        vmsTanpi(VLEN, (const float *)varg1, (float *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_float(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                 vres1[i], vref1[i], vref1[i], "Tanpi", acc);
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

int vTanpiAccuracyLiteTest_double() {
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
        vmdTanpi(VLEN, (const double *)varg1, (double *)vres1,
                 accuracy_mode[acc]);
      }
#pragma omp taskwait
    }

    for (int i = 0; i < VLEN; ++i) {
      errs += check_result_double(ARG1_RES1, varg1[i], varg1[i], vres1[i],
                                  vres1[i], vref1[i], vref1[i], "Tanpi", acc);
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

  printf("Running %s functions:\n", "Tanpi");

  printf("\tRunning %s with single precision real data type:\n", "Tanpi");
  errs = vTanpiAccuracyLiteTest_float();
  printf("\t%s single precision real result: %s\n\n", "Tanpi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;

  printf("\tRunning %s with double precision real data type:\n", "Tanpi");
  errs = vTanpiAccuracyLiteTest_double();
  printf("\t%s double precision real result: %s\n\n", "Tanpi",
         (errs == 0) ? "PASS" : "FAIL");
  total_errs += errs;
  printf("%s function result: %s\n\n", "Tanpi",
         (total_errs == 0) ? "PASS" : "FAIL");

  return total_errs;
}