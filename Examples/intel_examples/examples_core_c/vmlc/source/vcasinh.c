/*******************************************************************************
* Copyright 2001-2020 Intel Corporation.
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
!  Content:
!    vcAsinh  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  MKL_Complex8 cA[10],cB[10];
  MKL_Complex8 cBha0[10],cBha1[10],cBha2[10];
  MKL_Complex8           cBla1[10],cBla2[10];
  MKL_Complex8           cBep1[10],cBep2[10];
  float CurRMS,MaxRMS=0.0;

  MKL_Complex8 cA_I[10*INCA],cB_I[10*INCB];
  MKL_Complex8 cBha0_I[10*INCB],cBha1_I[10*INCB],cBha2_I[10*INCB];
  MKL_Complex8                  cBla1_I[10*INCB],cBla2_I[10*INCB];
  MKL_Complex8                  cBep1_I[10*INCB],cBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  cA[0].real=3.0000;cA[0].imag=10000.0000;
  cA[1].real=1113.7777;cA[1].imag=8889.2227;
  cA[2].real=2224.5554;cA[2].imag=7778.4443;
  cA[3].real=3335.3333;cA[3].imag=6667.6665;
  cA[4].real=4446.1113;cA[4].imag=5556.8887;
  cA[5].real=5556.8887;cA[5].imag=4446.1113;
  cA[6].real=6667.6665;cA[6].imag=3335.3333;
  cA[7].real=7778.4443;cA[7].imag=2224.5554;
  cA[8].real=8889.2227;cA[8].imag=1113.7777;
  cA[9].real=10000.0000;cA[9].imag=3.0000;
  cB[0].real=9.9034875950361272e+000;cB[0].imag=1.5704963207244873e+000;
  cB[1].real=9.7935305399540784e+000;cB[1].imag=1.4461505413055420e+000;
  cB[2].real=9.6915674216669974e+000;cB[2].imag=1.2922420501708984e+000;
  cB[3].real=9.6098341803762892e+000;cB[3].imag=1.1069687604904175e+000;
  cB[4].real=9.5633416449802731e+000;cB[4].imag=8.9598953723907471e-001;
  cB[5].real=9.5633416471460873e+000;cB[5].imag=6.7480677366256714e-001;
  cB[6].real=9.6098341857711080e+000;cB[6].imag=4.6382755041122437e-001;
  cB[7].real=9.6915674281509592e+000;cB[7].imag=2.7855432033538818e-001;
  cB[8].real=9.7935305459913469e+000;cB[8].imag=1.2464573234319687e-001;
  cB[9].real=9.9034876000361258e+000;cB[9].imag=2.9999998514540493e-004;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcAsinh(vec_len,cA,cBha0);
  vcAsinhI(vec_len,cA_I,inca,cBha0_I,incb);

  vmcAsinh(vec_len,cA,cBep1,VML_EP);
  vmcAsinhI(vec_len,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcAsinh(vec_len,cA,cBep2);
  vcAsinhI(vec_len,cA_I,inca,cBep2_I,incb);

  vmcAsinh(vec_len,cA,cBla1,VML_LA);
  vmcAsinhI(vec_len,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcAsinh(vec_len,cA,cBla2);
  vcAsinhI(vec_len,cA_I,inca,cBla2_I,incb);

  vmcAsinh(vec_len,cA,cBha1,VML_HA);
  vmcAsinhI(vec_len,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcAsinh(vec_len,cA,cBha2);
  vcAsinhI(vec_len,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcAsinh and vmcAsinh in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcAsinhI and vmcAsinhI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcAsinh and vmcAsinh in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcAsinhI and vmcAsinhI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcAsinh and vmcAsinh in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcAsinhI and vmcAsinhI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcAsinh test/example program\n\n");
  printf("           Argument                           vcAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA[i].real,cA[i].imag,cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcAsinhI test/example program\n\n");
  printf("           Argument                           vcAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA_I[i*inca].real,cA_I[i*inca].imag,cBha0_I[i*incb].real,
        cBha0_I[i*incb].imag);
    CurRMS_I=crelerr(cB_I[i*incb],cBha0_I[i*incb]);
    if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I;
  }
  printf("\n");
  if(MaxRMS>=1e-5) {
    printf("Error! Relative accuracy is %.16f\n",MaxRMS);
    return 1;
  }
  else {
    printf("Relative accuracy is %.16f\n",MaxRMS);
  }

  printf("\n");

  if(MaxRMS_I>=1e-5) {
    printf("Error! Relative strided accuracy is %.16f\n",MaxRMS_I);
    return 1;
  }
  else {
    printf("Relative strided accuracy is %.16f\n",MaxRMS_I);
  }
  return 0;
}
