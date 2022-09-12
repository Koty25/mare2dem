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
!    vcDiv  Example Program Text
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

  cA[0].real=-10000.0000;cA[0].imag=10000.0000;
  cA[1].real=-7777.7778;cA[1].imag=7777.7778;
  cA[2].real=-5555.5557;cA[2].imag=5555.5557;
  cA[3].real=-3333.3333;cA[3].imag=3333.3333;
  cA[4].real=-1111.1111;cA[4].imag=1111.1111;
  cA[5].real=1111.1111;cA[5].imag=-1111.1111;
  cA[6].real=3333.3333;cA[6].imag=-3333.3333;
  cA[7].real=5555.5557;cA[7].imag=-5555.5557;
  cA[8].real=7777.7778;cA[8].imag=-7777.7778;
  cA[9].real=10000.0000;cA[9].imag=-10000.0000;
  cB[0].real=1.0000000000000000e+000;cB[0].imag=0.0000000000000000e+000;
  cB[1].real=1.0000000000000000e+000;cB[1].imag=0.0000000000000000e+000;
  cB[2].real=1.0000000000000000e+000;cB[2].imag=0.0000000000000000e+000;
  cB[3].real=1.0000000000000000e+000;cB[3].imag=0.0000000000000000e+000;
  cB[4].real=1.0000000000000000e+000;cB[4].imag=0.0000000000000000e+000;
  cB[5].real=1.0000000000000000e+000;cB[5].imag=0.0000000000000000e+000;
  cB[6].real=1.0000000000000000e+000;cB[6].imag=0.0000000000000000e+000;
  cB[7].real=1.0000000000000000e+000;cB[7].imag=0.0000000000000000e+000;
  cB[8].real=1.0000000000000000e+000;cB[8].imag=0.0000000000000000e+000;
  cB[9].real=1.0000000000000000e+000;cB[9].imag=0.0000000000000000e+000;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcDiv(vec_len,cA,cA,cBha0);
  vcDivI(vec_len,cA_I,inca,cA_I,inca,cBha0_I,incb);

  vmcDiv(vec_len,cA,cA,cBep1,VML_EP);
  vmcDivI(vec_len,cA_I,inca,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcDiv(vec_len,cA,cA,cBep2);
  vcDivI(vec_len,cA_I,inca,cA_I,inca,cBep2_I,incb);

  vmcDiv(vec_len,cA,cA,cBla1,VML_LA);
  vmcDivI(vec_len,cA_I,inca,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcDiv(vec_len,cA,cA,cBla2);
  vcDivI(vec_len,cA_I,inca,cA_I,inca,cBla2_I,incb);

  vmcDiv(vec_len,cA,cA,cBha1,VML_HA);
  vmcDivI(vec_len,cA_I,inca,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcDiv(vec_len,cA,cA,cBha2);
  vcDivI(vec_len,cA_I,inca,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcDiv and vmcDiv in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcDivI and vmcDivI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcDiv and vmcDiv in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcDivI and vmcDivI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcDiv and vmcDiv in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcDivI and vmcDivI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcDiv test/example program\n\n");
  printf("           Arguments                           vcDiv\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .2f %+.2f*i\n",cA[i].real,cA[i].imag,cA[i].real,cA[i].imag,
        cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcDivI test/example program\n\n");
  printf("           Arguments                           vcDiv\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .2f %+.2f*i\n",cA_I[i*inca].real,cA_I[i*inca].imag,
        cA_I[i*inca].real,cA_I[i*inca].imag,cBha0_I[i*incb].real,cBha0_I[i*incb].imag);
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
