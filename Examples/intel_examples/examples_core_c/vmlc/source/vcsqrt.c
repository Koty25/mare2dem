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
!    vcSqrt  Example Program Text
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

  cA[0].real=0.0000;cA[0].imag=10000.0000;
  cA[1].real=1111.1111;cA[1].imag=8888.8887;
  cA[2].real=2222.2222;cA[2].imag=7777.7778;
  cA[3].real=3333.3333;cA[3].imag=6666.6665;
  cA[4].real=4444.4443;cA[4].imag=5555.5557;
  cA[5].real=5555.5557;cA[5].imag=4444.4443;
  cA[6].real=6666.6665;cA[6].imag=3333.3333;
  cA[7].real=7777.7778;cA[7].imag=2222.2222;
  cA[8].real=8888.8887;cA[8].imag=1111.1111;
  cA[9].real=10000.0000;cA[9].imag=0.0000;
  cB[0].real=7.0710678118654755e+001;cB[0].imag=7.0710678100585937e+001;
  cB[1].real=7.0954827284892332e+001;cB[1].imag=6.2637660980224609e+001;
  cB[2].real=7.1802622491520651e+001;cB[2].imag=5.4160820007324219e+001;
  cB[3].real=7.3440087809658650e+001;cB[3].imag=4.5388469696044922e+001;
  cB[4].real=7.6023111087747310e+001;cB[4].imag=3.6538597106933594e+001;
  cB[5].real=7.9593147214584178e+001;cB[5].imag=2.7919767379760742e+001;
  cB[6].real=8.4024479310894492e+001;cB[6].imag=1.9835489273071289e+001;
  cB[7].real=8.9069604404033797e+001;cB[7].imag=1.2474637985229492e+001;
  cB[8].real=9.4464153566176705e+001;cB[8].imag=5.8811254501342773e+000;
  cB[9].real=1.0000000000000000e+002;cB[9].imag=0.0000000000000000e+000;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcSqrt(vec_len,cA,cBha0);
  vcSqrtI(vec_len,cA_I,inca,cBha0_I,incb);

  vmcSqrt(vec_len,cA,cBep1,VML_EP);
  vmcSqrtI(vec_len,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcSqrt(vec_len,cA,cBep2);
  vcSqrtI(vec_len,cA_I,inca,cBep2_I,incb);

  vmcSqrt(vec_len,cA,cBla1,VML_LA);
  vmcSqrtI(vec_len,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcSqrt(vec_len,cA,cBla2);
  vcSqrtI(vec_len,cA_I,inca,cBla2_I,incb);

  vmcSqrt(vec_len,cA,cBha1,VML_HA);
  vmcSqrtI(vec_len,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcSqrt(vec_len,cA,cBha2);
  vcSqrtI(vec_len,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcSqrt and vmcSqrt in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcSqrtI and vmcSqrtI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcSqrt and vmcSqrt in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcSqrtI and vmcSqrtI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcSqrt and vmcSqrt in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcSqrtI and vmcSqrtI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcSqrt test/example program\n\n");
  printf("           Argument                           vcSqrt\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA[i].real,cA[i].imag,cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcSqrtI test/example program\n\n");
  printf("           Argument                           vcSqrt\n");
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
