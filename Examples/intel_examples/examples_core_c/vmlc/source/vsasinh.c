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
!    vsAsinh  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  float fA[10],fB[10];
  float fBha0[10],fBha1[10],fBha2[10];
  float           fBla1[10],fBla2[10];
  float           fBep1[10],fBep2[10];
  float CurRMS,MaxRMS=0.0;

  float fA_I[10*INCA],fB_I[10*INCB];
  float fBha0_I[10*INCB],fBha1_I[10*INCB],fBha2_I[10*INCB];
  float                  fBla1_I[10*INCB],fBla2_I[10*INCB];
  float                  fBep1_I[10*INCB],fBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  fA[0]=3.0000;
  fA[1]=1113.7777;
  fA[2]=2224.5554;
  fA[3]=3335.3333;
  fA[4]=4446.1113;
  fA[5]=5556.8887;
  fA[6]=6667.6665;
  fA[7]=7778.4443;
  fA[8]=8889.2227;
  fA[9]=10000.0000;
  fB[0]=1.8184464592320668e+000;
  fB[1]=7.7086602404420574e+000;
  fB[2]=8.4004595944431433e+000;
  fB[3]=8.8054750820135883e+000;
  fB[4]=9.0929323274814546e+000;
  fB[5]=9.3159408278816027e+000;
  fB[6]=9.4981724143920001e+000;
  fB[7]=9.6522588250499251e+000;
  fB[8]=9.7857420681667548e+000;
  fB[9]=9.9034875550361274e+000;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsAsinh(vec_len,fA,fBha0);
  vsAsinhI(vec_len,fA_I,inca,fBha0_I,incb);

  vmsAsinh(vec_len,fA,fBep1,VML_EP);
  vmsAsinhI(vec_len,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsAsinh(vec_len,fA,fBep2);
  vsAsinhI(vec_len,fA_I,inca,fBep2_I,incb);

  vmsAsinh(vec_len,fA,fBla1,VML_LA);
  vmsAsinhI(vec_len,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsAsinh(vec_len,fA,fBla2);
  vsAsinhI(vec_len,fA_I,inca,fBla2_I,incb);

  vmsAsinh(vec_len,fA,fBha1,VML_HA);
  vmsAsinhI(vec_len,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsAsinh(vec_len,fA,fBha2);
  vsAsinhI(vec_len,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsAsinh and vmsAsinh in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsAsinhI and vmsAsinhI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsAsinh and vmsAsinh in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsAsinhI and vmsAsinhI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsAsinh and vmsAsinh in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsAsinhI and vmsAsinhI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsAsinh test/example program\n\n");
  printf("           Argument                     vsAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsAsinhI test/example program\n\n");
  printf("           Argument                     vsAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA_I[i*inca],fBha0_I[i*incb]);
    CurRMS_I=srelerr(fB_I[i*incb],fBha0_I[i*incb]);
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
