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
!    vsPow  Example Program Text
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

  fA[0]=0.1000;
  fA[1]=0.8666;
  fA[2]=1.6333;
  fA[3]=2.4000;
  fA[4]=3.1666;
  fA[5]=3.9333;
  fA[6]=4.7000;
  fA[7]=5.4666;
  fA[8]=6.2333;
  fA[9]=7.0000;
  fB[0]=7.9432823318248802e-001;
  fB[1]=8.8331105044955371e-001;
  fB[2]=2.2284382351152026e+000;
  fB[3]=8.1753632374188427e+000;
  fB[4]=3.8474983327516213e+001;
  fB[5]=2.1845299866537277e+002;
  fB[6]=1.4416496181341763e+003;
  fB[7]=1.0784631419538207e+004;
  fB[8]=8.9890987095794902e+004;
  fB[9]=8.2354300000000000e+005;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsPow(vec_len,fA,fA,fBha0);
  vsPowI(vec_len,fA_I,inca,fA_I,inca,fBha0_I,incb);

  vmsPow(vec_len,fA,fA,fBep1,VML_EP);
  vmsPowI(vec_len,fA_I,inca,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsPow(vec_len,fA,fA,fBep2);
  vsPowI(vec_len,fA_I,inca,fA_I,inca,fBep2_I,incb);

  vmsPow(vec_len,fA,fA,fBla1,VML_LA);
  vmsPowI(vec_len,fA_I,inca,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsPow(vec_len,fA,fA,fBla2);
  vsPowI(vec_len,fA_I,inca,fA_I,inca,fBla2_I,incb);

  vmsPow(vec_len,fA,fA,fBha1,VML_HA);
  vmsPowI(vec_len,fA_I,inca,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsPow(vec_len,fA,fA,fBha2);
  vsPowI(vec_len,fA_I,inca,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsPow and vmsPow in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsPowI and vmsPowI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsPow and vmsPow in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsPowI and vmsPowI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsPow and vmsPow in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsPowI and vmsPowI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsPow test/example program\n\n");
  printf("                     Arguments                               vsPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14f % 25.14e\n",fA[i],fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsPowI test/example program\n\n");
  printf("                     Arguments                               vsPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14f % 25.14e\n",fA_I[i*inca],fA_I[i*inca],fBha0_I[i*incb]);
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
