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
!    vsCdfNorm  Example Program Text
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

  fA[0]=-8.0000;
  fA[1]=-6.5555;
  fA[2]=-5.1111;
  fA[3]=-3.6666;
  fA[4]=-2.2222;
  fA[5]=-0.7777;
  fA[6]=0.6666;
  fA[7]=2.1111;
  fA[8]=3.5555;
  fA[9]=5.0000;
  fB[0]=6.22096057427178e-16;
  fB[1]=2.77278689399088e-11;
  fB[2]=1.60144022698396e-07;
  fB[3]=1.22898417236556e-04;
  fB[4]=1.31348990439007e-02;
  fB[5]=2.18372944048033e-01;
  fB[6]=7.47486161912039e-01;
  fB[7]=9.82618139490682e-01;
  fB[8]=9.99811369549446e-01;
  fB[9]=9.99999713348428e-01;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsCdfNorm(vec_len,fA,fBha0);
  vsCdfNormI(vec_len,fA_I,inca,fBha0_I,incb);

  vmsCdfNorm(vec_len,fA,fBep1,VML_EP);
  vmsCdfNormI(vec_len,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsCdfNorm(vec_len,fA,fBep2);
  vsCdfNormI(vec_len,fA_I,inca,fBep2_I,incb);

  vmsCdfNorm(vec_len,fA,fBla1,VML_LA);
  vmsCdfNormI(vec_len,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsCdfNorm(vec_len,fA,fBla2);
  vsCdfNormI(vec_len,fA_I,inca,fBla2_I,incb);

  vmsCdfNorm(vec_len,fA,fBha1,VML_HA);
  vmsCdfNormI(vec_len,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsCdfNorm(vec_len,fA,fBha2);
  vsCdfNormI(vec_len,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsCdfNorm and vmsCdfNorm in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsCdfNormI and vmsCdfNormI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsCdfNorm and vmsCdfNorm in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsCdfNormI and vmsCdfNormI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsCdfNorm and vmsCdfNorm in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsCdfNormI and vmsCdfNormI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsCdfNorm test/example program\n\n");
  printf("           Argument                     vsCdfNorm\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsCdfNormI test/example program\n\n");
  printf("           Argument                     vsCdfNorm\n");
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