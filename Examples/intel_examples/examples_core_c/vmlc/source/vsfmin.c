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
!     vsFmin Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5
#define INCC 7

int main()
{
  float fA1[10],fA2[10],fB[10];
  float fBha0[10],fBha1[10],fBha2[10];
  float           fBla1[10],fBla2[10];
  float           fBep1[10],fBep2[10];
  float CurRMS,MaxRMS=0.0;

  float fA1_I[10*INCA],fA2_I[10*INCB],fB_I[10*INCC];
  float fBha0_I[10*INCC],fBha1_I[10*INCC],fBha2_I[10*INCC];
  float                  fBla1_I[10*INCC],fBla2_I[10*INCC];
  float                  fBep1_I[10*INCC],fBep2_I[10*INCC];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB,incc=INCC;

  MKL_INT i=0,vec_len=10;

  fA1[0]=0.1000;
  fA1[1]=0.8666;
  fA1[2]=1.6333;
  fA1[3]=2.4000;
  fA1[4]=3.1666;
  fA1[5]=3.9333;
  fA1[6]=4.7000;
  fA1[7]=5.4666;
  fA1[8]=6.2333;
  fA1[9]=7.0000;
  fA2[0]=-10.0000;
  fA2[1]=-7.7777;
  fA2[2]=-5.5555;
  fA2[3]=-3.3333;
  fA2[4]=-1.1111;
  fA2[5]=1.1111;
  fA2[6]=3.3333;
  fA2[7]=5.5555;
  fA2[8]=7.7777;
  fA2[9]=10.0000;
  fB[0]=-1.0000000000000000e+01;
  fB[1]=-7.7777000000000003e+00;
  fB[2]=-5.5555000000000003e+00;
  fB[3]=-3.3332999999999999e+00;
  fB[4]=-1.1111000000000000e+00;
  fB[5]=1.1111000000000000e+00;
  fB[6]=3.3332999999999999e+00;
  fB[7]=5.4665999999999997e+00;
  fB[8]=6.2332999999999998e+00;
  fB[9]=7.0000000000000000e+00;

  for(i=0;i<10;i++) {
    fA1_I[i*inca]=fA1[i];
    fA2_I[i*incb]=fA2[i];
    fB_I[i*incc]=fB[i];
  }

  vsFmin(vec_len,fA1,fA2,fBha0);
  vsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBha0_I,incc);

  vmsFmin(vec_len,fA1,fA2,fBep1,VML_EP);
  vmsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBep1_I,incc,VML_EP);

  vmlSetMode(VML_EP);
  vsFmin(vec_len,fA1,fA2,fBep2);
  vsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBep2_I,incc);

  vmsFmin(vec_len,fA1,fA2,fBla1,VML_LA);
  vmsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBla1_I,incc,VML_LA);

  vmlSetMode(VML_LA);
  vsFmin(vec_len,fA1,fA2,fBla2);
  vsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBla2_I,incc);

  vmsFmin(vec_len,fA1,fA2,fBha1,VML_HA);
  vmsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBha1_I,incc,VML_HA);

  vmlSetMode(VML_HA);
  vsFmin(vec_len,fA1,fA2,fBha2);
  vsFminI(vec_len,fA1_I,inca,fA2_I,incb,fBha2_I,incc);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsFmin and vmsFmin in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incc]!=fBha1_I[i*incc] || fBha1_I[i*incc]!=fBha2_I[i*incc]) {
      printf("Error! Difference between vsFminI and vmsFminI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsFmin and vmsFmin in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incc]!=fBla2_I[i*incc]) {
      printf("Error! Difference between vsFminI and vmsFminI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsFmin and vmsFmin in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incc]!=fBep2_I[i*incc]) {
      printf("Error! Difference between vsFminI and vmsFminI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsFmin test/example program\n\n");
  printf("           Argument                     vsFmin\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f %25.14f % 25.14e\n",fA1[i],fA2[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsFminI test/example program\n\n");
  printf("           Argument                     vsFmin\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f %25.14f % 25.14e\n",fA1_I[i*inca],fA2_I[i*incb],fBha0_I[i*incc]);
    CurRMS_I=srelerr(fB_I[i*incc],fBha0_I[i*incc]);
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
