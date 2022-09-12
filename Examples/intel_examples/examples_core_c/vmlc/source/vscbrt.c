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
!    vsCbrt  Example Program Text
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

  fA[0]=0.0000;
  fA[1]=1111.1111;
  fA[2]=2222.2222;
  fA[3]=3333.3333;
  fA[4]=4444.4443;
  fA[5]=5555.5557;
  fA[6]=6666.6665;
  fA[7]=7777.7778;
  fA[8]=8888.8887;
  fA[9]=10000.0000;
  fB[0]=0.0000000000000000e+000;
  fB[1]=1.0357441602223785e+001;
  fB[2]=1.3049558697698629e+001;
  fB[3]=1.4938015700291331e+001;
  fB[4]=1.6441413695069230e+001;
  fB[5]=1.7710976268349352e+001;
  fB[6]=1.8820720424457154e+001;
  fB[7]=1.9813073221945444e+001;
  fB[8]=2.0714883204447570e+001;
  fB[9]=2.1544346900318835e+001;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsCbrt(vec_len,fA,fBha0);
  vsCbrtI(vec_len,fA_I,inca,fBha0_I,incb);

  vmsCbrt(vec_len,fA,fBep1,VML_EP);
  vmsCbrtI(vec_len,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsCbrt(vec_len,fA,fBep2);
  vsCbrtI(vec_len,fA_I,inca,fBep2_I,incb);

  vmsCbrt(vec_len,fA,fBla1,VML_LA);
  vmsCbrtI(vec_len,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsCbrt(vec_len,fA,fBla2);
  vsCbrtI(vec_len,fA_I,inca,fBla2_I,incb);

  vmsCbrt(vec_len,fA,fBha1,VML_HA);
  vmsCbrtI(vec_len,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsCbrt(vec_len,fA,fBha2);
  vsCbrtI(vec_len,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsCbrt and vmsCbrt in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsCbrtI and vmsCbrtI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsCbrt and vmsCbrt in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsCbrtI and vmsCbrtI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsCbrt and vmsCbrt in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsCbrtI and vmsCbrtI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsCbrt test/example program\n\n");
  printf("           Argument                     vsCbrt\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsCbrtI test/example program\n\n");
  printf("           Argument                     vsCbrt\n");
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
