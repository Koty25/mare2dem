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
!    vsAcosh  Example Program Text
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

  fA[0]=3.5000;
  fA[1]=1114.2222;
  fA[2]=2224.9443;
  fA[3]=3335.6665;
  fA[4]=4446.3887;
  fA[5]=5557.1113;
  fA[6]=6667.8335;
  fA[7]=7778.5557;
  fA[8]=8889.2773;
  fA[9]=10000.0000;
  fB[0]=1.9248473002384139e+000;
  fB[1]=7.7090588124392951e+000;
  fB[2]=8.4006343067899572e+000;
  fB[3]=8.8055749477197853e+000;
  fB[4]=9.0929946791926213e+000;
  fB[5]=9.3159808793974239e+000;
  fB[6]=9.4981974479040083e+000;
  fB[7]=9.6522731290732260e+000;
  fB[8]=9.7857482139329655e+000;
  fB[9]=9.9034875500361288e+000;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsAcosh(vec_len,fA,fBha0);
  vsAcoshI(vec_len,fA_I,inca,fBha0_I,incb);

  vmsAcosh(vec_len,fA,fBep1,VML_EP);
  vmsAcoshI(vec_len,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsAcosh(vec_len,fA,fBep2);
  vsAcoshI(vec_len,fA_I,inca,fBep2_I,incb);

  vmsAcosh(vec_len,fA,fBla1,VML_LA);
  vmsAcoshI(vec_len,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsAcosh(vec_len,fA,fBla2);
  vsAcoshI(vec_len,fA_I,inca,fBla2_I,incb);

  vmsAcosh(vec_len,fA,fBha1,VML_HA);
  vmsAcoshI(vec_len,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsAcosh(vec_len,fA,fBha2);
  vsAcoshI(vec_len,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsAcosh and vmsAcosh in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsAcoshI and vmsAcoshI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsAcosh and vmsAcosh in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsAcoshI and vmsAcoshI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsAcosh and vmsAcosh in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsAcoshI and vmsAcoshI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsAcosh test/example program\n\n");
  printf("           Argument                     vsAcosh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsAcoshI test/example program\n\n");
  printf("           Argument                     vsAcosh\n");
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
