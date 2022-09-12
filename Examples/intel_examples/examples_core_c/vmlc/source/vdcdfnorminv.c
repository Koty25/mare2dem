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
!    vdCdfNormInv  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  double dA[10],dB[10];
  double dBha0[10],dBha1[10],dBha2[10];
  double           dBla1[10],dBla2[10];
  double           dBep1[10],dBep2[10];
  float CurRMS,MaxRMS=0.0;

  double dA_I[10*INCA],dB_I[10*INCB];
  double dBha0_I[10*INCB],dBha1_I[10*INCB],dBha2_I[10*INCB];
  double                  dBla1_I[10*INCB],dBla2_I[10*INCB];
  double                  dBep1_I[10*INCB],dBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  dA[0]=0.0001;
  dA[1]=0.1111;
  dA[2]=0.2222;
  dA[3]=0.3333;
  dA[4]=0.4444;
  dA[5]=0.5555;
  dA[6]=0.6666;
  dA[7]=0.7777;
  dA[8]=0.8888;
  dA[9]=0.9999;
  dB[0]=-3.71901648545568e+00;
  dB[1]=-1.22069901725486e+00;
  dB[2]=-7.64784296900383e-01;
  dB[3]=-4.30818976996713e-01;
  dB[4]=-1.39822798052277e-01;
  dB[5]= 1.39569677405137e-01;
  dB[6]= 4.30543954751993e-01;
  dB[7]= 7.64448526414231e-01;
  dB[8]= 1.22017115276020e+00;
  dB[9]= 3.71901648545568e+00;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
  }

  vdCdfNormInv(vec_len,dA,dBha0);
  vdCdfNormInvI(vec_len,dA_I,inca,dBha0_I,incb);

  vmdCdfNormInv(vec_len,dA,dBep1,VML_EP);
  vmdCdfNormInvI(vec_len,dA_I,inca,dBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vdCdfNormInv(vec_len,dA,dBep2);
  vdCdfNormInvI(vec_len,dA_I,inca,dBep2_I,incb);

  vmdCdfNormInv(vec_len,dA,dBla1,VML_LA);
  vmdCdfNormInvI(vec_len,dA_I,inca,dBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vdCdfNormInv(vec_len,dA,dBla2);
  vdCdfNormInvI(vec_len,dA_I,inca,dBla2_I,incb);

  vmdCdfNormInv(vec_len,dA,dBha1,VML_HA);
  vmdCdfNormInvI(vec_len,dA_I,inca,dBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vdCdfNormInv(vec_len,dA,dBha2);
  vdCdfNormInvI(vec_len,dA_I,inca,dBha2_I,incb);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdCdfNormInv and vmdCdfNormInv in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dBha1_I[i*incb]!=dBha2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormInvI and vmdCdfNormInvI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdCdfNormInv and vmdCdfNormInv in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormInvI and vmdCdfNormInvI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdCdfNormInv and vmdCdfNormInv in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormInvI and vmdCdfNormInvI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdCdfNormInv test/example program\n\n");
  printf("           Argument                     vdCdfNormInv\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdCdfNormInvI test/example program\n\n");
  printf("           Argument                     vdCdfNormInv\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA_I[i*inca],dBha0_I[i*incb]);
    CurRMS_I=drelerr(dB_I[i*incb],dBha0_I[i*incb]);
    if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I;
  }
  printf("\n");
  if(MaxRMS>=1e-14) {
    printf("Error! Relative accuracy is %.16f\n",MaxRMS);
    return 1;
  }
  else {
    printf("Relative accuracy is %.16f\n",MaxRMS);
  }

  printf("\n");

  if(MaxRMS_I>=1e-14) {
    printf("Error! Relative strided accuracy is %.16f\n",MaxRMS_I);
    return 1;
  }
  else {
    printf("Relative strided accuracy is %.16f\n",MaxRMS_I);
  }
  return 0;
}
