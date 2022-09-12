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
!    vdCdfNorm  Example Program Text
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

  dA[0]=-8.0000;
  dA[1]=-6.5555;
  dA[2]=-5.1111;
  dA[3]=-3.6666;
  dA[4]=-2.2222;
  dA[5]=-0.7777;
  dA[6]=0.6666;
  dA[7]=2.1111;
  dA[8]=3.5555;
  dA[9]=5.0000;
  dB[0]=6.22096057427178e-16;
  dB[1]=2.77278746107643e-11;
  dB[2]=1.60144189579716e-07;
  dB[3]=1.22898411922205e-04;
  dB[4]=1.31348962351314e-02;
  dB[5]=2.18372946100606e-01;
  dB[6]=7.47486165446128e-01;
  dB[7]=9.82618141277281e-01;
  dB[8]=9.99811369527550e-01;
  dB[9]=9.99999713348428e-01;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
  }

  vdCdfNorm(vec_len,dA,dBha0);
  vdCdfNormI(vec_len,dA_I,inca,dBha0_I,incb);

  vmdCdfNorm(vec_len,dA,dBep1,VML_EP);
  vmdCdfNormI(vec_len,dA_I,inca,dBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vdCdfNorm(vec_len,dA,dBep2);
  vdCdfNormI(vec_len,dA_I,inca,dBep2_I,incb);

  vmdCdfNorm(vec_len,dA,dBla1,VML_LA);
  vmdCdfNormI(vec_len,dA_I,inca,dBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vdCdfNorm(vec_len,dA,dBla2);
  vdCdfNormI(vec_len,dA_I,inca,dBla2_I,incb);

  vmdCdfNorm(vec_len,dA,dBha1,VML_HA);
  vmdCdfNormI(vec_len,dA_I,inca,dBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vdCdfNorm(vec_len,dA,dBha2);
  vdCdfNormI(vec_len,dA_I,inca,dBha2_I,incb);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdCdfNorm and vmdCdfNorm in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dBha1_I[i*incb]!=dBha2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormI and vmdCdfNormI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdCdfNorm and vmdCdfNorm in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormI and vmdCdfNormI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdCdfNorm and vmdCdfNorm in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb]) {
      printf("Error! Difference between vdCdfNormI and vmdCdfNormI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdCdfNorm test/example program\n\n");
  printf("           Argument                     vdCdfNorm\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdCdfNormI test/example program\n\n");
  printf("           Argument                     vdCdfNorm\n");
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
