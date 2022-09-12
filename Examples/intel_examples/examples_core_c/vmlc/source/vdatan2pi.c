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
!     vdAtan2pi Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5
#define INCC 7

int main()
{
  double dA1[10],dA2[10],dB[10];
  double dBha0[10],dBha1[10],dBha2[10];
  double           dBla1[10],dBla2[10];
  double           dBep1[10],dBep2[10];
  float CurRMS,MaxRMS=0.0;

  double dA1_I[10*INCA],dA2_I[10*INCB],dB_I[10*INCC];
  double dBha0_I[10*INCC],dBha1_I[10*INCC],dBha2_I[10*INCC];
  double                  dBla1_I[10*INCC],dBla2_I[10*INCC];
  double                  dBep1_I[10*INCC],dBep2_I[10*INCC];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB,incc=INCC;

  MKL_INT i=0,vec_len=10;

  dA1[0]=0.1000;
  dA1[1]=0.8666;
  dA1[2]=1.6333;
  dA1[3]=2.4000;
  dA1[4]=3.1666;
  dA1[5]=3.9333;
  dA1[6]=4.7000;
  dA1[7]=5.4666;
  dA1[8]=6.2333;
  dA1[9]=7.0000;
  dA2[0]=-10.0000;
  dA2[1]=-7.7777;
  dA2[2]=-5.5555;
  dA2[3]=-3.3333;
  dA2[4]=-1.1111;
  dA2[5]=1.1111;
  dA2[6]=3.3333;
  dA2[7]=5.5555;
  dA2[8]=7.7777;
  dA2[9]=10.0000;
  dB[0]=9.9681700723509170e-01;
  dB[1]=9.6467924216936096e-01;
  dB[2]=9.0898236246242259e-01;
  dB[3]=8.0136578364848177e-01;
  dB[4]=6.0741670213403143e-01;
  dB[5]=4.1236547979917204e-01;
  dB[6]=3.0364026366201680e-01;
  dB[7]=2.4743269111355592e-01;
  dB[8]=2.1505451546684570e-01;
  dB[9]=1.9440011221421480e-01;

  for(i=0;i<10;i++) {
    dA1_I[i*inca]=dA1[i];
    dA2_I[i*incb]=dA2[i];
    dB_I[i*incc]=dB[i];
  }

  vdAtan2pi(vec_len,dA1,dA2,dBha0);
  vdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBha0_I,incc);

  vmdAtan2pi(vec_len,dA1,dA2,dBep1,VML_EP);
  vmdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBep1_I,incc,VML_EP);

  vmlSetMode(VML_EP);
  vdAtan2pi(vec_len,dA1,dA2,dBep2);
  vdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBep2_I,incc);

  vmdAtan2pi(vec_len,dA1,dA2,dBla1,VML_LA);
  vmdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBla1_I,incc,VML_LA);

  vmlSetMode(VML_LA);
  vdAtan2pi(vec_len,dA1,dA2,dBla2);
  vdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBla2_I,incc);

  vmdAtan2pi(vec_len,dA1,dA2,dBha1,VML_HA);
  vmdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBha1_I,incc,VML_HA);

  vmlSetMode(VML_HA);
  vdAtan2pi(vec_len,dA1,dA2,dBha2);
  vdAtan2piI(vec_len,dA1_I,inca,dA2_I,incb,dBha2_I,incc);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdAtan2pi and vmdAtan2pi in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incc]!=dBha1_I[i*incc] || dBha1_I[i*incc]!=dBha2_I[i*incc]) {
      printf("Error! Difference between vdAtan2piI and vmdAtan2piI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdAtan2pi and vmdAtan2pi in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incc]!=dBla2_I[i*incc]) {
      printf("Error! Difference between vdAtan2piI and vmdAtan2piI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdAtan2pi and vmdAtan2pi in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incc]!=dBep2_I[i*incc]) {
      printf("Error! Difference between vdAtan2piI and vmdAtan2piI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdAtan2pi test/example program\n\n");
  printf("           Argument                     vdAtan2pi\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f %25.14f % 25.14e\n",dA1[i],dA2[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdAtan2piI test/example program\n\n");
  printf("           Argument                     vdAtan2pi\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f %25.14f % 25.14e\n",dA1_I[i*inca],dA2_I[i*incb],dBha0_I[i*incc]);
    CurRMS_I=drelerr(dB_I[i*incc],dBha0_I[i*incc]);
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