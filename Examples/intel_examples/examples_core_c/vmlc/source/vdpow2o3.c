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
!    vdPow2o3  Example Program Text
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

  dA[0]=1.0000;
  dA[1]=1.6666;
  dA[2]=2.3333;
  dA[3]=3.0000;
  dA[4]=3.6666;
  dA[5]=4.3333;
  dA[6]=5.0000;
  dA[7]=5.6666;
  dA[8]=6.3333;
  dA[9]=7.0000;
  dB[0]=1.0000000000000000e+000;
  dB[1]=1.4056836226901026e+000;
  dB[2]=1.7591939415405531e+000;
  dB[3]=2.0800838230519041e+000;
  dB[4]=2.3778020076302173e+000;
  dB[5]=2.6579440692151435e+000;
  dB[6]=2.9240177382128660e+000;
  dB[7]=3.1784474694102087e+000;
  dB[8]=3.4231035769933422e+000;
  dB[9]=3.6593057100229713e+000;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
  }

  vdPow2o3(vec_len,dA,dBha0);
  vdPow2o3I(vec_len,dA_I,inca,dBha0_I,incb);

  vmdPow2o3(vec_len,dA,dBep1,VML_EP);
  vmdPow2o3I(vec_len,dA_I,inca,dBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vdPow2o3(vec_len,dA,dBep2);
  vdPow2o3I(vec_len,dA_I,inca,dBep2_I,incb);

  vmdPow2o3(vec_len,dA,dBla1,VML_LA);
  vmdPow2o3I(vec_len,dA_I,inca,dBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vdPow2o3(vec_len,dA,dBla2);
  vdPow2o3I(vec_len,dA_I,inca,dBla2_I,incb);

  vmdPow2o3(vec_len,dA,dBha1,VML_HA);
  vmdPow2o3I(vec_len,dA_I,inca,dBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vdPow2o3(vec_len,dA,dBha2);
  vdPow2o3I(vec_len,dA_I,inca,dBha2_I,incb);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdPow2o3 and vmdPow2o3 in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dBha1_I[i*incb]!=dBha2_I[i*incb]) {
      printf("Error! Difference between vdPow2o3I and vmdPow2o3I in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdPow2o3 and vmdPow2o3 in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb]) {
      printf("Error! Difference between vdPow2o3I and vmdPow2o3I in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdPow2o3 and vmdPow2o3 in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb]) {
      printf("Error! Difference between vdPow2o3I and vmdPow2o3I in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdPow2o3 test/example program\n\n");
  printf("           Argument                     vdPow2o3\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdPow2o3I test/example program\n\n");
  printf("           Argument                     vdPow2o3\n");
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
