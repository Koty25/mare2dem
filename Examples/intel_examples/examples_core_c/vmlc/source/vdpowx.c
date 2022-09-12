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
!    vdPowx  Example Program Text
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

  dA[0]=0.1000;
  dA[1]=0.8666;
  dA[2]=1.6333;
  dA[3]=2.4000;
  dA[4]=3.1666;
  dA[5]=3.9333;
  dA[6]=4.7000;
  dA[7]=5.4666;
  dA[8]=6.2333;
  dA[9]=7.0000;
  dB[0]=1.1660038918514482e-004;
  dB[1]=5.6940625356680863e-001;
  dB[2]=6.8873530949113730e+000;
  dB[3]=3.1295714919202130e+001;
  dB[4]=9.3107199233060683e+001;
  dB[5]=2.1845298918747244e+002;
  dB[6]=4.4011128506374024e+002;
  dB[7]=7.9737652246254356e+002;
  dB[8]=1.3361783325501383e+003;
  dB[9]=2.1087452667453094e+003;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
  }

  vdPowx(vec_len,dA,dA[5],dBha0);
  vdPowxI(vec_len,dA_I,inca,dA[5],dBha0_I,incb);

  vmdPowx(vec_len,dA,dA[5],dBep1,VML_EP);
  vmdPowxI(vec_len,dA_I,inca,dA[5],dBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vdPowx(vec_len,dA,dA[5],dBep2);
  vdPowxI(vec_len,dA_I,inca,dA[5],dBep2_I,incb);

  vmdPowx(vec_len,dA,dA[5],dBla1,VML_LA);
  vmdPowxI(vec_len,dA_I,inca,dA[5],dBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vdPowx(vec_len,dA,dA[5],dBla2);
  vdPowxI(vec_len,dA_I,inca,dA[5],dBla2_I,incb);

  vmdPowx(vec_len,dA,dA[5],dBha1,VML_HA);
  vmdPowxI(vec_len,dA_I,inca,dA[5],dBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vdPowx(vec_len,dA,dA[5],dBha2);
  vdPowxI(vec_len,dA_I,inca,dA[5],dBha2_I,incb);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdPowx and vmdPowx in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dBha1_I[i*incb]!=dBha2_I[i*incb]) {
      printf("Error! Difference between vdPowxI and vmdPowxI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdPowx and vmdPowx in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb]) {
      printf("Error! Difference between vdPowxI and vmdPowxI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdPowx and vmdPowx in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb]) {
      printf("Error! Difference between vdPowxI and vmdPowxI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdPowx test/example program\n\n");
  printf("           Argument                     vdPowx\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdPowxI test/example program\n\n");
  printf("           Argument                     vdPowx\n");
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
