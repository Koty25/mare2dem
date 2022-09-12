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
!    vdLGamma  Example Program Text
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

  dA[0]=  .5;
  dA[1]= 3.5;
  dA[2]= 6.5;
  dA[3]= 9.5;
  dA[4]=12.5;
  dA[5]=15.5;
  dA[6]=18.5;
  dA[7]=21.5;
  dA[8]=24.5;
  dA[9]=27.5;
  dB[0]=5.7236494292470009e-01;
  dB[1]=1.2009736023470742e+00;
  dB[2]=5.6625620598571415e+00;
  dB[3]=1.1689333420797268e+01;
  dB[4]=1.8734347511936446e+01;
  dB[5]=2.6536914491115614e+01;
  dB[6]=3.4943315776876818e+01;
  dB[7]=4.3851925860675161e+01;
  dB[8]=5.3190494526169265e+01;
  dB[9]=6.2904990828876504e+01;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
  }

  vdLGamma(vec_len,dA,dBha0);
  vdLGammaI(vec_len,dA_I,inca,dBha0_I,incb);

  vmdLGamma(vec_len,dA,dBep1,VML_EP);
  vmdLGammaI(vec_len,dA_I,inca,dBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vdLGamma(vec_len,dA,dBep2);
  vdLGammaI(vec_len,dA_I,inca,dBep2_I,incb);

  vmdLGamma(vec_len,dA,dBla1,VML_LA);
  vmdLGammaI(vec_len,dA_I,inca,dBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vdLGamma(vec_len,dA,dBla2);
  vdLGammaI(vec_len,dA_I,inca,dBla2_I,incb);

  vmdLGamma(vec_len,dA,dBha1,VML_HA);
  vmdLGammaI(vec_len,dA_I,inca,dBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vdLGamma(vec_len,dA,dBha2);
  vdLGammaI(vec_len,dA_I,inca,dBha2_I,incb);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dBha1[i]!=dBha2[i]) {
      printf("Error! Difference between vdLGamma and vmdLGamma in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dBha1_I[i*incb]!=dBha2_I[i*incb]) {
      printf("Error! Difference between vdLGammaI and vmdLGammaI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i]) {
      printf("Error! Difference between vdLGamma and vmdLGamma in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb]) {
      printf("Error! Difference between vdLGammaI and vmdLGammaI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i]) {
      printf("Error! Difference between vdLGamma and vmdLGamma in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb]) {
      printf("Error! Difference between vdLGammaI and vmdLGammaI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdLGamma test/example program\n\n");
  printf("           Argument                     vdLGamma\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",dA[i],dBha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdLGammaI test/example program\n\n");
  printf("           Argument                     vdLGamma\n");
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
