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
!    vsExpInt1  Example Program Text
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

  fA[0]=0.0099;
  fA[1]=1.1188;
  fA[2]=2.2277;
  fA[3]=3.3366;
  fA[4]=4.4455;
  fA[5]=5.5544;
  fA[6]=6.6633;
  fA[7]=7.7722;
  fA[8]=8.8811;
  fA[9]=10.9900;
  fB[0]=4.0478804136108807e+000;
  fB[1]=1.8040243730827092e-001;
  fB[2]=3.5823720740759158e-002;
  fB[3]=8.5395671171776873e-003;
  fB[4]=2.2125289768999900e-003;
  fB[5]=6.0177553854344804e-004;
  fB[6]=1.6900054139095165e-004;
  fB[7]=4.8555131001862730e-005;
  fB[8]=1.4190391394966996e-005;
  fB[9]=1.4155671423752270e-006;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsExpInt1(vec_len,fA,fBha0);
  vsExpInt1I(vec_len,fA_I,inca,fBha0_I,incb);

  vmsExpInt1(vec_len,fA,fBep1,VML_EP);
  vmsExpInt1I(vec_len,fA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsExpInt1(vec_len,fA,fBep2);
  vsExpInt1I(vec_len,fA_I,inca,fBep2_I,incb);

  vmsExpInt1(vec_len,fA,fBla1,VML_LA);
  vmsExpInt1I(vec_len,fA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsExpInt1(vec_len,fA,fBla2);
  vsExpInt1I(vec_len,fA_I,inca,fBla2_I,incb);

  vmsExpInt1(vec_len,fA,fBha1,VML_HA);
  vmsExpInt1I(vec_len,fA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsExpInt1(vec_len,fA,fBha2);
  vsExpInt1I(vec_len,fA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsExpInt1 and vmsExpInt1 in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsExpInt1I and vmsExpInt1I in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsExpInt1 and vmsExpInt1 in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsExpInt1I and vmsExpInt1I in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsExpInt1 and vmsExpInt1 in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsExpInt1I and vmsExpInt1I in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsExpInt1 test/example program\n\n");
  printf("           Argument                     vsExpInt1\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsExpInt1I test/example program\n\n");
  printf("           Argument                     vsExpInt1\n");
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
