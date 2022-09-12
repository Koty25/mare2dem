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
!    vsPowx  Example Program Text
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

  fA[0]=0.1000;
  fA[1]=0.8666;
  fA[2]=1.6333;
  fA[3]=2.4000;
  fA[4]=3.1666;
  fA[5]=3.9333;
  fA[6]=4.7000;
  fA[7]=5.4666;
  fA[8]=6.2333;
  fA[9]=7.0000;
  fB[0]=1.1660039110312152e-004;
  fB[1]=5.6940619267520864e-001;
  fB[2]=6.8873522741537876e+000;
  fB[3]=3.1295720312255717e+001;
  fB[4]=9.3107199918766597e+001;
  fB[5]=2.1845299866537277e+002;
  fB[6]=4.4011122728416541e+002;
  fB[7]=7.9737651355941364e+002;
  fB[8]=1.3361785535769964e+003;
  fB[9]=2.1087453418813361e+003;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    fB_I[i*incb]=fB[i];
  }

  vsPowx(vec_len,fA,fA[5],fBha0);
  vsPowxI(vec_len,fA_I,inca,fA[5],fBha0_I,incb);

  vmsPowx(vec_len,fA,fA[5],fBep1,VML_EP);
  vmsPowxI(vec_len,fA_I,inca,fA[5],fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vsPowx(vec_len,fA,fA[5],fBep2);
  vsPowxI(vec_len,fA_I,inca,fA[5],fBep2_I,incb);

  vmsPowx(vec_len,fA,fA[5],fBla1,VML_LA);
  vmsPowxI(vec_len,fA_I,inca,fA[5],fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vsPowx(vec_len,fA,fA[5],fBla2);
  vsPowxI(vec_len,fA_I,inca,fA[5],fBla2_I,incb);

  vmsPowx(vec_len,fA,fA[5],fBha1,VML_HA);
  vmsPowxI(vec_len,fA_I,inca,fA[5],fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vsPowx(vec_len,fA,fA[5],fBha2);
  vsPowxI(vec_len,fA_I,inca,fA[5],fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vsPowx and vmsPowx in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vsPowxI and vmsPowxI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vsPowx and vmsPowx in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vsPowxI and vmsPowxI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vsPowx and vmsPowx in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vsPowxI and vmsPowxI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vsPowx test/example program\n\n");
  printf("           Argument                     vsPowx\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e\n",fA[i],fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vsPowxI test/example program\n\n");
  printf("           Argument                     vsPowx\n");
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
