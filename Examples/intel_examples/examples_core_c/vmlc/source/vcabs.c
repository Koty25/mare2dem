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
!    vcAbs  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  MKL_Complex8 cA[10];
  float fB[10];
  float fBha0[10],fBha1[10],fBha2[10];
  float           fBla1[10],fBla2[10];
  float           fBep1[10],fBep2[10];
  float CurRMS,MaxRMS=0.0;

  MKL_Complex8 cA_I[10*INCA];
  float fB_I[10*INCB];
  float fBha0_I[10*INCB],fBha1_I[10*INCB],fBha2_I[10*INCB];
  float                  fBla1_I[10*INCB],fBla2_I[10*INCB];
  float                  fBep1_I[10*INCB],fBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  cA[0].real=-1000.0000;cA[0].imag=1000.0000;
  cA[1].real=-777.7777;cA[1].imag=777.7777;
  cA[2].real=-555.5555;cA[2].imag=555.5555;
  cA[3].real=-333.3333;cA[3].imag=333.3333;
  cA[4].real=-111.1111;cA[4].imag=111.1111;
  cA[5].real=111.1111;cA[5].imag=-111.1111;
  cA[6].real=333.3333;cA[6].imag=-333.3333;
  cA[7].real=555.5555;cA[7].imag=-555.5555;
  cA[8].real=777.7777;cA[8].imag=-777.7777;
  cA[9].real=1000.0000;cA[9].imag=-1000.0000;
  fB[0]=1.4142135623730951e+003;
  fB[1]=1.0999437859382454e+003;
  fB[2]=7.8567409582014136e+002;
  fB[3]=4.7140449201878312e+002;
  fB[4]=1.5713482347986556e+002;
  fB[5]=1.5713482347986556e+002;
  fB[6]=4.7140449201878312e+002;
  fB[7]=7.8567409582014136e+002;
  fB[8]=1.0999437859382454e+003;
  fB[9]=1.4142135623730951e+003;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    fB_I[i*incb]=fB[i];
  }

  vcAbs(vec_len,cA,fBha0);
  vcAbsI(vec_len,cA_I,inca,fBha0_I,incb);

  vmcAbs(vec_len,cA,fBep1,VML_EP);
  vmcAbsI(vec_len,cA_I,inca,fBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcAbs(vec_len,cA,fBep2);
  vcAbsI(vec_len,cA_I,inca,fBep2_I,incb);

  vmcAbs(vec_len,cA,fBla1,VML_LA);
  vmcAbsI(vec_len,cA_I,inca,fBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcAbs(vec_len,cA,fBla2);
  vcAbsI(vec_len,cA_I,inca,fBla2_I,incb);

  vmcAbs(vec_len,cA,fBha1,VML_HA);
  vmcAbsI(vec_len,cA_I,inca,fBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcAbs(vec_len,cA,fBha2);
  vcAbsI(vec_len,cA_I,inca,fBha2_I,incb);

  for(i=0;i<10;i++) {
    if(fBha0[i]!=fBha1[i] || fBha1[i]!=fBha2[i]) {
      printf("Error! Difference between vcAbs and vmcAbs in VML_HA mode detected.\n");
      return 1;
    }

    if(fBha0_I[i*incb]!=fBha1_I[i*incb] || fBha1_I[i*incb]!=fBha2_I[i*incb]) {
      printf("Error! Difference between vcAbsI and vmcAbsI in VML_HA mode detected.\n");
      return 1;
    }

    if(fBla1[i]!=fBla2[i]) {
      printf("Error! Difference between vcAbs and vmcAbs in VML_LA mode detected.\n");
      return 1;
    }

    if(fBla1_I[i*incb]!=fBla2_I[i*incb]) {
      printf("Error! Difference between vcAbsI and vmcAbsI in VML_LA mode detected.\n");
      return 1;
    }

    if(fBep1[i]!=fBep2[i]) {
      printf("Error! Difference between vcAbs and vmcAbs in VML_EP mode detected.\n");
      return 1;
    }

    if(fBep1_I[i*incb]!=fBep2_I[i*incb]) {
      printf("Error! Difference between vcAbsI and vmcAbsI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcAbs test/example program\n\n");
  printf("           Argument                           vcAbs\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f\n",cA[i].real,cA[i].imag,fBha0[i]);
    CurRMS=srelerr(fB[i],fBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcAbsI test/example program\n\n");
  printf("           Argument                           vcAbs\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f\n",cA_I[i*inca].real,cA_I[i*inca].imag,fBha0_I[i*incb]);
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
