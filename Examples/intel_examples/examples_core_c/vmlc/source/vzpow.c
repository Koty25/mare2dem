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
!    vzPow  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  MKL_Complex16 zA[10],zB[10];
  MKL_Complex16 zBha0[10],zBha1[10],zBha2[10];
  MKL_Complex16           zBla1[10],zBla2[10];
  MKL_Complex16           zBep1[10],zBep2[10];
  float CurRMS,MaxRMS=0.0;

  MKL_Complex16 zA_I[10*INCA],zB_I[10*INCB];
  MKL_Complex16 zBha0_I[10*INCB],zBha1_I[10*INCB],zBha2_I[10*INCB];
  MKL_Complex16                  zBla1_I[10*INCB],zBla2_I[10*INCB];
  MKL_Complex16                  zBep1_I[10*INCB],zBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  zA[0].real=0.1000;zA[0].imag=7.0000;
  zA[1].real=0.8666;zA[1].imag=6.2333;
  zA[2].real=1.6333;zA[2].imag=5.4666;
  zA[3].real=2.4000;zA[3].imag=4.6999;
  zA[4].real=3.1666;zA[4].imag=3.9333;
  zA[5].real=3.9333;zA[5].imag=3.1666;
  zA[6].real=4.7000;zA[6].imag=2.3999;
  zA[7].real=5.4666;zA[7].imag=1.6333;
  zA[8].real=6.2333;zA[8].imag=0.8666;
  zA[9].real=7.0000;zA[9].imag=0.0999;
  zB[0].real=7.9222171309252660e-006;zB[0].imag=2.1083585790345553e-005;
  zB[1].real=6.4512100825773902e-004;zB[1].imag=9.1691995410974377e-005;
  zB[2].real=9.0503657535304545e-003;zB[2].imag=-1.2801470190500392e-002;
  zB[3].real=-1.5959437975218965e-001;zB[3].imag=-2.6566220521977135e-001;
  zB[4].real=-4.8997205920402100e+000;zB[4].imag=1.1363610159141824e+000;
  zB[5].real=4.1042328112091866e+000;zB[5].imag=6.8100661769434140e+001;
  zB[6].real=7.9823060676727789e+002;zB[6].imag=-5.7795808684580606e+001;
  zB[7].real=-2.3514338397471715e+003;zB[7].imag=-8.1467677119527843e+003;
  zB[8].real=-6.5479995707607908e+004;zB[8].imag=5.3650124957407497e+004;
  zB[9].real=7.8757355451202882e+005;zB[9].imag=2.3871476187431620e+005;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzPow(vec_len,zA,zA,zBha0);
  vzPowI(vec_len,zA_I,inca,zA_I,inca,zBha0_I,incb);

  vmzPow(vec_len,zA,zA,zBep1,VML_EP);
  vmzPowI(vec_len,zA_I,inca,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzPow(vec_len,zA,zA,zBep2);
  vzPowI(vec_len,zA_I,inca,zA_I,inca,zBep2_I,incb);

  vmzPow(vec_len,zA,zA,zBla1,VML_LA);
  vmzPowI(vec_len,zA_I,inca,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzPow(vec_len,zA,zA,zBla2);
  vzPowI(vec_len,zA_I,inca,zA_I,inca,zBla2_I,incb);

  vmzPow(vec_len,zA,zA,zBha1,VML_HA);
  vmzPowI(vec_len,zA_I,inca,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzPow(vec_len,zA,zA,zBha2);
  vzPowI(vec_len,zA_I,inca,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzPow and vmzPow in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzPowI and vmzPowI in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzPow and vmzPow in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzPowI and vmzPowI in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzPow and vmzPow in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzPowI and vmzPowI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzPow test/example program\n\n");
  printf("           Arguments                           vzPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .5f %+.5f*i\n",zA[i].real,zA[i].imag,zA[i].real,zA[i].imag,
        zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzPowI test/example program\n\n");
  printf("           Arguments                           vzPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .5f %+.5f*i\n",zA_I[i*inca].real,zA_I[i*inca].imag,
        zA_I[i*inca].real,zA_I[i*inca].imag,zBha0_I[i*incb].real,zBha0_I[i*incb].imag);
    CurRMS_I=zrelerr(zB_I[i*incb],zBha0_I[i*incb]);
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
