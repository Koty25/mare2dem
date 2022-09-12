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
!    vzAsinh  Example Program Text
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

  zA[0].real=3.0000;zA[0].imag=10000.0000;
  zA[1].real=1113.7777;zA[1].imag=8889.2222;
  zA[2].real=2224.5555;zA[2].imag=7778.4444;
  zA[3].real=3335.3333;zA[3].imag=6667.6666;
  zA[4].real=4446.1111;zA[4].imag=5556.8888;
  zA[5].real=5556.8888;zA[5].imag=4446.1111;
  zA[6].real=6667.6666;zA[6].imag=3335.3333;
  zA[7].real=7778.4444;zA[7].imag=2224.5555;
  zA[8].real=8889.2222;zA[8].imag=1113.7777;
  zA[9].real=10000.0000;zA[9].imag=3.0000;
  zB[0].real=9.9034875950361272e+000;zB[0].imag=1.5704963268023966e+000;
  zB[1].real=9.7935304892829596e+000;zB[1].imag=1.4461505920337059e+000;
  zB[2].real=9.6915674320018415e+000;zB[2].imag=1.2922419942501833e+000;
  zB[3].real=9.6098341947868935e+000;zB[3].imag=1.1069687573846247e+000;
  zB[4].real=9.5633416390116412e+000;zB[4].imag=8.9598953878147525e-001;
  zB[5].real=9.5633416411774572e+000;zB[5].imag=6.7480677838165393e-001;
  zB[6].real=9.6098342001817105e+000;zB[6].imag=4.6382756221178462e-001;
  zB[7].real=9.6915674384858015e+000;zB[7].imag=2.7855432850566603e-001;
  zB[8].real=9.7935304953202280e+000;zB[8].imag=1.2464573322417867e-001;
  zB[9].real=9.9034876000361258e+000;zB[9].imag=2.9999998950000077e-004;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzAsinh(vec_len,zA,zBha0);
  vzAsinhI(vec_len,zA_I,inca,zBha0_I,incb);

  vmzAsinh(vec_len,zA,zBep1,VML_EP);
  vmzAsinhI(vec_len,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzAsinh(vec_len,zA,zBep2);
  vzAsinhI(vec_len,zA_I,inca,zBep2_I,incb);

  vmzAsinh(vec_len,zA,zBla1,VML_LA);
  vmzAsinhI(vec_len,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzAsinh(vec_len,zA,zBla2);
  vzAsinhI(vec_len,zA_I,inca,zBla2_I,incb);

  vmzAsinh(vec_len,zA,zBha1,VML_HA);
  vmzAsinhI(vec_len,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzAsinh(vec_len,zA,zBha2);
  vzAsinhI(vec_len,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzAsinh and vmzAsinh in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzAsinhI and vmzAsinhI in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzAsinh and vmzAsinh in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzAsinhI and vmzAsinhI in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzAsinh and vmzAsinh in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzAsinhI and vmzAsinhI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzAsinh test/example program\n\n");
  printf("           Argument                           vzAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA[i].real,zA[i].imag,zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzAsinhI test/example program\n\n");
  printf("           Argument                           vzAsinh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA_I[i*inca].real,zA_I[i*inca].imag,zBha0_I[i*incb].real,
        zBha0_I[i*incb].imag);
    CurRMS_I=zrelerr(zB_I[i*incb],zBha0_I[i*incb]);
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
