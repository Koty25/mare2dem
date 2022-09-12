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
!    vzAtanh  Example Program Text
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

  zA[0].real=0.7500;zA[0].imag=0.9000;
  zA[1].real=0.7666;zA[1].imag=0.8833;
  zA[2].real=0.7833;zA[2].imag=0.8666;
  zA[3].real=0.8000;zA[3].imag=0.8500;
  zA[4].real=0.8166;zA[4].imag=0.8333;
  zA[5].real=0.8333;zA[5].imag=0.8166;
  zA[6].real=0.8500;zA[6].imag=0.8000;
  zA[7].real=0.8666;zA[7].imag=0.7833;
  zA[8].real=0.8833;zA[8].imag=0.7666;
  zA[9].real=0.9000;zA[9].imag=0.7500;
  zB[0].real=3.7257322955361039e-001;zB[0].imag=8.8743011636078650e-001;
  zB[1].real=3.8548667280375853e-001;zB[1].imag=8.8805592388369625e-001;
  zB[2].real=3.9865922482169486e-001;zB[2].imag=8.8905517994623651e-001;
  zB[3].real=4.1200697742599235e-001;zB[3].imag=8.9044257000699012e-001;
  zB[4].real=4.2558003685990275e-001;zB[4].imag=8.9211910143550066e-001;
  zB[5].real=4.3940258266500865e-001;zB[5].imag=8.9423397552554096e-001;
  zB[6].real=4.5338331719440389e-001;zB[6].imag=8.9679874027714557e-001;
  zB[7].real=4.6758174382061612e-001;zB[7].imag=8.9971577920437273e-001;
  zB[8].real=4.8200889769925359e-001;zB[8].imag=9.0314797406567959e-001;
  zB[9].real=4.9656448002689885e-001;zB[9].imag=9.0710287418476410e-001;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzAtanh(vec_len,zA,zBha0);
  vzAtanhI(vec_len,zA_I,inca,zBha0_I,incb);

  vmzAtanh(vec_len,zA,zBep1,VML_EP);
  vmzAtanhI(vec_len,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzAtanh(vec_len,zA,zBep2);
  vzAtanhI(vec_len,zA_I,inca,zBep2_I,incb);

  vmzAtanh(vec_len,zA,zBla1,VML_LA);
  vmzAtanhI(vec_len,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzAtanh(vec_len,zA,zBla2);
  vzAtanhI(vec_len,zA_I,inca,zBla2_I,incb);

  vmzAtanh(vec_len,zA,zBha1,VML_HA);
  vmzAtanhI(vec_len,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzAtanh(vec_len,zA,zBha2);
  vzAtanhI(vec_len,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzAtanh and vmzAtanh in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzAtanhI and vmzAtanhI in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzAtanh and vmzAtanh in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzAtanhI and vmzAtanhI in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzAtanh and vmzAtanh in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzAtanhI and vmzAtanhI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzAtanh test/example program\n\n");
  printf("           Argument                           vzAtanh\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA[i].real,zA[i].imag,zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzAtanhI test/example program\n\n");
  printf("           Argument                           vzAtanh\n");
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
