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
!    vzSqrt  Example Program Text
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

  zA[0].real=0.0000;zA[0].imag=10000.0000;
  zA[1].real=1111.1111;zA[1].imag=8888.8888;
  zA[2].real=2222.2222;zA[2].imag=7777.7777;
  zA[3].real=3333.3333;zA[3].imag=6666.6666;
  zA[4].real=4444.4444;zA[4].imag=5555.5555;
  zA[5].real=5555.5555;zA[5].imag=4444.4444;
  zA[6].real=6666.6666;zA[6].imag=3333.3333;
  zA[7].real=7777.7777;zA[7].imag=2222.2222;
  zA[8].real=8888.8888;zA[8].imag=1111.1111;
  zA[9].real=10000.0000;zA[9].imag=0.0000;
  zB[0].real=7.0710678118654755e+001;zB[0].imag=7.0710678118654755e+001;
  zB[1].real=7.0954827796266002e+001;zB[1].imag=6.2637660297921109e+001;
  zB[2].real=7.1802622191669840e+001;zB[2].imag=5.4160819358644105e+001;
  zB[3].real=7.3440088338943667e+001;zB[3].imag=4.5388470730262000e+001;
  zB[4].real=7.6023111008727710e+001;zB[4].imag=3.6538596134024317e+001;
  zB[5].real=7.9593146422574208e+001;zB[5].imag=2.7919768219763874e+001;
  zB[6].real=8.4024479916461544e+001;zB[6].imag=1.9835489034350779e+001;
  zB[7].real=8.9069603701822302e+001;zB[7].imag=1.2474638415588544e+001;
  zB[8].real=9.4464154246982829e+001;zB[8].imag=5.8811255383440271e+000;
  zB[9].real=1.0000000000000000e+002;zB[9].imag=0.0000000000000000e+000;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzSqrt(vec_len,zA,zBha0);
  vzSqrtI(vec_len,zA_I,inca,zBha0_I,incb);

  vmzSqrt(vec_len,zA,zBep1,VML_EP);
  vmzSqrtI(vec_len,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzSqrt(vec_len,zA,zBep2);
  vzSqrtI(vec_len,zA_I,inca,zBep2_I,incb);

  vmzSqrt(vec_len,zA,zBla1,VML_LA);
  vmzSqrtI(vec_len,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzSqrt(vec_len,zA,zBla2);
  vzSqrtI(vec_len,zA_I,inca,zBla2_I,incb);

  vmzSqrt(vec_len,zA,zBha1,VML_HA);
  vmzSqrtI(vec_len,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzSqrt(vec_len,zA,zBha2);
  vzSqrtI(vec_len,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzSqrt and vmzSqrt in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzSqrtI and vmzSqrtI in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzSqrt and vmzSqrt in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzSqrtI and vmzSqrtI in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzSqrt and vmzSqrt in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzSqrtI and vmzSqrtI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzSqrt test/example program\n\n");
  printf("           Argument                           vzSqrt\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA[i].real,zA[i].imag,zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzSqrtI test/example program\n\n");
  printf("           Argument                           vzSqrt\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA_I[i*inca].real,zA_I[i*inca].imag,zBha0_I[i*incb].real,
        zBha0_I[i*incb].imag);
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
