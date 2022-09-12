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
!    vzExp  Example Program Text
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

  zA[0].real=-17.0000;zA[0].imag=18.0000;
  zA[1].real=-13.1111;zA[1].imag=14.1111;
  zA[2].real=-9.2222;zA[2].imag=10.2222;
  zA[3].real=-5.3333;zA[3].imag=6.3333;
  zA[4].real=-1.4444;zA[4].imag=2.4444;
  zA[5].real=2.4444;zA[5].imag=-1.4444;
  zA[6].real=6.3333;zA[6].imag=-5.3333;
  zA[7].real=10.2222;zA[7].imag=-9.2222;
  zA[8].real=14.1111;zA[8].imag=-13.1111;
  zA[9].real=18.0000;zA[9].imag=-17.0000;
  zB[0].real=2.7336700468037275e-008;zB[0].imag=-3.1090404292366861e-008;
  zB[1].real=5.2718434590867303e-008;zB[1].imag=2.0219669530527123e-006;
  zB[2].real=-6.9031805334250157e-005;zB[2].imag=-7.0712150283650660e-005;
  zB[3].real=4.8220493390719000e-003;zB[3].imag=2.4185802936554069e-004;
  zB[4].real=-1.8084266353749434e-001;zB[4].imag=1.5145585212254198e-001;
  zB[5].real=1.4526697215961737e+000;zB[5].imag=-1.1431704865446497e+001;
  zB[6].real=3.2754677379631244e+002;zB[6].imag=4.5792469439319871e+002;
  zB[7].real=-2.6944627012817422e+004;zB[7].imag=-5.5343007569010069e+003;
  zB[8].real=1.1494097020286692e+006;zB[8].imag=-6.9640053097577486e+005;
  zB[9].real=-1.8067216284192696e+007;zB[9].imag=6.3125329645518661e+007;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzExp(vec_len,zA,zBha0);
  vzExpI(vec_len,zA_I,inca,zBha0_I,incb);

  vmzExp(vec_len,zA,zBep1,VML_EP);
  vmzExpI(vec_len,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzExp(vec_len,zA,zBep2);
  vzExpI(vec_len,zA_I,inca,zBep2_I,incb);

  vmzExp(vec_len,zA,zBla1,VML_LA);
  vmzExpI(vec_len,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzExp(vec_len,zA,zBla2);
  vzExpI(vec_len,zA_I,inca,zBla2_I,incb);

  vmzExp(vec_len,zA,zBha1,VML_HA);
  vmzExpI(vec_len,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzExp(vec_len,zA,zBha2);
  vzExpI(vec_len,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzExp and vmzExp in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzExpI and vmzExpI in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzExp and vmzExp in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzExpI and vmzExpI in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzExp and vmzExp in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzExpI and vmzExpI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzExp test/example program\n\n");
  printf("           Argument                           vzExp\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA[i].real,zA[i].imag,zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzExpI test/example program\n\n");
  printf("           Argument                           vzExp\n");
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
