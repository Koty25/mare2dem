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
!    vcCIS  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  float fA[10];
  MKL_Complex8 cB[10];
  MKL_Complex8 cBha0[10],cBha1[10],cBha2[10];
  MKL_Complex8           cBla1[10],cBla2[10];
  MKL_Complex8           cBep1[10],cBep2[10];
  float CurRMS,MaxRMS=0.0;

  float fA_I[10*INCA];
  MKL_Complex8 cB_I[10*INCB];
  MKL_Complex8 cBha0_I[10*INCB],cBha1_I[10*INCB],cBha2_I[10*INCB];
  MKL_Complex8                  cBla1_I[10*INCB],cBla2_I[10*INCB];
  MKL_Complex8                  cBep1_I[10*INCB],cBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  fA[0]=-10000.0000;
  fA[1]=-7777.7778;
  fA[2]=-5555.5557;
  fA[3]=-3333.3333;
  fA[4]=-1111.1111;
  fA[5]=1111.1111;
  fA[6]=3333.3333;
  fA[7]=5555.5557;
  fA[8]=7777.7778;
  fA[9]=10000.0000;
  cB[0].real=-9.5215536825901481e-001;cB[0].imag=3.0561438888825215e-001;
  cB[1].real=6.9269429374198155e-001;cB[1].imag=7.2123131893817349e-001;
  cB[2].real=3.4378424494653864e-001;cB[2].imag=-9.3904866376910323e-001;
  cB[3].real=-9.9465418116759940e-001;cB[3].imag=1.0326209316981849e-001;
  cB[4].real=5.2955928772685279e-001;cB[4].imag=8.4827292823844636e-001;
  cB[5].real=5.2955928772685279e-001;cB[5].imag=-8.4827292823844636e-001;
  cB[6].real=-9.9465418116759940e-001;cB[6].imag=-1.0326209316981849e-001;
  cB[7].real=3.4378424494653864e-001;cB[7].imag=9.3904866376910323e-001;
  cB[8].real=6.9269429374198155e-001;cB[8].imag=-7.2123131893817349e-001;
  cB[9].real=-9.5215536825901481e-001;cB[9].imag=-3.0561438888825215e-001;

  for(i=0;i<10;i++) {
    fA_I[i*inca]=fA[i];
    cB_I[i*incb]=cB[i];
  }

  vcCIS(vec_len,fA,cBha0);
  vcCISI(vec_len,fA_I,inca,cBha0_I,incb);

  vmcCIS(vec_len,fA,cBep1,VML_EP);
  vmcCISI(vec_len,fA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcCIS(vec_len,fA,cBep2);
  vcCISI(vec_len,fA_I,inca,cBep2_I,incb);

  vmcCIS(vec_len,fA,cBla1,VML_LA);
  vmcCISI(vec_len,fA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcCIS(vec_len,fA,cBla2);
  vcCISI(vec_len,fA_I,inca,cBla2_I,incb);

  vmcCIS(vec_len,fA,cBha1,VML_HA);
  vmcCISI(vec_len,fA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcCIS(vec_len,fA,cBha2);
  vcCISI(vec_len,fA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcCIS and vmcCIS in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcCISI and vmcCISI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcCIS and vmcCIS in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcCISI and vmcCISI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcCIS and vmcCIS in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcCISI and vmcCISI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcCIS test/example program\n\n");
  printf("           Argument                           vcCIS\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f      % .10f %+.10f*i\n",fA[i],cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcCISI test/example program\n\n");
  printf("           Argument                           vcCIS\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f      % .10f %+.10f*i\n",fA_I[i*inca],cBha0_I[i*incb].real,cBha0_I[i*incb].imag);
    CurRMS_I=crelerr(cB_I[i*incb],cBha0_I[i*incb]);
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
