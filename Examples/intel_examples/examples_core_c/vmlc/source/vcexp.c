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
!    vcExp  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5

int main()
{
  MKL_Complex8 cA[10],cB[10];
  MKL_Complex8 cBha0[10],cBha1[10],cBha2[10];
  MKL_Complex8           cBla1[10],cBla2[10];
  MKL_Complex8           cBep1[10],cBep2[10];
  float CurRMS,MaxRMS=0.0;

  MKL_Complex8 cA_I[10*INCA],cB_I[10*INCB];
  MKL_Complex8 cBha0_I[10*INCB],cBha1_I[10*INCB],cBha2_I[10*INCB];
  MKL_Complex8                  cBla1_I[10*INCB],cBla2_I[10*INCB];
  MKL_Complex8                  cBep1_I[10*INCB],cBep2_I[10*INCB];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB;

  MKL_INT i=0,vec_len=10;

  cA[0].real=-17.0000;cA[0].imag=18.0000;
  cA[1].real=-13.1111;cA[1].imag=14.1111;
  cA[2].real=-9.2222;cA[2].imag=10.2222;
  cA[3].real=-5.3333;cA[3].imag=6.3333;
  cA[4].real=-1.4444;cA[4].imag=2.4444;
  cA[5].real=2.4444;cA[5].imag=-1.4444;
  cA[6].real=6.3333;cA[6].imag=-5.3333;
  cA[7].real=10.2222;cA[7].imag=-9.2222;
  cA[8].real=14.1111;cA[8].imag=-13.1111;
  cA[9].real=18.0000;cA[9].imag=-17.0000;
  cB[0].real=2.7336700468037275e-008;cB[0].imag=-3.1090404917222259e-008;
  cB[1].real=5.2718026213238074e-008;cB[1].imag=2.0219670204824070e-006;
  cB[2].real=-6.9031750320313866e-005;cB[2].imag=-7.0712150773033500e-005;
  cB[3].real=4.8220487634171522e-003;cB[3].imag=2.4185802612919360e-004;
  cB[4].real=-1.8084268297689621e-001;cB[4].imag=1.5145584940910339e-001;
  cB[5].real=1.4526703648953678e+000;cB[5].imag=-1.1431704521179199e+001;
  cB[6].real=3.2754686308712672e+002;cB[6].imag=4.5792468261718750e+002;
  cB[7].real=-2.6944639799017234e+004;cB[7].imag=-5.5343007812500000e+003;
  cB[8].real=1.1494097911982329e+006;cB[8].imag=-6.9640050000000000e+005;
  cB[9].real=-1.8067216284192696e+007;cB[9].imag=6.3125328000000000e+007;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcExp(vec_len,cA,cBha0);
  vcExpI(vec_len,cA_I,inca,cBha0_I,incb);

  vmcExp(vec_len,cA,cBep1,VML_EP);
  vmcExpI(vec_len,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcExp(vec_len,cA,cBep2);
  vcExpI(vec_len,cA_I,inca,cBep2_I,incb);

  vmcExp(vec_len,cA,cBla1,VML_LA);
  vmcExpI(vec_len,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcExp(vec_len,cA,cBla2);
  vcExpI(vec_len,cA_I,inca,cBla2_I,incb);

  vmcExp(vec_len,cA,cBha1,VML_HA);
  vmcExpI(vec_len,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcExp(vec_len,cA,cBha2);
  vcExpI(vec_len,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcExp and vmcExp in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcExpI and vmcExpI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcExp and vmcExp in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcExpI and vmcExpI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcExp and vmcExp in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcExpI and vmcExpI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcExp test/example program\n\n");
  printf("           Argument                           vcExp\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA[i].real,cA[i].imag,cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcExpI test/example program\n\n");
  printf("           Argument                           vcExp\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA_I[i*inca].real,cA_I[i*inca].imag,cBha0_I[i*incb].real,
        cBha0_I[i*incb].imag);
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
