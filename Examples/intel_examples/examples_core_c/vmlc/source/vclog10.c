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
!    vcLog10  Example Program Text
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

  cA[0].real=0.1000;cA[0].imag=10000.0000;
  cA[1].real=1111.2000;cA[1].imag=8888.9004;
  cA[2].real=2222.2998;cA[2].imag=7777.7998;
  cA[3].real=3333.3999;cA[3].imag=6666.7002;
  cA[4].real=4444.5000;cA[4].imag=5555.6001;
  cA[5].real=5555.6001;cA[5].imag=4444.5000;
  cA[6].real=6666.7002;cA[6].imag=3333.3999;
  cA[7].real=7777.7998;cA[7].imag=2222.3000;
  cA[8].real=8888.9004;cA[8].imag=1111.2000;
  cA[9].real=10000.0000;cA[9].imag=0.1000;
  cB[0].real=4.0000000000217151e+000;cB[0].imag=6.8218386173248291e-001;
  cB[1].real=3.9522152564316348e+000;cB[1].imag=6.2817728519439697e-001;
  cB[2].real=3.9078977067930758e+000;cB[2].imag=5.6132048368453979e-001;
  cB[3].real=3.8723672294320828e+000;cB[3].imag=4.8082599043846130e-001;
  cB[4].real=3.8521536605751314e+000;cB[4].imag=3.8915094733238220e-001;
  cB[5].real=3.8521536605751314e+000;cB[5].imag=2.9303723573684692e-001;
  cB[6].real=3.8723672294320828e+000;cB[6].imag=2.0136220753192902e-001;
  cB[7].real=3.9078977103941459e+000;cB[7].imag=1.2086769193410873e-001;
  cB[8].real=3.9522152564316348e+000;cB[8].imag=5.4010897874832153e-002;
  cB[9].real=4.0000000000217151e+000;cB[9].imag=4.3429449760878924e-006;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcLog10(vec_len,cA,cBha0);
  vcLog10I(vec_len,cA_I,inca,cBha0_I,incb);

  vmcLog10(vec_len,cA,cBep1,VML_EP);
  vmcLog10I(vec_len,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcLog10(vec_len,cA,cBep2);
  vcLog10I(vec_len,cA_I,inca,cBep2_I,incb);

  vmcLog10(vec_len,cA,cBla1,VML_LA);
  vmcLog10I(vec_len,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcLog10(vec_len,cA,cBla2);
  vcLog10I(vec_len,cA_I,inca,cBla2_I,incb);

  vmcLog10(vec_len,cA,cBha1,VML_HA);
  vmcLog10I(vec_len,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcLog10(vec_len,cA,cBha2);
  vcLog10I(vec_len,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcLog10 and vmcLog10 in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcLog10I and vmcLog10I in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcLog10 and vmcLog10 in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcLog10I and vmcLog10I in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcLog10 and vmcLog10 in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcLog10I and vmcLog10I in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcLog10 test/example program\n\n");
  printf("           Argument                           vcLog10\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",cA[i].real,cA[i].imag,cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcLog10I test/example program\n\n");
  printf("           Argument                           vcLog10\n");
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
