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
!    vcPow  Example Program Text
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

  cA[0].real=0.1000;cA[0].imag=7.0000;
  cA[1].real=0.8666;cA[1].imag=6.2333;
  cA[2].real=1.6333;cA[2].imag=5.4666;
  cA[3].real=2.4000;cA[3].imag=4.6999;
  cA[4].real=3.1666;cA[4].imag=3.9333;
  cA[5].real=3.9333;cA[5].imag=3.1666;
  cA[6].real=4.7000;cA[6].imag=2.3999;
  cA[7].real=5.4666;cA[7].imag=1.6333;
  cA[8].real=6.2333;cA[8].imag=0.8666;
  cA[9].real=7.0000;cA[9].imag=0.0999;
  cB[0].real=7.9222171168020995e-006;cB[0].imag=2.1083585901275048e-005;
  cB[1].real=6.4512072154176666e-004;cB[1].imag=9.1692323654694654e-005;
  cB[2].real=9.0503621800775043e-003;cB[2].imag=-1.2801471360156119e-002;
  cB[3].real=-1.5959425975423699e-001;cB[3].imag=-2.6566230943845598e-001;
  cB[4].real=-4.8997204132252605e+000;cB[4].imag=1.1363607778135290e+000;
  cB[5].real=4.1042341669135736e+000;cB[5].imag=6.8100665478272219e+001;
  cB[6].real=7.9823020507184594e+002;cB[6].imag=-5.7795938272522804e+001;
  cB[7].real=-2.3514348248209863e+003;cB[7].imag=-8.1467661426596324e+003;
  cB[8].real=-6.5480032831151962e+004;cB[8].imag=5.3650159355771051e+004;
  cB[9].real=7.8757355465682549e+005;cB[9].imag=2.3871476140480436e+005;

  for(i=0;i<10;i++) {
    cA_I[i*inca]=cA[i];
    cB_I[i*incb]=cB[i];
  }

  vcPow(vec_len,cA,cA,cBha0);
  vcPowI(vec_len,cA_I,inca,cA_I,inca,cBha0_I,incb);

  vmcPow(vec_len,cA,cA,cBep1,VML_EP);
  vmcPowI(vec_len,cA_I,inca,cA_I,inca,cBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vcPow(vec_len,cA,cA,cBep2);
  vcPowI(vec_len,cA_I,inca,cA_I,inca,cBep2_I,incb);

  vmcPow(vec_len,cA,cA,cBla1,VML_LA);
  vmcPowI(vec_len,cA_I,inca,cA_I,inca,cBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vcPow(vec_len,cA,cA,cBla2);
  vcPowI(vec_len,cA_I,inca,cA_I,inca,cBla2_I,incb);

  vmcPow(vec_len,cA,cA,cBha1,VML_HA);
  vmcPowI(vec_len,cA_I,inca,cA_I,inca,cBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vcPow(vec_len,cA,cA,cBha2);
  vcPowI(vec_len,cA_I,inca,cA_I,inca,cBha2_I,incb);

  for(i=0;i<10;i++) {
    if(cBha0[i].real!=cBha1[i].real || cBha0[i].imag!=cBha1[i].imag || cBha1[i].real!=cBha2[i].real ||
        cBha1[i].imag!=cBha2[i].imag) {
      printf("Error! Difference between vcPow and vmcPow in VML_HA mode detected.\n");
      return 1;
    }

    if(cBha0_I[i*incb].real!=cBha1_I[i*incb].real || cBha0_I[i*incb].imag!=cBha1_I[i*incb].imag ||
        cBha1_I[i*incb].real!=cBha2_I[i*incb].real || cBha1_I[i*incb].imag!=cBha2_I[i*incb].imag) {
      printf("Error! Difference between vcPowI and vmcPowI in VML_HA mode detected.\n");
      return 1;
    }

    if(cBla1[i].real!=cBla2[i].real || cBla1[i].imag!=cBla2[i].imag) {
      printf("Error! Difference between vcPow and vmcPow in VML_LA mode detected.\n");
      return 1;
    }

    if(cBla1_I[i*incb].real!=cBla2_I[i*incb].real || cBla1_I[i*incb].imag!=cBla2_I[i*incb].imag) {
      printf("Error! Difference between vcPowI and vmcPowI in VML_LA mode detected.\n");
      return 1;
    }

    if(cBep1[i].real!=cBep2[i].real || cBep1[i].imag!=cBep2[i].imag) {
      printf("Error! Difference between vcPow and vmcPow in VML_EP mode detected.\n");
      return 1;
    }

    if(cBep1_I[i*incb].real!=cBep2_I[i*incb].real || cBep1_I[i*incb].imag!=cBep2_I[i*incb].imag) {
      printf("Error! Difference between vcPowI and vmcPowI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vcPow test/example program\n\n");
  printf("           Arguments                           vcPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .5f %+.5f*i\n",cA[i].real,cA[i].imag,cA[i].real,cA[i].imag,
        cBha0[i].real,cBha0[i].imag);
    CurRMS=crelerr(cB[i],cBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vcPowI test/example program\n\n");
  printf("           Arguments                           vcPow\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .2f %+.2f*i   % .2f %+.2f*i      % .5f %+.5f*i\n",cA_I[i*inca].real,cA_I[i*inca].imag,
        cA_I[i*inca].real,cA_I[i*inca].imag,cBha0_I[i*incb].real,cBha0_I[i*incb].imag);
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
