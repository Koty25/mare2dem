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
!    vzLog10  Example Program Text
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

  zA[0].real=0.1000;zA[0].imag=10000.0000;
  zA[1].real=1111.1999;zA[1].imag=8888.9000;
  zA[2].real=2222.2999;zA[2].imag=7777.8000;
  zA[3].real=3333.3999;zA[3].imag=6666.7000;
  zA[4].real=4444.5000;zA[4].imag=5555.6000;
  zA[5].real=5555.6000;zA[5].imag=4444.5000;
  zA[6].real=6666.7000;zA[6].imag=3333.4000;
  zA[7].real=7777.8000;zA[7].imag=2222.3000;
  zA[8].real=8888.9000;zA[8].imag=1111.2000;
  zA[9].real=10000.0000;zA[9].imag=0.1000;
  zB[0].real=4.0000000000217151e+000;zB[0].imag=6.8218383397610183e-001;
  zB[1].real=3.9522152373323811e+000;zB[1].imag=6.2817728449739596e-001;
  zB[2].real=3.9078977182816068e+000;zB[2].imag=5.6132048847527494e-001;
  zB[3].real=3.8723672191923448e+000;zB[3].imag=4.8082597824966516e-001;
  zB[4].real=3.8521536559202567e+000;zB[4].imag=3.8915095567898106e-001;
  zB[5].real=3.8521536559202567e+000;zB[5].imag=2.9303722124193959e-001;
  zB[6].real=3.8723672217981222e+000;zB[6].imag=2.0136220388273199e-001;
  zB[7].real=3.9078977197566052e+000;zB[7].imag=1.2086769360797429e-001;
  zB[8].real=3.9522152379337570e+000;zB[8].imag=5.4010897234153650e-002;
  zB[9].real=4.0000000000217151e+000;zB[9].imag=4.3429448188877533e-006;

  for(i=0;i<10;i++) {
    zA_I[i*inca]=zA[i];
    zB_I[i*incb]=zB[i];
  }

  vzLog10(vec_len,zA,zBha0);
  vzLog10I(vec_len,zA_I,inca,zBha0_I,incb);

  vmzLog10(vec_len,zA,zBep1,VML_EP);
  vmzLog10I(vec_len,zA_I,inca,zBep1_I,incb,VML_EP);

  vmlSetMode(VML_EP);
  vzLog10(vec_len,zA,zBep2);
  vzLog10I(vec_len,zA_I,inca,zBep2_I,incb);

  vmzLog10(vec_len,zA,zBla1,VML_LA);
  vmzLog10I(vec_len,zA_I,inca,zBla1_I,incb,VML_LA);

  vmlSetMode(VML_LA);
  vzLog10(vec_len,zA,zBla2);
  vzLog10I(vec_len,zA_I,inca,zBla2_I,incb);

  vmzLog10(vec_len,zA,zBha1,VML_HA);
  vmzLog10I(vec_len,zA_I,inca,zBha1_I,incb,VML_HA);

  vmlSetMode(VML_HA);
  vzLog10(vec_len,zA,zBha2);
  vzLog10I(vec_len,zA_I,inca,zBha2_I,incb);

  for(i=0;i<10;i++) {
    if(zBha0[i].real!=zBha1[i].real || zBha0[i].imag!=zBha1[i].imag || zBha1[i].real!=zBha2[i].real ||
        zBha1[i].imag!=zBha2[i].imag) {
      printf("Error! Difference between vzLog10 and vmzLog10 in VML_HA mode detected.\n");
      return 1;
    }

    if(zBha0_I[i*incb].real!=zBha1_I[i*incb].real || zBha0_I[i*incb].imag!=zBha1_I[i*incb].imag ||
        zBha1_I[i*incb].real!=zBha2_I[i*incb].real || zBha1_I[i*incb].imag!=zBha2_I[i*incb].imag) {
      printf("Error! Difference between vzLog10I and vmzLog10I in VML_HA mode detected.\n");
      return 1;
    }

    if(zBla1[i].real!=zBla2[i].real || zBla1[i].imag!=zBla2[i].imag) {
      printf("Error! Difference between vzLog10 and vmzLog10 in VML_LA mode detected.\n");
      return 1;
    }

    if(zBla1_I[i*incb].real!=zBla2_I[i*incb].real || zBla1_I[i*incb].imag!=zBla2_I[i*incb].imag) {
      printf("Error! Difference between vzLog10I and vmzLog10I in VML_LA mode detected.\n");
      return 1;
    }

    if(zBep1[i].real!=zBep2[i].real || zBep1[i].imag!=zBep2[i].imag) {
      printf("Error! Difference between vzLog10 and vmzLog10 in VML_EP mode detected.\n");
      return 1;
    }

    if(zBep1_I[i*incb].real!=zBep2_I[i*incb].real || zBep1_I[i*incb].imag!=zBep2_I[i*incb].imag) {
      printf("Error! Difference between vzLog10I and vmzLog10I in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vzLog10 test/example program\n\n");
  printf("           Argument                           vzLog10\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   % .4f %+.4f*i      % .10f %+.10f*i\n",zA[i].real,zA[i].imag,zBha0[i].real,zBha0[i].imag);
    CurRMS=zrelerr(zB[i],zBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vzLog10I test/example program\n\n");
  printf("           Argument                           vzLog10\n");
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
