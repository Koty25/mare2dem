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
!    vdModf  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

#include "_rms.h"

#define INCA 3
#define INCB 5
#define INCC 7

int main()
{
  double dA[10],dB[10],dC[10];
  double dBha0[10],dBha1[10],dBha2[10];
  double           dBla1[10],dBla2[10];
  double           dBep1[10],dBep2[10];
  double dCha0[10],dCha1[10],dCha2[10];
  double           dCla1[10],dCla2[10];
  double           dCep1[10],dCep2[10];
  float CurRMS,MaxRMS=0.0;

  double dA_I[10*INCA],dB_I[10*INCB],dC_I[10*INCC];
  double dBha0_I[10*INCB],dBha1_I[10*INCB],dBha2_I[10*INCB];
  double                  dBla1_I[10*INCB],dBla2_I[10*INCB];
  double                  dBep1_I[10*INCB],dBep2_I[10*INCB];
  double dCha0_I[10*INCC],dCha1_I[10*INCC],dCha2_I[10*INCC];
  double                  dCla1_I[10*INCC],dCla2_I[10*INCC];
  double                  dCep1_I[10*INCC],dCep2_I[10*INCC];
  float CurRMS_I,MaxRMS_I=0.0;
  MKL_INT inca=INCA,incb=INCB,incc=INCC;

  MKL_INT i=0,vec_len=10;

  dA[0]=-2.0000;
  dA[1]=-1.5555;
  dA[2]=-1.1111;
  dA[3]=-0.6666;
  dA[4]=-0.2222;
  dA[5]=0.2222;
  dA[6]=0.6666;
  dA[7]=1.1111;
  dA[8]=1.5555;
  dA[9]=2.0000;
  dB[0]=-2.0000000000000000e+000;
  dB[1]=-1.0000000000000000e+000;
  dB[2]=-1.0000000000000000e+000;
  dB[3]=0.0000000000000000e+000;
  dB[4]=0.0000000000000000e+000;
  dB[5]=0.0000000000000000e+000;
  dB[6]=0.0000000000000000e+000;
  dB[7]=1.0000000000000000e+000;
  dB[8]=1.0000000000000000e+000;
  dB[9]=2.0000000000000000e+000;
  dC[0]=0.0000000000000000e+000;
  dC[1]=-5.5550000000000010e-001;
  dC[2]=-1.1109999999999998e-001;
  dC[3]=-6.6659999999999997e-001;
  dC[4]=-2.2220000000000001e-001;
  dC[5]=2.2220000000000001e-001;
  dC[6]=6.6659999999999997e-001;
  dC[7]=1.1109999999999998e-001;
  dC[8]=5.5550000000000010e-001;
  dC[9]=0.0000000000000000e+000;

  for(i=0;i<10;i++) {
    dA_I[i*inca]=dA[i];
    dB_I[i*incb]=dB[i];
    dC_I[i*incc]=dC[i];
  }

  vdModf(vec_len,dA,dBha0,dCha0);
  vdModfI(vec_len,dA_I,inca,dBha0_I,incb,dCha0_I,incc);

  vmdModf(vec_len,dA,dBep1,dCep1,VML_EP);
  vmdModfI(vec_len,dA_I,inca,dBep1_I,incb,dCep1_I,incc,VML_EP);

  vmlSetMode(VML_EP);
  vdModf(vec_len,dA,dBep2,dCep2);
  vdModfI(vec_len,dA_I,inca,dBep2_I,incb,dCep2_I,incc);

  vmdModf(vec_len,dA,dBla1,dCla1,VML_LA);
  vmdModfI(vec_len,dA_I,inca,dBla1_I,incb,dCla1_I,incc,VML_LA);

  vmlSetMode(VML_LA);
  vdModf(vec_len,dA,dBla2,dCla2);
  vdModfI(vec_len,dA_I,inca,dBla2_I,incb,dCla2_I,incc);

  vmdModf(vec_len,dA,dBha1,dCha1,VML_HA);
  vmdModfI(vec_len,dA_I,inca,dBha1_I,incb,dCha1_I,incc,VML_HA);

  vmlSetMode(VML_HA);
  vdModf(vec_len,dA,dBha2,dCha2);
  vdModfI(vec_len,dA_I,inca,dBha2_I,incb,dCha2_I,incc);

  for(i=0;i<10;i++) {
    if(dBha0[i]!=dBha1[i] || dCha0[i]!=dCha1[i] || dBha1[i]!=dBha2[i] || dCha1[i]!=dCha2[i]) {
      printf("Error! Difference between vdModf and vmdModf in VML_HA mode detected.\n");
      return 1;
    }

    if(dBha0_I[i*incb]!=dBha1_I[i*incb] || dCha0_I[i*incc]!=dCha1_I[i*incc] || dBha1_I[i*incb]!=dBha2_I[i*incb] ||
        dCha1_I[i*incc]!=dCha2_I[i*incc]) {
      printf("Error! Difference between vdModfI and vmdModfI in VML_HA mode detected.\n");
      return 1;
    }

    if(dBla1[i]!=dBla2[i] || dCla1[i]!=dCla2[i]) {
      printf("Error! Difference between vdModf and vmdModf in VML_LA mode detected.\n");
      return 1;
    }

    if(dBla1_I[i*incb]!=dBla2_I[i*incb] || dCla1_I[i*incc]!=dCla2_I[i*incc]) {
      printf("Error! Difference between vdModfI and vmdModfI in VML_LA mode detected.\n");
      return 1;
    }

    if(dBep1[i]!=dBep2[i] || dCep1[i]!=dCep2[i]) {
      printf("Error! Difference between vdModf and vmdModf in VML_EP mode detected.\n");
      return 1;
    }

    if(dBep1_I[i*incb]!=dBep2_I[i*incb] || dCep1_I[i*incc]!=dCep2_I[i*incc]) {
      printf("Error! Difference between vdModfI and vmdModfI in VML_EP mode detected.\n");
      return 1;
    }
  }

  printf("vdModf test/example program\n\n");
  printf("           Argument                            vdModf\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e % 25.14e\n",dA[i],dBha0[i],dCha0[i]);
    CurRMS=drelerr(dB[i],dBha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
    CurRMS=drelerr(dC[i],dCha0[i]);
    if(CurRMS>MaxRMS) MaxRMS=CurRMS;
  }
  printf("\n");
  printf("vdModfI test/example program\n\n");
  printf("           Argument                            vdModf\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14e % 25.14e\n",dA_I[i*inca],dBha0_I[i*incb],dCha0_I[i*incc]);
    CurRMS_I=drelerr(dB_I[i*incb],dBha0_I[i*incb]);
    if(CurRMS_I>MaxRMS_I) MaxRMS_I=CurRMS_I;
    CurRMS_I=drelerr(dC_I[i*incc],dCha0_I[i*incc]);
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
