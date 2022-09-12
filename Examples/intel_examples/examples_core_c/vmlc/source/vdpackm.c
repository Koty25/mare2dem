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
!    vdPackM  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  double dA[10],dB1[10];

  MKL_INT i=0,vec_len=10,ma[10];

  for(i=0;i<vec_len;i++) {
    dA[i]=(double)i+1.0;
    dB1[i]=0.0;
    ma[i]=i&1;
  }

  vdPackM(vec_len,dA,ma,dB1);

  printf("vdPackM test/example program\n\n");
  printf("           Before packing             After packing\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f ",dA[i]);
    printf("% 25.14f",dB1[i]);
    printf("\n");
  }

  return 0;
}
