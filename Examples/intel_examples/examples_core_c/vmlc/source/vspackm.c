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
!    vsPackM  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  float fA[10],fB1[10];

  MKL_INT i=0,vec_len=10,ma[10];

  for(i=0;i<vec_len;i++) {
    fA[i]=(float)i+1.0;
    fB1[i]=0.0;
    ma[i]=i&1;
  }

  vsPackM(vec_len,fA,ma,fB1);

  printf("vsPackM test/example program\n\n");
  printf("           Before packing             After packing\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f ",fA[i]);
    printf("% 25.14f",fB1[i]);
    printf("\n");
  }

  return 0;
}
