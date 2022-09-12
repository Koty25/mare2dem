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
!    vsUnpackI  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  float fA[10],fB1[10],fB2[10];

  MKL_INT i=0,vec_len=10,incra=3;

  for(i=0;i<vec_len;i++) {
    fA[i]=(float)i+1.0;
    fB1[i]=0.0;
    fB2[i]=0.0;
  }

  vec_len=vec_len/incra+1;
  vsPackI(vec_len,fA,incra,fB1);
  vsUnpackI(vec_len,fB1,fB2,incra);

  printf("vsUnpackI test/example program\n\n");
  printf("           Before packing             After packing          After Unpacking\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f ",fA[i]);
    printf("% 25.14f ",fB1[i]);
    printf("% 25.14f",fB2[i]);
    printf("\n");
  }

  return 0;
}
