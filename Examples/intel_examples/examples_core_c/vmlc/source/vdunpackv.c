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
!    vdUnpackV  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  double dA[10],dB1[10],dB2[10];

  MKL_INT i=0,vec_len=10,ia[10];

  for(i=0;i<vec_len;i++) {
    dA[i]=(double)i+1.0;
    dB1[i]=0.0;
    dB2[i]=0.0;
    ia[i]=vec_len-i-1;
  }

  vdPackV(vec_len,dA,ia,dB1);
  vdUnpackV(vec_len,dB1,dB2,ia);

  printf("vdUnpackV test/example program\n\n");
  printf("           Before packing             After packing          After Unpacking\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("% 25.14f % 25.14f % 25.14f\n",dA[i],dB1[i],dB2[i]);
  }

  return 0;
}
