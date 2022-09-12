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
!    vcPackV  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  MKL_Complex8 cA[10],cB1[10];

  MKL_INT i=0,vec_len=10,ia[10];

  for(i=0;i<vec_len;i++) {
    cA[i].real=(float)i+1.0;cA[i].imag=(float)i+1.0;
    cB1[i].real=0.0;cB1[i].imag=0.0;
    ia[i]=vec_len-i-1;
  }

  vcPackV(vec_len,cA,ia,cB1);

  printf("vcPackV test/example program\n\n");
  printf("           Before packing             After packing\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("         %.4f + %.4f*i             %.4f + %.4f*i\n",cA[i].real,cA[i].imag,cB1[i].real,cB1[i].imag);
  }

  return 0;
}
