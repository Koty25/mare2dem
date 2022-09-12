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
!    vcUnpackM  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  MKL_Complex8 cA[10],cB1[10],cB2[10];

  MKL_INT i=0,vec_len=10,ma[10];

  for(i=0;i<vec_len;i++) {
    cA[i].real=(double)i+1.0;cA[i].imag=(double)i+1.0;
    cB1[i].real=0.0;cB1[i].imag=0.0;
    cB2[i].real=0.0;cB2[i].imag=0.0;
    ma[i]=i&1;
  }

  vcPackM(vec_len,cA,ma,cB1);
  vcUnpackM(vec_len,cB1,cB2,ma);

  printf("vcUnpackM test/example program\n\n");
  printf("     Before packing             After packing          After Unpacking\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   %.4f + %.4f*i         %.4f + %.4f*i         %.4f + %.4f*i\n",cA[i].real,cA[i].imag,cB1[i].real,cB1[i].imag,cB2[i].real,cB2[i].imag);
  }

  return 0;
}
