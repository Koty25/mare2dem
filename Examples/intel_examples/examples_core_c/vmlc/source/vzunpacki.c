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
!    vzUnpackI  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int main()
{
  MKL_Complex16 zA[10],zB1[10],zB2[10];

  MKL_INT i=0,vec_len=10,incra=3;

  for(i=0;i<vec_len;i++) {
    zA[i].real=(double)i+1.0;zA[i].imag=(double)i+1.0;
    zB1[i].real=0.0;zB1[i].imag=0.0;
    zB2[i].real=0.0;zB2[i].imag=0.0;
  }

  vec_len=vec_len/incra+1;
  vzPackI(vec_len,zA,incra,zB1);
  vzUnpackI(vec_len,zB1,zB2,incra);

  printf("vzUnpackI test/example program\n\n");
  printf("     Before packing             After packing          After Unpacking\n");
  printf("===============================================================================\n");
  for(i=0;i<10;i++) {
    printf("   %.4f + %.4f*i         %.4f + %.4f*i         %.4f + %.4f*i\n",zA[i].real,zA[i].imag,zB1[i].real,zB1[i].imag,zB2[i].real,zB2[i].imag);
  }

  return 0;
}
