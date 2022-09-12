/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
!
!******************************************************************************/
#include <mkl_types.h>
#include "common_func.h"
/* 
**  Print matrix of any type
*/                          

size_t print_number(void* data, char data_type)
{ 
  switch(data_type)
  {
    case 's':
      printf("%4f\t", *((float*)data));  
      return sizeof(float);
    case 'd':
      printf("%4lf\t", *((double*)data));  
      return sizeof(double);
    case 'c':
      printf("(%4f,%4fi)\t", ((MKL_Complex8*)data)->real, ((MKL_Complex8*)data)->imag);  
      return sizeof(MKL_Complex8);
    case 'z':
      printf("(%4lf,%4lfi)\t", ((MKL_Complex16*)data)->real, ((MKL_Complex16*)data)->imag);  
      return sizeof(MKL_Complex16);
  }
  return 0;
}

void print_matrix(size_t rows, size_t cols, char data_type, void *src)
{
  size_t i=0;
  size_t j=0;
  char* p_char=((char*) src);

  for(i=0; i<rows; i++)
  {
    for(j=0; j<cols; j++)
      p_char += print_number( p_char , data_type);
    printf("\n");
  }
}
