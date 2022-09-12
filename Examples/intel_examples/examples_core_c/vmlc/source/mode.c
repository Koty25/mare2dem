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
!    vmlSetMode/vmlGetMode  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

void PrintTextVmlMode(unsigned int );

int main()
{
  unsigned MKL_INT mode;

  printf("vmlSetMode/vmlGetMode example program\n\n");

  mode=(unsigned MKL_INT)vmlGetMode();
  printf("Default value of vmlMode: ");
  PrintTextVmlMode((unsigned int)mode);
  printf(" (0x%x)\n",(unsigned int)mode);

  mode=(unsigned MKL_INT)(VML_LA|VML_FLOAT_CONSISTENT|VML_ERRMODE_IGNORE);
  vmlSetMode(mode);
  mode=(unsigned MKL_INT)vmlGetMode();
  printf("Value of vmlMode after using vmlSetMode: ");
  PrintTextVmlMode((unsigned int)mode);
  printf(" (0x%x)\n",(unsigned int)mode);

  return 0;
}

void PrintTextVmlMode(unsigned int mode)
{
  switch(mode&VML_ACCURACY_MASK) {
    case 0x001: {
                printf("VML_LA");
                break;
    }
    case 0x002: {
                printf("VML_HA");
                break;
    }
  }
  switch(mode&VML_FPUMODE_MASK) {
    case 0x000: {
                printf(" | VML_DEFAULT_PRECISION");
                break;
    }
    case 0x010: {
                printf(" | VML_FLOAT_CONSISTENT");
                break;
    }
    case 0x020: {
                printf(" | VML_DOUBLE_CONSISTENT");
                break;
    }
    case 0x030: {
                printf(" | VML_RESTORE");
                break;
    }
  }
  if(mode&VML_ERRMODE_IGNORE) {
    printf(" | VML_ERRMODE_IGNORE");
  }
  else if (mode&VML_ERRMODE_NOERR)
  {
    printf(" | VML_ERRMODE_NOERR");
  }
  else {
    if(mode&VML_ERRMODE_ERRNO) {
      printf(" | VML_ERRMODE_ERRNO");
    }
    if(mode&VML_ERRMODE_STDERR) {
      printf(" | VML_ERRMODE_STDERR");
    }
    if(mode&VML_ERRMODE_EXCEPT) {
      printf(" | VML_ERRMODE_EXCEPT");
    }
    if(mode&VML_ERRMODE_CALLBACK) {
      printf(" | VML_ERRMODE_CALLBACK");
    }
  }
}
