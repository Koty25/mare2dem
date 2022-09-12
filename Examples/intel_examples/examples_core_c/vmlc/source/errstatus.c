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
!    vmlSetErrStatus/vmlGetErrStatus/vmlClearErrStatus example program  Example
!    Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

void PrintTextVmlErrStatus(int );

int main()
{
  MKL_INT errst;

  printf("vmlSetErrStatus/vmlGetErrStatus/vmlClearErrStatus example program\n\n");

  errst=(MKL_INT)vmlGetErrStatus();

  printf("Initial value of vmlErrStatus: ");
  PrintTextVmlErrStatus((int)errst);
  printf(" (0x%x)\n",(int)errst);

  errst=(MKL_INT)(VML_STATUS_BADMEM);
  vmlSetErrStatus(errst);
  errst=(int)vmlGetErrStatus();
  printf("Value of vmlErrStatus after using vmlSetErrStatus: ");
  PrintTextVmlErrStatus((int)errst);
  printf(" (0x%x)\n",(int)errst);

  vmlClearErrStatus();
  errst=(MKL_INT)vmlGetErrStatus();
  printf("Value of vmlErrStatus after using vmlClearErrStatus: ");
  PrintTextVmlErrStatus((int)errst);
  printf(" (0x%x)\n",(int)errst);

  return 0;
}

void PrintTextVmlErrStatus(int errst)
{
  switch(errst) {
    case VML_STATUS_OK: {
                printf("VML_STATUS_OK");
                break;
    }
    case VML_STATUS_BADSIZE: {
                printf("VML_STATUS_BADSIZE");
                break;
    }
    case VML_STATUS_BADMEM: {
                printf("VML_STATUS_BADMEM");
                break;
    }
    case VML_STATUS_ERRDOM: {
                printf("VML_STATUS_ERRDOM");
                break;
    }
    case VML_STATUS_SING: {
                printf("VML_STATUS_SING");
                break;
    }
    case VML_STATUS_OVERFLOW: {
                printf("VML_STATUS_OVERFLOW");
                break;
    }
    case VML_STATUS_UNDERFLOW: {
                printf("VML_STATUS_UNDERFLOW");
                break;
    }
  }
}
