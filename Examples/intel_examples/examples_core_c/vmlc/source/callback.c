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
!    vmlSetErrorCallBack/vmlGetErrorCallBack/vmlClearErrorCallBack  Example
!    Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vml.h"

int UserCallBack(DefVmlErrorContext* );

int main()
{
  VMLErrorCallBack errcb;
  double dbA = 0.0;
  double dbR;

  printf("Set/Get/Clear CallBack example program\n\n");

  errcb=vmlGetErrorCallBack();
  printf("Initial adress of CallBack function: 0x%p\n",errcb);

  errcb=UserCallBack;
  vmlSetErrorCallBack(errcb);
  errcb=vmlGetErrorCallBack();
  printf("Adress of CallBack function after using Set CallBack: 0x%p\n",errcb);
  printf("Test user callback on vdLn function\n");
  vdLn(1, &dbA, &dbR);


  vmlClearErrorCallBack();
  errcb=vmlGetErrorCallBack();
  printf("Adress of CallBack function after using Clear CallBack: 0x%p\n",errcb);

  return 0;
}

int UserCallBack(DefVmlErrorContext* pdefVmlErrorContext)
{
    printf("In function %s argument a[%d]=%f is wrong.\n",
      pdefVmlErrorContext->cFuncName,
      pdefVmlErrorContext->iIndex,
      pdefVmlErrorContext->dbA1);
    return 0;
}
