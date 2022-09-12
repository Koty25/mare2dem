/*******************************************************************************
* Copyright 1999-2020 Intel Corporation.
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
!      MKLGetVersion example program to obtain an MKLVersion structure that
!      contains the version information.
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl_service.h"

int main(void)
{
   MKLVersion ver;
   int len=198;
   char buf[198];

   MKL_Get_Version_String(buf, len);
   printf("\n%s\n",buf);
   printf("\n");

   MKL_Get_Version(&ver);
   printf("Major version:          %d\n",ver.MajorVersion);
   printf("Minor version:          %d\n",ver.MinorVersion);
   printf("Update version:         %d\n",ver.UpdateVersion);
   printf("Product status:         %s\n",ver.ProductStatus);
   printf("Build:                  %s\n",ver.Build);
   printf("Processor optimization: %s\n",ver.Processor);
   printf("================================================================\n");
   printf("\n");
   return 0;
}
