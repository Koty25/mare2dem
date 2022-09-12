/*******************************************************************************
* Copyright 2011-2020 Intel Corporation.
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
! Content:
! An example of error processing when using DFTI functions.
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mkl_dfti.h"

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

static void report_error(MKL_LONG status);

int main(void)
{
    /* Execution status */
    MKL_LONG status = 0;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    printf("Example error_processing\n");

    printf("Try to create a DFTI descriptor with misplaced arguments\n");
    status = DftiCreateDescriptor(&hand, DFTI_COMPLEX, DFTI_SINGLE, 1, 1024);
    if (status != DFTI_NO_ERROR) report_error(status);

    if (status == DFTI_NO_ERROR) {
        printf("TEST FAILED\n");
        return 1;
    } else {
        printf("TEST PASSED\n");
        return 0;
    }
}

static void report_error(MKL_LONG status)
{
    printf(" Nonzero status = "LI"\n", status);
    printf(" Check if the status indicates of error\n");
    if (!DftiErrorClass(status, DFTI_NO_ERROR)) {
        printf("   Error: %s\n", DftiErrorMessage(status));
    } else {
        printf("   Not an error: %s\n", DftiErrorMessage(status));
    }
}
