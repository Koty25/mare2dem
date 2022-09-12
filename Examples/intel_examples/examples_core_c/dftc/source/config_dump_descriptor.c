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
! Example of printing DFTI descriptor's configuration
!
!****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl_dfti.h"

static void dump_descriptor(DFTI_DESCRIPTOR_HANDLE hand);

/* Define the format to printf MKL_LONG values */
#if !defined(MKL_ILP64)
#define LI "%li"
#else
#define LI "%lli"
#endif

int main(void)
{
    /* Sizes of a 4D FFT */
    MKL_LONG N[4] = {10, 20, 30, 40};

    /* Execution status */
    MKL_LONG status = 0;

    DFTI_DESCRIPTOR_HANDLE hand = NULL;

    printf("Create a DFTI descriptor\n");
    status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 4, N);
    if (status != DFTI_NO_ERROR) goto failed;

    printf("Dump the descriptor\n");
    dump_descriptor(hand);

 cleanup:

    printf("Free DFTI descriptor\n");
    DftiFreeDescriptor(&hand);

    printf("TEST %s\n", (status == 0) ? "PASSED" : "FAILED");
    return status;

 failed:
    printf(" ERROR, status = "LI"\n", status);
    status = 1;
    goto cleanup;
}

/* Maximum supported rank of multidimensional FFTs */
#define MAX_RANK 7

static void dump_descriptor(DFTI_DESCRIPTOR_HANDLE hand)
{
    /* Execution status */
    MKL_LONG status = 0;

    char version[DFTI_VERSION_LENGTH];

    enum DFTI_CONFIG_VALUE placement, precision, domain, storage, packfmt,
                           wspace, cmtstatus;
    MKL_LONG rank, lengths[MAX_RANK];
    double fwd_scale, bwd_scale;
    MKL_LONG nut, is[1+MAX_RANK], os[1+MAX_RANK], ntr, idist, odist, tlimit;

    DftiGetValue(0, DFTI_VERSION, version);
    printf("%s\n", version);

    printf("  PRECISION = "); fflush(0);
    status = DftiGetValue(hand, DFTI_PRECISION, &precision);
    if (status != DFTI_NO_ERROR) goto failed;
    if      (precision == DFTI_SINGLE) printf("DFTI_SINGLE\n");
    else if (precision == DFTI_DOUBLE) printf("DFTI_DOUBLE\n");
    else { printf("unknown (%i)\n", precision); goto failed; }

    printf("  FORWARD_DOMAIN = "); fflush(0);
    status = DftiGetValue(hand, DFTI_FORWARD_DOMAIN, &domain);
    if (status != DFTI_NO_ERROR) goto failed;
    if      (domain == DFTI_COMPLEX) printf("DFTI_COMPLEX\n");
    else if (domain == DFTI_REAL)    printf("DFTI_REAL\n");
    else { printf("unknown (%i)\n", domain); goto failed; }

    printf("  DIMENSION = "); fflush(0);
    status = DftiGetValue(hand, DFTI_DIMENSION, &rank);
    if (status != DFTI_NO_ERROR) goto failed;
    printf(LI"\n", rank);

    printf("  LENGTHS = "); fflush(0);
    status = DftiGetValue(hand, DFTI_LENGTHS, lengths);
    if (status != DFTI_NO_ERROR) goto failed;

    {
        int r = 0;
        printf(LI, lengths[0]);
        for (r = 1; r < rank; ++r) printf(", "LI, lengths[r]);
        printf("\n");
    }

    printf("  PLACEMENT = "); fflush(0);
    status = DftiGetValue(hand, DFTI_PLACEMENT, &placement);
    if (status != DFTI_NO_ERROR) goto failed;
    if      (placement == DFTI_NOT_INPLACE) printf("DFTI_NOT_INPLACE\n");
    else if (placement == DFTI_INPLACE)     printf("DFTI_INPLACE\n");
    else { printf("unknown (%i)\n", placement); goto failed; }

    printf("  F/B SCALES = "); fflush(0);
    if (precision == DFTI_DOUBLE) {
        status = DftiGetValue(hand, DFTI_FORWARD_SCALE, &fwd_scale);
        if (status != DFTI_NO_ERROR) goto failed;
        status = DftiGetValue(hand, DFTI_BACKWARD_SCALE, &bwd_scale);
        if (status != DFTI_NO_ERROR) goto failed;
    } else {
        float fs, bs;
        status = DftiGetValue(hand, DFTI_FORWARD_SCALE, &fs);
        if (status != DFTI_NO_ERROR) goto failed;
        status = DftiGetValue(hand, DFTI_BACKWARD_SCALE, &bs);
        if (status != DFTI_NO_ERROR) goto failed;
        fwd_scale = (double)fs;
        bwd_scale = (double)bs;
    }
    printf(" %lg, %lg\n", fwd_scale, bwd_scale);

    printf("  NO OF USER THREADS = "); fflush(0);
    status = DftiGetValue(hand, DFTI_NUMBER_OF_USER_THREADS, &nut);
    if (status != DFTI_NO_ERROR) goto failed;
    printf(LI"\n", nut);

    printf("  INPUT  STRIDES = "); fflush(0);
    status = DftiGetValue(hand, DFTI_INPUT_STRIDES, is);
    if (status != DFTI_NO_ERROR) goto failed;

    {
        int r = 0;
        printf(LI, is[0]);
        for (r = 1; r <= rank; ++r) printf(", "LI, is[r]);
        printf("\n");
    }

    printf("  OUTPUT STRIDES = "); fflush(0);
    status = DftiGetValue(hand, DFTI_OUTPUT_STRIDES, os);
    if (status != DFTI_NO_ERROR) goto failed;

    {
        int r = 0;
        printf(LI, os[0]);
        for (r = 1; r <= rank; ++r) printf(", "LI, os[r]);
        printf("\n");
    }

    printf("  NO OF TRANSFORMS = "); fflush(0);
    status = DftiGetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, &ntr);
    if (status != DFTI_NO_ERROR) goto failed;
    printf(LI"\n", ntr);

    printf("  I/O DISTANCES = "); fflush(0);
    status = DftiGetValue(hand, DFTI_INPUT_DISTANCE, &idist);
    if (status != DFTI_NO_ERROR) goto failed;
    status = DftiGetValue(hand, DFTI_OUTPUT_DISTANCE, &odist);
    if (status != DFTI_NO_ERROR) goto failed;
    printf(LI", "LI"\n", idist, odist);

    if (domain == DFTI_COMPLEX) {
        printf("  COMPLEX STORAGE = "); fflush(0);
        status = DftiGetValue(hand, DFTI_COMPLEX_STORAGE, &storage);
        if (status != DFTI_NO_ERROR) goto failed;
        if (storage == DFTI_COMPLEX_COMPLEX) printf("DFTI_COMPLEX_COMPLEX\n");
        else if (storage == DFTI_REAL_REAL)  printf("DFTI_REAL_REAL\n");
        else { printf("wrong (%i)\n", storage); goto failed; }
    } else {
        printf("  CONJUGATE EVEN STORAGE = "); fflush(0);
        status = DftiGetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, &storage);
        if (status != DFTI_NO_ERROR) goto failed;
        if (storage == DFTI_COMPLEX_COMPLEX)   printf("DFTI_COMPLEX_COMPLEX\n");
        else if (storage == DFTI_COMPLEX_REAL) printf("DFTI_COMPLEX_REAL\n");
        else { printf("wrong (%i)\n", storage); goto failed; }
        if (storage == DFTI_COMPLEX_REAL) {
            printf("     PACKED FORMAT = "); fflush(0);
            status = DftiGetValue(hand, DFTI_PACKED_FORMAT, &packfmt);
            if (status != DFTI_NO_ERROR) goto failed;
            if      (packfmt == DFTI_CCS_FORMAT) printf("DFTI_CCS_FORMAT\n");
            else if (packfmt == DFTI_PACK_FORMAT) printf("DFTI_PACK_FORMAT\n");
            else if (packfmt == DFTI_PERM_FORMAT) printf("DFTI_PERM_FORMAT\n");
            else { printf("wrong (%i)\n", packfmt); goto failed; }
        }
    }

    printf("  WORKSPACE = "); fflush(0);
    status = DftiGetValue(hand, DFTI_WORKSPACE, &wspace);
    if (status != DFTI_NO_ERROR) goto failed;
    if      (wspace == DFTI_ALLOW) printf("DFTI_ALLOW\n");
    else if (wspace == DFTI_AVOID) printf("DFTI_AVOID\n");
    else if (wspace == DFTI_NONE)  printf("DFTI_NONE\n");
    else { printf("wrong (%i)\n", wspace); goto failed; }

    printf("  COMMIT STATUS = "); fflush(0);
    status = DftiGetValue(hand, DFTI_COMMIT_STATUS, &cmtstatus);
    if (status != DFTI_NO_ERROR) goto failed;
    if      (cmtstatus == DFTI_COMMITTED) printf("DFTI_COMMITTED\n");
    else if (cmtstatus == DFTI_UNCOMMITTED) printf("DFTI_UNCOMMITTED\n");
    else { printf("wrong (%i)\n", cmtstatus); goto failed; }

    printf("  THREAD LIMIT = "); fflush(0);
    status = DftiGetValue(hand, DFTI_THREAD_LIMIT, &tlimit);
    if (status != DFTI_NO_ERROR) goto failed;
    printf(LI"\n", tlimit); fflush(0);

    return;

 failed:
    printf("Error, status = "LI"\n", status);
    exit(1);
}
