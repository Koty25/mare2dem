/*******************************************************************************
* Copyright 2004-2020 Intel Corporation.
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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) PARDISO C example
*
********************************************************************************
*/
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "PARDISO" routine */
/* with low-rank update on unsymmetric linear systems */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

MKL_INT main (void)
{
    /* Matrix data. */
    MKL_INT n = 5;
    MKL_INT ia[6] = { 1, 4, 6, 9, 12, 14};
    MKL_INT ja[13] = 
    { 1, 2,    4,
      1, 2,
            3, 4, 5,
      1,    3, 4,
         2,       5
    };
    double a1[13] = 
    { 1.0,-1.0,     -3.0,
     -2.0, 5.0,
                4.0, 6.0, 4.0,
     -4.0,      2.0, 7.0,
           8.0,          -5.0
    };
    /* Matrix A2 is structurally identical to A1, but some values have changed */
    double a2[13] = 
    { 1.0,-1.0,     -3.3,
     -2.0, 5.0,
                4.0, 6.0, 4.0,
     -4.0,      2.0, 7.0,
           8.8,          -5.0
    };
    MKL_INT lowrank[5] = { 2, 1, 4, 5, 2 }; /* List of changed values from matrix A1 to A2 */
    MKL_INT mtype = 11;       /* Real unsymmetric matrix */
    /* RHS and solution vectors. */
    double b1[5] = {-3, 3, 14, 5, 3.0};
    double b2[5] = {-3.3, 3, 14, 5, 3.8};
    double x[5];
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i, j;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    char transa = 'N';
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for (i = 1; i < sizeof(lowrank) / sizeof(MKL_INT); ++i) lowrank[i] -= 1; /* low-rank update only accepts 0-based indices for list of changed values */
    for ( i = 0; i < 64; i++ ) iparm[i] = 0;

    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    iparm[23] = 10;       /* Set improved version of two-level factorization (required for low-rank update) */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 0;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a1, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d\n", iparm[18]);
/* -------------------------------------------------------------------- */
/* .. Full numerical factorization and solution for A1. */
/* -------------------------------------------------------------------- */
    phase = 23;
    printf ("\nStarting phase 23 for matrix A1 (full update) ... \n");
    fflush(0);
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a1, ia, ja, &idum, &nrhs, iparm, &msglvl, b1, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during phase 23 for matrix A1: %d", error);
        exit (2);
    }
    printf ("\nPhase 23 completed for matrix A1 (full update) ... \n");
    fflush(0);
    for ( j = 1; j <= n; j++ )
    {
        printf("x1[%d] = %1.2f\n",j-1,x[j-1]);
    }
/* -------------------------------------------------------------------- */
/* .. Low-rank update factorization and solution for A2. */
/* -------------------------------------------------------------------- */
    phase = 23;
    iparm[38] = 1; // turn on low rank update
    printf ("\nStarting phase 23 for matrix A2 (low-rank update) ... \n");
    fflush(0);
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a2, ia, ja, lowrank, &nrhs, iparm, &msglvl, b2, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during phase 23 for matrix A2: %d", error);
        exit (2);
    }
    printf ("\nPhase 23 completed for matrix A2 (low-rank update) ... \n");
    for ( j = 1; j <= n; j++ )
    {
        printf("x2[%d] = %1.2f\n",j-1,x[j-1]);
    }

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
    return 0;
}
