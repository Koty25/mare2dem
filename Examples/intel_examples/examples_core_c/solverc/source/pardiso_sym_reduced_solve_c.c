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
/* -------------------------------------------------------------------------- */
/* Example program to show the use of the partial solve for sparse right-hand */
/* sides and sparse solution in the case of symmetric linear systems.         */
/* The feature can used either a few components of the solution vector are    */
/* needed   or the user wants to reduce computation cost at solver step.      */
/* -------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mkl_pardiso.h"
#include "mkl_types.h"

MKL_INT main (void)
{
    /* Matrix data. */
    MKL_INT n = 8;
    MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19};
    MKL_INT ja[18] = 
    { 1,   3,       6, 7,
        2, 3,    5,
           3,            8,
              4,       7,
                 5, 6, 7,
                    6,   8,
                       7,
                         8
    };
    double a[18] = 
    { 7.0,      1.0,           2.0, 7.0,
          -4.0, 8.0,      2.0,
                1.0,                     5.0,
                     7.0,           9.0,
                          5.0, 1.0, 5.0,
                              -1.0,      5.0,
                                   11.0,
                                         5.0
    };
    MKL_INT mtype = -2;       /* Real symmetric matrix */
    /* RHS and solution vectors. */
    double b[8], x[8], x_reduced[8];
    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64], perm[8];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i, j;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    double eps, t;        /* Auxiliary parameters */
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off */
    /* (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 0;           /* Supress printing statistical information */
    error = 0;            /* Initialize error flag */
/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    /* Set all components except the first and last one of right hand side to zero. */
    for ( i = 1; i < n - 1; i++ )
    {
        b[i] = 0.0;
    }
    b[0] = 1.0;
    b[n - 1] = 1.0;

/* -------------------------------------------------------------------- */
/* .. Compute the solution vector of the full system in regular way */
/* -------------------------------------------------------------------- */
    phase = 13;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during the solution of the system: %d", error);
        exit (1);
    }
    printf ("\n Solution vector is computed ... ");
    printf ("\nThe solution of the system is: ");
    for ( i = 0; i < n; i++ )
    {
        printf ("\n x [%d] = % f", i, x[i]);
    }
    printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
/* -------------------------------------------------------------------- */
/*  ..  Compute the first and last components of the solution only      */
/*    by setting iparm[30]=1.                                           */
/*    To use this option the user has to define  the input vector PERM  */
/*    so that PERM[i]=1 means that the i-th component in the right-hand */
/*    side is nonzero.                                                  */
/*    In this case, PERM[i]=1 also means that  the i-th component in    */
/*    the solution vector should be computed.                           */
/* -------------------------------------------------------------------- */
    iparm[30] = 1;
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    for ( i = 1; i < n - 1; i++ )
    {
        perm[i] = 0;
    }
    perm[0] = 1;
    perm[n - 1] = 1;
    phase = 13;
    /* compute the first and last components of the solution */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x_reduced, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during computing the first and last components: %d",
                error);
        exit (1);
    }
    printf
    ("\n The first and last components of the solution vector are computed. ");
    printf
    ("\n Compare the computed components with the solution vector found in  ");
    printf ("\n the regular way: ");
    eps = 0.0;
    for ( i = n - 2; i < n; i++ )
    {
        j = perm[i] - 1;
        t = x[j] - x_reduced[j];
        eps += t * t;
        printf ("\n x [%d] = % f    %f", j, x[j], x_reduced[j]);
    }
    eps = sqrt (eps);
    printf ("\n Accuracy of computing the first and last component ");
    printf ("\n of the solution vector %f", eps);
    printf ("\n");
/* --------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* --------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

/* ---------------------------------------------------------------------- */
/*  ..Setting IPARM[30]=2 allows reducing computational cost              */
/*    at the solver step. All components of the solution vector are       */
/*    computed.                                                           */
/*    To use IPARM[30]=2, define the components of the permutation vector */
/*    PERM so that PERM[i]=1 means that the i-th component in the         */
/*    right-hand side is nonzero. Please note that if PERM[i] is not     */
/*    equal to 1,  the i-th component of the right hand side must         */
/*    be set to zero explicitly.                                          */
/* ---------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    iparm[30] = 2;
    for ( i = 1; i < n - 1; i++ )
    {
        perm[i] = 0;
        b[i] = 0.0;
    }
    perm[0] = 1;
    perm[n - 1] = 1;
    phase = 13;
    /* compute the first and last components of the solution */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x_reduced, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during the solution of the system : %d", error);
        exit (1);
    }
    printf
    ("\n Compare the computed components with the solution vector found in  ");
    printf ("\n the regular way: ");
    eps = 0.0;
    for ( i = 0; i < n; i++ )
    {
        t = x[i] - x_reduced[i];
        eps += t * t;
        printf ("\n x [%d] = % f    %f", i, x[i], x_reduced[i]);
    }
    eps = sqrt (eps);
    printf ("\n Accuracy of the solution vector found with the help  ");
    printf ("\n of reduced forward step is %f. All components of the solution ",
            eps);
    printf ("\n vector are computed. ");
    printf ("\n");
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);
    return 0;
}
