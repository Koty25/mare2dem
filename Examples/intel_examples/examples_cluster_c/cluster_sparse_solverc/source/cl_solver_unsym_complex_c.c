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
*
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Cluster Sparse Solver
*             C example for complex, double precision, unsymmetric matrix
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "mkl.h"
#include "mkl_cluster_sparse_solver.h"

int main (void)
{
    /* Matrix data. */
    MKL_INT n = 8;
    MKL_INT ia[9] = { 1, 5, 8, 10, 12, 13, 16, 18, 21 };
    MKL_INT ja[20] = { 1, 3, 6, 7,         /* index of non-zeros in 1 row*/
                       2, 3, 5,            /* index of non-zeros in 2 row*/
                       3, 8,               /* index of non-zeros in 3 row*/
                       4, 7,               /* index of non-zeros in 4 row*/
                       2,                  /* index of non-zeros in 5 row*/
                       3, 6, 8,            /* index of non-zeros in 6 row*/
                       2, 7,               /* index of non-zeros in 7 row*/
                       3, 7, 8             /* index of non-zeros in 8 row*/
                    };

    MKL_Complex16 a[20] = {
        a[0].real = 7.0,   a[0].imag = 1.0,
        a[1].real = 1.0,   a[1].imag = 1.0,
        a[2].real = 2.0,   a[2].imag = 1.0,
        a[3].real = 7.0,   a[3].imag = 1.0,
        a[4].real = -4.0,  a[4].imag = 0.0,
        a[5].real = 8.0,   a[5].imag = 1.0,
        a[6].real = 2.0,   a[6].imag = 1.0,
        a[7].real = 1.0,   a[7].imag = 1.0,
        a[8].real = 5.0,   a[8].imag = 1.0,
        a[9].real = 7.0,   a[9].imag = 0.0,
        a[10].real = 9.0,  a[10].imag = 1.0,
        a[11].real = -4.0, a[11].imag = 1.0,
        a[12].real = 7.0,  a[12].imag = 1.0,
        a[13].real = 3.0,  a[13].imag = 1.0,
        a[14].real = 8.0,  a[14].imag = 0.0,
        a[15].real = 1.0,  a[15].imag = 1.0,
        a[16].real = 11.0, a[16].imag = 1.0,
        a[17].real = -3.0, a[17].imag = 1.0,
        a[18].real = 2.0,  a[18].imag = 1.0,
        a[19].real = 5.0,  a[19].imag = 0.0
    };

    MKL_INT mtype = 13; /* set matrix type to "complex unsymmetric matrix" */

    /* RHS and solution vectors. */
    MKL_Complex16 b[8], x[8], bs[8];
    double        res, res0;
    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64] = { 0 };
    /* Cluster Sparse Solver control parameters. */
    MKL_INT iparm[64] = { 0 };
    MKL_INT maxfct, mnum, phase, msglvl, error;

    /* Auxiliary variables. */
    double  ddum; /* Double dummy   */
    MKL_INT idum; /* Integer dummy. */
    MKL_INT i, j;
    int     mpi_stat = 0;
    int     argc = 0;
    int     comm, rank;
    char*   uplo;
    char**  argv;

/* -------------------------------------------------------------------- */
/* .. Init MPI.                                                         */
/* -------------------------------------------------------------------- */
    mpi_stat = MPI_Init( &argc, &argv );
    mpi_stat = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    comm =  MPI_Comm_c2f( MPI_COMM_WORLD );

/* -------------------------------------------------------------------- */
/* .. Setup Cluster Sparse Solver control parameters.                                 */
/* -------------------------------------------------------------------- */
    iparm[ 0] =  1; /* Solver default parameters overriden with provided by iparm */
    iparm[ 1] =  2; /* Use METIS for fill-in reordering */
    iparm[ 5] =  0; /* Write solution into x */
    iparm[ 7] =  2; /* Max number of iterative refinement steps */
    iparm[ 9] = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] =  1; /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] =  1; /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[26] =  1; /* Check input data for correctness */
    iparm[39] =  0; /* Input: matrix/rhs/solution stored on master */
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 1; /* Print statistical information in file */
    error  = 0; /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates   */
/* all memory that is necessary for the factorization.                  */
/* -------------------------------------------------------------------- */
    phase = 11;
    cluster_sparse_solver ( pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &comm, &error );
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during symbolic factorization: %lli", (long long int)error);
        goto final;
    }
    if ( rank == 0 ) printf ("\nReordering completed ... ");

/* -------------------------------------------------------------------- */
/* .. Numerical factorization.                                          */
/* -------------------------------------------------------------------- */
    phase = 22;
    cluster_sparse_solver ( pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &comm, &error );
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during numerical factorization: %lli", (long long int)error);
        goto final;
    }
    if ( rank == 0 ) printf ("\nFactorization completed ... ");

/* -------------------------------------------------------------------- */
/* .. Back substitution and iterative refinement.                       */
/* -------------------------------------------------------------------- */
    phase = 33;
    /* Set right hand side to one. */
    if ( rank == 0 )
    {
        for (i = 0; i < n; i++)
        {
            b[i].real = 1.;
            b[i].imag = 1.;
            x[i].real = 0.;
            x[i].imag = 0.;
        }
    }
    if ( rank == 0 ) printf ("\nSolving system...");
    cluster_sparse_solver ( pt, &maxfct, &mnum, &mtype, &phase,
                &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &comm, &error );
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during solution: %lli", (long long int)error);
        goto final;
    }
    if ( rank == 0 )
    {
        printf ("\nThe solution of the system is: ");
        for ( j = 0; j < n; j++ )
        {
            printf ("\n x [%lli] = % f % f", (long long int)j, x[j].real, x[j].imag);
        }
        /*  Compute residual */
        uplo = "non-transposed";
        mkl_zcsrgemv (uplo, &n, a, ia, ja, x, bs);
        res = 0.0;
        res0 = 0.0;
        for(i=0;i<n;i++){
            res=res + (bs[i].real-b[i].real)*(bs[i].real-b[i].real);
            res=res + (bs[i].imag-b[i].imag)*(bs[i].imag-b[i].imag);
            res0=res0 + b[i].real*b[i].real;
            res0=res0 + b[i].imag*b[i].imag;
        }
        res = sqrt ( res ) / sqrt ( res0 );
        printf ( "\nRelative residual = %e\n", res );
    }

/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    cluster_sparse_solver ( pt, &maxfct, &mnum, &mtype, &phase,
                &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &comm, &error );
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during release memory: %lli", (long long int)error);
        goto final;
    }
    /* Check residual */
	if(rank == 0)
	{
        if ( res > 1e-10 )
        {
            printf ("\nError: residual is too high!\n");
            error = 5;
            goto final;
        }
	}
final:
    if ( rank == 0 )
    {
        if ( error != 0 )
        {
            printf("\n TEST FAILED\n");
        } else {
            printf("\n TEST PASSED\n");
        }
    }
    mpi_stat = MPI_Finalize();
    return error;
}
