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
*   Intel(R) Math Kernel Library (Intel(R) MKL) Cluster Sparse Solver example
*   demonstrating the case when initial data (matrix and rhs) distributed
*   between several MPI processes, final solution is distributed between
*   MPI processes in the same way as they hold initial data.
*
********************************************************************************
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "mkl.h"
#include "mkl_cluster_sparse_solver.h"

#ifdef MKL_ILP64
#define MPI_DT MPI_LONG
#else
#define MPI_DT MPI_INT
#endif
#define MPI_REDUCE_AND_BCAST \
        MPI_Reduce(&err_mem, &error, 1, MPI_DT, MPI_SUM, 0, MPI_COMM_WORLD); \
        MPI_Bcast(&error, 1, MPI_DT, 0, MPI_COMM_WORLD);
int main(void)
{
    /* Matrix data. */
    MKL_INT n = 5;

    MKL_INT mtype = 11; /* Real unsymmetric matrix */
    MKL_INT *ia = NULL;
    MKL_INT *ja = NULL;
    double  *a = NULL;
    /* RHS and solution vectors. */
    double  *b = NULL;
    double  *x = NULL;

    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64] = { 0 };
    /* Cluster Sparse Solver control parameters. */
    MKL_INT iparm[64] = { 0 };
    MKL_INT maxfct, mnum, phase, msglvl, error, err_mem;;

    /* Auxiliary variables. */
    double  ddum; /* Double dummy   */
    MKL_INT idum; /* Integer dummy. */
    MKL_INT j;
    int     mpi_stat = 0;
    int     argc = 0;
    int     comm, rank, size;
    char**  argv;

    /* -------------------------------------------------------------------- */
    /* .. Init MPI.                                                         */
    /* -------------------------------------------------------------------- */
    mpi_stat = MPI_Init( &argc, &argv );
    mpi_stat = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	mpi_stat = MPI_Comm_size( MPI_COMM_WORLD, &size );
    comm =  MPI_Comm_c2f( MPI_COMM_WORLD );

	if( size < 2 )
	{
		printf("\nERROR: this example doesn't work on number of MPI less than 2");
        error = 1;
        goto final;
	}

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
    iparm[39] =  2; /* Input: matrix/rhs/solution are distributed between MPI processes  */
    /* If iparm[39]=2, the matrix is provided in distributed assembled matrix input          
       format. In this case, each MPI process stores only a part (or domain) of the matrix A 
       data. The bounds of the domain should be set via iparm(41) and iparm(42). Solution    
       vector is distributed between process in same manner with rhs. */		

    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 1; /* Print statistical information in file */
    error  = 0; /* Initialize error flag */
    err_mem = 0; /* Initialize error flag for memory allocation */

    /* Initialize matrix and rhs components on each process:
       In this example initial matrix is distributed between 2 processes
       so for MPI processes with rank > 1 input domains are empty */
    if ( rank == 0 )
    {
        iparm[40] = 1; /* The number of row in global matrix, rhs element and solution vector
                          that begins the input domain belonging to this MPI process */
        iparm[41] = 3; /* The number of row in global matrix, rhs element and solution vector
                          that ends the input domain belonging to this MPI process   */
        ia = (MKL_INT*) MKL_malloc (sizeof (MKL_INT) * 4, 64);
        ja = (MKL_INT*) MKL_malloc (sizeof (MKL_INT) * 6, 64);
        a = (double*) MKL_malloc (sizeof (double ) * 6, 64);
        x = (double*) MKL_malloc (sizeof (double) * 3, 64);
        b = (double*) MKL_malloc (sizeof (double) * 3, 64);
        if ( ia == NULL || ja == NULL || a == NULL || x == NULL || b == NULL )
        {
            printf ("\nERROR during memory allocation on 0 rank");
            err_mem = 1;
        }
        MPI_REDUCE_AND_BCAST;
        if (error) goto final;
        ia[0] = 1;
        ia[1] = 4;
        ia[2] = 6;
        ia[3] = 7;

        ja[0] = 1;
        ja[1] = 2;
        ja[2] = 3;
        ja[3] = 1;
        ja[4] = 2;
        ja[5] = 4;

        a[0] = 1.;
        a[1] = -1.;
        a[2] = -3.;
        a[3] = -2.;
        a[4] = 5.;
        a[5] = 6.;

        b[0] = 1.;
        b[1] = 1.;
        b[2] = 0.25;
    }
    else
    {
        if ( rank == 1 )
        {
            iparm[40] = 3; /* The number of row in global matrix, rhs element and solution vector
                              that begins the input domain belonging to this MPI process*/
            iparm[41] = 5; /* The number of row in global matrix, rhs element and solution vector
                              that ends the input domain belonging to this MPI process*/
            ia = (MKL_INT*) MKL_malloc (sizeof (MKL_INT) * 4, 64);
            ja = (MKL_INT*) MKL_malloc (sizeof (MKL_INT) * 7, 64);
            a = (double*) MKL_malloc (sizeof (double ) * 7, 64);
            x = (double*) MKL_malloc (sizeof (double) * 3, 64);
            b = (double*) MKL_malloc (sizeof (double) * 3, 64);
            if ( ia == NULL || ja == NULL || a == NULL || x == NULL || b == NULL )
            {
                err_mem = 1;
                printf ("\nERROR during memory allocation on 1 rank");
            }
            MPI_REDUCE_AND_BCAST;
            if (error) goto final;
            ia[0] = 1;
            ia[1] = 3;
            ia[2] = 6;
            ia[3] = 8;

            ja[0] = 3;
            ja[1] = 5;
            ja[2] = 1;
            ja[3] = 2;
            ja[4] = 4;
            ja[5] = 2;
            ja[6] = 5;

            a[0] = 4.;
            a[1] = 4.;
            a[2] = -4.;
            a[3] = 2.;
            a[4] = 7.;
            a[5] = 8.;
            a[6] = -5.;

            b[0] = 0.75;
            b[1] = 1.;
            b[2] = 1.;
        }
        else
        {
            MPI_REDUCE_AND_BCAST;
            /* In this example MPI processes with rank > 1 doesn't have input domain
               so iparm[40] need to be greater then iparm[41] */
            iparm[40] = 2;
            iparm[41] = 1;
        }
    }

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

    if ( rank == 0 ) printf ("\nSolving system...");
    cluster_sparse_solver ( pt, &maxfct, &mnum, &mtype, &phase,
        &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &comm, &error );
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during solution: %lli", (long long int)error);
        goto final;
    }
    /* The solution of the system is distributed between MPI processes like as input matrix
       so MPI processes with rank 0 and 1 keep only part of solution */
    if ( rank == 0 )
    {
        printf ("\nThe solution of the system is: ");
        for ( j = 0; j < 3; j++ )
        {
            printf ("\n on zero process x [%lli] = % f", (long long int)j, x[j]);
        }
        printf ("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if ( rank == 1 )
    {
        printf ("\nThe solution of the system is: ");
        for ( j = 0; j < 3; j++ )
        {
            printf ("\n on first process x [%lli] = % f", (long long int)j, x[j]);
        }
        printf ("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);

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
final:
    MPI_Barrier(MPI_COMM_WORLD);
    if ( rank < 2 )
    {
        MKL_free(ia);
        MKL_free(ja);
        MKL_free(a);
        MKL_free(x);
        MKL_free(b);
    }
    mpi_stat = MPI_Finalize();
    if ( rank == 0 ) 
    {
        if ( error != 0 )
        {
            printf("\n TEST FAILED\n");
        } else {
            printf("\n TEST PASSED\n");
        }
    }
    return error;
}
