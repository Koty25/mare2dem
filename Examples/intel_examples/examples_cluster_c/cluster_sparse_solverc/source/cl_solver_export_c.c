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
*             C example which demonstrates usage of export functionality which
*             provides the factors L, U and permutations P and Q such that
*             P * A * Q = L * U.
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
    MKL_INT ia[9] = { 0, 4, 7, 9, 11, 14, 16, 17, 18 };
    MKL_INT ja[18] = { 0,   2,       5, 6,      /* index of non-zeros in 0 row*/
                         1, 2,    4,            /* index of non-zeros in 1 row*/
                            2,             7,   /* index of non-zeros in 2 row*/
                               3,       6,      /* index of non-zeros in 3 row*/
                                  4, 5, 6,      /* index of non-zeros in 4 row*/
                                     5,    7,   /* index of non-zeros in 5 row*/
                                        6,      /* index of non-zeros in 6 row*/
                                           7    /* index of non-zeros in 7 row*/
    };
    float a[18] = { 7.0, /*0*/ 1.0, /*0*/ /*0*/  2.0,  7.0, /*0*/
                         -4.0, 8.0, /*0*/ 2.0,  /*0*/ /*0*/ /*0*/
                               1.0, /*0*/ /*0*/ /*0*/ /*0*/ 5.0,
                                    7.0,  /*0*/ /*0*/ 9.0,  /*0*/
                                          5.0,  1.0,  5.0,  /*0*/
                                                -1.0, /*0*/ 5.0,
                                                      11.0, /*0*/
                                                            5.0
    };

    MKL_INT mtype = -2;  /* set matrix type to "real symmetric indefinite matrix" */
    MKL_INT nrhs  =  1;  /* number of right hand sides. */
    float b[8], x[8], bs[8], res, res0; /* RHS and solution vectors. */

    /* Internal solver memory pointer pt
     *       32-bit:      int pt[64] or void *pt[64];
     *       64-bit: long int pt[64] or void *pt[64]; */
    void *pt[64] = { 0 };

    /* Cluster Sparse Solver control parameters. */
    MKL_INT iparm[64] = { 0 };
    MKL_INT maxfct, mnum, phase, msglvl, error;

    /* Variables and pointers required for the export functionality */
    int ierror = 0;
    float *l_a = NULL, *u_a = NULL;
    MKL_INT *l_ia = NULL, *l_ja = NULL, *u_ia = NULL, *u_ja = NULL;
    MKL_INT *row_perm = NULL, *col_perm = NULL;
    MKL_INT l_nrows, l_nnz, u_nrows, u_nnz;

    /* Auxiliary variables. */
    float   ddum; /* float dummy   */
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
    iparm[10] =  0; /* Don't use nonsymmetric permutation and scaling MPS */
    iparm[12] =  0; /* Switch off Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[27] =  1; /* Single precision mode of Cluster Sparse Solver */
    iparm[34] =  1; /* Cluster Sparse Solver use C-style indexing for ia and ja arrays */
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
    for ( i = 0; i < n; i++ )
    {
        b[i] = 1.0;
        x[i] = 0.0;
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
            printf ( "\n x [%lli] = % f", (long long int)j, x[j] );
        }
        /*  Compute residual */
        uplo = "Upper-triangle";
        mkl_cspblas_scsrsymv ( uplo, &n, a, ia, ja, x, bs );
        res  = 0.0;
        res0 = 0.0;
        for ( j = 0; j < n; j++ )
        {
            res  += (bs[j] - b[j]) * (bs[j] - b[j]);
            res0 += b[j] * b[j];
        }
        res = sqrt ( res ) / sqrt ( res0 );
        printf ( "\nRelative residual = %e\n", res );
    }
/* -------------------------------------------------------------------- */
/* .. Calling the export functionality.                                 */
/* -------------------------------------------------------------------- */
    /* Getting the local number of rows and local number of nnz in each of the factor */
    if ( rank == 0 ) printf ("\nGetting size and numbers of nonzeros for the factors L and U...");
    cluster_sparse_solver_get_csr_size(pt, SPARSE_PTLUQT_L, &l_nrows, &l_nnz, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during getting sizes for L data: %lli", (long long int)error);
        goto final;
    }
    printf("rank %d: l_nrows = %lld l_nnz = %lld \n",
                                         rank, (long long int)l_nrows, (long long int)l_nnz);

    cluster_sparse_solver_get_csr_size(pt, SPARSE_PTLUQT_U, &u_nrows, &u_nnz, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during getting sizes for U data: %lli", (long long int)error);
        goto final;
    }
    printf("rank %d: u_nrows = %lld u_nnz = %lld \n",
                                         rank, (long long int)u_nrows, (long long int)u_nnz);

    if ( rank == 0 )
    {
        printf ("\nAllocating data for the factors L and U and permutations P and Q...");
        l_ia = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*(l_nrows+1)   , 128);
        if (l_ia == NULL) ierror = -1;
        l_ja = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*l_nnz  , 128);
        if (l_ja == NULL) ierror = -1;
        l_a  = (float*)mkl_malloc(sizeof(float)*l_nnz, 128);
        if (l_a == NULL) ierror = -1;

        u_ia = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*(u_nrows+1)   , 128);
        if (u_ia == NULL) ierror = -1;
        u_ja = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*u_nnz  , 128);
        if (u_ja == NULL) ierror = -1;
        u_a  = (float*)mkl_malloc(sizeof(float)*u_nnz, 128);
        if (u_a == NULL) ierror = -1;

        row_perm = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*(n) , 128);
        if (row_perm == NULL) ierror = -1;
        col_perm = (MKL_INT*)mkl_malloc(sizeof(MKL_INT)*(n) , 128);
        if (col_perm == NULL) ierror = -1;
    }

    MPI_Bcast(&ierror, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ( ierror != 0 )
    {
        error = ierror;
        if ( rank == 0 ) printf ("\nERROR during memory allocation for export: %lli",
                                                                           (long long int)error);
        goto final;
    }

    /* Saving the allocated pointers to csr arrays for L and U and arrays for P and Q into the solver*/
    if ( rank == 0 ) printf ("\nSaving allocated pointers inside the cluster sparse solver...");
    cluster_sparse_solver_set_csr_ptrs(pt, SPARSE_PTLUQT_L, l_ia, l_ja, l_a, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during setting ptrs to L data: %lli", (long long int)error);
        goto final;
    }
    cluster_sparse_solver_set_csr_ptrs(pt, SPARSE_PTLUQT_U, u_ia, u_ja, u_a, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during setting ptrs to U data: %lli", (long long int)error);
        goto final;
    }
    cluster_sparse_solver_set_ptr(pt, SPARSE_PTLUQT_P, row_perm, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during setting ptr to P data: %lli", (long long int)error);
        goto final;
    }
    cluster_sparse_solver_set_ptr(pt, SPARSE_PTLUQT_Q, col_perm, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during setting ptr to Q data: %lli", (long long int)error);
        goto final;
    }

    /* Calling the main export functionality which computes the values 
       for the user's provided pointers */
    if ( rank == 0 ) printf ("\nCalling the main export functionality...");
    cluster_sparse_solver_export(pt, SPARSE_PTLUQT, &comm, &error);
    if ( error != 0 )
    {
        if ( rank == 0 ) printf ("\nERROR during exporting data: %lli", (long long int)error);
        goto final;
    }

    if ( rank == 0 ) printf ("\nMain export functionality has been caled successfully...");

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
    if( rank == 0 )
    {
        if ( res > 1e-5 )
        {
            printf ("\nError: residual is too high!\n");
            error = 5;
            goto final;
        }
    }
final:
    if (rank == 0 )
    {
        mkl_free(l_ia);
        mkl_free(l_ja);
        mkl_free(u_ia);
        mkl_free(u_ja);
    }
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
