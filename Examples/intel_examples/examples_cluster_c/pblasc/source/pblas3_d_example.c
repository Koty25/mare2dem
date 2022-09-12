/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
!      Intel(R) Math Kernel Library PBLAS C example source file
!
!******************************************************************************/

/*=========================================================================
 *==== PBLAS Level 3 example program ======================================
 *=========================================================================
 *
 * The exapmle performs multiplication  of  orthonormal matrix  by  random
 * matrix with  check of result.  It also shows the use of BLACS routines,
 * work with 2D process grid and block-ciclyc distribution.
 *
 * List of PBLAS Level 3, ScaLAPACK routines demonstrated in the example:
 *     PDGEMM,
 *
 *     PDLANGE (ScaLAPACK)
 *
 * List of BLACS routines demonstrated in the example:
 *     BLACS_GET,
 *     BLACS_PINFO,
 *     BLACS_GRIDINIT,
 *     BLACS_BARRIER,
 *     BLACS_GRIDEXIT,
 *     BLACS_EXIT,
 *     IGEBS2D,
 *     IGEBR2D,
 *     DGEBS2D,
 *     DGEBR2D
 *
 *
 *
 *
 *
 * PBLAS Level 3 routine used for block-cyclic distribution of matrices:
 *     PDGEADD
 *     *******
 *     This routine is originally designed for addition of matrices
 *         C := beta * op( A ) + gamma * C
 *         op( A ) = A  if  trans = 'N',
 *         op( A ) = A' if  trans = 'T' or trans = 'C',
 *     Note that matrices may have  different distribution,  and result is
 *     distributed as matrix C. When gamma = 0 the routine just copies mat-
 *     rix A to  matrix C, if  A is contained on one process ( case mb = m,
 *     nb = n ) the routine performs block-cyclic distribution of matrix A
 *     over all process grid.
 *
 *
 * ** Description of the example **
 *
 * User  specifies  dimension of task,  size of blocks,  dimensions of the
 * process grid and treshold for residual check in the 'pblas3ex.in' file.
 * In the beginning of the program  temporary  1D process grid  is created
 * by  means  of  BLACS  routines.  Data  is  read  from  the  input  file
 * by #0 process  and sent to all other processes.  Then temporary process
 * grid is destroyed and working 2D process grid is initialized.
 * Orthonormal matrix A and random matrix B are generated on (0,0) process
 * and distributed over the process grid.
 * 
 *     /  1/q_1 ........   1/q_n-1     1/q_n  \
 *     |        .                             |
 *     |         `.           :         :     |
 *     | -1/q_1    `.         :         :     |
 *     |        .    `.       :         :     |  =  A
 *     |   0     `.    `                      |
 *     |   : `.    `.      1/q_n-1     1/q_n  |
 *     |   :   `.    `.                       |
 *     \   0 .... 0     -(n-1)/q_n-1   1/q_n  /
 *
 *     where q_i = sqrt( i^2 + i ), i=1..n-1, q_n = sqrt( n )
 *
 *     A  -  n*n real matrix, orthonormal (thus inv_A = transposed(A))
 *     B  -  random n*n real matrix
 *
 * Product C=A*B is computed by means of p?gemm,  difference  B-inv_A*C  is
 * also computed by means of p?gemm (but with transa='T'). Norm of the dif-
 * ference and norms of matrices A and B are computed using p?lange.
 * Sheme of 2D block-ciclyc distribution of n*n matrix:
 *
 *        0     1     0     1    0             (0,0)              (0,1)
 *      __________________________         ______________      ___________
 *     |1  1 |2  2 |3  3 |4  4 |5 |       |1  1 |3  3 |5 |    |2  2 |4  4 |
 * 0   |1  1 |2  2 |3  3 |4  4 |5 |       |1  1 |3  3 |5 |    |2  2 |4  4 |
 *     |_____|_____|_____|_____|__|       |_____|_____|__|    |_____|_____|
 *     |6  6 |7  7 |8  8 |9  9 |10|       |11 11|13 13|15|    |12 12|14 14|
 * 1   |6  6 |7  7 |8  8 |9  9 |10|       |11 11|13 13|15|    |12 12|14 14|
 *     |_____|_____|_____|_____|__|       |_____|_____|__|    |_____|_____|
 *     |11 11|12 12|13 13|14 14|15| --->  |21 21|23 23|25|    |22 22|24 24|
 * 0   |11 11|12 12|13 13|14 14|15|       |_____|_____|__|    |_____|_____|
 *     |_____|_____|_____|_____|__|
 *     |16 16|17 17|18 18|19 19|20|        
 * 1   |16 16|17 17|18 18|19 19|20|            (1,0)              (1,1)
 *     |_____|_____|_____|_____|__|        ______________      ___________
 * 0   |21 21|22 22|23 23|24 24|25|       |6  6 |8  8 |10|    |7  7 |9  9 |
 *     |_____|_____|_____|_____|__|       |6  6 |8  8 |10|    |7  7 |9  9 |
 *                                        |_____|_____|__|    |_____|_____|
 *                                        |16 16|18 18|20|    |17 17|19 19|
 *                                        |16 16|18 18|20|    |17 17|19 19|
 *                                        |_____|_____|__|    |_____|_____|
 *
 *========================================================================*/

/* Header files*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_scalapack.h>
#include "mkl.h"
#include "pblas_examples.h"

/* Parameters */
const double  zero = 0.0E+0, one = 1.0E+0, two = 2.0E+0, negone = -1.0E+0;
const MKL_INT i_zero = 0, i_one = 1, i_four = 4, i_negone = -1;
const char    trans = 'N';


/*==== MAIN FUNCTION =================================================*/
int main( int argc, char *argv[] ){

/*  ==== Declarations =================================================== */

/*  File variables */
    FILE    *fin;

/*  Matrix descriptors */
    MDESC   descA, descB, descC, descA_local, descB_local;

/*  Local scalars */
    MKL_INT iam, nprocs, ictxt, myrow, mycol, nprow, npcol;
    MKL_INT n, nb, mp, nq, lld, lld_local;
    MKL_INT i, j, info;
    int     n_int, nb_int, nprow_int, npcol_int;
    double  thresh, diffnorm, anorm, bnorm, residual, eps;

/*  Local arrays */
    double  *A_local, *B_local, *A, *B, *C, *work;
    MKL_INT iwork[ 4 ];


/*  ==== Executable statements ========================================== */

/*  Get information about how many processes are used for program execution
    and number of current process */
    blacs_pinfo_( &iam, &nprocs );

/*  Init temporary 1D process grid */
    blacs_get_( &i_negone, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "C", &nprocs, &i_one );

/*  Open input file */
    if ( iam == 0 ) {
        fin = fopen( "../in/pblas3ex.in", "r" );
        if ( fin == NULL ) {
            printf( "Error while open input file." );
            return 2;
        }
    }

/*  Read data and send it to all processes */
    if ( iam == 0 ) {

/*      Read parameters */
        fscanf( fin, "%d n, dimension of vectors, 0 < n <= 1000 ", &n_int );
        fscanf( fin, "%d nb, size of blocks, must be > 0 ", &nb_int );
        fscanf( fin, "%d p, number of rows in the process grid, must be > 0", &nprow_int );
        fscanf( fin, "%d q, number of columns in the process grid, must be > 0, p*q = number of processes", &npcol_int );
        fscanf( fin, "%lf threshold for residual check (to switch off check set it < 0.0) ", &thresh );
        fclose( fin );

        n = (MKL_INT) n_int;
        nb = (MKL_INT) nb_int;
        nprow = (MKL_INT) nprow_int;
        npcol = (MKL_INT) npcol_int;

/*      Check if all parameters are correct */
        if( ( n<=0 )||( n>1000 )||( nb<=0 )||( nprow<=0 )||( npcol<=0 )||( nprow*npcol != nprocs ) ) {
            printf( "One or several input parameters has incorrect value. Limitations:\n" );
            printf( "1000 >= n > 0, nb > 0, p > 0, q > 0 - integer\n" );
            printf( "p*q = number of processes\n" );
            printf( "threshold - double (set to negative to swicth off check)\n");
            return 2;
        }

/*      Pack data into array and send it to other processes */
        iwork[ 0 ] = n;
        iwork[ 1 ] = nb;
        iwork[ 2 ] = nprow;
        iwork[ 3 ] = npcol;
        igebs2d_( &ictxt, "All", " ", &i_four, &i_one, iwork, &i_four );
        dgebs2d_( &ictxt, "All", " ", &i_one, &i_one, &thresh, &i_one );
    } else {

/*      Recieve and unpack data */
        igebr2d_( &ictxt, "All", " ", &i_four, &i_one, iwork, &i_four, &i_zero, &i_zero );
        dgebr2d_( &ictxt, "All", " ", &i_one, &i_one, &thresh, &i_one, &i_zero, &i_zero );
        n = iwork[ 0 ];
        nb = iwork[ 1 ];
        nprow = iwork[ 2 ];
        npcol = iwork[ 3 ];
    }

/*  Destroy temporary process grid */
    blacs_gridexit_( &ictxt );

/*  Init workind 2D process grid */
    blacs_get_( &i_negone, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "R", &nprow, &npcol );
    blacs_gridinfo_( &ictxt, &nprow, &npcol, &myrow, &mycol );

/*  Create on process 0 two matrices: A - orthonormal, B -random */
    if ( ( myrow == 0 ) && ( mycol == 0 ) ){

/*      Allocate arrays */
        A_local = (double*) mkl_calloc(n*n, sizeof( double ), 64);
        B_local = (double*) mkl_calloc(n*n, sizeof( double ), 64);

/*      Set arrays */
        for ( i=0; i<n; i++ ){
            for ( j=0; j<n; j++ ){
                B_local[ i+n*j ] = one*rand()/RAND_MAX;
            }
            B_local[ i+n*i ] += two;
        }
        for ( j=0; j<n; j++ ){
            for ( i=0; i<n; i++ ){
                if ( j < n-1 ){
                    if ( i <= j ){
                        A_local[ i+n*j ] = one / sqrt( ( double )( (j+1)*(j+2) ) );
                    } else if ( i == j+1 ) {
                        A_local[ i+n*j ] = -one / sqrt( one + one/( double )(j+1) );
                    } else {
                        A_local[ i+n*j ] = zero;
                    }
                } else {
                    A_local[ i+n*(n-1) ] = one / sqrt( ( double )n );
                }
            }
        }

/*      Print information of task */
        printf( "=== START OF EXAMPLE ===================\n" );
        printf( "Matrix-matrix multiplication: A*B = C\n\n" );
        printf( "/  1/q_1 ........   1/q_n-1     1/q_n  \\ \n" );
        printf( "|        .                             | \n" );
        printf( "|         `.           :         :     | \n" );
        printf( "| -1/q_1    `.         :         :     | \n" );
        printf( "|        .    `.       :         :     |  =  A \n" );
        printf( "|   0     `.    `                      | \n" );
        printf( "|   : `.    `.      1/q_n-1     1/q_n  | \n" );
        printf( "|   :   `.    `.                       | \n" );
        printf( "\\   0 .... 0     -(n-1)/q_n-1   1/q_n  / \n\n" );
        printf( "q_i = sqrt( i^2 + i ), i=1..n-1, q_n = sqrt( n )\n\n" );
        printf( "A  -  n*n real matrix (orthonormal) \n" );
        printf( "B  -  random n*n real matrix\n\n" );
        printf( "n = %d, nb = %d; %dx%d - process grid\n\n", n, nb, nprow, npcol );
        printf( "=== PROGRESS ===========================\n" );
    } else {

/*      Other processes don't contain parts of initial arrays */
        A_local = NULL;
        B_local = NULL;
    }

/*  Compute precise length of local pieces and allocate array on
    each process for parts of distributed vectors */
    mp = numroc_( &n, &nb, &myrow, &i_zero, &nprow );
    nq = numroc_( &n, &nb, &mycol, &i_zero, &npcol );
    A = (double*) mkl_calloc(mp*nq, sizeof( double ), 64);
    B = (double*) mkl_calloc(mp*nq, sizeof( double ), 64);
    C = (double*) mkl_calloc(mp*nq, sizeof( double ), 64);

/*  Compute leading dimensions */
    lld_local = MAX( numroc_( &n, &n, &myrow, &i_zero, &nprow ), 1 );
    lld = MAX( mp, 1 );

/*  Initialize descriptors for initial arrays located on 0 process */
    descinit_( descA_local, &n, &n, &n, &n, &i_zero, &i_zero, &ictxt, &lld_local, &info );
    descinit_( descB_local, &n, &n, &n, &n, &i_zero, &i_zero, &ictxt, &lld_local, &info );

/*  Initialize descriptors for distributed arrays */
    descinit_( descA, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descB, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descC, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info );

/*  Distribute matrices from 0 process over process grid */
    pdgeadd_( &trans, &n, &n, &one, A_local, &i_one, &i_one, descA_local, &zero, A, &i_one, &i_one, descA );
    pdgeadd_( &trans, &n, &n, &one, B_local, &i_one, &i_one, descB_local, &zero, B, &i_one, &i_one, descB );
    if( iam == 0 ){ printf( ".. Arrays are distributed ( p?geadd ) ..\n" ); }

/*  Destroy arrays on 0 process - they are not necessary anymore */
    if( ( myrow == 0 ) && ( mycol == 0 ) ){
        mkl_free( A_local );
        mkl_free( B_local );
    }

/*  Compute norm of A and B */
    work = (double*) mkl_calloc(mp, sizeof( double ), 64);
    anorm = pdlange_( "I", &n, &n, A, &i_one, &i_one, descA, work );
    bnorm = pdlange_( "I", &n, &n, B, &i_one, &i_one, descB, work );
    if( iam == 0 ){ printf( ".. Norms of A and B are computed ( p?lange ) ..\n" ); }

/*  Compute product C = A*B */
    pdgemm_( "N", "N", &n, &n, &n, &one, A, &i_one, &i_one, descA, B, &i_one, &i_one, descB,
             &zero, C, &i_one, &i_one, descC );
    if( iam == 0 ){ printf( ".. Multiplication A*B=C is done ( p?gemm ) ..\n" ); }

/*  Compute difference  B - inv_A*C (inv_A = transpose(A) because A is orthonormal) */
    pdgemm_( "T", "N", &n, &n, &n, &one, A, &i_one, &i_one, descA, C, &i_one, &i_one, descC,
             &negone, B, &i_one, &i_one, descB );
    if( iam == 0 ){ printf( ".. Difference is computed ( p?gemm ) ..\n" ); }

/*  Compute norm of B - inv_A*C (which is contained in B) */
    diffnorm = pdlange_( "I", &n, &n, B, &i_one, &i_one, descB, work );
    mkl_free( work );
    if( iam == 0 ){ printf( ".. Norms of the difference B-inv_A*C is computed ( p?lange ) ..\n" ); }

/*  Print results */
    if( iam == 0 ){
        printf( ".. Solutions are compared ..\n" );
        printf( "== Results ==\n" );
        printf( "||A|| = %03.11f\n", anorm );
        printf( "||B|| = %03.11f\n", bnorm );
        printf( "=== END OF EXAMPLE =====================\n" );
    }

/*  Compute machine epsilon */
    eps = pdlamch_( &ictxt, "e" );

/*  Compute residual */
    residual = diffnorm /( two*anorm*bnorm*eps );

/*  Destroy arrays */
    mkl_free( A );
    mkl_free( B );
    mkl_free( C );

/*  Destroy process grid */    
    blacs_gridexit_( &ictxt );
    blacs_exit_( &i_zero );
    
/*  Check if residual passed or failed the threshold */
    if ( ( iam == 0 ) && ( thresh >= zero ) && !( residual <= thresh ) ){
        printf( "FAILED. Residual = %05.16f\n", residual );
        return 1;
    } else {
        return 0;
    }

/*========================================================================
  ====== End of PBLAS Level 3 example ====================================
  ======================================================================*/
}


/* Stub for proper work on Windows */
#ifdef _WIN_
void PXERBLA( MKL_INT *ICTXT, char *SRNAME, MKL_INT *INFO, MKL_INT lenSRNAME ) {
    MKL_INT myrow, mycol, nprow, npcol, i;
    int     INFO_int, myrow_int, mycol_int;

    Cblacs_gridinfo( *ICTXT, &nprow, &npcol, &myrow, &mycol );
    
    INFO_int = (int) *INFO;
    myrow_int = (int) myrow;
    mycol_int = (int) mycol;
    printf( "{%5i,%5i}:  On entry to ", myrow_int, mycol_int );
    for( i = 0; i < lenSRNAME; i++ ) printf( "%c", SRNAME[ i ] );
    printf( " parameter number %4i had an illegal value\n", INFO_int );

    return;
}
#endif
