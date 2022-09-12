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
 *==== PBLAS Level 2 example program ======================================
 *=========================================================================
 *
 * The  example solves system  of linear algebraic equations with triangle
 * matrix  which can be simply inversed. It shows the use of PBLAS Level 2
 * routines,  BLACS routines,  work with  2D process grid and block-ciclyc
 * distribution.
 *
 * List of PBLAS Level 2, ScaLAPACK routines used in the example:
 *     PSTRSV,
 *     PSTRMV,
 *
 *     PSCOPY  (PBLAS Level 1),
 *     PSNRM2  (PBLAS Level 1),
 *     PSLANGE (ScaLAPACK)
 *
 * List of BLACS routines used in the example:
 *     BLACS_GET,
 *     BLACS_PINFO,
 *     BLACS_GRIDINIT,
 *     BLACS_BARRIER,
 *     BLACS_GRIDEXIT,
 *     BLACS_EXIT,
 *     IGEBS2D,
 *     IGEBR2D,
 *     SGEBS2D,
 *     SGEBR2D
 *
 * PBLAS Level 3 routine used for block-cyclic distribution of matrices:
 *     PSGEADD
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
 * process grid and treshold for residual check in the 'pblas2ex.in' file.
 * In the beginning of the program  temporary  1D process grid  is created
 * by  means  of  BLACS  routines.  Data  is  read  from  the  input  file
 * by #0 process  and sent to all other processes.  Then temporary process
 * grid is destroyed and working 2D process grid is initialized.
 * Matrix A, inversed matrix A  and random vector b of right-hand side are
 * generated on (0,0) process and distributed over the process grid.
 *
 *         / 1 ......... 1 \
 *         | 0 .         : |
 *         | : .`.       : |
 *     A = | :  `.`.     : |
 *         | :    `.`.   : |
 *         | :      `.`. : |
 *         \ 0 ....... 0 1 /
 *
 *         / 1 -1 0 .... 0 \
 *         | 0 . .  .    : |
 *         | : .`.`. `.  : |
 * inv_A = | :  `.`.`. ` : |
 *         | :    `.`.`. 0 |
 *         | :      `.`.-1 |
 *         \ 0 ....... 0 1 /
 *
 * System A*x=b and right-hand side b is solved by means of p?trsv. Then the
 * system is solved using another way -  x'=inv_A*b  -  and both results are
 * compared.
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
 * Distribution of right-hand side vector:
 *
 *       0          (0,0)   (0,1)
 *      __           __
 *     |26|         |26|     null
 * 0   |26|         |26|
 *     |__|         |__|
 *     |27|         |28|
 * 1   |27|         |28|
 *     |__|         |__|
 *     |28|  --->   |30|
 * 0   |28|         |__|
 *     |__|
 *     |29|
 * 1   |29|         (1,0)   (1,1)
 *     |__|          __
 * 0   |30|         |27|     null
 *     |__|         |27|
 *                  |__|
 *                  |29|
 *                  |29|
 *                  |__|
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
const float   zero = 0.0e+0, one = 1.0e+0, two = 2.0e+0;
const MKL_INT i_zero = 0, i_one = 1, i_four = 4, i_negone = -1;
const char    trans = 'N';

/*==== MAIN FUNCTION =================================================*/
int main( int argc, char *argv[] ){

/*  ==== Declarations =================================================== */

/*  File variables */
    FILE    *fin;

/*  Matrix descriptors */
    MDESC   descA, descinvA, descb, descx, descA_local, descinvA_local, descb_local;

/*  Local scalars */
    MKL_INT iam, nprocs, ictxt, myrow, mycol, nprow, npcol;
    MKL_INT n, nb, mp, nq, lld, lld_local, nq_rhs;
    MKL_INT i, j, info;
    int     n_int, nb_int, nprow_int, npcol_int;
    float   thresh, diffnorm, anorm, residual, eps;

/*  Local arrays */
    float   *A_local, *inv_A_local, *b_local, *A, *inv_A, *x, *b, *work;
    MKL_INT iwork[ 4 ];

/*  ==== Executable statements ========================================== */

/*  Get information about how many processes areis used for program execution
    and number of current process */
    blacs_pinfo_( &iam, &nprocs );

/*  Init temporary 1D process grid */
    blacs_get_( &i_negone, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "C", &nprocs, &i_one );

/*  Open input file */
    if ( iam == 0 ) {
        fin = fopen( "../in/pblas2ex.in", "r" );
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
        fscanf( fin, "%f threshold for residual check (to switch off check set it < 0.0) ", &thresh );
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
            printf( "threshold - float (set to negative to swicth off check)\n");
            return 2;
        }

/*      Pack data into array and send it to other processes */
        iwork[ 0 ] = n;
        iwork[ 1 ] = nb;
        iwork[ 2 ] = nprow;
        iwork[ 3 ] = npcol;
        igebs2d_( &ictxt, "All", " ", &i_four, &i_one, iwork, &i_four );
        sgebs2d_( &ictxt, "All", " ", &i_one, &i_one, &thresh, &i_one );
    } else {

/*      Recieve and unpack data */
        igebr2d_( &ictxt, "All", " ", &i_four, &i_one, iwork, &i_four, &i_zero, &i_zero );
        sgebr2d_( &ictxt, "All", " ", &i_one, &i_one, &thresh, &i_one, &i_zero, &i_zero );
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

/*  Create on process 0 vector of the right-hand side (random),
    matrix A and inversed matrix A */
    if ( ( myrow == 0 ) && ( mycol == 0 ) ){

/*      Allocate arrays */
        A_local = (float*) mkl_calloc(n*n, sizeof( float ), 64);
        b_local = (float*) mkl_calloc(n, sizeof( float ), 64);
        inv_A_local = (float*) mkl_calloc(n*n, sizeof( float ), 64);

/*      Set arrays */
        for( i=0; i<n; i++ ){
            b_local[ i ] = one*rand()/RAND_MAX;
        }
        slaset_( "U", &n, &n, &one, &one, A_local, &n );
        slaset_( "U", &n, &n, &zero, &one, inv_A_local, &n );
        for( i=0; i<n-1; i++ ){
            inv_A_local[ i+n*(i+1) ] = -one;
        }

/*      Print information of task */
        printf( "=== START OF EXAMPLE ===================\n" );
        printf( "Task to be solved - A*x=b, where:\n" );
        printf( "/ 1 ......... 1 \\ \n" );
        printf( "| 0 .         : | \n" );
        printf( "| : .`.       : | \n" );
        printf( "| :  `.`.     : |  =  A  -  n*n real matrix \n" );
        printf( "| :    `.`.   : | \n" );
        printf( "| :      `.`. : | \n" );
        printf( "\\ 0 ....... 0 1 / \n\n" );
        printf( "b - random n-dimensional real vector" );
        printf( "Matrix A can be inversed very simply:\n" );
        printf( "/ 1 -1 0 .... 0 \\ \n" );
        printf( "| 0 . .  .    : | \n" );
        printf( "| : .`.`. `.  : | \n" );
        printf( "| :  `.`.`. ` : |  =  inv_A \n" );
        printf( "| :    `.`.`. 0 | \n" );
        printf( "| :      `.`.-1 | \n" );
        printf( "\\ 0 ....... 0 1 / \n\n" );
        printf( "n = %d, nb = %d; %dx%d - process grid\n\n", n, nb, nprow, npcol );
        printf( "=== PROGRESS ===========================\n" );
    } else {

/*      Other processes don't contain parts of initial arrays */
        A_local = NULL;
        inv_A_local = NULL;
        b_local = NULL;
    }

/*  Compute precise length of local pieces and allocate array on
    each process for parts of distributed vectors */
    mp = numroc_( &n, &nb, &myrow, &i_zero, &nprow );
    nq = numroc_( &n, &nb, &mycol, &i_zero, &npcol );
    nq_rhs = numroc_( &i_one, &i_one, &mycol, &i_zero, &npcol );
    A = (float*) mkl_calloc(mp*nq, sizeof( float ), 64);
    inv_A = (float*) mkl_calloc(mp*nq, sizeof( float ), 64);

/*  Only #0 process column contains distributed right-hand side vector */
    if ( nq_rhs != 0 ){
        x = (float*) mkl_calloc(mp, sizeof( float ), 64);
        b = (float*) mkl_calloc(mp, sizeof( float ), 64);
    } else {
        x = NULL;
        b = NULL;
    }

/*  Compute leading dimensions */
    lld_local = MAX( numroc_( &n, &n, &myrow, &i_zero, &nprow ), 1 );
    lld = MAX( mp, 1 );

/*  Initialize descriptors for initial arrays located on 0 process */
    descinit_( descA_local, &n, &n, &n, &n, &i_zero, &i_zero, &ictxt, &lld_local, &info );
    descinit_( descinvA_local, &n, &n, &n, &n, &i_zero, &i_zero, &ictxt, &lld_local, &info );
    descinit_( descb_local, &n, &i_one, &n, &i_one, &i_zero, &i_zero, &ictxt, &lld_local, &info );

/*  Initialize descriptors for distributed arrays */
    descinit_( descA, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descinvA, &n, &n, &nb, &nb, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descb, &n, &i_one, &nb, &i_one, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descx, &n, &i_one, &nb, &i_one, &i_zero, &i_zero, &ictxt, &lld, &info );

/*  Distribute vectors and matrices from 0 process over process grid */
    psgeadd_( &trans, &n, &n, &one, A_local, &i_one, &i_one, descA_local, &zero, A, &i_one, &i_one, descA );
    psgeadd_( &trans, &n, &n, &one, inv_A_local, &i_one, &i_one, descinvA_local, &zero, inv_A, &i_one, &i_one, descinvA );
    psgeadd_( &trans, &n, &i_one, &one, b_local, &i_one, &i_one, descb_local, &zero, b, &i_one, &i_one, descb );
    if( iam == 0 ){ printf( ".. Arrays are distributed ( p?geadd ) ..\n" ); }

/*  Destroy arrays on 0 process - they are not necessary anymore */
    if( ( myrow == 0 ) && ( mycol == 0 ) ){
        mkl_free( inv_A_local );
        mkl_free( A_local );
        mkl_free( b_local );
    }

/*  Copy vector b to x */
    pscopy_( &n, b, &i_one, &i_one, descb, &i_one, x, &i_one, &i_one, descx, &i_one );


/*  Solve system of equations A*x=b ( on entry array x contains vector b,
    on exit - solution of the system ) */
    pstrsv_( "U", "N", "U", &n, A, &i_one, &i_one, descA, x, &i_one, &i_one, descx, &i_one );
    if( iam == 0 ){ printf( ".. System of equations A*x=b is solved ( p?trsv ) ..\n" ); }

/*  Compute vector inv_A*b ( on entry array b contains vector b,  on exit
    - vector inv_A*b ) */
    pstrmv_( "U", "N", "U", &n, inv_A, &i_one, &i_one, descA, b, &i_one, &i_one, descb, &i_one );
    if( iam == 0 ){ printf( ".. Vector inv_A*b ( = x' ) is computed ( p?trmv ) ..\n" ); }


/*  Compute norm of matrix A */
    work = (float*) mkl_calloc(mp, sizeof( float ), 64);
    anorm = pslange_( "I", &n, &n, A, &i_one, &i_one, descA, work );
    mkl_free( work );

/*  Compute difference  x - inv_A*b  ( it can be done in parallel because
    both vectors have the same distribution ) */
    if ( nq_rhs != 0 ){
        for( i=0; i<mp; i++ ){
            x[ i ] = x[ i ] - b[ i ];
        }
    }
    psnrm2_( &n, &diffnorm, x, &i_one, &i_one, descx, &i_one );

/*  Compute machine epsilon */
    eps = pslamch_( &ictxt, "e" );

/*  Compute residual */
    residual = diffnorm / ( two*eps*anorm );

    if( iam == 0 ){
        printf( ".. Solutions are compared ..\n" );
        printf( "== Results ==\n" );
        printf( "|| x - inv_A*b ||/(2*eps*||A||) = %03.6f\n", residual );
        printf( "=== END OF EXAMPLE =====================\n" );
    }

/*  Destroy arrays */
    if ( nq_rhs != 0 ){
        mkl_free( x );
        mkl_free( b );
    }
    mkl_free( A );
    mkl_free( inv_A );

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
  ====== End of PBLAS Level 2 example ====================================
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
