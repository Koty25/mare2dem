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
 *==== PBLAS Level 1 example program ======================================
 *=========================================================================
 *
 * The example performs computation of angle between two multi-dimensional
 * vectors by means of   PBLAS Level 1   routines.  It is designed to show
 * the use of PBLAS Level 1 routines, BLACS routines, work with 1D process
 * grid and block-cyclic distribution.
 *
 * List of PBLAS Level 1 routines used in the example:
 *     PSNRM2,
 *     PSSCAL,
 *     PSDOT
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
 *
 *
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
 * User  can specify dimension of vectors,  size of blocks,  length of each
 * vector, angle between vectors and threshold for residual checking in the
 * input file 'pblas1ex.in'.
 * In  the  beginning of program 1D process grid is initialized by means of
 * BLACS routines. Data is read from the input file by  #0 process and sent
 * to all other processes.
 * A  pair  of vectors with given angle between them and random orientation
 * are created on  #0 process and distributed over process grid.  Scheme of
 * block-cyclic distribution of vector:
 *
 *         x on process P0                    P1                   P2
 *          _               _
 *      x0 |   ------->  x0|
 *      x1 |_0   copy    x1|_0
 *      x2 |   --\  /->  x6|
 *      x3 |_1    \/     x7|_3                  _
 *      x4 |   --\/\---------------------->  x2|
 *      x5 |_2   /\          send/recieve    x3|_1
 *      x6 |   -/  \/--------------------->  x8|_4                   _
 *      x7 |_3     /\------------------------------------------>  x4|
 *      x8 |_4 ---/                               send/recieve    x5|_2
 *
 *      This scheme demonstrates distribution of vector with
 *      n=9, nb=2 and 1D process grid with 3 processes.
 *
 * Cosine of angle between vectors is computed by means of PBLAS Level 1
 * and formula:
 *
 *                   dot( X, Y )
 *      cos( X^Y ) = -----------
 *                   ||X||*||Y||
 *
 * Finally,  computed  angle  and  length of each vectors are compared with
 * input ones.
 *
 *========================================================================*/

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_scalapack.h>
#include "mkl.h"
#include "pblas_examples.h"

/* Parameters */
const float   zero = 0.0e+0, one = 1.0e+0, two = 2.0e+0, oneeighty = 180.0e+0;
const MKL_INT i_zero = 0, i_one = 1, i_two = 2, i_four = 4, i_negone = -1;
const char    trans = 'N';


/*==== MAIN FUNCTION =================================================*/
int main( int argc, char *argv[] ){

/*  ==== Declarations ==================================================*/

/*  File variables */
    FILE    *fin;

/*  Matrix descriptors */
    MDESC   descx, descy, descx_local, descy_local;

/*  Local scalars */
    MKL_INT iam, nprocs, ictxt, myrow, mycol, nprow, npcol;
    MKL_INT lld, lld_local, info, n, nb, nq, i, j;
    int     n_int, nb_int;
    float   norm_x, norm_y, rec_norm_x, rec_norm_y, cosine, new_val_i;
    float   length_x, length_y, alpha_deg, alpha_rad, thresh, eps, residual;

/*  Local arrays */
    float   *x_local, *y_local, *x, *y;
    float   swork[ 4 ];
    MKL_INT iwork[ 2 ];

/*  ==== Create process grid ===========================================

    Get information about how many processes are used for program execu-
    tion and number of current process */
    blacs_pinfo_( &iam, &nprocs );

/*  Init working process grid */
    blacs_get_( &i_negone, &i_zero, &ictxt );
    blacs_gridinit_( &ictxt, "R", &nprocs, &i_one );

/*  ==== Read data =====================================================

    Open input file */
    if ( iam == 0 ) {
        fin = fopen( "../in/pblas1ex.in", "r" );
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
        fscanf( fin, "%f alpha, angle between vectors (measured in degrees), 0.0 <= alpha <= 180.0 ", &alpha_deg );
        fscanf( fin, "%f length of vector x, must be > 0.0 ", &length_x );
        fscanf( fin, "%f length of vector y, must be > 0.0 ", &length_y );
        fscanf( fin, "%f threshold for residual check (to switch off check set it < 0.0) ", &thresh );
        fclose( fin );

        n = (MKL_INT) n_int;
        nb = (MKL_INT) nb_int;

/*      Check if all parameters are correct */
        if( ( n<=0 )||( n>1000 )||( nb<=0 )||( alpha_deg<zero )||( alpha_deg>oneeighty )||( length_x<=zero )||(length_y<=zero) ) {
            printf( "One or several input parameters has incorrect value. Limitations:\n" );
            printf( "1000 >= n > 0, nb > 0 - integer\n" );
            printf( "0.0 <= alpha <= 180.0 - float (degrees)\n" );
            printf( "length_x > 0.0, length_y > 0.0 - float\n" );
            printf( "threshold - float (set to negative to swicth off check)\n");
            return 2;
        }

/*      Pack data into array and send it to other processes */
        iwork[ 0 ] = n;
        iwork[ 1 ] = nb;
        swork[ 0 ] = alpha_deg;
        swork[ 1 ] = length_x;
        swork[ 2 ] = length_y;
        swork[ 3 ] = thresh;
        sgebs2d_( &ictxt, "All", " ", &i_four, &i_one, swork, &i_four );
        igebs2d_( &ictxt, "All", " ", &i_two, &i_one, iwork, &i_two );
    } else {

/*      Recieve and unpack data */
        sgebr2d_( &ictxt, "All", " ", &i_four, &i_one, swork, &i_four, &i_zero, &i_zero );
        igebr2d_( &ictxt, "All", " ", &i_two, &i_one, iwork, &i_two, &i_zero, &i_zero );
        n = iwork[ 0 ];
        nb = iwork[ 1 ];
        alpha_deg = swork[ 0 ];
        length_x = swork[ 1 ];
        length_y = swork[ 2 ];
        thresh = swork[ 3 ];
    }

/*  Create on process 0 two random vectors with angle 'alpha' degrees
    between them */
    if ( iam == 0 ){

/*      Allocate arrays for vectors */
        x_local = (float*) mkl_calloc(n, sizeof( float ), 64);
        y_local = (float*) mkl_calloc(n, sizeof( float ), 64);

/*      Initialize vectors */
        for( i=0; i<n; i++ ){
            x_local[ i ] = zero;
            y_local[ i ] = zero;
        }

/*      Angle between vectors computed in radians */
        alpha_rad = ( alpha_deg * M_PI ) / oneeighty;

/*      Set components to obtain two vectors with given length and angle
        between them */
        x_local[ 0 ] =  length_x;
        y_local[ 0 ] =  length_y * cos(alpha_rad);
        y_local[ 1 ] = -length_y * sin(alpha_rad);

/*      Rotate pair of vectors in a random way. Angle between vectors
        stands to be a given constant */
        srand(1);
        for( j=1; j<n; j++ ){
            for( i=j-1; i>=0; i-- ){
                alpha_rad = two*M_PI*rand()/RAND_MAX;
                new_val_i    =   x_local[ i ]*cos( alpha_rad ) + x_local[ j ]*sin( alpha_rad );
                x_local[ j ] = - x_local[ i ]*sin( alpha_rad ) + x_local[ j ]*cos( alpha_rad );
                x_local[ i ] = new_val_i;
                new_val_i    =   y_local[ i ]*cos( alpha_rad ) + y_local[ j ]*sin( alpha_rad );
                y_local[ j ] = - y_local[ i ]*sin( alpha_rad ) + y_local[ j ]*cos( alpha_rad );
                y_local[ i ] = new_val_i;
            }
        }

/*      Print information of task */
        printf( "=== START OF EXAMPLE ===================\n" );
        printf( "Two vectors in %d-dimensional space\n", n );
        printf( "Block size for distribution: %d\n\n", nb );
        printf( "=== PROGRESS ===========================\n" );
    } else {

/*      Other processes don't contain parts of initial arrays */
        x_local = NULL;
        y_local = NULL;
    }

/*  Compute precise length of local pieces and allocate array on
    each process for parts of distributed vectors */
    nq = numroc_( &n, &nb, &iam, &i_zero, &nprocs );
    x = (float*) mkl_calloc(nq, sizeof( float ), 64);
    y = (float*) mkl_calloc(nq, sizeof( float ), 64);

/*  Compute leading dimensions */
    lld_local = MAX( numroc_( &n, &n, &iam, &i_zero, &nprocs ), 1 );
    lld = MAX( nq, 1 );

/*  Initialize descriptors for initial arrays located on 0 process */
    descinit_( descx_local, &n, &i_one, &n, &i_one, &i_zero, &i_zero, &ictxt, &lld_local, &info );
    descinit_( descy_local, &n, &i_one, &n, &i_one, &i_zero, &i_zero, &ictxt, &lld_local, &info );

/*  Initialize descriptors for distributed arrays */
    descinit_( descx, &n, &i_one, &nb, &i_one, &i_zero, &i_zero, &ictxt, &lld, &info );
    descinit_( descy, &n, &i_one, &nb, &i_one, &i_zero, &i_zero, &ictxt, &lld, &info );

/*  Distribute vectors from 0 process over process grid */
    psgeadd_( &trans, &n, &i_one, &one, x_local, &i_one, &i_one, descx_local, &zero, x, &i_one, &i_one, descx );
    psgeadd_( &trans, &n, &i_one, &one, y_local, &i_one, &i_one, descy_local, &zero, y, &i_one, &i_one, descy );
    if( iam == 0 ){ printf( ".. Vectors are distributed ( p?geadd ) ..\n" ); }

/*  Destroy arrays from 0 process - they are not necessary anymore */
    if ( iam == 0 ){
        mkl_free( x_local );
        mkl_free( y_local );
    }

/*  ====================================================================

    Computation of angle between vectors using PBLAS level 1 routines
    The formula is:

                 dot( X, Y )          X      Y
    cos( X^Y ) = ----------- = dot( -----, ----- )
                 ||X||*||Y||        ||X||  ||Y||

    where:
    X^Y         - angle between vectors X and Y,
    dot( X, Y ) - scalar product of two vectors,
    ||.||       - Euclidian norm of vector

    --------------------------------------------------------------------

    Call PBLAS level 1 routine for computation of Euclidian norm of
    vector to compute norms of distributed vectors, ||X|| and ||Y||   */
    psnrm2_( &n, &norm_x, x, &i_one, &i_one, descx, &i_one );
    psnrm2_( &n, &norm_y, y, &i_one, &i_one, descy, &i_one );
    if( iam == 0 ){ printf( ".. Value of norm for each vector is computed ( p?nrm2 ) ..\n" ); }

/*  Compute reciprocal to norms */
    rec_norm_x = 1 / norm_x;
    rec_norm_y = 1 / norm_y;

/*  Call PBLAS level 1 routine for multiplication vector by scalar
                  X           Y
    to compute  -----  and  -----
                ||X||       ||Y||                               */
    psscal_( &n, &rec_norm_x, x, &i_one, &i_one, descx, &i_one );
    psscal_( &n, &rec_norm_y, y, &i_one, &i_one, descy, &i_one );
    if( iam == 0 ){ printf( ".. Vectors are normalized ( p?scal ) ..\n" ); }

/*  Call PBLAS level 1 routine to compute scalar product of vectors
      X           Y
    -----  and  -----
    ||X||       ||Y||                                            */
    psdot_( &n, &cosine, x, &i_one, &i_one, descx, &i_one, y, &i_one, &i_one, descy, &i_one );
    if ( cosine > 1.0 ) { cosine = 1.0; }
    if ( cosine < -1.0 ) { cosine = -1.0; }
    if( iam == 0 ){ printf( ".. Cosine between vectors is computed ( p?dot ) ..\n" ); }

/*  ====================================================================

    Destroy arrays */
    mkl_free( x );
    mkl_free( y );

/*  Angle can be computed using acos(.) function */
    alpha_rad = acos( cosine );

/*  Print results */
    if ( iam == 0 ){
        printf( "== Results ==\n" );
        printf( "||x||=%03.6f - length of vector x\n", norm_x );
        printf( "||y||=%03.6f - length of vector y\n", norm_y );
        printf( "cos(x^y)=%01.6f\n", cosine );
        printf( "angle(x^y)=%03.6f degrees\n", oneeighty*alpha_rad/M_PI );
        printf( "=== END OF EXAMPLE =====================\n" );
    };

/*  Compute residual */
    eps = pslamch_( &ictxt, "e" );
    residual = fabs( oneeighty*alpha_rad/M_PI-alpha_deg )*( oneeighty*eps )
               + fabs( length_x - norm_x )/( length_x*eps )
               + fabs( length_y - norm_y )/( length_y*eps );

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
  ====== End of PBLAS Level 1 example ====================================
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
