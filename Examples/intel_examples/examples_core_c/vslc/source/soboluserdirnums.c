/*******************************************************************************
* Copyright 2003-2020 Intel Corporation.
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
!    User-defined direction numbers for Sobol QRNG  Example Program Text
!    Owen-scrambling strategy
!******************************************************************************/
#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"


/*
   This sample demonsrates registration of user-defined parameters
   in Sobol QRNG. Owen scrambling strategy with slight changes
   is also implemented, [1]. Set of initial Sobol parameters used here
   follows to [2].

   Bibliography
   1. H.S. Hong, and F.J. Hickernell. Algorithm 823: Implementing
   Scrambled Digital Sequences. ACM Transactions on Mathematical Software,
   Vol.29, No. 2, 95-109, 2003
   2. Sobol, I.M., and Levitan, Yu.L. The production of points uniformly
   distributed in a multidimensional cube. Preprint 40, Institute of
   Applied Mathematics, USSR Academy of Sciences, 1976 (In Russian)
*/

/* Dimension of quasi-random vectors */
#define DIM      40
/* Number of bits used for storing unsigned integer */
#define NBITS    32
/* Number of primitive polynomials */
#define NPOLY    39
/* Maximum degree of primitive polynomials */
#define MAXDEG    8

/* Number of digits to scramble */
#define S_DIGITS 32

/* Table of initial direction numbers */
static unsigned int uSobolMInit[NPOLY][MAXDEG];
/* Array of primitive polynomials */
static unsigned int uSobolPoly[NPOLY];
/* Array of degrees of primitive polynomials */
static int iSobolPolyDeg[NPOLY];
/* Generator matrix, table of direction numbers */
static unsigned int uSobolV[DIM][NBITS];
/* Generator matrix, table of modified direction numbers */
static unsigned int uSobolSV[DIM][NBITS];

/* Routine to initialize matrix of Sobol */
static void InitSobolMatrix( int dimen );
/* Routine to initialize Sobol matrix for Owen-like scrambled sequence */
static void InitSobolMatrixOs( int dimen );
/* Routine to generate sequence of low triangular matrices for
   Owen scrambling */
static void GenerateLSM( int dimen, unsigned int lsm[DIM][S_DIGITS] );

/* Dimension of vectors to generate */
#define DM 20

#define NN  10

/* Size of buffer for quasi-random sequence */
#define N   1000*DM
/* Size of params array */
#define NU  DIM*NBITS + 3

/* Buffer for quasi-random sequence */
double  r[N];
double rt[N];

int main()
{
    VSLStreamStatePtr stream;

    unsigned int params[NU];
    MKL_INT brng;
    int dim;
    int errcode;
    double a, b;
    int i, j, k;
    MKL_INT nparams;
    int errnum;


    /***** Initialize *****/
    brng = VSL_BRNG_SOBOL;
    a = 0.0;
    b = 1.0;
    dim = DM;

    /* GENERATION OF SOBOL SEQUENCE USING DIRECTION
       NUMBERS REGISTERED IN VSL */
    /* Initialize */
    errcode = vslNewStream( &stream, brng, (MKL_INT)dim );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, N, rt, a, b );
    CheckVslError( errcode );

    /***** Printing results *****/
    printf("Sobol sequence with default direction numbers:\n");
    printf("--------------------------------------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);

    k = 12;
    printf("Results: %d component of quasi-random sequence (first %d of %d):\n", k, NN, N/dim);

    for(i=0;i<NN;i++) {
        printf("%.4f ",rt[k+i*dim]);
    }
    printf("\n\n");

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );


    /* GENERATION OF SOBOL SEQUENCE USING EXTERNAL DIRECTION NUMBERS */
    /* For a moment they are the same as those registered in the library */

    /* Calculate Sobol matrix */
    InitSobolMatrix( dim );
    /* Pack Sobol matrix to pass into the generator */
    params[0] = dim;
    params[1] = VSL_USER_QRNG_INITIAL_VALUES;
    params[2] = VSL_USER_DIRECTION_NUMBERS;
    k = 3;
    for ( i = 0; i < dim; i++ )
    for ( j = 0; j < NBITS; j++ ) params[k++] = uSobolV[i][j];
    nparams = 3 + dim * NBITS;

    /***** Initialize ****/
    errcode = vslNewStreamEx( &stream, brng, nparams, params );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, N, r, a, b );
    CheckVslError( errcode );

    /***** Printing results *****/
    printf("Sobol sequence with user-defined direction numbers:\n");
    printf("--------------------------------------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);

    k = 12;
    printf("Results: %d component of quasi-random sequence (first %d of %d):\n", k, NN, N/dim);

    for(i=0;i<NN;i++) {
        printf("%.4f ",r[k+i*dim]);
    }
    printf("\n\n");

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    errnum = 0;
    for ( i = 0; i < N; i++ )
    {
        if (rt[i] != r[i] ) errnum++;
    }

    if ( errnum == 0 )
        printf("Sobol sequences obtained using two methods are similar\n\n");
     else
        printf("Sobol sequences obtained using two methods are not similar\n\n");

    /* GENERATION OF RANDOMIZED SEQUENCE */
    /* Owen-like scrambling strategy with slight changes
       is used for demonstration */

    /* Modify matrix of the generator */
    InitSobolMatrixOs( dim );
    /* pack Sobol generator matrix to pass into the generator */
    params[0] = dim;
    params[1] = VSL_USER_QRNG_INITIAL_VALUES;
    params[2] = VSL_USER_DIRECTION_NUMBERS;
    k = 3;
    for ( i = 0; i < dim; i++ )
    for ( j = 0; j < NBITS; j++ ) params[k++] = uSobolSV[i][j];
    nparams = 3 + dim * NBITS;

    /***** Initialize *****/
    errcode = vslNewStreamEx( &stream, brng, nparams, params );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, N, r, a, b );
    CheckVslError( errcode );

    /***** Printing results *****/
    printf("Scrambled Sobol sequence:\n");
    printf("--------------------------------------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);

    k = 12;
    printf("Results: %d component of quasi-random sequence (first %d of %d):\n", k, NN, N/dim);

    for(i=0;i<NN;i++) {
        printf("%.4f ",r[k+i*dim]);
    }
    printf("\n\n");

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}

/* Routine to initialize matrix of Sobol direction numbers */
static void InitSobolMatrix( int dimen )
{
    int  i, j;
    int  d;
    int  iPolyDeg;
    unsigned int uPoly;
    unsigned int uMNew;

    for (j = 0; j < NBITS; j++)  uSobolV[0][j] = 1;

    for (i = 1; i < dimen; i++)
    {
        iPolyDeg = iSobolPolyDeg[i-1];

        for (j = 0; j < iPolyDeg; j++)
        {
            uSobolV[i][j] = uSobolMInit[i-1][j];
        }

         for ( j = iPolyDeg ; j < NBITS; j++)
        {
            uPoly = uSobolPoly[i-1];
            uMNew = uSobolV[i][j-iPolyDeg];
            for (d = iPolyDeg; d > 0 ; d--)
            {
                if ( uPoly & 0x1 )
                {
                    uMNew ^= uSobolV[i][j-d] << d;
                }
                uPoly >>= 1;
            }

            uSobolV[i][j] = uMNew;
        }
    }

    for (j = 0; j < NBITS; j++)
    {
        for (i = 0; i < dimen; i++)
        {
            uSobolV[i][j] <<= NBITS-j-1;
        }
    }

}


/* Routines for Owen-like scrambling */
/* Generator matrix uSobolV should be initialized before
   its transformation for scrambling */
static void InitSobolMatrixOs( int dimen )
{
    unsigned int lsm[DIM][S_DIGITS];
    int i,j,k,l;
    unsigned int shift;
    unsigned int column, sum, a, b, mask, sDigits, nBits;

    sDigits = S_DIGITS;
    nBits = NBITS;

    if (sDigits >= 32)
      mask = 0x0;
    else
      mask  = ( sDigits < nBits )?( 0xffffffff << sDigits ) : 0x0;

    /* Generate dimen low triangular matrices */
    GenerateLSM( dimen, lsm );

    /* Modify generator matrix to support randomized sequence */
    for ( i = 0; i < dimen; i++ )
    {
        for ( j = 0; j < nBits; j++ )
        {
            shift = 0;
            column = 0;
            for ( l = sDigits-1; l >= 0; l-- )
            {
                sum = 0;
                a = lsm[i][l];
                b = uSobolV[i][j];
                for ( k = 0; k < nBits; k++ )
                {
                    sum += ((a & 0x1) & (b & 0x1));
                    a >>= 0x1;
                    b >>= 0x1;
                }
                sum = sum & 0x1;
                sum = sum << shift;
                column  = column | sum;
                shift = shift + 1;
            }

            if ( sDigits  < nBits )
               uSobolSV[i][j] = ( mask & uSobolV[i][j] ) | column;
            else
               uSobolSV[i][j] = column;
        }
    }
}


static void GenerateLSM( int dimen, unsigned int lsm[DIM][S_DIGITS] )
{
   int i,j,k,l;
   VSLStreamStatePtr stream;
   int errcode;
   int method, seed;
   MKL_INT brng;
   double p;
   int r[DIM*S_DIGITS*(S_DIGITS-1)/2];

   unsigned int shift, elem;

   method = VSL_RNG_METHOD_BERNOULLI_ICDF;
   brng = VSL_BRNG_MCG31;
   seed = 77987391;
   p = 0.5;
   errcode = vslNewStream( &stream, brng, (MKL_INT)seed );
   CheckVslError( errcode );
   /* Generate sequence of Bernoulli random variates */
   errcode = viRngBernoulli
       ( method, stream, dimen*0.5*S_DIGITS*(S_DIGITS-1), r, p );
   CheckVslError( errcode );
   errcode = vslDeleteStream( &stream );
   CheckVslError( errcode );

   l = 0;
   for ( i = 0; i < dimen; i++ )
   {
        for ( j = S_DIGITS-1;  j >= 0; j-- )
        {
            shift = NBITS-j-1;
            lsm[i][j] = (1 << shift);
            for ( k = j-1; k>=0; k-- )
            {
                shift = shift + 1;
                elem = r[l] << shift;
                lsm[i][j] = lsm[i][j] | elem;
                l = l + 1;
            }
        }
   }

}


/*
// Set of initial Sobol parameters used here
// follows to [2]
*/
/* Table of initial direction numbers */
static unsigned int uSobolMInit[NPOLY][MAXDEG] = {
   1,0,0, 0, 0, 0,  0, 0,
   1,1,0, 0, 0, 0,  0, 0,
   1,3,7, 0, 0, 0,  0, 0,
   1,1,5, 0, 0, 0,  0, 0,
   1,3,1, 1, 0, 0,  0, 0,
   1,1,3, 7, 0, 0,  0, 0,
   1,3,3, 9, 9, 0,  0, 0,
   1,3,7,13, 3, 0,  0, 0,
   1,1,5,11,27, 0,  0, 0,
   1,3,5, 1,15, 0,  0, 0,
   1,1,7, 3,29, 0,  0, 0,
   1,3,7, 7,21, 0,  0, 0,
   1,1,1, 9,23,37,  0, 0,
   1,3,3, 5,19,33,  0, 0,
   1,1,3,13,11, 7,  0, 0,
   1,1,7,13,25, 5,  0, 0,
   1,3,5,11, 7,11,  0, 0,
   1,1,1, 3,13,39,  0, 0,
   1,3,1,15,17,63, 13, 0,
   1,1,5, 5, 1,27, 33, 0,
   1,3,3, 3,25,17,115, 0,
   1,1,3,15,29,15, 41, 0,
   1,3,1, 7, 3,23, 79, 0,
   1,3,7, 9,31,29, 17, 0,
   1,1,5,13,11, 3, 29, 0,
   1,3,1, 9, 5,21,119, 0,
   1,1,3, 1,23,13, 75, 0,
   1,3,3,11,27,31, 73, 0,
   1,1,7, 7,19,25,105, 0,
   1,3,5, 5,21, 9,  7, 0,
   1,1,1,15, 5,49, 59, 0,
   1,1,1, 1, 1,33, 65, 0,
   1,3,5,15,17,19, 21, 0,
   1,1,7,11,13,29,  3, 0,
   1,3,7, 5, 7,11,113, 0,
   1,1,5, 3,15,19, 61, 0,
   1,3,1, 1, 9,27, 89, 7,
   1,1,3, 7,31,15, 45,23,
   1,3,3, 9, 9,25,107,39
};


/* Array of primitive polynomials */
static unsigned int uSobolPoly[NPOLY] = {
        3,   7,  11,  13,  19,  25,  37,  59,  47,
  61,  55,  41,  67,  97,  91, 109, 103, 115, 131,
 193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
 213, 191, 253, 203, 211, 239, 247, 285, 369, 299
};

/* Array of degrees of primitive polynomials */
static int iSobolPolyDeg[NPOLY] = {
      1,  2,  3,  3,  4,  4,  5,  5,  5,
  5,  5,  5,  6,  6,  6,  6,  6,  6,  7,
  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
  7,  7,  7,  7,  7,  7,  7,  8,  8,  8
};
