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
!    Use of Intel(R) MKL Sobol QRNG with direction numbers designed
!    by S. Joe,  F. Y. Kuo and supporting dimension up to 21201
!    Example Program Text
!******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl_service.h"
#include "mkl_vml.h"
#include "mkl_vsl.h"
#include "errcheck.inc"


/*
   This sample demonsrates registration of user-defined parameters
   in Sobol QRNG. Set of initial Sobol parameters used here
   follows to [1],[2].

   Bibliography
   1. S. Joe,  F. Y. Kuo. Remark on Algorithm 659: Implementing Sobol_s
   Quasirandom Sequence Generator. ACM Transactions on Mathematical Software,
   Vol. 29, No. 1, 50-57, March 2003
   2. http://web.maths.unsw.edu.au/~fkuo/sobol/
*/

/* Number of bits used for storing unsigned integer */
#define NBITS       32
/* Number of primitive polynomials */
#define NPOLY       110
/* Maximum degree of primitive polynomials for dimensions up to 111 */
#define MAXDEG      10

/* Dimension of vectors to generate */
#define DIM 111
/* Number of vectors to generate in one block */
#define M       10000
/* Number of blocks of vectors to generate */
#define NBLOCKS 100

#define NN  10

/* Size of buffer for quasi-random sequence */
#define N   M*DIM

#define NU (DIM*NBITS + 3)

#define ALIGN   32

/* Routine to check memory allocation status */
static void CheckMalloc(void *buf);

/* Routine to read initial Sobol direction numbers, polynomials
   and polynomial degrees from file */
static int ReadSobolData( int npoly, int maxdeg,
    unsigned int **uSobolMInit, unsigned int *uSobolPoly, int *iSobolPolyDeg );

/* Routine to initialize matrix of Sobol */
static void InitSobolMatrix( int dimen, int nbits, unsigned int **uSobolV,
    unsigned int **uSobolMInit, unsigned int *uSobolPoly, int *iSobolPolyDeg );

/* Routine to test generated quasi-random values by calculating an integral
   of the model function */
static void Integrate( int dimen, int n, int nblocks, MKL_INT method, VSLStreamStatePtr stream,
                       double a, double b, double *r, double *res, double *err );

int main(int argc, char *argv[])
{
    /* Table of initial direction numbers */
    int uMInit_sz;
    unsigned int **uSobolMInit;
    unsigned int *uSobolMInit_buf;
    /* Generator matrix, table of direction numbers */
    int uSobolV_sz;
    unsigned int **uSobolV;
    unsigned int *uSobolV_buf;
    /* Array of primitive polynomials */
    unsigned int *uSobolPoly;
    /* Array of degrees of primitive polynomials */
    int *iSobolPolyDeg;

    VSLStreamStatePtr stream;

    unsigned int *params;
    double *r;
    double res, err;
    MKL_INT brng;
    int dim;
    int errcode, status = 0;
    double a, b;
    int i, j, k;
    MKL_INT nparams;

    /***** Allocate memory *****/
    uMInit_sz = NPOLY*MAXDEG*sizeof(unsigned int);
    uSobolMInit_buf = (unsigned int *)mkl_malloc(uMInit_sz, ALIGN);
    CheckMalloc(uSobolMInit_buf);

    uSobolMInit = (unsigned int **)mkl_malloc(NPOLY*sizeof(unsigned int*),
        ALIGN);
    CheckMalloc(uSobolMInit);

    uSobolV_sz = DIM*NBITS*sizeof(unsigned int);
    uSobolV_buf = (unsigned int *)mkl_malloc(uSobolV_sz, ALIGN);
    CheckMalloc(uSobolV_buf);

    uSobolV = (unsigned int **)mkl_malloc(DIM*sizeof(unsigned int*), ALIGN);
    CheckMalloc(uSobolV);

    uSobolPoly = (unsigned int *)mkl_malloc(NPOLY*sizeof(unsigned int), ALIGN);
    CheckMalloc(uSobolPoly);

    iSobolPolyDeg = (int *)mkl_malloc(NPOLY*sizeof(int), ALIGN);
    CheckMalloc(iSobolPolyDeg);

    params = (unsigned int *)mkl_malloc(NU*sizeof(unsigned int), ALIGN);
    CheckMalloc(params);

    r = (double *)mkl_malloc(N*sizeof(double), ALIGN);
    CheckMalloc(r);

    for (i = 0; i < NPOLY; i++)
    {
        uSobolMInit[i] = uSobolMInit_buf + i*MAXDEG;
    }

    for (i = 0; i < DIM; i++)
    {
        uSobolV[i] = uSobolV_buf + i*NBITS;
    }

    /***** Initialize *****/
    brng = VSL_BRNG_SOBOL;
    a = 0.0;
    b = 1.0;
    dim = DIM;

    /***** Read initial Sobol direction numbers, polynomials and
           polynomial degrees from file *****/
    errcode = ReadSobolData( NPOLY, MAXDEG, uSobolMInit, uSobolPoly,
        iSobolPolyDeg );
    if (errcode != 0) return errcode;

    /* Calculate table of Sobol direction numbers using Joe and Kuo
       initial direction numbers and primitive polynomials  */
    InitSobolMatrix( dim, NBITS, uSobolV, uSobolMInit, uSobolPoly,
        iSobolPolyDeg );

    /* Pack Sobol direction numbers to pass into the generator */
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
    printf("Sobol sequence with Joe and Kuo direction numbers:\n");
    printf("--------------------------------------------------\n\n");
    printf("Parameters:\n");
    printf("    dim = %d\n", dim);
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);
    k = 12;
    printf("Results: %d component of quasi-random sequence (first %d of %d):\n",
        k, NN, N/dim);

    for ( i = 0; i < NN; i++ ) printf("%.4f ",r[k+i*dim]);
    printf("\n\n");

    /***** Test generated values by calculating intagral of the model
           function f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj)).
           Model function is taken from [1] *****/
    Integrate( dim, M, NBLOCKS, VSL_RNG_METHOD_UNIFORM_STD, stream, a, b, r, &res, &err );

    printf("Calculate integral for the function:\n");
    printf("  f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj)),  j = 1,...,dim\n");
    printf("  Cj = j^(1/3)\n");
    printf("over unit hypercube [0,1]^dim\n\n");
    printf("Computed result: %10.7lf\n", res);
    printf("Expected result: %10.7lf\n", 1.0);
    printf("Calculation error: %10.7lf\n", err);

    if (  1.0 - 3.0*err <= res && res <= 1.0 + 3.0*err )
    {
        printf("\nNumerical integration results agree with theory\n");
    }
    else
    {
        printf("\nError: Numerical integration results do not agree with theory\n");
        status = 1;
    }

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    /***** Release memory *****/
    mkl_free(uSobolMInit_buf);
    mkl_free(uSobolMInit);
    mkl_free(uSobolV_buf);
    mkl_free(uSobolV);
    mkl_free(uSobolPoly);
    mkl_free(iSobolPolyDeg);
    mkl_free(params);
    mkl_free(r);

    return status;
}

static void CheckMalloc(void *buf)
{
    if (!buf)
    {
        printf("Error: Memory allocation failed\n");
        exit(1);
    }
}

/* Routine to read initial Sobol direction numbers, polynomials
   and polynomial degrees from file */
static int ReadSobolData( int npoly, int maxdeg,
    unsigned int **uSobolMInit, unsigned int *uSobolPoly, int *iSobolPolyDeg )
{
    int i, j, dim;
    FILE *f;
    char filename[100];
    char line[255];

    scanf("%99s", filename);

    f = fopen(filename, "r");
    if (f == 0)
    {
        printf("Error while openning file with dimension values\n");
        return 1;
    }

    fgets(line, 255, f);

    for (i = 0; i < npoly; i++)
    {
        fscanf(f, "%d %d %u ", &dim, iSobolPolyDeg+i, uSobolPoly+i);

        /* Modify the coding of polynomials to the format expected
           by the library */
        uSobolPoly[i] <<= 1;
        uSobolPoly[i] |= (1 | (1 << iSobolPolyDeg[i]));

        for (j = 0; j < iSobolPolyDeg[i]; j++)
        {
            fscanf(f, "%u ", &uSobolMInit[i][j]);
        }
        fscanf(f, "\n");
    }

    fclose(f);
    return 0;
}

/* Routine to initialize matrix of Sobol direction numbers */
static void InitSobolMatrix( int dimen, int nbits, unsigned int **uSobolV,
    unsigned int **uSobolMInit, unsigned int *uSobolPoly, int *iSobolPolyDeg )
{
    int  i, j;
    int  d;
    int  iPolyDeg;
    unsigned int uPoly;
    unsigned int uMNew;

    for (j = 0; j < nbits; j++)  uSobolV[0][j] = 1;

    for (i = 1; i < dimen; i++)
    {
        iPolyDeg = iSobolPolyDeg[i-1];

        for (j = 0; j < iPolyDeg; j++)
        {
            uSobolV[i][j] = uSobolMInit[i-1][j];
        }

        for ( j = iPolyDeg ; j < nbits; j++)
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

/* Routine to calculate integral of the function
   f(X) = Prod((|4*Xj - 2| + Cj)/(1+Cj)), j = 1,...,dimen
   over unit hypercube [0,1]^dimen. Integration result should be equal to 1.
   Here Cj = j^(1/3).
   Model function is taken from [1] */
static void Integrate( int dimen, int m, int nblocks, MKL_INT method, VSLStreamStatePtr stream,
                       double a, double b, double *r, double *res, double *err )
{
    int i, j, iblock, errcode;
    double cj[1], dj[1], prod;
    double f, f2;
    double inv_n = 1.0 / ((double)m * (double)nblocks);

    f = 0.0;
    f2 = 0.0;

    for (iblock = 0; iblock < NBLOCKS; iblock++)
    {
        errcode = vdRngUniform( method, stream, m*dimen, r, a, b );
        CheckVslError( errcode );
        for (i = 0; i < m; i++)
        {
            prod = 1.0;
            for (j = 0, dj[0] = 1.0; j < dimen; j++, dj[0]+=1.0)
            {
                vdCbrt(1, dj, cj);
                prod *= (fabs(4.0*r[i*dimen + j] - 2.0) + cj[0]);
                prod /= (1.0 + cj[0]);
            }
            f += prod;
            f2 += prod*prod;
        }
    }
    f *= inv_n;
    f2 *= inv_n;

    *res = f;
    *err = sqrt((f2 - f*f)*inv_n);
}
