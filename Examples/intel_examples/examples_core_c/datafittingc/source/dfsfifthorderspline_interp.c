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
!    Interpolation, and computation of the 2nd derivative for the spline
!  of fifth order Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              4 // number of breakpoints
#define NY             1 // number of datasets to interpolate
#define NSITE_BLOCK   10 // number of sites for spline-based interpolation
                         // in one block
#define NBLOCKS        4 // number of blocks
#define NDORDER        3 // size of array describing derivative orders
                         // to compute
#define NDER           2 // number of derivatives to compute

#define SPLINE_ORDER   5
#define NSCOEFF        NY*(N-1)*SPLINE_ORDER   // total number of spline
                                               // coefficients
#define NSITE          (NSITE_BLOCK*NBLOCKS)   // total number of interpolation
                                               // sites

#define LLIM_X       0.0 // left limit of interpolation interval
#define RLIM_X       1.5 // right limit of interpolation interval

int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type
    MKL_INT nx;                         // number of break points
    MKL_INT xhint;                      // additional info about break points
    MKL_INT ny;                         // number of functions
    MKL_INT yhint;                      // additional info about function
    MKL_INT nscoeff;                    // number of spline coefficients
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT nsite_bl;                   // number of interpolaton sites
                                        // in one block
    MKL_INT nblocks;                    // number of blocks
    MKL_INT nsite;                      // total number of interpolation sites
    MKL_INT sitehint;                   // additional info about interpolation
                                        // sites
    MKL_INT ndorder;                    // size of array describing derivative
                                        // orders
    MKL_INT dorder[] = { 1, 0, 1 };     // values and 2nd derivatives
                                        // will be computed
    MKL_INT rhint;                      // interpolation results storage format
    float left = LLIM_X;                // left limit of the interpolation
                                        // interval
    float right = RLIM_X;               // right limit of the interpolation
                                        // interval
    float x[N] = {0.0, 0.5, 1.0, 1.5};  // limits of the interpolation interval
    float y[NY*N];                      // function values
    // array of spline coefficients
    float scoeff[NSCOEFF] = {
        /* Func 1 */
        0.0,   1.0, 0.0, 0.0, 0.0,
        0.5,  -1.0, 0.0, 0.0, 0.0,
        0.0,   1.0, 0.0, 0.0, 0.0
    };

    float site[NSITE];                 // array of interpolation sites
    float r[NDER*NSITE];               // spline evaluation results

    float *site_ptr, *r_ptr;
    float ref_r[NDER*NSITE];

    float xx[N];

    int i, j, errcode = 0;
    int nder;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder     = SPLINE_ORDER;
    stype      = DF_PP_NATURAL;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_QUASI_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    nscoeff    = NSCOEFF;
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing interpolation sites *****/
    nsite_bl   = NSITE_BLOCK;
    nsite      = NSITE;
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_MATRIX_STORAGE_ROWS;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /***** Number of derivatives to compute *****/
    nder = NDER;

    /***** Generate interpolation sites *****/
    errcode = sUniformRandData( site, left, right, nsite );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfsNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/

    errcode = dfsEditPPSpline1D( task, sorder, stype, 0, 0, 0, 0,
                                 scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Interpolate, and compute 2nd order derivatives using STD method *****/

    nblocks = NBLOCKS;

    site_ptr = site;
    r_ptr    = r;
    for ( i = 0; i < nblocks; i++ )
    {
        errcode = dfsInterpolate1D( task, DF_INTERP, DF_METHOD_PP,
                                    nsite_bl, site_ptr, sitehint, ndorder,
                                    dorder, 0, r_ptr, rhint, 0 );
        CheckDfError(errcode);

        site_ptr += nsite_bl;
        r_ptr    += nder*nsite_bl;
    }

    /***** Check results of interpolation *****/
    errcode = sUniformData( xx, left, right, nx );
    CheckDfError(errcode);

    errcode = sCheckFifthOrderInterpRes( nx, xx, ny, scoeff,
                                         nsite, site, ndorder, dorder,
                                         r, ref_r );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print function *****/

    printf("\n i  x(i)\n");
    for( j = 0; j < nx; j++ )
    {
        printf(" %1d %+lf \n", j, xx[j]);
    }

    /***** Print spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3 + Ei*(x - x(i))^4\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci            Di            Ei\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f\n", j,
            scoeff[sorder*j], scoeff[sorder*j + 1], scoeff[sorder*j + 2],
            scoeff[sorder*j + 3], scoeff[sorder*j + 4]);
    }

    /***** Print interpolation results ******/
    printf("\nResults of interpolation:\n");
    printf("                    Function value");
    printf("             Second derivative\n");
    printf("    Sites         Obtained     Expected");
    printf("      Obtained    Expected\n");
    for ( i = 0; i < nsite; i++ )
    {
        printf(" %+11.6lf", site[i]);

        for ( j = 0; j < nder; j++ )
        {
            printf("   %+11.6lf  %+11.6lf", r[i*nder + j],
                    ref_r[i*nder + j]);
        }
        printf("\n");
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed interpolation results");
        printf(" are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed interpolation results");
        printf(" are correct\n");
    }

    return 0;
}
