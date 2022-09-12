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
!    Construction of natural cubic spline, interpolation, and computation
!    of the 2nd derivative Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              9 // number of break points
#define NY             1 // number of functions
#define NSITE_BLOCK   10 // number of sites for spline-based interpolation
                         // in one block
#define NBLOCKS        4 // number of blocks
#define NDORDER        3 // size of array describing derivative orders
                         // to compute
#define NDER           2 // number of derivatives to compute

#define NSCOEFF     NY*(N-1)*DF_PP_CUBIC    // total number of spline
                                            // coefficients
#define NSITE       (NSITE_BLOCK*NBLOCKS)   // total number of interpolation
                                            // sites

#define LLIM_X      -1.0 // left limit of interpolation interval
#define RLIM_X       2.0 // right limit of interpolation interval
#define FREQ         1.7


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
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT nbc;                        // number of boundary conditions
    MKL_INT ic_type;                    // internal conditions type
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
    MKL_INT *cell_idx;                  // indices of cells containing
                                        // interpolation sites
    double left = LLIM_X;               // left limit of the interpolation
                                        // interval
    double right = RLIM_X;              // right limit of the interpolation
                                        // interval
    double x[2];                        // limits of the interpolation interval
    double y[NY*N];                     // function values
    double *ic;                         // array of internal conditions
    double *bc;                         // array of boundary conditions
    double scoeff[NSCOEFF];             // array of spline coefficients
    double site[NSITE];                 // array of interpolation sites
    double r[NDER*NSITE];               // spline evaluation results
    double *datahint;                   // additional info about structure
                                        // of arrays x and y

    double *site_ptr, *r_ptr;
    double ref_r[NDER*NSITE];

    double freq;

    double xx[N];

    int i, j, errcode = 0;
    int nder;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder     = DF_PP_CUBIC;
    stype      = DF_PP_NATURAL;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_UNIFORM_PARTITION;
    /* Limits of interpolation interval are provided in case
       of uniform partition */
    x[0]       = left;
    x[1]       = right;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    nscoeff    = NSCOEFF;
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type    = DF_BC_FREE_END;
    /* No additional parameters are provided for Free-End boundary conditions */
    bc         = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing interpolation sites *****/
    nsite_bl   = NSITE_BLOCK;
    nsite      = NSITE;
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_MATRIX_STORAGE_ROWS;

    /**** Parameter describing array of cell indices *****/
    /* No cell indices are computed */
    cell_idx   = 0;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /***** Number of derivatives to compute *****/
    nder = NDER;

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    freq = FREQ;
    errcode = dSinDataUniformGrid( y, freq, left, right, nx );
    CheckDfError(errcode);

    /***** Generate interpolation sites *****/
    errcode = dUniformRandData( site, left, right, nsite );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/

    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, 0, ic_type, 0,
                                 scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct natural cubic spline using STD method *****/
    errcode = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Interpolate, and compute 2nd order derivatives using STD method *****/

    nblocks = NBLOCKS;

    site_ptr = site;
    r_ptr    = r;
    for ( i = 0; i < nblocks; i++ )
    {
        errcode = dfdInterpolate1D( task, DF_INTERP, DF_METHOD_PP,
                                    nsite_bl, site_ptr, sitehint, ndorder,
                                    dorder, datahint, r_ptr, rhint, cell_idx );
        CheckDfError(errcode);

        site_ptr += nsite_bl;
        r_ptr    += nder*nsite_bl;
    }

    /***** Check results of interpolation *****/
    errcode = dUniformData( xx, left, right, nx );
    CheckDfError(errcode);

    errcode = dCheckCubInterpRes( nx, xx, ny, scoeff,
                                  nsite, site, ndorder, dorder,
                                  r, ref_r );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print function *****/

    printf("\n i  x(i)        y(i)\n");
    for( j = 0; j < nx; j++ )
    {
        printf(" %1d %+lf   %+lf\n", j, xx[j], y[j]);
    }

    /***** Print spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci            Di\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f\n", j,
            scoeff[sorder*j], scoeff[sorder*j + 1], scoeff[sorder*j + 2],
            scoeff[sorder*j + 3]);
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
