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
!    Construction and integration of linear spline Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N             10 // number of break points
#define NY             2 // number of functions
#define NLIM           3 // number of pairs of integration limits

#define NSCOEFF        (NY*(N-1)*DF_PP_LINEAR)  // total number of spline
                                                // coefficients

#define LLIM_X      -2.0 // left  limit of interpolation interval
#define RLIM_X       2.0 // right limit of interpolation interval
#define LEFT_LLIM   -1.5
#define RIGHT_LLIM  -0.5
#define LEFT_RLIM    0.5
#define RIGHT_RLIM   1.5

#define FREQ         1.3

int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type
    MKL_INT nx;                         // number of break points
    MKL_INT xhint;                      // additional info about break points
    MKL_INT ny;                         // number of functions
    MKL_INT yhint;                      // additional info about function
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT ic_type;                    // internal conditions type
    MKL_INT nlim;                       // number of pairs of integration limits
    MKL_INT llimhint;                   // additional info about the structure
                                        // of left integration limits
    MKL_INT rlimhint;                   // additional info about the structure
                                        // of right integration limits
    MKL_INT rhint;                      // integration results storage format
    float x[N];                         // array of break points
    float y[NY*N];                      // function values
    float *ic;                          // internal conditions
    float *bc;                          // boundary conditions
    float scoeff[NSCOEFF];              // array of spline coefficients
    float llim[2], rlim[2];             // left and right limits
                                        // of the intervals are provided
                                        // if integration limits are uniform

    float *ldatahint, *rdatahint;       // additional info about the structure
                                        // of x and integration limits
    float r[NY*NLIM];                   // integration results
    float *scoeff_ptr;

    float llim_full[NLIM];              // left integration limits
                                        // in full format
    float rlim_full[NLIM];              // right integration limits
                                        // in full format
    float left = LLIM_X, right = RLIM_X;
    float freq;

    int i, j, errcode = DF_STATUS_OK;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder = DF_PP_LINEAR;
    stype  = DF_PP_DEFAULT;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    /* No boundary conditions are provided for linear spline */
    bc_type    = DF_NO_BC;
    bc         = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for linear spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters decsribing integration limits *****/
    nlim = NLIM;
    llimhint = DF_UNIFORM_PARTITION;
    rlimhint = DF_UNIFORM_PARTITION;
    /* Left and right points of interval are provided in case of uniform
       partition */
    llim[0] = LEFT_LLIM;
    llim[1] = RIGHT_LLIM;
    rlim[0] = LEFT_RLIM;
    rlim[1] = RIGHT_RLIM;

    /***** Additional information about the structure of integration
           limits *****/
    /* No additional info is provided */
    ldatahint = 0;
    rdatahint = 0;

    /***** Parameter dascribing integration results storage format *****/
    rhint = DF_NO_HINT;

    /***** Generate array of uniformly distributed break points *****/
    errcode = sUniformRandSortedData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Generate functions y = sin(2 * Pi * freq * x) *****/
    for( i = 0; i < ny; i++ )
    {
        freq = (i + 1) * FREQ;
        errcode = sSinDataNotUniformGrid( y + i*nx, x, freq, nx );
        CheckDfError(errcode);
    }

    /***** Create Data Fitting task *****/
    errcode = dfsNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for linear spline construction *****/
    errcode = dfsEditPPSpline1D( task, sorder, stype, bc_type, bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct linear spline using STD method *****/
    errcode =  dfsConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Compute integrals *****/
    errcode = dfsIntegrate1D( task, DF_METHOD_PP, nlim, llim, llimhint,
                              rlim, rlimhint, ldatahint, rdatahint,
                              r, rhint );

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/
    errcode = sCheckLinSplineBreakPoints( nx, x, ny, y, sorder,
                                          scoeff );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given function *****/
    printf("\n  X           Y1          Y2\n");

    for( j = 0; j < nx; j++ )
    {
        printf(" %+lf   %+lf   %+lf\n", x[j], y[j], y[j+nx]);
    }

    /***** Print computed spline coefficients *****/

    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i))\n");
    printf("    where x(i) <= x < x(i+1)\n");

    scoeff_ptr = scoeff;
    for( i = 0; i < ny; i++ )
    {
        printf("\nSpline coefficients for Y%d :\n", i+1 );
        printf(" i    Ai            Bi\n");
        for( j = 0; j < nx-1; j++ )
        {
            printf(" %1d %+11.6f   %+11.6f\n", j,
                    scoeff_ptr[sorder*j],
                    scoeff_ptr[sorder*j + 1]);
        }
        scoeff_ptr += sorder*(nx-1);
    }

    /***** Print integration results *****/

    sUniformData( llim_full, llim[0], llim[1], nlim );
    CheckDfError(errcode);
    sUniformData( rlim_full, rlim[0], rlim[1], nlim );
    CheckDfError(errcode);

    for( i = 0; i < ny; i++ )
    {
        printf("\nIntegration results for Y%d :\n", i+1 );
        printf("Integration interval      Result\n");
        for ( j = 0; j < nlim; j++ )
        {
            printf(" ( %+4.1lf, %+4.1lf )        %+11.6lf\n", llim_full[j],
                rlim_full[j], r[i*nlim + j]);
        }
    }

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed linear spline coefficients or integrals");
        printf(" are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed linear spline coefficients and integrals");
        printf(" are correct\n");
    }

    return 0;
}
