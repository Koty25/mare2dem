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
!    Calculation of Akima cubic spline coefficients and integral
!    computation with callback function  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N             10 /* number of breakpoints */
#define NY             2 /* number of datasets to interpolate */
#define NLIM           2 /* number of pairs of integration limits */

#define NSCOEFF        (N-1)*DF_PP_CUBIC*NY

#define LLIM_X         0.0  // left  limit of interpolation interval
#define RLIM_X         3.0  // right limit of interpolation interval
#define FREQ    0.75

MKL_INT* llim_cell_p;
MKL_INT* rlim_cell_p;
float * llim_p;
float * rlim_p;

/*******************************************************************************
!   Definition of the integration call back for integral calculations
! on the interval (-inf, x[0]).
!
! API
!   int  left_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
!                           MKL_INT64 rcell[], float  rlim[], float  r[],
!                           void *x0,
!                           dfIntegrCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  lcell  - array of size n with indices of the cells that contain
!           the left-side integration limits
!  llim   - array of size n that holds the left-side integration limits
!  rcell  - array of size n with indices of the cells that contain
!           the right-side integration limits
!  rlim   - array of size n that holds the right-side integration limits
!  x0     - left  limit of interpolation interval
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int  left_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
                        MKL_INT64 rcell[], float  rlim[], float  r[], void *x0,
                        dfIntegrCallBackLibraryParams* library_params );

/*******************************************************************************
!   Definition of the integration call back for integral calculations
! on the interval [x[N-1], +inf).
!
! API
!   int right_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
!                           MKL_INT64 rcell[], float  rlim[], float  r[],
!                           void *xN,
!                           dfIntegrCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  lcell  - array of size n with indices of the cells that contain
!           the left-side integration limits
!  llim   - array of size n that holds the left-side integration limits
!  rcell  - array of size n with indices of the cells that contain
!           the right-side integration limits
!  rlim   - array of size n that holds the right-side integration limits
!  xN     - right limit of interpolation interval
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int right_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
                        MKL_INT64 rcell[], float  rlim[], float  r[], void *xN,
                        dfIntegrCallBackLibraryParams* library_params );

/*******************************************************************************
!   Definition of the search call back, used by integration routine
! for search of left and right integration limits.
!
! API
!   int integr_search_cb( MKL_INT64* n, float  site[], MKL_INT64 cell[],
!                     int flag[], void* user_params,
!                     dfSearchCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  site   - array of size n with integration limits (left or right,
!           depending on value of library_params->limit_type_flag)
!  user_params - pointer to user-defined parameters
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  cell   - array of size n with indices of the cells that contain
!           corresponding integration limits in site array
!           (calculated by this search call back)
!  flag   - array of size n filled by ones here to indicate
!           that each corresponding integration limit was found
!
! RETURN VALUE:
!  The status returned by the callback function:
!   DF_STATUS_EXACT_RESULT - indicates successful completion of the callback
!   <0 - error
*******************************************************************************/
int integr_search_cb( MKL_INT64* n, float  site[], MKL_INT64 cell[],
                      int flag[], void* user_params,
                      dfSearchCallBackLibraryParams* library_params );

/*******************************************************************************
!   Definition of the structure containing user-defined parameters,
! to be passed to search call back, used during integration.
*******************************************************************************/
typedef struct _dfUserParams
{
    float  v0;
    float  v1;
    float  v2;
} dfUserParams;

int main()
{
    DFTaskPtr task;                       // Data Fitting task descriptor
    MKL_INT sorder;                       // spline order
    MKL_INT stype;                        // spline type
    MKL_INT nx;                           // number of break points
    MKL_INT xhint;                        // additional info about break points
    MKL_INT ny;                           // number of functions
    MKL_INT yhint;                        // functions storage format
    MKL_INT scoeffhint;                   // additional info about spline
                                          // coefficients
    MKL_INT bc_type;                      // boundary conditions type
    MKL_INT ic_type;                      // internal conditions type
    MKL_INT nlim;                         // number of pairs of integration
                                          // limits
    float   llim[NLIM] = {-0.5,0.0};      // left  limits of integration intervals
    float   rlim[NLIM] = { 3.5,4.0};      // right limits of integration intervals
    MKL_INT llim_cell[NLIM];  // search results for left  integration limits
    MKL_INT rlim_cell[NLIM];  // search results for right integration limits
    MKL_INT llimhint, rlimhint;           // integration limits storage formats
    MKL_INT rhint;                        // integration results storage format
    float  left = LLIM_X, right = RLIM_X; // limits of interpolation interval
    float  x[N];                          // break points
    float  y[N*NY];                       // function values
    float  bc[] = { 0.0, -0.5 };          // boundary conditions
    float  *ic;                           // internal conditions
    float  scoeff[NSCOEFF];               // Akima spline coefficients
    float  *ldatahint, *rdatahint;        // additional info about the
                                          // integration limits
    float  r[NLIM*NY];                    // integration results

    dfsIntegrCallBack le_cb, re_cb;       // integration call backs
    float  le_params, re_params;          // integration call backs parameters

    dfsSearchCellsCallBack  search_cb;    // search call back
    dfUserParams user_params = {1.0, 2.0, 3.0}; // Parameters passed to seach
                                          // callback, used during integration

    float  freq;
    float  r_ref;
    float  left_val[NY*(N-1)], right_val[NY*(N-1)];
    float  left_der1[NY*(N-1)], right_der1[NY*(N-1)];

    int i, j, errcode = 0;
    int yi;
    int errnums = 0;

    llim_p = llim;
    rlim_p = rlim;
    llim_cell_p = llim_cell;
    rlim_cell_p = rlim_cell;

    /***** Initializing parameters for Data Fitting task *****/

    /***** Parameters describing order and type of the spline *****/
    sorder       = DF_PP_CUBIC;
    stype        = DF_PP_AKIMA;

    /***** Parameters describing interpolation interval *****/
    nx           = N;
    xhint        = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny           = NY;
    yhint        = DF_NO_HINT;

    /***** Parameter describing additional info about spline
           coefficients *****/
    scoeffhint   = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type      = DF_BC_2ND_LEFT_DER | DF_BC_1ST_RIGHT_DER;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for Akima cubic spline */
    ic_type      = DF_NO_IC;
    ic = 0;

    /***** Parameters decsribing integration limits *****/
    nlim         = NLIM;
    llimhint     = DF_NO_HINT;
    rlimhint     = DF_NO_HINT;

    /***** Additional information about the structure of integration
           limits *****/
    /* No additional info is provided */
    ldatahint = 0;
    rdatahint = 0;

    /***** Parameter dascribing integration results storage format *****/
    rhint = DF_NO_HINT;

    /***** Generate partition with uniformly distributed break points *****/
    errcode = sUniformRandSortedData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Call backs for integration on the outer intervals *****/
    le_cb =  left_akima_integr;
    re_cb = right_akima_integr;
    /* Use limits of interpolation interval as call back parameters */
    le_params = left;
    re_params = right;

    /***** Search call back *****/
    search_cb = integr_search_cb;

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    for( yi=0; yi<ny; yi++ )
    {
        freq = (float )(yi+1)*FREQ;
        errcode = sSinDataNotUniformGrid( y+yi*nx, x, freq, nx );
        CheckDfError(errcode);
    }

    /***** Create Data Fitting task *****/
    errcode = dfsNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Search for left integration limits, save result to llim_cell *****/
    errcode = dfsSearchCells1D( task, DF_METHOD_STD, NLIM,llim,
        DF_NON_UNIFORM_PARTITION,NULL, llim_cell );
    CheckDfError(errcode);

    /***** Search for right integration limits, save result to rlim_cell *****/
    errcode = dfsSearchCells1D( task, DF_METHOD_STD, NLIM,rlim,
        DF_NON_UNIFORM_PARTITION,NULL, rlim_cell );
    CheckDfError(errcode);

    /***** Edit task parameters for Akima cubic spline construction *****/
    errcode = dfsEditPPSpline1D( task, sorder, stype, bc_type, bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct Akima cubic spline using STD method *****/
    errcode =  dfsConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Print results of integration limit searches *****/
    printf("Print results of integration limit searches:\n");
    for(i=0;i<NLIM;i++)
    {
        printf("i=%2d:",i);
        printf("    llim[i]=%9f",llim[i]);
        printf("    llim_cell[i]=%d",(int)llim_cell[i]);
        printf("    rlim[i]=%9f",rlim[i]);
        printf("    rlim_cell[i]=%d",(int)rlim_cell[i]);
        printf("\n");
    }


    /***** Compute integral for the spline on the interval (llim, rlim) *****/
    /***** and using search results *****/
    errcode = dfsIntegrateEx1D( task, DF_METHOD_PP, nlim, llim, llimhint,
                                rlim, rlimhint, ldatahint, rdatahint,
                                r, rhint, le_cb, &le_params, re_cb, &re_params,
                                0, 0, search_cb, &user_params );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/

    /***** Check spline values in break points *****/
    errcode = sCheckCubBreakPoints( nx, x, ny, y, scoeff, left_val, right_val );
    if ( errcode < 0 ) errnums++;

    /***** Check that spline 1st derivatives are equal for left
           and right piece of the spline for each break point *****/
    errcode = sCheckCub1stDerConsistency( nx, x, ny, scoeff,
                                          left_der1, right_der1 );
    if ( errcode < 0 ) errnums++;


    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given function *****/
    printf("\n i   x(i)");
    for( i = 0; i < ny; i++ )
    {
        printf("        y%d(i)",i);
    }
    printf("\n");

    for( j = 0; j < nx; j++ )
    {
        printf(" %1d %+10lf", j, x[j]);
        for( yi = 0; yi < ny; yi++ )
        {
            printf("   %+10lf", y[j+yi*nx]);
        }
        printf("\n");
    }


    /***** Print computed spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i      Ai              Bi              Ci              Di      ");
    printf("      Pi(x(i))      Pi(x(i+1))    Pi'(x(i))     Pi'(x(i+1))\n");
    for( yi = 0; yi < ny; yi++ )
    {
        printf("   for function y%d:\n",yi);
        for( j = 0; j < nx-1; j++ )
        {
            printf(" %1d %+13.6f   %+13.6f   %+13.6f   %+13.6f   %+11.6f   %+11.6f",
                    j,
                    scoeff[yi*(nx-1)*sorder + sorder*j + 0],
                    scoeff[yi*(nx-1)*sorder + sorder*j + 1],
                    scoeff[yi*(nx-1)*sorder + sorder*j + 2],
                    scoeff[yi*(nx-1)*sorder + sorder*j + 3],
                    right_val[yi*(nx-1) + j],
                    left_val [yi*(nx-1) + j]);
            printf("   %+11.6f   %+11.6f\n", right_der1[yi*(nx-1) + j], left_der1[yi*(nx-1) + j]);
        }
    }

    /***** Print computed integration results *****/
    printf("\nSpline-based integrals for %d functions on %d intervals:\n",(int)ny,NLIM);
    for( j = 0; j < NLIM; j++ )
    {
        printf("interval_%d=[ %4.1lf, %4.1lf ): ",j, llim[j], rlim[j] );
        for( yi = 0; yi < ny; yi++ )
        {
            printf("  integral_of_y%d=%lf",yi,r[j+yi*NLIM]);
        }
        printf("\n");
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0) {
        printf("\n\nError: Computed Akima cubic spline coefficients");
        printf(" or integartion results are incorrect\n");
        return 1;
    }
    else {
        printf("\n\nComputed Akima cubic spline coefficients");
        printf(" and integration results are correct\n");
    }

    return 0;
}

/*******************************************************************************
!   Integration call back for integral calculations on the interval
!   (-inf, x[0]).
!
! API
!   int  left_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
!                           MKL_INT64 rcell[], float  rlim[], float  r[],
!                           void *x0,
!                           dfIntegrCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  lcell  - array of size n with indices of the cells that contain
!           the left-side integration limits
!  llim   - array of size n that holds the left-side integration limits
!  rcell  - array of size n with indices of the cells that contain
!           the right-side integration limits
!  rlim   - array of size n that holds the right-side integration limits
!  x0     - left  limit of interpolation interval
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int left_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
                       MKL_INT64 rcell[], float  rlim[], float  r[], void *x0,
                       dfIntegrCallBackLibraryParams* library_params )
{
    MKL_INT64 i;
    float  *x = (float *)x0;

    for ( i = 0; i < n[0]; i++ )
    {
        r[i] = x[0] * ( x[0] - llim[i] );
    }

    return 0;
}

/*******************************************************************************
!   Integration call back for integral calculations on the interval
!   [x[N-1], +inf).
!
! API
!   int right_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
!                           MKL_INT64 rcell[], float  rlim[], float  r[],
!                           void *xN,
!                           dfIntegrCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  lcell  - array of size n with indices of the cells that contain
!           the left-side integration limits
!  llim   - array of size n that holds the left-side integration limits
!  rcell  - array of size n with indices of the cells that contain
!           the right-side integration limits
!  rlim   - array of size n that holds the right-side integration limits
!  xN     - right limit of interpolation interval
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int right_akima_integr( MKL_INT64* n, MKL_INT64 lcell[], float  llim[],
                        MKL_INT64 rcell[], float  rlim[], float  r[], void *xN,
                        dfIntegrCallBackLibraryParams* library_params )
{
    MKL_INT64 i;
    float  *x = (float *)xN;
    for ( i = 0; i < n[0]; i++ )
    {
        r[i] = x[0] * ( rlim[i] - x[0] );
    }

    return 0;
}

/*******************************************************************************
!   Definition of the search call back, used by integration routine
! for search of left and right integration limits.
!
! API
!   int integr_search_cb( MKL_INT64* n, float  site[], MKL_INT64 cell[],
!                     int flag[], void* user_params,
!                     dfSearchCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of pairs of integration limits
!  site   - array of size n with integration limits (left or right,
!           depending on value of library_params->limit_type_flag)
!  user_params - pointer to user-defined parameters
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  cell   - array of size n with indices of the cells that contain
!           corresponding integration limits in site array
!           (calculated by this search call back)
!  flag   - array of size n filled by ones here to indicate
!           that each corresponding integration limit was found
!
! RETURN VALUE:
!  The status returned by the callback function:
!   DF_STATUS_EXACT_RESULT - indicates successful completion of the callback
!   <0 - error
*******************************************************************************/
int integr_search_cb( MKL_INT64* n, float  site[], MKL_INT64 cell[],
                      int flag[], void* user_params,
                      dfSearchCallBackLibraryParams* library_params )
{
    int lim_type;
    int xxi,j;
    int nlim1=*n;
    int res;

    /* Check if user_params passed here correctly */
    if( (((dfUserParams*)user_params)->v0!=1.0) ||
        (((dfUserParams*)user_params)->v1!=2.0) ||
        (((dfUserParams*)user_params)->v2!=3.0) )
    {
        printf("Error: one of user_params is incorrect.\n");
        res = -1;
    }
    else
    {
        if(library_params==NULL)
        {
            printf("Error: library_params is NULL.\n");
            res = -1;
        }
        else
        {
            lim_type = library_params->limit_type_flag;

            if(lim_type==DF_INTEGR_SEARCH_CB_LLIM_FLAG)
            {
                for (xxi = 0; xxi < nlim1; xxi++) {
                    for (j = 0; j < NLIM; j++) {
                        if ( site[xxi] == llim_p[j] ) { break;}
                    }
                    cell[xxi] = llim_cell_p[j];
                    flag[xxi] = 1;
                }
            }
            else /*Here if(lim_type==DF_INTEGR_SEARCH_CB_RLIM_FLAG)*/
            {
                for (xxi = 0; xxi < nlim1; xxi++) {
                    for (j = 0; j < NLIM; j++) {
                        if ( site[xxi] == rlim_p[j] ) { break;}
                    }
                    cell[xxi] = rlim_cell_p[j];
                    flag[xxi] = 1;
                }
            }
            res = DF_STATUS_EXACT_RESULT;
        }
    }

    return res;
}
