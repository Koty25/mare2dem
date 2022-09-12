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
!    Construction of co-monotone cubic spline, accordingly Hyman algorithm
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N       11  // number of break points
#define NY      1   // number of functions
#define NDORDER 2   // size of array describing derivative orders to compute
#define NDER    2   // number of derivatives to compute
#define NSCOEFF NY*(N-1)*DF_PP_CUBIC // total number of spline coefficients
#define NSITE   10  // total number of interpolation sites

int main()
{
    DFTaskPtr interpTask = 0;   // Data Fitting task descriptor
    MKL_INT sorder;             // spline order
    MKL_INT stype;              // spline type
    MKL_INT nx;                 // number of break points
    MKL_INT xhint;              // additional info about break points
    MKL_INT ny;                 // number of functions
    MKL_INT yhint;              // additional info about functions
    MKL_INT nscoeff;            // number of spline coefficients
    MKL_INT scoeffhint;         // additional info about spline
                                // coefficients
    MKL_INT bc_type;            // boundary conditions type
    MKL_INT ic_type;            // internal conditions type
    MKL_INT nsite;              // total number of interpolation sites
    MKL_INT sitehint;           // additional info about interpolation
                                // sites
    MKL_INT ndorder;            // size of array describing derivative
                                // orders
    MKL_INT dorder[] = {1,1};   // function values and 1st derivatives
                                // will be computed
    MKL_INT rhint;              // interpolation results storage format
    MKL_INT *cell_idx;          // indices of cells containing
                                // interpolation sites
    double y[N*NY]={10,10,10,10,10,10,10.5,15,50,60,85}; // function values
    double x[N   ]={0 ,2 ,3 ,5 ,6 ,8 ,9   ,11,12,14,15}; // breakpoints
    double *ic;                 // array of internal conditions
    double bc[2];               // array of boundary conditions
    double scoeff[NSCOEFF];     // array of spline coefficients
    double site[NSITE];         // array of interpolation sites
    double r[NDER*NSITE*NY];    // spline evaluation results
    double *datahint;           // additional info about structure
                                // of arrays x and y
    int xi,xxi,errcode;
    int comonot_status,comonot_status_1interval;
    double slope,di,dip1,alpha,beta;

    /***** Initializing parameters for Data Fitting task *****/

    sorder     = DF_PP_CUBIC;
    stype      = DF_PP_HYMAN;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_MATRIX_STORAGE_ROWS;

    /***** Parameters describing spline coefficients storage *****/
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type    = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;
    bc[0]      = 0.0;   /* 1st derivative on left boundary */
    bc[1]      = 0.0;   /* 1st derivative on right boundary */

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing interpolation sites *****/
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


    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D(&interpTask, nx,x,xhint, ny,y,yhint);
    CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/
    errcode = dfdEditPPSpline1D(interpTask, sorder, stype, bc_type, bc,
        ic_type, ic, scoeff, scoeffhint);
    CheckDfError(errcode);

    /***** Construct co-monotone cubic spline *****/
    errcode = dfdConstruct1D(interpTask, DF_PP_SPLINE, DF_METHOD_STD);
    CheckDfError(errcode);

    /***** Co-monotonicity test *****/
    /*
    //  Test accordingly Fritch-Carlson criteria described in
    //  Fritsch, F. N.; Carlson, R. E. (1980). "Monotone Piecewise Cubic
    //  Interpolation". SIAM Journal on Numerical Analysis
    //  (SIAM) 17 (2): 238-246.
    */
    comonot_status=1;
    for(xi=0;xi<=nx-2;xi++) /* For each sub-interval */
    {
        printf("Interval %d, x=[%f;%f]: ",xi,x[xi],x[xi+1]);

        /* Interpolate on current interval (and get 1st derivatives on its ends) */
        for(xxi=0;xxi<nsite;xxi++) {
            site[xxi] = x[xi] + (x[xi+1]-x[xi])/(double)(NSITE-1)*(double)xxi;
        }
        errcode = dfdInterpolate1D(interpTask, DF_INTERP, DF_METHOD_PP,
            nsite,site,sitehint, ndorder,dorder,datahint, r,rhint, cell_idx);
        CheckDfError(errcode);


        comonot_status_1interval=0;

        /***** Co-monotonicity test on current sub-interval started *****/
        /* Slope of the input data on current sub-interval */
        slope = (y[xi+1]-y[xi]) / (x[xi+1]-x[xi]);
        /* 1st derivative on left  side of current sub-interval */
        di    = r[        0*NDORDER+1];
        /* 1st derivative on right side of current sub-interval */
        dip1  = r[(nsite-1)*NDORDER+1];

        /* Using Fritsch-Carlson criteria */
        if(slope==0.0)
        {
            /* Here, when input function is constant */
            if( (di==0.0) && (dip1==0.0) )
            {
                comonot_status_1interval=1; /* Ok, co-monotone. */
            }
        }
        else
        {
            /* Here, when input function is increasing or decreasing */
            alpha = di   / slope;
            beta  = dip1 / slope;

            if( ((di>=0.0) && (dip1>=0.0) && (slope>=0.0))   ||
                ((di<=0.0) && (dip1<=0.0) && (slope<=0.0)) )
            {
                if( ((alpha>=0.0) && (alpha<=3.0)    &&
                    (beta >=0.0) && (beta<=3.0))     ||
                    ((2.0*alpha+beta-3.0)*(2.0*alpha+beta-3.0)
                     <= 3.0*alpha*(alpha+beta-2.0)) )
                {
                    comonot_status_1interval=1; /* Ok, co-monotone. */
                }
            }
        }
        /***** Co-monotonicity test on current sub-interval finished *****/

        if(comonot_status_1interval)
        {
            printf(" spline is co-monotone here.\n");
        }
        else
        {
            printf(" spline is not co-monotone here.\n");
            printf("  slope=%e, di=%e, dip1=%e\n",slope,di,dip1);
            printf("  alpha=%e, beta=%e\n",alpha,beta);

            for(xxi=0;xxi<nsite;xxi++) {
                printf("  site[%d]=%e,",xxi,site[xxi]);
                printf(" func=%e," ,r[xxi*NDORDER+0]);
                printf(" deriv=%e,",r[xxi*NDORDER+1]);
                printf("\n");
            }
        }
        comonot_status &= comonot_status_1interval;
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask(&interpTask);
    CheckDfError(errcode);

    /***** Print summary of co-monotonicity test *****/
    if (comonot_status == 1)
    {
        printf("\n\nConstructed spline is co-monotone with input data.\n");
        return 0;
    }
    else
    {
        printf("\n\nError: constructed spline is not co-monotone with input data.\n");
        return 1;
    }
}
