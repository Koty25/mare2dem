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

#include <stdlib.h>

void CheckDfError(int num)
{
    switch(num)
    {
    case DF_STATUS_OK: {
        break;
    }
    case DF_ERROR_NULL_TASK_DESCRIPTOR:
    {
        printf( "Error: null task descriptor (code %d).\n", num );
        break;
    }
    case DF_ERROR_MEM_FAILURE:
    {
        printf( "Error: memory allocation failure in DF functionality (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_NX:
    {
        printf( "Error: the number of breakpoints is invalid (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_X:
    {
        printf( "Error: the array which contains the breakpoints is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_X_HINT:
    {
        printf( "Error: invalid flag describing structure of partition (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_NY:
    {
        printf( "Error: invalid dimension of vector-valued function y (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_Y:
    {
        printf( "Error: the array which contains function values is invalid (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_Y_HINT:
    {
        printf( "Error: invalid flag describing structure of function y (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_SPLINE_ORDER:
    {
        printf( "Error: invalid spline order (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_SPLINE_TYPE:
    {
        printf( "Error: invalid type of the spline (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_IC_TYPE:
    {
        printf( "Error: invalid type of internal conditions used in the spline construction (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_IC:
    {
        printf( "Error: array of internal conditions for spline construction is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_BC_TYPE:
    {
        printf( "Error: invalid type of boundary conditions used in the spline construction (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_BC:
    {
        printf( "Error: array which presents boundary conditions for spline construction is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_PP_COEFF:
    {
        printf( "Error: array of piece-wise polynomial spline coefficients is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_PP_COEFF_HINT:
    {
        printf( "Error: invalid flag describing structure of the piece-wise polynomial spline coefficients (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_PERIODIC_VAL:
    {
        printf( "Error: function values at the end points of the interpolation interval are not equal as required in periodic boundary conditions (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_DATA_ATTR:
    {
        printf( "Error: invalid attribute of the pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_DATA_IDX:
    {
        printf( "Error: index of pointer to be set or modified in Data Fitting task descriptor with EditIdxPtr editor is out of range (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_NSITE:
    {
        printf( "Error: invalid number of interpolation sites (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_SITE:
    {
        printf( "Error: array of interpolation sites is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_SITE_HINT:
    {
        printf( "Error: invalid flag describing structure of interpolation sites (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_NDORDER:
    {
        printf( "Error: invalid size of array that defines order of the derivatives to be computed at the interpolation sites (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_DORDER:
    {
        printf( "Error: array defining derivative orders to be computed at interpolation sites is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_DATA_HINT:
    {
        printf( "Error: invalid flag providing a-priori information about partition and/or interpolation sites (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_INTERP:
    {
        printf( "Error: array of spline based interpolation results is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_INTERP_HINT:
    {
        printf( "Error: invalid flag defining structure of spline based interpolation results (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_CELL_IDX:
    {
        printf( "Error: array of indices of partition cells containing interpolation sites is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_NLIM:
    {
        printf( "Error: invalid size of arrays containing integration limits (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_LLIM:
    {
        printf( "Error: array of left integration limits is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_RLIM:
    {
        printf( "Error: array of right integration limits is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_INTEGR:
    {
        printf( "Error: array of spline based integration results is not defined (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_INTEGR_HINT:
    {
        printf( "Error: invalid flag defining structure of spline based integration results (code %d).\n", num );
        break;
    }
    case DF_ERROR_BAD_LOOKUP_INTERP_SITE:
    {
        printf( "Error: bad site provided for interpolation with look-up interpolator (code %d).\n", num );
        break;
    }
    case DF_ERROR_NULL_PTR:
    {
        printf( "Error: bad pointer provided in DF function (code %d).\n", num );
        break;
    }
    default: break;
    }

    if(num < 0) {
       exit(1);
    }
}
