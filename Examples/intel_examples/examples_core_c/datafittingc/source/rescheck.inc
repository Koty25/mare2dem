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

#define EPSILON_SINGLE      1e-02
#define EPSILON_DOUBLE      1e-13

#define DF_ABS(x)   (((x) > 0) ? (x) : -(x))


/*
//  Compute single precision piecewise polynomial spline derivative
//  in a given point x
*/
static float sComputePPSplineDer( MKL_INT der, MKL_INT n, MKL_INT piece_idx,
                                  float x, float x_i, MKL_INT data_set,
                                  MKL_INT sorder, float sppcoeff[] )
{
    int i, k;
    int mult;
    float *coeff;
    double der_val;
    double t, xpow;

    coeff = sppcoeff + data_set * sorder * ( n - 1 ) + sorder*piece_idx;
    t = (double)x - (double)x_i;
    xpow = 1.0;
    der_val = 0.0;
    for ( i = der; i < sorder; i++ )
    {
        mult = 1;
        for ( k = i; k > i-der; k-- ) { mult *= k; }

        der_val += (double)mult * (double)coeff[i] * xpow;
        xpow *= t;
    }

    return (float)der_val;
}

/*
//  Compute double precision piecewise polynomial spline derivative
//  in a given point x
*/
static double dComputePPSplineDer( MKL_INT der, MKL_INT n, MKL_INT piece_idx,
                                  double x, double x_i, MKL_INT data_set,
                                  MKL_INT sorder, double sppcoeff[] )
{
    int i, k;
    int mult;
    double *coeff;
    double der_val;
    double t, xpow;

    coeff = sppcoeff + data_set * sorder * ( n - 1 ) + sorder*piece_idx;
    t = x - x_i;
    xpow = 1.0;
    der_val = 0.0;
    for ( i = der; i < sorder; i++ )
    {
        mult = 1;
        for ( k = i; k > i-der; k-- ) { mult *= k; }

        der_val += (double)mult * coeff[i] * xpow;
        xpow *= t;
    }

    return der_val;
}

/*
//++
//  LINEAR SPLINE RELATED CHECKS
//--
*/

/*
// Check single precision linear spline coefficients
// by testing spline values in break points
*/
static int sCheckLinSplineBreakPoints( MKL_INT n, float x[],
                                       MKL_INT ny, float y[], MKL_INT sorder,
                                       float sppcoeff[] )
{
    int i, j, errcode = 0;
    float *coeff;
    float spline_val_right, spline_val_left, dif;

    for ( i = 0; i < ny; i++ )
    {
        coeff = sppcoeff + i * sorder * ( n - 1 );
        spline_val_left  = ( coeff[0] + coeff[1] * ( x[0] - x[0] ) );
        dif        = DF_ABS( spline_val_left - y[i*n + 0] );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
        for ( j = 0; j < n-2; j++ )
        {

            spline_val_right = ( coeff[sorder*j    +0] + coeff[sorder*j    +1] * ( x[j+1] - x[j] ) );
            spline_val_left  = ( coeff[sorder*(j+1)+0] + coeff[sorder*(j+1)+1] * ( x[j+1] - x[j+1] ) );
            dif        = DF_ABS( spline_val_right - spline_val_left );
            if ( dif > EPSILON_SINGLE )
            {
                errcode = -1;
                return errcode;
            }
        }
        spline_val_right = coeff[sorder*(n-2)+0] + coeff[sorder*(n-2)+1] * ( x[n-1] - x[n-2] );
        dif        = DF_ABS( spline_val_right - y[i*n + n-1] );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check double precision linear spline coefficients
// by testing spline values in break points
*/
static int dCheckLinSplineBreakPoints( MKL_INT n, double x[],
                                       MKL_INT ny, double y[], MKL_INT sorder,
                                       double sppcoeff[] )
{
    int i, j, errcode = 0;
    double *coeff;
    double spline_val_right, spline_val_left, dif;

    for ( i = 0; i < ny; i++ )
    {
        coeff = sppcoeff + i * sorder * ( n - 1 );
        spline_val_left  = ( coeff[0] + coeff[1] * ( x[0] - x[0] ) );
        dif        = DF_ABS( spline_val_left - y[i*n + 0] );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
        for ( j = 0; j < n-2; j++ )
        {

            spline_val_right = ( coeff[sorder*j    +0] + coeff[sorder*j    +1] * ( x[j+1] - x[j] ) );
            spline_val_left  = ( coeff[sorder*(j+1)+0] + coeff[sorder*(j+1)+1] * ( x[j+1] - x[j+1] ) );
            dif        = DF_ABS( spline_val_right - spline_val_left );
            if ( dif > EPSILON_SINGLE )
            {
                errcode = -1;
                return errcode;
            }
        }
        spline_val_right = coeff[sorder*(n-2)+0] + coeff[sorder*(n-2)+1] * ( x[n-1] - x[n-2] );
        dif        = DF_ABS( spline_val_right - y[i*n + n-1] );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
//++
//  CUBIC SPLINE RELATED CHECKS
//--
*/

/*
// Definition of auxiliary functions for cubic spline coefficients testing
*/
static int sCheckCub1stDerBC ( MKL_INT n, float  x[], MKL_INT ny, MKL_INT sorder, float  sppcoeff[], float  bc, int isLeft );
static int sCheckCub2ndDerBC ( MKL_INT n, float  x[], MKL_INT ny, MKL_INT sorder, float  sppcoeff[], float  bc, int isLeft );

static int dCheckCub1stDerBC ( MKL_INT n, double x[], MKL_INT ny, MKL_INT sorder, double sppcoeff[], double bc, int isLeft );
static int dCheckCub2ndDerBC ( MKL_INT n, double x[], MKL_INT ny, MKL_INT sorder, double sppcoeff[], double bc, int isLeft );

/*
// Check that spline value consist with given tabular function
// for each break-point.
*/
static int sCheckCubBreakPoints( MKL_INT n, float x[], MKL_INT ny, float y[],
                                 float sppcoeff[], float *left, float *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n-1; j++ )
        {
            right[i*(n-1) + j] = sComputePPSplineDer( 0, n, j, x[j], x[j], i,
                                                      sorder, sppcoeff );
            dif = DF_ABS( right[i*(n-1) + j]  - y[i * n + j+0] );
            if ( dif > EPSILON_SINGLE ) errnums++;

            left[i*(n-1) + j] = sComputePPSplineDer( 0, n, j, x[j+1], x[j], i,
                                                      sorder, sppcoeff );
            dif = DF_ABS( left[i*(n-1) + j] - y[i * n + j+1] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that spline value consist with given tabular function
// for each break-point.
*/
static int dCheckCubBreakPoints( MKL_INT n, double x[], MKL_INT ny, double y[],
                                 double sppcoeff[], double *left, double *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n-1; j++ )
        {
            right[i*(n-1) + j] = dComputePPSplineDer( 0, n, j, x[j], x[j], i,
                                                      sorder, sppcoeff );
            dif = DF_ABS( right[i*(n-1) + j]  - y[i * n + j+0] );
            if ( dif > EPSILON_SINGLE ) errnums++;

            left[i*(n-1) + j] = dComputePPSplineDer( 0, n, j, x[j+1], x[j], i,
                                                     sorder, sppcoeff );
            dif = DF_ABS( left[i*(n-1) + j]  - y[i * n + j+1] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }
    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
//  Check that spline 1st derivatives are the same for the left piece
//  of spline and the right piece of spline for each break-point
*/
static int sCheckCub1stDerConsistency( MKL_INT n, float x[], MKL_INT ny,
                                       float sppcoeff[],
                                       float *left, float *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = sComputePPSplineDer(1, n, 0, x[0], x[0], i, sorder,
                                                 sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = sComputePPSplineDer(1, n, j-1, x[j], x[j-1],
                                                      i, sorder, sppcoeff);
            right[i*(n-1) + j] = sComputePPSplineDer(1, n, j, x[j], x[j], i,
                                                     sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = sComputePPSplineDer(1, n, n-2, x[n-1], x[n-2],
                                                  i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
//  Check that spline 1st derivatives are the same for the left piece
//  of spline and the right piece of spline for each break-point
*/
static int dCheckCub1stDerConsistency( MKL_INT n, double x[], MKL_INT ny,
                                       double sppcoeff[],
                                       double *left, double *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = dComputePPSplineDer(1, n, 0, x[0], x[0], i, sorder,
                                               sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = dComputePPSplineDer(1, n, j-1, x[j], x[j-1],
                                                      i, sorder, sppcoeff);
            right[i*(n-1) + j] = dComputePPSplineDer(1, n, j, x[j], x[j],
                                                     i, sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = dComputePPSplineDer(1, n, n-2, x[n-1], x[n-2],
                                                  i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
//  Check that spline 2nd derivatives are the same for the left piece
//  of spline and the right piece of spline for each break-point
*/
static int sCheckCub2ndDerConsistency( MKL_INT n, float x[], MKL_INT ny,
                                       float sppcoeff[],
                                       float *left, float *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = sComputePPSplineDer(2, n, 0, x[0], x[0], i, sorder,
                                                 sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = sComputePPSplineDer(2, n, j-1, x[j], x[j-1],
                                                      i, sorder, sppcoeff);
            right[i*(n-1) + j] = sComputePPSplineDer(2, n, j, x[j], x[j],
                                                     i, sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = sComputePPSplineDer(2, n, n-2, x[n-1], x[n-2],
                                                  i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
//  Check that spline 2nd derivatives are the same for the left piece
//  of spline and the right piece of spline for each break-point
*/
static int dCheckCub2ndDerConsistency( MKL_INT n, double x[], MKL_INT ny,
                                       double sppcoeff[],
                                       double *left, double *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = dComputePPSplineDer(2, n, 0, x[0], x[0], i, sorder,
                                                 sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = dComputePPSplineDer(2, n, j-1, x[j], x[j-1],
                                                      i, sorder, sppcoeff);
            right[i*(n-1) + j] = dComputePPSplineDer(2, n, j, x[j], x[j],
                                                     i, sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = dComputePPSplineDer(2, n, n-2, x[n-1], x[n-2],
                                                  i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that spline coefficienst satisfy given boundary conditions
// in single precision case
*/
static int sCheckCubBC( MKL_INT n, float x[], MKL_INT ny,
                        float sppcoeff[], MKL_INT bc_type, float bc[] )
{
    int errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;

    if ( bc_type == DF_BC_NOT_A_KNOT ) return errcode;

    if ( (bc_type & DF_BC_FREE_END) == DF_BC_FREE_END )
    {
        errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0, 1 );
        if ( errcode < 0 ) return errcode;
        errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0, 0 );
    }
    else
    {
        if ( (bc_type & DF_BC_1ST_LEFT_DER)  == DF_BC_1ST_LEFT_DER )
        {
            errcode = sCheckCub1stDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_1ST_RIGHT_DER) == DF_BC_1ST_RIGHT_DER )
        {
            errcode = sCheckCub1stDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_LEFT_DER)  == DF_BC_2ND_LEFT_DER )
        {
            errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_RIGHT_DER) == DF_BC_2ND_RIGHT_DER )
        {
            errcode = sCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
    }

    return errcode;
}

/*
// Check that spline coefficienst satisfy given boundary conditions
// in double precision case
*/
static int dCheckCubBC( MKL_INT n, double x[], MKL_INT ny,
                        double sppcoeff[], MKL_INT bc_type, double bc[] )
{
    int errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;

    if ( bc_type == DF_BC_NOT_A_KNOT ) return errcode;

    if ( (bc_type & DF_BC_FREE_END) == DF_BC_FREE_END )
    {
        errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0, 1 );
        if ( errcode < 0 ) return errcode;
        errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, 0.0, 0 );
    }
    else
    {
        if ( (bc_type & DF_BC_1ST_LEFT_DER)  == DF_BC_1ST_LEFT_DER )
        {
            errcode = dCheckCub1stDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_1ST_RIGHT_DER) == DF_BC_1ST_RIGHT_DER )
        {
            errcode = dCheckCub1stDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_LEFT_DER)  == DF_BC_2ND_LEFT_DER )
        {
            errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_RIGHT_DER) == DF_BC_2ND_RIGHT_DER )
        {
            errcode = dCheckCub2ndDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
    }

    return errcode;
}


/*
// Compute single precision integral for a given spline piece
*/
static float sComputeCubSplineIntegr( MKL_INT n, MKL_INT piece_idx,
                                      float llim, float rlim,
                                      float x_i, MKL_INT data_set,
                                      MKL_INT sorder, float sppcoeff[] )
{
    float *coeff;
    float coeff0, coeff1, coeff2, coeff3;
    float integr_val;
    float t;

    coeff = sppcoeff + data_set * sorder * ( n - 1 );

    /* Get spline coefficients */
    coeff0 = coeff[sorder*piece_idx + 0];
    coeff1 = coeff[sorder*piece_idx + 1];
    coeff2 = coeff[sorder*piece_idx + 2];
    coeff3 = coeff[sorder*piece_idx + 3];

    /* Compute integartion coefficients */
    coeff1 *= 0.5;
    coeff2 /= 3.0;
    coeff3 *= 0.25;

    /* Spline value computation */
    t = rlim - x_i;
    integr_val  = t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)));

    t = llim - x_i;
    integr_val -= t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)));

    return integr_val;
}


/*
// Compute double precision integral for a given spline piece
*/
static double dComputeCubSplineIntegr( MKL_INT n, MKL_INT piece_idx,
                                       double llim, double rlim,
                                       double x_i, MKL_INT data_set,
                                       MKL_INT sorder, double sppcoeff[] )
{
    double *coeff;
    double coeff0, coeff1, coeff2, coeff3;
    double integr_val;
    double t;

    coeff = sppcoeff + data_set * sorder * ( n - 1 );

    /* Get spline coefficients */
    coeff0 = coeff[sorder*piece_idx + 0];
    coeff1 = coeff[sorder*piece_idx + 1];
    coeff2 = coeff[sorder*piece_idx + 2];
    coeff3 = coeff[sorder*piece_idx + 3];

    /* Compute integartion coefficients */
    coeff1 *= 0.5;
    coeff2 /= 3.0;
    coeff3 *= 0.25;

    /* Spline value computation */
    t = rlim - x_i;
    integr_val  = t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)));

    t = llim - x_i;
    integr_val -= t*(coeff0 + t*(coeff1 + t*(coeff2 + t*coeff3)));

    return integr_val;
}


/*
// Check 1st derivatrives boundary conditions of cubic splines
// in single precision case
*/
static int sCheckCub1stDerBC( MKL_INT n, float x[], MKL_INT ny,
                              MKL_INT sorder, float sppcoeff[],
                              float bc, int isLeft )
{
    int i, piece_idx, x_idx;
    float spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-2) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = sComputePPSplineDer(1, n, piece_idx, x[x_idx],
            x[piece_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 1st derivatrives boundary conditions of cubic splines
// in double precision case
*/
static int dCheckCub1stDerBC( MKL_INT n, double x[], MKL_INT ny,
                              MKL_INT sorder, double sppcoeff[],
                              double bc, int isLeft )
{
    int i, piece_idx, x_idx;
    double spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-2) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = dComputePPSplineDer(1, n, piece_idx, x[x_idx],
            x[piece_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_DOUBLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 2nd derivatrives boundary conditions in single precision case
*/
static int sCheckCub2ndDerBC( MKL_INT n, float x[], MKL_INT ny,
                              MKL_INT sorder, float sppcoeff[],
                              float bc, int isLeft )
{
    int i, piece_idx, x_idx;
    float spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-2) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = sComputePPSplineDer(2, n, piece_idx, x[x_idx],
             x[piece_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 2nd derivatrives boundary conditions in double precision case
*/
static int dCheckCub2ndDerBC( MKL_INT n, double x[], MKL_INT ny,
                              MKL_INT sorder, double sppcoeff[],
                              double bc, int isLeft )
{
    int i, piece_idx, x_idx;
    double spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-2) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = dComputePPSplineDer(2, n, piece_idx, x[x_idx],
             x[piece_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_DOUBLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
//++
//  QUADRATIC SPLINE RELATED CHECKS
//--
*/

/*
// Definition of auxiliary functions for quadratic spline coefficients testing
*/
static int    sCheckQuad1stDerBC    ( MKL_INT n, float  x[], MKL_INT ny, MKL_INT sorder, float  sppcoeff[], float  bc, int isLeft );
static int    sCheckQuad2ndDerBC    ( MKL_INT n, float  x[], MKL_INT ny, MKL_INT sorder, float  sppcoeff[], float  bc, int isLeft );

static int    dCheckQuad1stDerBC    ( MKL_INT n, double x[], MKL_INT ny, MKL_INT sorder, double sppcoeff[], double bc, int isLeft );
static int    dCheckQuad2ndDerBC    ( MKL_INT n, double x[], MKL_INT ny, MKL_INT sorder, double sppcoeff[], double bc, int isLeft );

/*
// Check that quadratic spline value consist with given tabular function
// for each break-point.
*/
static int sCheckQuadBreakPoints( MKL_INT n, float x[], MKL_INT ny, float y[],
                                  float sppcoeff[], float *left, float *right)
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n-1; j++ )
        {
            right[i*(n-1) + j] = sComputePPSplineDer( 0, n, j, x[j], x[j],
                                                      i, sorder, sppcoeff );
            dif = DF_ABS( right[i*(n-1) + j]  - y[i * n + j+0] );
            if ( dif > EPSILON_SINGLE ) errnums++;

            left[i*(n-1) + j] = sComputePPSplineDer( 0, n, j, x[j+1], x[j],
                                                     i, sorder, sppcoeff );
            dif = DF_ABS( left[i*(n-1) + j]  - y[i * n + j+1] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic spline value consist with given tabular function
// for each break-point.
*/
static int dCheckQuadBreakPoints( MKL_INT n, double x[], MKL_INT ny, double y[],
                                  double sppcoeff[], double *left, double *right)
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n-1; j++ )
        {
            right[i*(n-1) + j] = dComputePPSplineDer( 0, n, j, x[j], x[j],
                                                      i, sorder, sppcoeff );
            dif = DF_ABS( right[i*(n-1) + j]  - y[i * n + j+0] );
            if ( dif > EPSILON_SINGLE ) errnums++;

            left[i*(n-1) + j] = dComputePPSplineDer( 0, n, j, x[j+1], x[j],
                                                     i, sorder, sppcoeff );
            dif = DF_ABS( left[i*(n-1) + j]  - y[i * n + j+1] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic spline values are consistent in each q-node.
*/
static int sCheckSubbQNodes( MKL_INT nx, float x[], MKL_INT nqnode,
                             float qnode[], MKL_INT ny,
                             float sppcoeff[],
                             float *left, float *right)
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(nqnode-1) + 0] = sComputePPSplineDer(0, nqnode, 0,
                                qnode[0], x[0], i, sorder, sppcoeff);

        for ( j = 1; j < nqnode-1; j++ )
        {
            left[i*(nqnode-1) + j-1] = sComputePPSplineDer(0, nqnode, j-1,
                                qnode[j], x[j-1], i, sorder, sppcoeff);
            right[i*(nqnode-1) + j]  = sComputePPSplineDer(0, nqnode, j,
                                qnode[j], x[j], i, sorder, sppcoeff);
            dif = DF_ABS(left[i*(nqnode-1) + j-1] - right[i*(nqnode-1) + j]);
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(nqnode-1) + nqnode-2] = sComputePPSplineDer(0, nqnode,
                nqnode-2, qnode[nqnode-1],  x[nqnode-2], i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic spline values are consistent in each q-node.
*/
static int dCheckSubbQNodes( MKL_INT nx, double x[], MKL_INT nqnode,
                             double qnode[], MKL_INT ny,
                             double sppcoeff[],
                             double *left, double *right)
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(nqnode-1) + 0] = dComputePPSplineDer(0, nqnode, 0,
                                qnode[0], x[0], i, sorder, sppcoeff);

        for ( j = 1; j < nqnode-1; j++ )
        {
                        left[i*(nqnode-1) + j-1] = dComputePPSplineDer(0, nqnode,
                                 j-1, qnode[j], x[j-1], i, sorder, sppcoeff);
            right[i*(nqnode-1) + j]  = dComputePPSplineDer(0, nqnode, j,
                                qnode[j], x[j], i, sorder, sppcoeff);
            dif = DF_ABS(left[i*(nqnode-1) + j-1] - right[i*(nqnode-1) + j]);
                        if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(nqnode-1) + nqnode-2] = dComputePPSplineDer(0, nqnode,
                nqnode-2, qnode[nqnode-1], x[nqnode-2], i, sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}


/*
// Check that Subbotin spline value consist with given tabular function
// for each break-point.
*/
static int sCheckSubbBreakPoints( MKL_INT n, float x[], MKL_INT ny, float y[],
                                  float sppcoeff[], float *spline_val )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            spline_val[i*n + j] = sComputePPSplineDer( 0, n+1, j, x[j], x[j], i,
                                                       sorder, sppcoeff );
            dif = DF_ABS( spline_val[i*n + j]  - y[i*n + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that Subbotin spline value consist with given tabular function
// for each break-point.
*/
static int dCheckSubbBreakPoints( MKL_INT n, double x[], MKL_INT ny, double y[],
                                  double sppcoeff[], double *spline_val )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            spline_val[i*n + j] = dComputePPSplineDer( 0, n+1, j, x[j],  x[j], i,
                                                       sorder, sppcoeff );
            dif = DF_ABS( spline_val[i*n + j]  - y[i*n + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic spline 1st derivatives are the same
// for the left piece of spline and the right piece of spline
// for each break-point
*/
static int sCheckQuad1stDerConsistency( MKL_INT n, float x[], MKL_INT ny,
                                        float sppcoeff[],
                                        float *left, float *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = sComputePPSplineDer(1, n, 0, x[0], x[0], i, sorder,
                                                sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = sComputePPSplineDer(1, n, j-1, x[j], x[j-1],
                                             i, sorder, sppcoeff);
            right[i*(n-1) + j]  = sComputePPSplineDer(1, n, j, x[j], x[j],
                                             i, sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = sComputePPSplineDer(1, n, n-2, x[n-1], x[n-2], i,
                                            sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic spline 1st derivatives are the same
// for the left piece of spline and the right piece of spline
// for each break-point
*/
static int dCheckQuad1stDerConsistency( MKL_INT n, double x[], MKL_INT ny,
                                        double sppcoeff[],
                                        double *left, double *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = dComputePPSplineDer(1, n, 0, x[0], x[0], i, sorder,
                                                sppcoeff);

        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = dComputePPSplineDer(1, n, j-1, x[j], x[j-1],
                                            i, sorder, sppcoeff);
            right[i*(n-1) + j]  = dComputePPSplineDer(1, n, j, x[j], x[j],
                                            i, sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }

        left[i*(n-1) + n-2] = dComputePPSplineDer(1, n, n-2, x[n-1], x[n-2], i,
                                            sorder, sppcoeff);
    }

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Check that quadratic Subbotin spline 1st derivatives are the same
// for the left piece of spline and the right piece of spline
// for each break-point
*/
static int sCheckQuadSubb1stDerConsistency( MKL_INT n, float x[], float qnode[],
                                            MKL_INT ny, float sppcoeff[],
                                            float *left, float *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    float dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = sComputePPSplineDer(1, n, 0, qnode[0], x[0], i, sorder,
                                                sppcoeff);
        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = sComputePPSplineDer(1, n, j-1, qnode[j], x[j-1], i,
                                                     sorder, sppcoeff);
            right[i*(n-1) + j]  = sComputePPSplineDer(1, n, j, qnode[j], x[j], i,
                                                     sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
        left[i*(n-1) + n-2] = sComputePPSplineDer(1, n, n-2, qnode[n-1], x[n-2], i,
                                                 sorder, sppcoeff);
    }
    if (errnums > 0) errcode = -1;
    return errcode;
 }
 /*
 // Check that quadratic subbotin spline 1st derivatives are the same
 // for the left piece of spline and the right piece of spline
 // for each break-point
 */
 static int dCheckQuadSubb1stDerConsistency( MKL_INT n, double x[], double qnode[],
                                             MKL_INT ny, double sppcoeff[],
                                             double *left, double *right )
{
    int i, j, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    int errnums = 0;
    double dif;

    for ( i = 0; i < ny; i++ )
    {
        right[i*(n-1) + 0] = dComputePPSplineDer(1, n, 0, qnode[0], x[0], i, sorder,
                                                sppcoeff);
        for ( j = 1; j < n-1; j++ )
        {
            left[i*(n-1) + j-1] = dComputePPSplineDer(1, n, j-1, qnode[j], x[j-1], i,
                                                     sorder, sppcoeff);
            right[i*(n-1) + j]  = dComputePPSplineDer(1, n, j, qnode[j], x[j], i,
                                                     sorder, sppcoeff);
            dif = DF_ABS( left[i*(n-1) + j-1] - right[i*(n-1) + j] );
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
        left[i*(n-1) + n-2] = dComputePPSplineDer(1, n, n-2, qnode[n-1], x[n-2], i,
                                                 sorder, sppcoeff);
    }
    if (errnums > 0) errcode = -1;
    return errcode;
}

/*
// Check that spline coefficienst satisfy given boundary conditions
// in single precision case
*/
static int sCheckQuadBC( MKL_INT n, float x[], MKL_INT ny,
                         float sppcoeff[], MKL_INT bc_type, float bc[] )
{
    int i, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    float left_val, spline_left_val, dif;

    if ( (bc_type & DF_BC_Q_VAL) == DF_BC_Q_VAL )
    {
        left_val = ( x[0] + x[1] ) / 2.0;
        for ( i = 0; i < ny; i++ )
        {
            spline_left_val = sComputePPSplineDer(0, n, 0, left_val, x[0], i,
                                                    sorder, sppcoeff);
            dif = DF_ABS(spline_left_val - bc[0]);
            if (dif > EPSILON_SINGLE)
            {
                errcode = -1;
                return errcode;
            }
        }
    }
    else
    {
        if ( (bc_type & DF_BC_1ST_LEFT_DER)  == DF_BC_1ST_LEFT_DER )
        {
            errcode = sCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_1ST_RIGHT_DER) == DF_BC_1ST_RIGHT_DER )
        {
            errcode = sCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_LEFT_DER)  == DF_BC_2ND_LEFT_DER )
        {
            errcode = sCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_RIGHT_DER) == DF_BC_2ND_RIGHT_DER )
        {
            errcode = sCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
    }

    return errcode;
}

/*
// Check that spline coefficienst satisfy given boundary conditions
// in single precision case
*/
static int dCheckQuadBC( MKL_INT n, double x[], MKL_INT ny,
                         double sppcoeff[], MKL_INT bc_type, double bc[] )
{
    int i, errcode = 0;
    MKL_INT sorder = DF_PP_QUADRATIC;
    double left_val, spline_left_val, dif;

    if ( (bc_type & DF_BC_Q_VAL) == DF_BC_Q_VAL )
    {
        left_val = ( x[0] + x[1] ) / 2.0;
        for ( i = 0; i < ny; i++ )
        {
            spline_left_val = dComputePPSplineDer(0, n, 0, left_val, x[0], i,
                                                    sorder, sppcoeff);
            dif = DF_ABS(spline_left_val - bc[0]);
            if (dif > EPSILON_SINGLE)
            {
                errcode = -1;
                return errcode;
            }
        }
    }
    else
    {
        if ( (bc_type & DF_BC_1ST_LEFT_DER)  == DF_BC_1ST_LEFT_DER )
        {
            errcode = dCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_1ST_RIGHT_DER) == DF_BC_1ST_RIGHT_DER )
        {
            errcode = dCheckQuad1stDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_LEFT_DER)  == DF_BC_2ND_LEFT_DER )
        {
            errcode = dCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff, bc[0], 1 );
            if ( errcode < 0 ) return errcode;
        }
        if ( (bc_type & DF_BC_2ND_RIGHT_DER) == DF_BC_2ND_RIGHT_DER )
        {
            errcode = dCheckQuad2ndDerBC( n, x, ny, sorder, sppcoeff, bc[1], 0 );
            if ( errcode < 0 ) return errcode;
        }
    }

    return errcode;
}

/*
// Check 1st derivatrives boundary conditions of quadratic splines
// in single precision case
*/
static int sCheckQuad1stDerBC( MKL_INT n, float x[], MKL_INT ny,
                               MKL_INT sorder, float sppcoeff[],
                               float bc, int isLeft )
{
    int i, piece_idx, x_idx;
    float spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-1) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = sComputePPSplineDer(1, n, piece_idx, x[x_idx],
                                             x[piece_idx], i, sorder,
                                             sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 1st derivatrives boundary conditions of cubic splines
// in double precision case
*/
static int dCheckQuad1stDerBC( MKL_INT n, double x[], MKL_INT ny,
                               MKL_INT sorder, double sppcoeff[],
                               double bc, int isLeft )
{
    int i, piece_idx, x_idx;
    double spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-1) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = dComputePPSplineDer(1, n, piece_idx, x[x_idx],
                                             x[piece_idx], i, sorder,
                                             sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_DOUBLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 2nd derivatrives boundary conditions of quadratic splines
// in single precision case
*/
static int sCheckQuad2ndDerBC( MKL_INT n, float x[], MKL_INT ny,
                               MKL_INT sorder, float sppcoeff[],
                               float bc, int isLeft )
{
    int i, piece_idx, x_idx;
    float spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-1) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = sComputePPSplineDer(2, n, piece_idx, x[x_idx],
                                             x[x_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}

/*
// Check 2nd derivatrives boundary conditions of quadratic splines
// in double precision case
*/
static int dCheckQuad2ndDerBC( MKL_INT n, double x[], MKL_INT ny,
                               MKL_INT sorder, double sppcoeff[],
                               double bc, int isLeft )
{
    int i, piece_idx, x_idx;
    double spline_der_val, dif;
    int errcode = 0;

    piece_idx = ( isLeft ? 0 : (n-1) );
    x_idx     = ( isLeft ? 0 : (n-1) );

    for ( i = 0; i < ny; i++ )
    {
        spline_der_val = dComputePPSplineDer(2, n, piece_idx, x[x_idx],
                                             x[x_idx], i, sorder, sppcoeff);
        dif = fabs( spline_der_val - bc );
        if ( dif > EPSILON_SINGLE )
        {
            errcode = -1;
            return errcode;
        }
    }

    return errcode;
}


/*
//++
//  SEARCH RELATED CHECKS
//--
*/

/*
// Single precision cell search
*/
static int sFindCells( MKL_INT n, float x[], MKL_INT nsite, float site[],
                       MKL_INT bin[] )
{
    int i, j;

    for ( j = 0; j < nsite; j++ )
    {
        if( site[j] == x[n-1])
        {
            bin[j] = n-1;
        }
        else
        {
        for ( i = 0; i < n; i++ )
        {
                bin[j] = n; /* default result */
            if ( site[j] < x[i] )
            {
                    bin[j] = i; /* found result */
                break;
            }
        }
        }
    }

    return 1;
}

/*
// Double precision cell search
*/
static int dFindCells( MKL_INT n, double x[], MKL_INT nsite, double site[],
                       MKL_INT bin[] )
{
    int i, j;

    for ( j = 0; j < nsite; j++ )
    {
        if( site[j] == x[n-1])
        {
            bin[j] = n-1;
        }
        else
        {
        for ( i = 0; i < n; i++ )
        {
                bin[j] = n; /* default result */
            if ( site[j] < x[i] )
            {
                    bin[j] = i; /* found result */
                break;
            }
        }
        }
    }

    return 1;
}

/*
//++
//  INTERPOLATION RELATED CHECKS
//--
*/


/*
// Perform single precision piecewise polynomial spline evaluation
// and compare results with previously computed interpolation results
*/
static int sCheckPPSplineInterpRes( MKL_INT n, float x[], MKL_INT ny,
                               MKL_INT sorder, float scoeff[], MKL_INT nsite, float site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               float r[], float ref_r[] )
{
    int i, j, j1, jder, errcode = 0;
    int errnums = 0;
    int nderorder = 0;
    MKL_INT *bin;
    float dif;
    float *ref_r_ptr, *r_ptr;

    bin = malloc( nsite * sizeof(MKL_INT) );
    if(bin==NULL) exit(1);

    errcode = sFindCells(n, x, nsite, site, bin);
    for ( i = 0; i < nsite; i++ )
    {
        bin[i]--;
        if (bin[i] > n-2) bin[i] = n-2;
    }

    if (dorder)
    {
        for ( i = 0; i < ndorder; i++ )
        {
            if ( dorder[i] ) nderorder++;
        }
    }
    else
    {
        nderorder = 1;
    }

    for ( i = 0; i < ny*nderorder*nsite; i++ )
        ref_r[i] = 0.0;

    for ( i = 0; i < ny; i++ )
    {
        ref_r_ptr = ref_r + i*nderorder*nsite;
        r_ptr     = r     + i*nderorder*nsite;

        for ( j = 0; j < nsite; j++ )
        {
            j1 = 0;

            if ((site[j] < x[0]) || (site[j] > x[n-1]))
                continue;

            if ( !dorder || dorder[0] )
            {
                ref_r_ptr[j*nderorder+j1] = sComputePPSplineDer(0, n, bin[j],
                                                site[j], x[bin[j]], i, sorder, scoeff);
                dif = DF_ABS(ref_r_ptr[j*nderorder+j1]-r_ptr[j*nderorder+j1]);
                if ( dif > EPSILON_SINGLE ) errnums++;
                j1++;
            }

            if ( dorder )
            {
                jder = 1;
                while (jder <= nderorder)
                {
                    if ( (ndorder > jder) && dorder[jder] )
                    {
                        ref_r_ptr[j*nderorder+j1] = sComputePPSplineDer(jder, n, bin[j],
                                                site[j], x[bin[j]], i, sorder, scoeff);
                        dif = DF_ABS(ref_r_ptr[j*nderorder+j1]-r_ptr[j*nderorder+j1]);
                        if ( dif > EPSILON_SINGLE ) errnums++;
                        j1++;
                    }
                    jder++;
                }
            }
        }
    }

    free(bin);

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Perform double precision piecewise polynomial spline evaluation
// and compare results with previously computed interpolation results
*/
static int dCheckPPSplineInterpRes( MKL_INT n, double x[], MKL_INT ny,
                               MKL_INT sorder, double scoeff[], MKL_INT nsite, double site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               double r[], double ref_r[] )
{
    int i, j, j1, jder, errcode = 0;
    int errnums = 0;
    int nderorder = 0;
    MKL_INT *bin;
    double dif;
    double *ref_r_ptr, *r_ptr;

    bin = malloc( nsite * sizeof(MKL_INT) );
    if(bin==NULL) exit(1);

    errcode = dFindCells(n, x, nsite, site, bin);
    for ( i = 0; i < nsite; i++ )
    {
        bin[i]--;
        if (bin[i] > n-2) bin[i] = n-2;
    }

    if (dorder)
    {
        for ( i = 0; i < ndorder; i++ )
        {
            if ( dorder[i] ) nderorder++;
        }
    }
    else
    {
        nderorder = 1;
    }

    for ( i = 0; i < ny*nderorder*nsite; i++ )
        ref_r[i] = 0.0;

    for ( i = 0; i < ny; i++ )
    {
        ref_r_ptr = ref_r + i*nderorder*nsite;
        r_ptr     = r     + i*nderorder*nsite;

        for ( j = 0; j < nsite; j++ )
        {
            j1 = 0;

            if ((site[j] < x[0]) || (site[j] > x[n-1]))
                continue;

            if ( !dorder || dorder[0] )
            {
                ref_r_ptr[j*nderorder+j1] = dComputePPSplineDer(0, n, bin[j],
                                                site[j], x[bin[j]], i, sorder, scoeff);
                dif = DF_ABS(ref_r_ptr[j*nderorder+j1]-r_ptr[j*nderorder+j1]);
                if ( dif > EPSILON_SINGLE ) errnums++;
                j1++;
            }

            if ( dorder )
            {
                jder = 1;
                while (jder <= nderorder)
                {
                    if ( (ndorder > jder) && dorder[jder] )
                    {
                        ref_r_ptr[j*nderorder+j1] = dComputePPSplineDer(jder, n, bin[j],
                                                site[j], x[bin[j]], i, sorder, scoeff);
                        dif = DF_ABS(ref_r_ptr[j*nderorder+j1]-r_ptr[j*nderorder+j1]);
                        if ( dif > EPSILON_SINGLE ) errnums++;
                        j1++;
                    }
                    jder++;
                }
            }
        }
    }

    free(bin);

    if (errnums > 0) errcode = -1;

    return errcode;
}


/*
// Perform single precision linear spline evaluation and compare results with
// previously computed interpolation results
*/
static int sCheckLinInterpRes( MKL_INT n, float x[], MKL_INT ny,
                               float scoeff[], MKL_INT nsite, float site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               float r[], float ref_r[] )
{
    return sCheckPPSplineInterpRes(n, x, ny, DF_PP_LINEAR, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}

/*
// Perform double precision linear spline evaluation and compare results with
// previously computed interpolation results
*/
static int dCheckLinInterpRes( MKL_INT n, double x[], MKL_INT ny,
                               double scoeff[], MKL_INT nsite, double site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               double r[], double ref_r[] )
{
    return dCheckPPSplineInterpRes(n, x, ny, DF_PP_LINEAR, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}


/*
// Perform single precision cubic spline evaluation and compare results with
// previously computed interpolation results
*/
static int sCheckCubInterpRes( MKL_INT n, float x[], MKL_INT ny,
                               float scoeff[], MKL_INT nsite, float site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               float r[], float ref_r[] )
{
    return sCheckPPSplineInterpRes(n, x, ny, DF_PP_CUBIC, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}

/*
// Perform double precision cubic spline evaluation and compare results with
// previously computed interpolation results
*/
static int dCheckCubInterpRes( MKL_INT n, double x[], MKL_INT ny,
                               double scoeff[], MKL_INT nsite, double site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               double r[], double ref_r[] )
{
    return dCheckPPSplineInterpRes(n, x, ny, DF_PP_CUBIC, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}

/*
// Perform single precision cubic spline integration and compare results with
// previously computed integration results
*/
static int sCheckCubIntegrRes( MKL_INT n, float *x, MKL_INT ny, float *scoeff,
                               MKL_INT nlim, float *llim, float *rlim,
                               float *r, float *ref_r )
{
    int i, j, jbin, tmp_bin, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    MKL_INT *lbin, *rbin;
    float sign, tmp_lim;
    float dif;
    float *ref_r_ptr, *r_ptr;

    lbin = malloc( 2 * nlim * sizeof(MKL_INT) );
    if(lbin==NULL) exit(1);
    rbin = lbin + nlim;

    errcode = sFindCells(n, x, nlim, llim, lbin);
    errcode = sFindCells(n, x, nlim, rlim, rbin);
    for ( i = 0; i < nlim; i++ )
    {
        lbin[i]--;
        rbin[i]--;
    }

    for ( i = 0; i < ny*nlim; i++ )
        ref_r[i] = 0.0;

    for ( i = 0; i < ny; i++ )
    {
        ref_r_ptr = ref_r + i*nlim;
        r_ptr     = r     + i*nlim;

        for ( j = 0; j < nlim; j++ )
        {
            sign = 1.0;
            if ( rlim[j] < llim[j] )
            {
                /* Here if right integration limit is lower than left */
                tmp_lim = llim[j];
                llim[j] = rlim[j];
                rlim[j] = tmp_lim;

                tmp_bin = lbin[j];
                lbin[j] = rbin[j];
                rbin[j] = tmp_bin;

                sign = -1.0;
            }

            if ( lbin[j] == rbin[j] )
            {
                /* Here if there are no break points between integration
                   limits */

                ref_r_ptr[j] = sComputeCubSplineIntegr( n, lbin[j], llim[j],
                    rlim[j], x[lbin[j]], i, sorder, scoeff );
            }
            else
            {
                /* Here if there is one or more break points between integration
                   limits */

                ref_r_ptr[j] = sComputeCubSplineIntegr( n, lbin[j], llim[j],
                    x[lbin[j]+1], x[lbin[j]], i, sorder, scoeff );

                for ( jbin = lbin[j] + 1; jbin < rbin[j]; jbin++ )
                {
                    ref_r_ptr[j] += sComputeCubSplineIntegr( n, jbin,
                        x[jbin], x[jbin+1], x[jbin], i, sorder, scoeff );
                }

                ref_r_ptr[j] += sComputeCubSplineIntegr( n, rbin[j],
                    x[rbin[j]], rlim[j], x[rbin[j]], i, sorder, scoeff );
            }

            ref_r_ptr[j] *= sign;

            dif = DF_ABS(ref_r_ptr[j]-r_ptr[j]);
            if ( dif > EPSILON_SINGLE ) errnums++;
        }
    }

    free(lbin);

    if (errnums > 0) errcode = -1;

    return errcode;
}

/*
// Perform double precision cubic spline integration and compare results with
// previously computed integration results
*/
static int dCheckCubIntegrRes( MKL_INT n, double *x, MKL_INT ny, double *scoeff,
                               MKL_INT nlim, double *llim, double *rlim,
                               double *r, double *ref_r )
{
    int i, j, jbin, tmp_bin, errcode = 0;
    MKL_INT sorder = DF_PP_CUBIC;
    int errnums = 0;
    MKL_INT *lbin, *rbin;
    double sign, tmp_lim;
    double dif;
    double *ref_r_ptr, *r_ptr;

    lbin = malloc( 2 * nlim * sizeof(MKL_INT) );
    if(lbin==NULL) exit(1);
    rbin = lbin + nlim;

    errcode = dFindCells(n, x, nlim, llim, lbin);
    errcode = dFindCells(n, x, nlim, rlim, rbin);
    for ( i = 0; i < nlim; i++ )
    {
        lbin[i]--;
        rbin[i]--;
    }

    for ( i = 0; i < ny*nlim; i++ )
        ref_r[i] = 0.0;

    for ( i = 0; i < ny; i++ )
    {
        ref_r_ptr = ref_r + i*nlim;
        r_ptr     = r     + i*nlim;

        for ( j = 0; j < nlim; j++ )
        {
            sign = 1.0;
            if ( rlim[j] < llim[j] )
            {
                /* Here if right integration limit is lower than left */
                tmp_lim = llim[j];
                llim[j] = rlim[j];
                rlim[j] = tmp_lim;

                tmp_bin = lbin[j];
                lbin[j] = rbin[j];
                rbin[j] = tmp_bin;

                sign = -1.0;
            }

            if ( lbin[j] == rbin[j] )
            {
                /* Here if there are no break points between integration
                   limits */

                ref_r_ptr[j] = dComputeCubSplineIntegr( n, lbin[j], llim[j],
                    rlim[j], x[lbin[j]], i, sorder, scoeff );
            }
            else
            {
                /* Here if there is one or more break points between integration
                   limits */

                ref_r_ptr[j] = dComputeCubSplineIntegr( n, lbin[j], llim[j],
                    x[lbin[j]+1], x[lbin[j]], i, sorder, scoeff );

                for ( jbin = lbin[j] + 1; jbin < rbin[j]; jbin++ )
                {
                    ref_r_ptr[j] += dComputeCubSplineIntegr( n, jbin,
                        x[jbin], x[jbin+1], x[jbin], i, sorder, scoeff );
                }

                ref_r_ptr[j] += dComputeCubSplineIntegr( n, rbin[j],
                    x[rbin[j]], rlim[j], x[rbin[j]], i, sorder, scoeff );
            }

            ref_r_ptr[j] *= sign;

            dif = DF_ABS(ref_r_ptr[j]-r_ptr[j]);
            if ( dif > EPSILON_DOUBLE ) errnums++;
        }
    }

    free(lbin);

    if (errnums > 0) errcode = -1;

    return errcode;
}


/*
// Perform single precision quadratic spline evaluation and compare results
// with previously computed interpolation results
*/
static int sCheckQuadInterpRes( MKL_INT n, float x[], MKL_INT ny,
                                float scoeff[], MKL_INT nsite, float site[],
                                MKL_INT ndorder, MKL_INT dorder[],
                                float r[], float ref_r[] )
{
    return sCheckPPSplineInterpRes(n, x, ny, DF_PP_QUADRATIC, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}


/*
// Perform double precision quadratic spline evaluation and compare results
// with previously computed interpolation results
*/
static int dCheckQuadInterpRes( MKL_INT n, double x[], MKL_INT ny,
                                double scoeff[], MKL_INT nsite, double site[],
                                MKL_INT ndorder, MKL_INT dorder[],
                                double r[], double ref_r[] )
{
    return dCheckPPSplineInterpRes(n, x, ny, DF_PP_QUADRATIC, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}


/*
// Perform single precision fifth order spline evaluation and compare results with
// previously computed interpolation results
*/
static int sCheckFifthOrderInterpRes( MKL_INT n, float x[], MKL_INT ny,
                               float scoeff[], MKL_INT nsite, float site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               float r[], float ref_r[] )
{
    MKL_INT sorder = 5;
    return sCheckPPSplineInterpRes(n, x, ny, sorder, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}

/*
// Perform double precision fifth order spline evaluation and compare results with
// previously computed interpolation results
*/
static int dCheckFifthOrderInterpRes( MKL_INT n, double x[], MKL_INT ny,
                               double scoeff[], MKL_INT nsite, double site[],
                               MKL_INT ndorder, MKL_INT dorder[],
                               double r[], double ref_r[] )
{
    MKL_INT sorder = 5;
    return dCheckPPSplineInterpRes(n, x, ny, sorder, scoeff, nsite, site,
                                   ndorder, dorder, r, ref_r);
}
