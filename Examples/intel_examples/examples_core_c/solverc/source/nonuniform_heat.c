/*******************************************************************************
* Copyright 2013-2020 Intel Corporation.
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
*   Content : Example of solving nonlinear nonuniform heat problem  
*             by Intel(R) Math Kernel Library (Intel(R) MKL) routines.
*             We will use Intel MKL routines from different
*             sub-domains: PDE Poisson, ISS, Sparse BLAS, BLAS
*
********************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "memory.h"

#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"
#include "mkl_dfti.h"
#include "mkl_poisson.h"

void construct_matrix(MKL_INT n, MKL_INT** ia_ptr, MKL_INT** ja_ptr, double** a_ptr, MKL_INT *error);
int finalize(double *a1, double *a2, double *a3, double *a4, double *a5, double *a6, 
             double *a7, double *a8, double *a9, double *a10, 
             DFTI_DESCRIPTOR_HANDLE xhandle, DFTI_DESCRIPTOR_HANDLE yhandle, MKL_INT *ipar_helm,
             int error);

MKL_INT main (void)
{
    /*---------------------------------------------------------------------------
    ** Parameters of problem
    **-------------------------------------------------------------------------*/
    MKL_INT one = 1;
    MKL_INT n = 64;
    MKL_INT number_unknowns;   /* number of unknowns */
    double *solution = NULL, *rhs = NULL;
    MKL_INT error;

    /*---------------------------------------------------------------------------
    ** Sparse matrix in Compressed Sparse Row (CSR format) and parameters
    **-------------------------------------------------------------------------*/
    MKL_INT *ia = NULL, *ja = NULL;
    double *a = NULL;
    // Descriptor of main sparse matrix properties
    struct matrix_descr descrA;
    // Structure with sparse matrix stored in CSR format
    sparse_matrix_t       csrA;

    /*---------------------------------------------------------------------------
    ** RSI interface parameters and arrays
    **-------------------------------------------------------------------------*/
    MKL_INT rci_request, itercount;
    MKL_INT ipar[128];
    double  dpar[128];
    double *tmp = NULL;

    /*---------------------------------------------------------------------------
    ** PDE Poisson parameters and arrays
    **-------------------------------------------------------------------------*/
    char *BCtype;
    double *bd_ax = NULL, *bd_bx = NULL, *bd_ay = NULL, *bd_by = NULL, *bd_az = NULL, *bd_bz = NULL;
    double ax, bx, ay, by, az, bz;
    DFTI_DESCRIPTOR_HANDLE xhandle = 0;
    DFTI_DESCRIPTOR_HANDLE yhandle = 0;
    MKL_INT ipar_helm[128];
    double *dpar_helm = NULL;
    /* Set the coefficient q to 0 */
    double q  = 0.0E0;

    /*---------------------------------------------------------------------------
    ** Some additional variables to use in program
    **-------------------------------------------------------------------------*/
    MKL_INT i, j, k, kk;
    number_unknowns = (n+1)*(n+1)*(n+1);
    BCtype = "DDDDDD";

    /* Allocate arrays */
    if (((tmp = (double*) mkl_malloc( sizeof (double) * number_unknowns * 4, 64))   == NULL ) ||    
        ((solution  = (double*) mkl_malloc( sizeof (double) * number_unknowns, 64)) == NULL ) ||
        ((rhs      = (double*) mkl_malloc( sizeof (double) * number_unknowns, 64 )) == NULL ) ||
        ((dpar_helm = (double*) mkl_malloc( sizeof (double) * (13*n + 9), 64))      == NULL ) ||
        ((bd_ax    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) ||
        ((bd_bx    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) ||
        ((bd_ay    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) ||
        ((bd_by    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) ||
        ((bd_az    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) ||
        ((bd_bz    = (double*) mkl_malloc((n+1)*(n+1)*sizeof(double), 64))   == NULL ) )
        return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                        bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                        xhandle, yhandle, ipar_helm, 1);

    /* Initialize rhs */
    for (kk=0;kk<number_unknowns;kk++)
        rhs[kk] = 0.;

    for(k=1;k<n;k++)
        for(j=1;j<n;j++)
            for(i=1;i<n;i++)
                rhs[k*(n+1)*(n+1)+j*(n+1)+i]=k*k+i*i+j*j;

    /* Initialize zero dirichlet boundary condition */
    for (k=0;k<(n+1)*(n+1);k++) 
    {
        bd_ax[k] = 0.0E0; bd_bx[k] = 0.0E0;
        bd_ay[k] = 0.0E0; bd_by[k] = 0.0E0;
        bd_az[k] = 0.0E0; bd_bz[k] = 0.0E0;
    }

    /* Define the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Helmholtz Solver */
    ax = 0.0E0; bx = 1.0E0;
    ay = 0.0E0; by = 1.0E0;
    az = 0.0E0; bz = 1.0E0;

    /* Initialize ipar_helm array to for garbage clean up */
    for( i=0; i<128; i++ )
        ipar_helm[i]=0;

    /* Initialize simple data structures of Poisson Library for 3D Helmholtz Solver */
    d_init_Helmholtz_3D( &ax, &bx, &ay, &by, &az, &bz, &n,
                         &n, &n, BCtype, &q, 
                         ipar_helm, dpar_helm, &error );
    if (error != 0)
        return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                        bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                        xhandle, yhandle, ipar_helm, 1);

    /* Only internal data for Helmholtz solver needs to be initialized in this step */
    ipar_helm[0] = 99;
    d_commit_Helmholtz_3D(&tmp[3*number_unknowns], bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                            &xhandle, &yhandle, ipar_helm, dpar_helm, &error);
    if (error != 0)
        return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                        bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                        xhandle, yhandle, ipar_helm, error);

    /*  set up matrix in 3-arrays csr format */
    construct_matrix(n, &ia, &ja, &a, &error);
    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ZERO, number_unknowns, number_unknowns, ia, ia+1, ja, a );

    if (error != 0) 
        return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                        bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                        xhandle, yhandle, ipar_helm, error);

    /* Initialize the initial guess                                              */
    for (k = 0; k < number_unknowns; k++)
        solution[k] = 1.E0;

    /* Initialize parameters for RCI solver                                      */
    for (k = 0; k < 128; k++)
    {
        ipar[k] = 0;
        dpar[k] = 0.E0;
    }

    dcg_init (&number_unknowns, solution, rhs, &rci_request, ipar, dpar, tmp);
    if (rci_request != 0)
        return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                        bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                        xhandle, yhandle, ipar_helm, 1);

    /*---------------------------------------------------------------------------
    ** Set the desired parameters:                                               
    ** Integer parameters:                                                       
    ** set the maximal number of iterations to 100                               
    ** LOGICAL parameters:                                                       
    ** run the preconditioned version of RCI (P)CG with C_inverse preconditioner 
    ** DOUBLE parameters                                                         
    ** set stop criteria to 10^-10                                               */
    ipar[ 4] = 100;
    ipar[ 8] = 1;
    ipar[ 9] = 0;
    ipar[10] = 1;
    dpar[ 0] = 1.E-10;
    
/*---------------------------------------------------------------------------
** Compute the solution by RCI (P)CG solver without preconditioning          
** Reverse Communication starts here                                         
**-------------------------------------------------------------------------*/
    do {
        /* Do the next step */
        dcg (&number_unknowns, solution, rhs, &rci_request, ipar, dpar, tmp);

        switch(rci_request)
        {
        /* If rci_request=0, then the solution was found with the required precision */
        case 0:
            break;
        /* If rci_request=1, then compute the vector A*tmp[0]and put the result in vector tmp[n] */
        case 1:
	        mkl_sparse_d_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrA, descrA, tmp, 0.0, &tmp[number_unknowns]);
            break;
        /*---------------------------------------------------------------------------
        ** If rci_request=3, then compute step applies the simplest preconditioning  
        ** on vector tmp[2*n] and puts the result in vector tmp[3*n]                 
        ** We use Helmholtz functionality as preconditioner of main functionality    
        ** Using such kind of preconditioner reduces condition number of the system  
        ** to mu_max/mu_min, where mu_max and mu_min are respectively the minimum    
        ** and maximum values of density (1 and 2 in this case)                      
        **-------------------------------------------------------------------------*/
        case 3:
            dcopy(&number_unknowns, &(tmp[2*number_unknowns]), &one, &(tmp[3*number_unknowns]), &one );

            /* In this step we need to multiply rhs by mesh parameters. Because boundary  */
            /* condition is equal to zero we can skip this step also, so iparm_helm[0]    */
            /* is set to 990                                                              */
            ipar_helm[0] = 990;
            d_commit_Helmholtz_3D(&tmp[3*number_unknowns], bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                                    &xhandle,&yhandle, ipar_helm, dpar_helm, &error);
            if (error != 0) 
                return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                                bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                                xhandle, yhandle, ipar_helm, 1);

            d_Helmholtz_3D(&tmp[3*number_unknowns], bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz,
                            &xhandle,&yhandle, ipar_helm, dpar_helm, &error);
            if (error != 0)
                return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                                bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                                xhandle, yhandle, ipar_helm, 1);
            break;            

        /* If rci_request=anything else, then dcg subroutine failed                  
        ** to compute the solution vector: solution[n]                               */
        default:
            return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                            bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                            xhandle, yhandle, ipar_helm, 1);
            break;
        }
    } while (rci_request != 0);

/*---------------------------------------------------------------------------
** Reverse Communication ends here                                           
** Get the current iteration number into itercount                           
**-------------------------------------------------------------------------*/

    dcg_get (&number_unknowns, solution, rhs, &rci_request, ipar, dpar, tmp, &itercount);
    printf ("Number of iteration equal to %d\n", (int)itercount);

    /*----------------------------------------------------------------------------------
    ** Release internal Intel MKL memory that might have been used for computations  
    ** NOTE: It is important to call the routine below to avoid memory leaks            
    ** unless you disable Intel MKL Memory Manager                                      
    **--------------------------------------------------------------------------------*/
    return finalize(tmp, solution, dpar_helm, rhs, bd_ax, 
                    bd_bx, bd_ay, bd_by, bd_az, bd_bz, 
                    xhandle, yhandle, ipar_helm, error);
    mkl_sparse_destroy(csrA);
}

/* This function sets matrix value of array a composed of 3 arrays in csr format. Matrix structure is based
   on the following  stencil:

    -a_1*u(i,j,k-1) - a_2*u(i,j,k+1) - 
     a_3*u(i,j-1,k) - a_4*u(i,j+1,k) - 
     a_5*u(i-1,j,k) - a_6*u(i+1,j,k) - 
     a_7*u(i,j,k)   = f(i,j,k)

    where a_1 - a_7 parameters are computed from mu(z) from the differential equation div(mu(z)grad(u(x, y, z))) = f(x, y, z) with 
              mu(z) equal to 1.0 if z is less than 0.5 and to 2.0 otherwise. 
    Dirichlet boundary condition is used, where

    Nx, Ny, Nz are mesh sizes of computational domain. */

void construct_matrix( MKL_INT n, MKL_INT** ia_ptr, MKL_INT** ja_ptr, double** a_ptr, MKL_INT *error )
{
    MKL_INT i, j, k, jj, kk, nonzero_elements, temp1, temp2, n_local, num_el;
    MKL_INT *ia, *ja;
    double *a;
    double alpha=2.0E0;

/* Allocate memory for ia array */
    if (( ia = (MKL_INT*) mkl_malloc (sizeof (MKL_INT) * ((n+1)*(n+1)*(n+1)+1), 64)) == NULL )
    {
        *error = -1;
        return;
    }

/* Calculate number of elements in each string of matrix; required to set ia array */
    for ( k = 0; k < n+1; k++ )
    {
        for ( j = 0; j < n+1; j++ )
        {
            for ( i = 0; i < n+1; i++ )
            {
                nonzero_elements = 7;
                if ((k==0)||(k==n)) nonzero_elements--;
                if ((j==0)||(j==n)) nonzero_elements--;
                if ((i==0)||(i==n)) nonzero_elements--;
                ia[k*(n+1)*(n+1)+j*(n+1)+i] = nonzero_elements;
            }
        }
    }

/* Set ia array */
    temp1 = ia[0];
    ia[0] = 0;

    for ( kk = 1; kk <= (n+1)*(n+1)*(n+1); kk++ )
    {
        temp2  = ia[kk];
        ia[kk] = ia[kk-1] + temp1;
        temp1  = temp2;
    }

/* Allocate memory for ja array */
    if (( ja = (MKL_INT*) mkl_malloc (sizeof (MKL_INT) * ia[(n+1)*(n+1)*(n+1)], 64)) == NULL )
    {
        free( ia );
        ia = NULL;
        *error = -1;
        return;
    }

/* Set ja array */
    n_local = -1;
    for ( k = 0; k < n+1; k++ )
    {
        for ( j = 0; j < n+1; j++ )
        {
            for ( i = 0; i < n+1; i++ )
            {
                if ( k > 0 ) {n_local++; ja[n_local]=(k-1)*(n+1)*(n+1)+(j  )*(n+1)+i  ;}
                if ( j > 0 ) {n_local++; ja[n_local]=(k  )*(n+1)*(n+1)+(j-1)*(n+1)+i  ;}
                if ( i > 0 ) {n_local++; ja[n_local]=(k  )*(n+1)*(n+1)+(j  )*(n+1)+i-1;}

                              n_local++; ja[n_local]=(k  )*(n+1)*(n+1)+(j  )*(n+1)+i;

                if ( i < n ) {n_local++; ja[n_local]=(k  )*(n+1)*(n+1)+(j  )*(n+1)+i+1;}
                if ( j < n ) {n_local++; ja[n_local]=(k  )*(n+1)*(n+1)+(j+1)*(n+1)+i  ;}
                if ( k < n ) {n_local++; ja[n_local]=(k+1)*(n+1)*(n+1)+(j  )*(n+1)+i  ;}
            }
        }
    }
    
/* Allocate memory for a array */
    if (( a = (double*) mkl_malloc (sizeof (double) * ia[(n+1)*(n+1)*(n+1)], 64)) == NULL )
    {
        free( ia );
        ia = NULL;
        free( ja );
        ja = NULL;
        *error = -1;
        return;
    }


    for ( kk = 0; kk < ((n+1)*(n+1)*(n+1)); kk++ )
    {
        for ( jj = ia[kk]; jj < ia[kk+1]; jj++ )
        {
            if (ja[jj] == kk)
                a[jj] = 1.0E0;
            else
                a[jj] = 0.0E0;
        }
    }

/* Set correct value in matrix */
    for ( k = 1; k < n; k++ )
    {
        for ( j = 1; j < n; j++ )
        {
            for ( i = 1; i < n; i++ )
            {
                num_el = ia[k*(n+1)*(n+1)+j*(n+1)+i] + 3;
                if (k < n/2)
                {
                  if ( k > 1   ) { a[num_el-3] = -1.0/(n*n); } else { a[num_el-3] = 0.0E0; }
                  if ( j > 1   ) { a[num_el-2] = -1.0/(n*n); } else { a[num_el-2] = 0.0E0; }
                  if ( i > 1   ) { a[num_el-1] = -1.0/(n*n); } else { a[num_el-1] = 0.0E0; }
                                   a[num_el  ] =  6.0/(n*n);
                  if ( i < n-1 ) { a[num_el+1] = -1.0/(n*n); } else { a[num_el+1] = 0.0E0; }
                  if ( j < n-1 ) { a[num_el+2] = -1.0/(n*n); } else { a[num_el+2] = 0.0E0; }
                  if ( k < n-1 ) { a[num_el+3] = -1.0/(n*n); } else { a[num_el+3] = 0.0E0; }
                }
                else
                {
                    if (k > n/2)
                    {
                      if ( k > 1   ) { a[num_el-3] = -alpha/(n*n); } else { a[num_el-3] = 0.0E0; }
                      if ( j > 1   ) { a[num_el-2] = -alpha/(n*n); } else { a[num_el-2] = 0.0E0; }
                      if ( i > 1   ) { a[num_el-1] = -alpha/(n*n); } else { a[num_el-1] = 0.0E0; }
                                       a[num_el  ] =6*alpha/(n*n);
                      if ( i < n-1 ) { a[num_el+1] = -alpha/(n*n); } else { a[num_el+1] = 0.0E0; }
                      if ( j < n-1 ) { a[num_el+2] = -alpha/(n*n); } else { a[num_el+2] = 0.0E0; }
                      if ( k < n-1 ) { a[num_el+3] = -alpha/(n*n); } else { a[num_el+3] = 0.0E0; }
                    }
                    else /* (k==n/2)*/
                    {
                      if ( k > 1   ) { a[num_el-3] =         -1.0/(n*n); } else { a[num_el-3] = 0.0E0; }
                      if ( j > 1   ) { a[num_el-2] =  -(1+alpha)/(2*n*n); } else { a[num_el-2] = 0.0E0; }
                      if ( i > 1   ) { a[num_el-1] =  -(1+alpha)/(2*n*n); } else { a[num_el-1] = 0.0E0; }
                                       a[num_el  ] = 6*(1+alpha)/(2*n*n);
                      if ( i < n-1 ) { a[num_el+1] =  -(1+alpha)/(2*n*n); } else { a[num_el+1] = 0.0E0; }
                      if ( j < n-1 ) { a[num_el+2] =  -(1+alpha)/(2*n*n); } else { a[num_el+2] = 0.0E0; }
                      if ( k < n-1 ) { a[num_el+3] =        -alpha/(n*n); } else { a[num_el+3] = 0.0E0; }
                    }
                }
            }
        }
    } 

    *ia_ptr = ia;
    *ja_ptr = ja;
    *a_ptr  = a;
}
/* This function clean-up memory, print status and return actual error */
int finalize(double *a1, double *a2, double *a3, double *a4, double *a5, double *a6, 
             double *a7, double *a8, double *a9, double *a10, 
             DFTI_DESCRIPTOR_HANDLE xhandle, DFTI_DESCRIPTOR_HANDLE yhandle, MKL_INT *ipar_helm,
             int error)
{
    if (a1 ) mkl_free(a1 );
    if (a2 ) mkl_free(a2 );
    if (a3 ) mkl_free(a3 );
    if (a4 ) mkl_free(a4 );
    if (a5 ) mkl_free(a5 );
    if (a6 ) mkl_free(a6 );
    if (a7 ) mkl_free(a7 );
    if (a8 ) mkl_free(a8 );
    if (a9 ) mkl_free(a9 );
    if (a10) mkl_free(a10);
    if ((xhandle != NULL) && (yhandle != NULL)) 
        free_Helmholtz_3D(&xhandle, &yhandle, ipar_helm, &error);
    
    MKL_Free_Buffers ();

    if (error != 0) 
    {
        printf ("This example has FAILED!\n");
        return 1;
    }
    else
    {
        printf ("This example has successfully PASSED through all steps of computation!\n");
        return 0;
    }
}
