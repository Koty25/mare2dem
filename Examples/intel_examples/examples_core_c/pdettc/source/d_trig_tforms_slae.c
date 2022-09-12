/*******************************************************************************
* Copyright 2006-2020 Intel Corporation.
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
!  Double precision C test example for trigonometric transforms
!*******************************************************************************
!
! This example gives the solution of the linear algebraic system,
! corresponding to the three-point approximation of the 1D differential
! problems with the equation  -u"+u=f(x), 0<x<1, and with 3 types of
! boundary conditions:
! u(0)=u(1)=0 (DD case), or u'(0)=u'(1)=0 (NN case), or u'(0)=u(1)=0 (ND case)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_service.h"
#include "mkl_trig_transforms.h"

#if !defined(MKL_ILP64)
#define LI "%d"
#else
#define LI "%lli"
#endif

int main(void)
{
  MKL_INT n=8, i, k, tt_type;
  MKL_INT ir, ipar[128];
  double pi=3.14159265358979324, xi, c;
  double maxErr = 0.0, maxSolErr = 0.0, maxExpErr = 0.0;
  double *u = 0, *f = 0, *dpar = 0, *lambda = 0;
  DFTI_DESCRIPTOR_HANDLE handle = 0;

  /* Printing the header for the example */
  printf("\n Example of use of Intel(R) Math Kernel Library (Intel(R) MKL) Trigonometric Transforms\n");
  printf(" **********************************************\n\n");
  printf(" This example gives the the exact solutions of the linear algebraic \n");
  printf(" systems, corresponding to the three-point approximations of the 1D \n");
  printf(" differential problems with the equation -u\"+u=f(x), 0<x<1, \n");
  printf(" and with 3 types of boundary conditions:\n");
  printf(" DD case: u(0)=u(1)=0,\n");
  printf(" NN case: u'(0)=u'(1)=0,\n");
  printf(" ND case: u'(0)=u(1)=0.\n");
  printf(" -----------------------------------------------------------------\n");
  printf(" The error expected to be close to the machine zero in all 3 cases\n");
  printf(" if everything is OK\n");
  printf(" -----------------------------------------------------------------\n");
  printf("                   DOUBLE PRECISION COMPUTATIONS                  \n");
  printf(" =================================================================\n\n");

  u=(double*)mkl_malloc((n+1)*sizeof(double), 64);
  if (u == 0) goto FAILURE;
  f=(double*)mkl_malloc((n+1)*sizeof(double), 64);
  if (f == 0) goto FAILURE;
  /* NOTE: This example uses shorter dpar array of size 3n/2+2 instead of 5n/2+2
  as only sine, cosine, and staggered cosine transforms are used. More details
  can be found in Chapter 13 of Intel MKL Manual. */
  dpar=(double*)mkl_malloc((3*n/2+2)*sizeof(double), 64);
  if (dpar == 0) goto FAILURE;
  lambda=(double*)mkl_malloc((n+1)*sizeof(double), 64);
  if (lambda == 0) goto FAILURE;

  for(i=0;i<=2;i++)
  {
    /* Varying the type of the transform */
    tt_type=i;

    /* Computing test solutions u(x) */
    for(k=0;k<=n;k++)
    {
      xi=1.0E0*k/n;
      /* The exact solution with Dirichlet-Dirichlet boundary conditions */
      if (tt_type==0) c=sin(pi*xi)*exp(xi);
      /* The exact solution with Neumann-Neumann boundary conditions */
      if (tt_type==1) c=(pi*cos(pi*xi)-sin(pi*xi))*exp(xi);
      /* The exact solution with Neumann-Dirichlet boundary conditions */
      if (tt_type==2) c=cos(0.5E0*pi*xi)*(exp(xi)-xi);
      u[k]=c;
    }
    /* Computing the right-hand side f(x) that corresponds to the finite element
    approximation and provides the exact solution computed above of the resulting
    algebraic system */
    for(k=1;k<=n-1;k++)
    {
      f[k]=(2.0E0+1.0E0/(n*n))*u[k]-u[k-1]-u[k+1];
    }
    if (tt_type==0)
    {
      /* The Dirichlet boundary conditions */
      f[0]=0.0E0;
      f[n]=0.0E0;
    }
    if (tt_type==1)
    {
      /* The Neumann boundary conditions */
      f[0]=(2.0E0+1.0E0/(n*n))*u[0]-2.0E0*u[1];
      f[n]=(2.0E0+1.0E0/(n*n))*u[n]-2.0E0*u[n-1];
    }
    if (tt_type==2)
    {
      /* The mixed Neumann-Dirichlet boundary conditions */
      f[0]=(2.0E0+1.0E0/(n*n))*u[0]-2.0E0*u[1];
      f[n]=0.0E0;
    }

    /* Computing the eigenvalues for the three-point finite-difference problem */
    if (tt_type==0||tt_type==1)
    {
      for(k=0;k<=n;k++)
      {
        lambda[k]=pow(2.0E0*sin(0.5E0*pi*k/n),2.0E0)+1.0E0/(n*n);
      }
    }
    if (tt_type==2)
    {
      for(k=0;k<=n;k++)
      {
        lambda[k]=pow(2.0E0*sin(0.25E0*pi*(2*k+1)/n),2.0E0)+1.0E0/(n*n);
      }
    }

    /* Computing the solution of 1D problem using trigonometric transforms
       First we initialize the transform */
    d_init_trig_transform(&n,&tt_type,ipar,dpar,&ir);
    if (ir!=0) goto FAILURE;
    /* Then we commit the transform. Note that the data in f will be changed at
    this stage !
    If you want to keep them, save them in some other array before the call to
    the routine */
    d_commit_trig_transform(f,&handle,ipar,dpar,&ir);
    if (ir!=0) goto FAILURE;
    /* Now we can apply trigonometric transform */
    d_forward_trig_transform(f,&handle,ipar,dpar,&ir);
    if (ir!=0) goto FAILURE;

    /* Scaling the solution by the eigenvalues */
    for(k=0;k<=n;k++)
    {
      f[k]=f[k]/lambda[k];
    }

    /* Now we can apply trigonometric transform once again as ONLY input vector f
    has changed */
    d_backward_trig_transform(f,&handle,ipar,dpar,&ir);
    if (ir!=0) goto FAILURE;
    /* Cleaning the memory used by handle
       Now we can use handle for other kind of trigonometric transform */
    free_trig_transform(&handle,ipar,&ir);
    if (ir!=0) goto FAILURE;

    /* Performing the error analysis */
    for(k=0;k<=n;k++)
    {
        maxExpErr = (fabs(u[k]) > maxExpErr) ? fabs(u[k]) : maxExpErr;
        maxSolErr = (fabs(f[k]-u[k]) > maxSolErr) ? fabs(f[k]-u[k]) : maxSolErr;
    }
    maxErr = maxExpErr ? maxSolErr / maxExpErr : maxErr;

    /* Printing the results */
    if (tt_type==0)
    {
      printf("The computed solution of DD problem is\n\n");
    }
    if (tt_type==1)
    {
      printf("The computed solution of NN problem is\n\n");
    }
    if (tt_type==2)
    {
      printf("The computed solution of ND problem is\n\n");
    }
    for(k=0;k<=n;k++)
    {
      printf("u["LI"]=",k);
      printf("%6.3f\n",f[k]);
    }
    printf("\nError=%6.3e\n\n", maxErr);
    if (maxErr > 1.0e-14)
    {
      printf("The computed solution seems to be inaccurate. ");
      goto FAILURE;
    }
    /* End of the loop over the different kind of transforms and problems */
  }

  /* Free Intel MKL memory if any was allocated */
  MKL_Free_Buffers();
  /* Success message to print if everything is OK */
  printf("This example has successfully PASSED through all steps of computation!\n");
  mkl_free(u); 
  mkl_free(f); 
  mkl_free(dpar); 
  mkl_free(lambda); 
  return 0;

  /* Failure message to print if something went wrong */
FAILURE: printf("Failed to compute the solution(s)...\n");
        /* Free Intel MKL memory if any was allocated */
        MKL_Free_Buffers();
        if (u != 0) mkl_free(u); 
        if (f != 0) mkl_free(f); 
        if (dpar != 0) mkl_free(dpar); 
        if (lambda != 0) mkl_free(lambda); 
        return 1;

  /* End of the example code */
}
