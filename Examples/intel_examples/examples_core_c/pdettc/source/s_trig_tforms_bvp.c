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
!  Single precision C test example for trigonometric transforms
!*******************************************************************************
!
! This example gives the solution of the 1D differential problems
! with the equation  -u"+u=f(x), 0<x<1, and with 3 types of boundary conditions:
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
  float pi=3.14159265358979324, xi;
  float maxErr = 0.0f, maxSolErr = 0.0f, maxExpErr = 0.0f;
  float *u = 0, *f = 0, *spar = 0, *lambda = 0;
  DFTI_DESCRIPTOR_HANDLE handle = 0;

  /* Printing the header for the example */
  printf("\n Example of use of Intel(R) Math Kernel Library (Intel(R) MKL) Trigonometric Transforms\n");
  printf(" **********************************************\n\n");
  printf(" This example gives the the solutions of the 1D differential problems\n");
  printf(" with the equation -u\"+u=f(x), 0<x<1, \n");
  printf(" and with 3 types of boundary conditions:\n");
  printf(" DD case: u(0)=u(1)=0,\n");
  printf(" NN case: u'(0)=u'(1)=0,\n");
  printf(" ND case: u'(0)=u(1)=0.\n");
  printf(" ----------------------------------------------------------------------\n");
  printf(" In general, the error should be of order O(1.0/n**2)\n");
  printf(" For this example, the value of n is "LI"\n", n);
  printf(" The approximation error should be of order 5.0e-002, ");
  printf(" if everything is OK\n");
  printf(" ----------------------------------------------------------------------\n");
  printf("                     SINGLE PRECISION COMPUTATIONS                     \n");
  printf(" ======================================================================\n");
  printf("\n");

  u=(float*)mkl_malloc((n+1)*sizeof(float), 64);
  if (u == 0) goto FAILURE;
  f=(float*)mkl_malloc((n+1)*sizeof(float), 64);
  if (f == 0) goto FAILURE;
  /* NOTE: This example uses shorter spar array of size 3n/2+2 instead of 5n/2+2
  as only sine, cosine, and staggered cosine transforms are used. More details
  can be found in Chapter 13 of Intel MKL Manual. */
  spar=(float*)mkl_malloc((3*n/2+2)*sizeof(float), 64);
  if (spar == 0) goto FAILURE;
  lambda=(float*)mkl_malloc((n+1)*sizeof(float), 64);
  if (lambda == 0) goto FAILURE;

  for(i=0;i<=2;i++)
  {
    /* Varying the type of the transform */
    tt_type=i;

    /* Computing test solutions u(x) */
    for(k=0;k<=n;k++)
    {
      xi=(float)1.0E0*(float)k/(float)n;
      u[k]=pow(sin(pi*xi),(float)2.0E0);
    }
    /* Computing the right-hand side f(x) */
    for(k=0;k<=n;k++)
    {
      f[k]=((float)4.0E0*(pi*pi)+(float)1.0E0)*u[k]-(float)2.0E0*(pi*pi);
    }
    /* Computing the right-hand side for the algebraic system */
    for(k=0;k<=n;k++)
    {
      f[k]=f[k]/(n*n);
    }
    if (tt_type==0)
    {
      /* The Dirichlet boundary conditions */
      f[0]=(float)0.0E0;
      f[n]=(float)0.0E0;
    }
    if (tt_type==2)
    {
      /* The mixed Neumann-Dirichlet boundary conditions */
      f[n]=(float)0.0E0;
    }

    /* Computing the eigenvalues for the three-point finite-difference problem */
    if (tt_type==0||tt_type==1)
    {
      for(k=0;k<=n;k++)
      {
        lambda[k]=pow((float)2.0E0*sin((float)0.5E0*pi*(float)k/(float)n),
                                              (float)2.0E0)+(float)1.0E0/(float)(n*n);
      }
    }
    if (tt_type==2)
    {
      for(k=0;k<=n;k++)
      {
        lambda[k]=pow((float)2.0E0*sin((float)0.25E0*pi*(float)(2*k+1)/(float)n),
                                              (float)2.0E0)+(float)1.0E0/(float)(n*n);
      }
    }

    /* Computing the solution of 1D problem using trigonometric transforms
    First we initialize the transform */
    s_init_trig_transform(&n,&tt_type,ipar,spar,&ir);
    if (ir!=0) goto FAILURE;
    /* Then we commit the transform. Note that the data in f will be changed at
    this stage !
    If you want to keep them, save them in some other array before the call to
    the routine */
    s_commit_trig_transform(f,&handle,ipar,spar,&ir);
    if (ir!=0) goto FAILURE;
    /* Now we can apply trigonometric transform */
    s_forward_trig_transform(f,&handle,ipar,spar,&ir);
    if (ir!=0) goto FAILURE;

    /* Scaling the solution by the eigenvalues */
    for(k=0;k<=n;k++)
    {
      f[k]=f[k]/lambda[k];
    }

    /* Now we can apply trigonometric transform once again as ONLY input vector f
    has changed */
    s_backward_trig_transform(f,&handle,ipar,spar,&ir);
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
      printf("u["LI"]=", k);
      printf("%6.3f\n",f[k]);
    }
    printf("\nError=%6.3e\n\n", maxErr);
    if (maxErr > 0.05)
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
  mkl_free(spar); 
  mkl_free(lambda); 
  return 0;

  /* Failure message to print if something went wrong */
FAILURE: printf("Failed to compute the solution(s)...\n");
        if (u != 0) mkl_free(u); 
        if (f != 0) mkl_free(f); 
        if (spar != 0) mkl_free(spar); 
        if (lambda != 0) mkl_free(lambda); 
        /* Free Intel MKL memory if any was allocated */
        MKL_Free_Buffers();
        return 1;

  /* End of the example code */
}
