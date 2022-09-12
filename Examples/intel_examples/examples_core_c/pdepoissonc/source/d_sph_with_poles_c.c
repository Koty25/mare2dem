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
*  Content:
*  C double precision example of solving Helmholtz problem on a whole sphere
*  using Intel(R) Math Kernel Library (Intel(R) MKL) Poisson Library
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl_service.h"
/* Include Poisson Library header files */
#include "mkl_poisson.h"

int main(void)
{
  /* Note that the size of the transform np must be even !!! */
  MKL_INT np=8, nt=8;
  double pi=3.14159265358979324;

  MKL_INT ip, it, i, stat;
  MKL_INT ipar[128];
  double ap, bp, at, bt, lp, lt, hp, ht, theta_i, ct, c1;
  double *dpar=NULL, *f=NULL, *u=NULL;
  double q;
  DFTI_DESCRIPTOR_HANDLE handle_s = 0;
  DFTI_DESCRIPTOR_HANDLE handle_c = 0;
  MKL_INT mem_error, error;

  /* Printing the header for the example */
  printf("\n Example of use of Intel MKL Poisson Library\n");
  printf(" **********************************************\n\n");
  printf(" This example gives the solution of Helmholtz problem on a whole sphere\n");
  printf(" 0<p<2*pi, 0<t<pi, with Helmholtz coefficient q=1 and right-hand side\n");
  printf(" f(p,t)=3*cos(t)\n");
  printf(" ----------------------------------------------------------------------\n");
  printf(" In general, the error should be of order O(1.0/np^2+1.0/nt^2)\n");
  printf(" For this example, the value of np=nt is %d\n", np);
  printf(" The approximation error should be of order 1.5e-01,");
  printf(" if everything is OK\n");
  printf(" ----------------------------------------------------------------------\n");
  printf("   Note that np should be even to solve the PERIODIC problem!\n");
  printf(" ----------------------------------------------------------------------\n");
  printf("                      DOUBLE PRECISION COMPUTATIONS                    \n");
  printf(" ======================================================================\n");
  printf("\n");

  error = 0;
  /* memory allocation */
  mem_error = 1;
  dpar=(double*)mkl_malloc((5*np/2+nt+10)*sizeof(double), 64);
  if(dpar == NULL) goto end;
  f=(double*)mkl_malloc((np+1)*(nt+1)*sizeof(double), 64);
  if(f == NULL) goto end;
  u=(double*)mkl_malloc((np+1)*(nt+1)*sizeof(double), 64);
  if(u == NULL) goto end;
  /* memory allocated correctly */
  mem_error = 0;

  /* Defining the rectangular domain on a sphere 0<p<2*pi, 0<t<pi
  for Helmholtz Solver on a sphere
  Poisson Library will automatically detect that this problem is on a whole sphere! */
  ap=0.0E0;
  bp=2*pi;
  at=0.0E0;
  bt=pi;

  /* Setting the coefficient q to 1.0E0 for Helmholtz problem
  If you like to solve Poisson problem, please set q to 0.0E0 */
  q=1.0E0;

  /* Computing the mesh size hp in phi-direction */
  lp=bp-ap;
  hp=lp/np;
  /* Computing the mesh size ht in theta-direction */
  lt=bt-at;
  ht=lt/nt;

  /* Filling in the values of the TRUE solution u(p,t)=cos(t)
  in the mesh points into the array u
  Filling in the right-hand side f(p,t)=(2+q)*cos(t)
  in the mesh points into the array f.
  We choose the right-hand side to correspond to the TRUE solution
  of Helmholtz equation on a sphere.
  Here we are using the mesh sizes hp and ht computed before to compute
  the coordinates (phi_i,theta_i) of the mesh points */
  for(it=0;it<=nt;it++)
  {
    for(ip=0;ip<=np;ip++)
    {
      theta_i=ht*it;
      ct=cos(theta_i);
      u[ip+it*(np+1)]=ct;
      f[ip+it*(np+1)]=ct*(2.+q);
    }
  }


  /* Initializing ipar array to make it free from garbage */
  for(i=0;i<128;i++)
  {
    ipar[i]=0;
  }

  /* Initializing simple data structures of Poisson Library
  for Helmholtz Solver on a sphere
  As we are looking for the solution on a whole sphere, this is a PERIODIC problem
  Therefore, the routines ending with "_p" are used to find the solution */
  d_init_sph_p(&ap,&bp,&at,&bt,&np,&nt,&q,ipar,dpar,&stat);
    if (stat!=0){ 
		error = 1;
		goto end;
	}

  /* Initializing complex data structures of Poisson Library
  for Helmholtz Solver on a sphere
  NOTE: Right-hand side f may be altered after the Commit step. If you want
  to keep it, you should save it in another memory location! */
  d_commit_sph_p(f,&handle_s,&handle_c,ipar,dpar,&stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Computing the approximate solution of Helmholtz problem on a whole sphere */
  d_sph_p(f,&handle_s,&handle_c,ipar,dpar,&stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Cleaning the memory used by handle_s and handle_c */
  free_sph_p(&handle_s,&handle_c,ipar,&stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Now we can use handle_s and handle_c to solve another Helmholtz problem */
  /* after a proper initialization */

  /* Printing the results */
  printf("The number of mesh intervals in phi-direction is np=%d\n", np);
  printf("The number of mesh intervals in theta-direction is nt=%d\n\n", nt);

  /* Watching the error along the line phi=hp */
  ip=1;
  c1=0.0;
  for(it=0;it<=nt;it++)
  {
    printf("In the mesh point (%5.3f,%5.3f) ", ip*hp, it*ht);
    printf("the error between the computed and the true solution ");
    printf("is equal to %10.3e\n", f[ip+it*(np+1)]-u[ip+it*(np+1)]);
    if (c1<fabs(f[ip+it*(np+1)]-u[ip+it*(np+1)]))
        c1=fabs(f[ip+it*(np+1)]-u[ip+it*(np+1)]);
  }
  if (c1>1.5e-01)
  {
    printf("The computed solution seems to be inaccurate. ");
	error = 1;
    goto end;
  }

end:
    /* Free Intel MKL memory if any was allocated */
	mkl_free (dpar);
	mkl_free (f);
	mkl_free (u);
  MKL_Free_Buffers();
/* Failure message to print if something went wrong */
	if (mem_error == 1) 
	{
		printf ("| insufficient memory \n");
		return 1;
	}
	if(error != 0)
	{
		 printf("\nDouble precision Helmholtz example on a whole sphere has ");
         printf("FAILED to compute the solution...\n");
         return 1;
	}
  /* Success message to print if everything is OK */
  printf("\n Double precision Helmholtz example on a whole sphere has ");
  printf("successfully PASSED\n through all steps of computation!\n");
  return 0;
  /* End of the example code */
}
