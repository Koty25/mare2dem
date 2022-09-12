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
*  C double precision example of solving 3D Laplace problem in a
*  parallelepiped domain using Intel(R) Math Kernel Library (Intel(R) MKL)
*  Poisson Library
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

  MKL_INT nx=4, ny=4, nz=4;
  double pi=3.14159265358979324;

  MKL_INT ix, iy, iz, i, stat;
  MKL_INT ipar[128];
  double ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, cx, cy, cz, c1;
  double *dpar=NULL, *f=NULL, *u=NULL, *bd_ax=NULL, *bd_bx=NULL, *bd_ay=NULL, *bd_by=NULL, *bd_az=NULL, *bd_bz=NULL;
  double q;
  DFTI_DESCRIPTOR_HANDLE xhandle = 0;
  DFTI_DESCRIPTOR_HANDLE yhandle = 0;
  char *BCtype;
  MKL_INT mem_error, error;

  /* Printing the header for the example */
  printf("\n Example of use of Intel MKL Poisson Library\n");
  printf(" **********************************************\n\n");
  printf(" This example gives the solution of 3D Laplace equation\n");
  printf(" u_xx+u_yy+u_zz=0, 0<x<1, 0<y<1, 0<z<1,\n");
  printf(" with the following Dirichlet boundary conditions:\n");
  printf(" u(0,y,z)=  y^4+z^4-3* (y^2)*(z^2)\n");
  printf(" u(1,y,z)=1+y^4+z^4-3*((y^2)*(z^2)+y^2+z^2)\n");
  printf(" u(x,0,z)=x^4+  z^4-3* (x^2)*(z^2)\n");
  printf(" u(x,1,z)=x^4+1+z^4-3*((x^2)*(z^2)+x^2+z^2)\n");
  printf(" u(x,y,0)=x^4+y^4-  3* (x^2)*(y^2)\n");
  printf(" u(x,y,1)=x^4+y^4+1-3*((x^2)*(y^2)+x^2+y^2)\n");
  printf(" ----------------------------------------------------------------------\n");
  printf(" In general, the error should be of order O(1.0/nx^2+1.0/ny^2+1.0/nz^2)\n");
  printf(" For this example, the value of nx=ny=nz is %d\n", nx);
  printf(" The approximation error should be of order 2.0e-02,");
  printf(" if everything is OK\n");
  printf(" ----------------------------------------------------------------------\n");
  printf("                      DOUBLE PRECISION COMPUTATIONS                    \n");
  printf(" ======================================================================\n");
  printf("\n");

  error = 0;
  /* memory allocation */
  mem_error = 1;
  dpar=(double*)mkl_malloc((13*(nx+ny)/2+9)*sizeof(double), 64);
  if(dpar == NULL) goto end;
  f=(double*)mkl_malloc((nx+1)*(ny+1)*(nz+1)*sizeof(double), 64);
  if(f == NULL) goto end;
  u=(double*)mkl_malloc((nx+1)*(ny+1)*(nz+1)*sizeof(double), 64);
  if(u == NULL) goto end;
  bd_ax=(double*)mkl_malloc((ny+1)*(nz+1)*sizeof(double), 64);
  if(bd_ax == NULL) goto end;
  bd_bx=(double*)mkl_malloc((ny+1)*(nz+1)*sizeof(double), 64);
  if(bd_bx == NULL) goto end;
  bd_ay=(double*)mkl_malloc((nx+1)*(nz+1)*sizeof(double), 64);
  if(bd_ay == NULL) goto end;
  bd_by=(double*)mkl_malloc((nx+1)*(nz+1)*sizeof(double), 64);
  if(bd_by == NULL) goto end;
  bd_az=(double*)mkl_malloc((nx+1)*(ny+1)*sizeof(double), 64);
  if(bd_az == NULL) goto end;
  bd_bz=(double*)mkl_malloc((nx+1)*(ny+1)*sizeof(double), 64);
  if(bd_bz == NULL) goto end;
  /* memory allocated correctly */
  mem_error = 0;

  /* Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Laplace Solver */
  ax=0.0E0;
  bx=1.0E0;
  ay=0.0E0;
  by=1.0E0;
  az=0.0E0;
  bz=1.0E0;

  /*******************************************************************************
  Setting the coefficient q to 0.
  Note that this is the first of two steps on the way to use Helmholtz Solver
  to solve Laplace problem!
  *******************************************************************************/
  q=0.0E0;

  /* Computing the mesh size hx in x-direction */
  lx=bx-ax;
  hx=lx/nx;
  /* Computing the mesh size hy in y-direction */
  ly=by-ay;
  hy=ly/ny;
  /* Computing the mesh size hx in z-direction */
  lz=bz-az;
  hz=lz/nz;

  /* Filling in the values of the TRUE solution
  u(x,y,z)=x^4+y^4+z^4-3*[(x^2)*(y^2)+(x^2)*(z^2)+(y^2)*(z^2)]
  in the mesh points into the array u
  *******************************************************************************
  Filling in the right-hand side f(x,y,z)=0.0 in the mesh points into the array f.
  Note that this the second of two steps on the way to use Helmholtz Solver
  to solve Laplace problem! Current implementation of Poisson Library requires
  the array f for the solution of Laplace problem
  *******************************************************************************
  Here we are using the mesh sizes hx and hy computed before to compute
  the coordinates (xi,yi,zi) of the mesh points */
  for(iz=0;iz<=nz;iz++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      for(ix=0;ix<=nx;ix++)
      {
        xi=hx*ix/lx;
        yi=hy*iy/ly;
        zi=hz*iz/lz;

        u[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]=xi*xi*xi*xi+yi*yi*yi*yi+zi*zi*zi*zi-
                                  3.0E0*(xi*xi*yi*yi+xi*xi*zi*zi+yi*yi*zi*zi);
        f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]=0.0E0;
      }
    }
  }

  /* Setting the type of the boundary conditions on each surface
  of the parallelepiped domain:
  On the boundary laying on the plane x=0(=ax) Dirichlet boundary condition
  will be used
  On the boundary laying on the plane x=1(=bx) Dirichlet boundary condition
  will be used
  On the boundary laying on the plane y=0(=ay) Dirichlet boundary condition
  will be used
  On the boundary laying on the plane y=1(=by) Dirichlet boundary condition
  will be used
  On the boundary laying on the plane z=0(=az) Dirichlet boundary condition
  will be used
  On the boundary laying on the plane z=1(=bz) Dirichlet boundary condition
  will be used */
  BCtype = "DDDDDD";

  /* Setting the values of the boundary function G(x,y,z) that is equal to
  the TRUE solution in the mesh points laying on Dirichlet boundaries */
  for(iy=0;iy<=ny;iy++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      yi=hy*iy/ly;
      zi=hz*iz/lz;
      bd_ax[iy+iz*(ny+1)]=yi*yi*yi*yi+zi*zi*zi*zi-3.0E0*yi*yi*zi*zi;
      bd_bx[iy+iz*(ny+1)]=1.0E0+yi*yi*yi*yi+zi*zi*zi*zi-
                         3.0E0*(yi*yi+zi*zi+yi*yi*zi*zi);
    }
  }
  for(ix=0;ix<=nx;ix++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      xi=hx*ix/lx;
      zi=hz*iz/lz;
      bd_ay[ix+iz*(nx+1)]=xi*xi*xi*xi+zi*zi*zi*zi-3.0E0*xi*xi*zi*zi;
      bd_by[ix+iz*(nx+1)]=xi*xi*xi*xi+1.0E0+zi*zi*zi*zi-
                   3.0E0*(xi*xi+xi*xi*zi*zi+zi*zi);
    }
  }
  for(ix=0;ix<=nx;ix++)
  {
    for(iy=0;iy<=ny;iy++)
    {
      xi=hx*ix/lx;
      yi=hy*iy/ly;
      bd_az[ix+iy*(nx+1)]=xi*xi*xi*xi+yi*yi*yi*yi-3.0E0*xi*xi*yi*yi;
      bd_bz[ix+iy*(nx+1)]=xi*xi*xi*xi+yi*yi*yi*yi+1.0E0-
                   3.0E0*(xi*xi*yi*yi+xi*xi+yi*yi);
    }
  }

  /* Initializing ipar array to make it free from garbage */
  for(i=0;i<128;i++)
  {
    ipar[i]=0;
  }

  /* Initializing simple data structures of Poisson Library for 3D Laplace Solver */
  d_init_Helmholtz_3D(&ax, &bx, &ay, &by, &az, &bz, &nx, &ny, &nz, BCtype, &q,
                                                                   ipar, dpar, &stat);
    if (stat!=0){ 
		error = 1;
		goto end;
	}

  /* Initializing complex data structures of Poisson Library for 3D Laplace Solver
  NOTE: Right-hand side f may be altered after the Commit step. If you want
  to keep it, you should save it in another memory location! */
  d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, &xhandle,
                                                         &yhandle, ipar, dpar, &stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Computing the approximate solution of 3D Laplace problem
  NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz
  should not be changed between the Commit step and the subsequent call to
  the Solver routine! Otherwise the results may be wrong. */
  d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, &xhandle,
                                                         &yhandle, ipar, dpar, &stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Cleaning the memory used by xhandle and yhandle */
  free_Helmholtz_3D(&xhandle, &yhandle, ipar, &stat);
    if (stat!=0){
		error = 1;
		goto end;
	}
  /* Now we can use xhandle and yhandle to solve another 3D Laplace problem  */

  /* Printing the results */
  printf("The number of mesh intervals in x-direction is nx=%d\n", nx);
  printf("The number of mesh intervals in y-direction is ny=%d\n", ny);
  printf("The number of mesh intervals in z-direction is nz=%d\n\n", nz);

  /* Watching the error in the plane x=hx */
  ix=1;
  c1=0.0;
  for(iy=0;iy<=ny;iy++)
  {
    for(iz=0;iz<=nz;iz++)
    {
      printf("In the mesh point (%4.2f,%4.2f,%4.2f) ", ix*hx, iy*hy, iz*hz);
      printf("the error between the computed and the true solution ");
      printf("is equal to %10.3e\n",
                  f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]-u[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]);
      if (c1<fabs(f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]-u[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]))
          c1=fabs(f[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]-u[ix+iy*(nx+1)+iz*(nx+1)*(ny+1)]);
    }
  }
  if (c1>2.0e-02)
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
	mkl_free (bd_ax);
	mkl_free (bd_bx);
	mkl_free (bd_ay);
	mkl_free (bd_by);
	mkl_free (bd_az);
	mkl_free (bd_bz);
	MKL_Free_Buffers();
/* Failure message to print if something went wrong */
	if (mem_error == 1) 
	{
		printf ("| insufficient memory \n");
		return 1;
	}
	if(error != 0)
	{
		printf("\nDouble precision 3D Laplace example FAILED to compute ");
		printf("the solution...\n");
		return 1;
	}
  /* Success message to print if everything is OK */
  printf("\n Double precision 3D Laplace example has successfully PASSED\n");
  printf(" through all steps of computation!\n");
  return 0;
  /* End of the example code */
}
