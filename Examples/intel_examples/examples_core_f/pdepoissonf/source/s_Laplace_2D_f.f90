!===============================================================================
! Copyright 2006-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!  Fortran-90 single precision example of solving 2D Laplace problem in a
!  rectangular domain using Intel(R) Math Kernel Library (Intel(R) MKL)
!  Poisson Library
!
!*******************************************************************************

program Laplace_2D_single_precision

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer nx,ny
parameter(nx=6, ny=6)

integer ix, iy, i, stat
integer ipar(128)
real ax, bx, ay, by, lx, ly, hx, hy, xi, yi, c1
real spar(13*nx/2+7)
! Note that proper packing of data in right-hand side array f is
! automatically provided by the following declarations of the array
real f(nx+1,ny+1), u(nx+1,ny+1)
real bd_ax(ny+1), bd_bx(ny+1), bd_ay(nx+1), bd_by(nx+1)
real q
type(DFTI_DESCRIPTOR), pointer :: xhandle
character(4) BCtype

! Printing the header for the example
print *, ''
print *, ' Example of use of Intel MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of the 2D Laplace problem'
print *, ' with the equation u_xx+u_yy=0, 0<x<1, 0<y<1,'
print *, ' with the following boundary conditions:'
print *, ' u(0,y)=  y^4'
print *, ' u(1,y)=1+y^4-6*(y^2)'
print *, ' u(x,0)=x^4'
print *, ' u(x,1)=x^4+1-6*(x^2)'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/nx^2+1.0/ny^2)'
print '(1x,a,I1)', ' For this example, the value of nx=ny is ', nx
print *, ' The approximation error should be of order 0.5E-2, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, '                      SINGLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the rectangular domain 0<x<1, 0<y<1 for 2D Laplace Solver
ax=0.0E0
bx=1.0E0
ay=0.0E0
by=1.0E0

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the first of two steps on the way to use Helmholtz Solver
! to solve Laplace problem!
!*******************************************************************************
q=0.0E0

! Computing the mesh size hx in x-direction
lx=bx-ax
hx=lx/nx
! Computing the mesh size hy in y-direction
ly=by-ay
hy=ly/ny

! Filling in the values of the TRUE solution u(x,y)=x^4+y^4-6*(y^2)*(x^2)
! in the mesh points into the array u
!*******************************************************************************
! Filling in the right-hand side f(x,y)=0.0 in the mesh points into the array f.
! Note that this the second of two steps on the way to use Helmholtz Solver
! to solve Laplace problem! Current implementation of Poisson Library requires
! the array f for the solution of Laplace problem
!*******************************************************************************
! Here we are using the mesh sizes hx and hy computed before to compute
! the coordinates (xi,yi) of the mesh points
do iy=1,ny+1
   do ix=1,nx+1
      xi=hx*(ix-1)/lx
      yi=hy*(iy-1)/ly

      u(ix,iy)=xi**4+yi**4-6.0E0*(yi**2)*(xi**2)
      f(ix,iy)=0.0E0
   enddo
enddo

! Setting the type of the boundary conditions on each side of the rectangular domain:
! On the boundary laying on the line x=0(=ax) Dirichlet boundary condition
! will be used
! On the boundary laying on the line x=1(=bx) Dirichlet boundary condition
! will be used
! On the boundary laying on the line y=0(=ay) Dirichlet boundary condition
! will be used
! On the boundary laying on the line y=1(=by) Dirichlet boundary condition
! will be used
BCtype = 'DDDD'

! Setting the values of the boundary function G(x,y) that is equal to
! the TRUE solution in the mesh points laying on Dirichlet boundaries
do iy = 1,ny+1
   yi=hy*(iy-1)/ly
   bd_ax(iy) = yi**4
   bd_bx(iy) = 1.0E0-6.0E0*yi**2+yi**4
enddo
do ix = 1,nx+1
   xi=hx*(ix-1)/lx
   bd_ay(ix) = xi**4
   bd_by(ix) = 1.0E0-6.0E0*xi**2+xi**4
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for 2D Laplace Solver
call s_init_Helmholtz_2D(ax, bx, ay, by, nx, ny, BCtype, q, ipar, spar, stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library for 2D Laplace Solver
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call s_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of 2D Laplace problem
! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by
! should not be changed between the Commit step and the subsequent call to
! the Solver routine! Otherwise the results may be wrong.
call s_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, spar, stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by xhandle
call free_Helmholtz_2D(xhandle, ipar, stat)
if (stat.ne.0) goto 999
! Now we can use xhandle to solve another 2D Laplace problem

! Printing the results
write(*,10) nx
write(*,11) ny
print *, ''
! Watching the error along the line x=hx
ix=2
c1 = 0.0
do iy=1,ny+1
   write(*,12) (ix-1)*hx, (iy-1)*hy, f(ix,iy)-u(ix,iy)
   if (abs(f(ix,iy)-u(ix,iy)).ge.c1) c1 = abs(f(ix,iy)-u(ix,iy))
enddo
print *, ''

if (c1.ge.0.5E-2) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free Intel MKL memory if any was allocated
call mkl_free_buffers
! Success message to print if everything is OK
print *, ' Single precision 2D Laplace example has successfully PASSED'
print *, ' through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Single precision 2D Laplace example FAILED to compute the solution...'
    ! Free Intel MKL memory if any was allocated
    call mkl_free_buffers
    stop 1

10    format(1x,'The number of mesh intervals in x-direction is nx=',I1)
11    format(1x,'The number of mesh intervals in y-direction is ny=',I1)
12    format(1x,'In the mesh point (',F5.3,',',F5.3,') the error between the ',&
                'computed and the true solution is equal to ', E10.3)

! End of the example code
end