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
!  Fortran-90 double precision example of solving 2D Helmholtz problem in a
!  rectangular domain using Intel(R) Math Kernel Library (Intel(R) MKL)
!  Poisson Library
!
!*******************************************************************************

program Helmholtz_2D_double_precision

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer nx,ny

parameter(nx=6, ny=6)
double precision pi
parameter(pi=3.14159265358979324D0)

integer ix, iy, i, stat
integer ipar(128)
double precision ax, bx, ay, by, lx, ly, hx, hy, xi, yi, cx, cy, c1
double precision dpar(13*nx/2+7)
! Note that proper packing of data in right-hand side array f is
! automatically provided by the following declarations of the array
double precision f(nx+1,ny+1), u(nx+1,ny+1)
double precision bd_ax(ny+1), bd_bx(ny+1), bd_ay(nx+1), bd_by(nx+1)
double precision q
type(DFTI_DESCRIPTOR), pointer :: xhandle
character(4) BCtype

! Printing the header for the example
print *, ''
print *, ' Example of use of Intel MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of 2D Helmholtz problem'
print *, ' with the equation -u_xx-u_yy+u=f(x,y), 0<x<1, 0<y<1,'
print *, ' f(x,y)=(8*pi*pi+1)*sin(2*pi*x)*sin(2*pi*y)+1,'
print *, ' and with the following boundary conditions:'
print *, '  u(0,y)=u(1,y)=1 (Dirichlet boundary conditions),'
print *, ' -u_y(x,0)=-2.0*pi*sin(2*pi*x) (Neumann boundary condition),'
print *, '  u_y(x,1)= 2.0*pi*sin(2*pi*x) (Neumann boundary condition).'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/nx^2+1.0/ny^2)'
print '(1x,a,I1)', ' For this example, the value of nx=ny is ', nx
print *, ' The approximation error should be of order 0.1D+0, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, '                      DOUBLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the rectangular domain 0<x<1, 0<y<1 for 2D Helmholtz Solver
ax=0.0D0
bx=1.0D0
ay=0.0D0
by=1.0D0

! Setting the coefficient q to 1
q=1.0D0

! Computing the mesh size hx in x-direction
lx=bx-ax
hx=lx/nx
! Computing the mesh size hy in y-direction
ly=by-ay
hy=ly/ny

! Filling in the values of the TRUE solution u(x,y)=sin(2*pi*x)*sin(2*pi*y)+1
! in the mesh points into the array u
! Filling in the right-hand side f(x,y)=(8*pi*pi+q)*sin(2*pi*x)*sin(2*pi*y)+q
! in the mesh points into the array f.
! We choose the right-hand side to correspond to the TRUE solution of
! Helmholtz equation.
! Here we are using the mesh sizes hx and hy computed before to compute
! the coordinates (xi,yi) of the mesh points
do iy=1,ny+1
   do ix=1,nx+1
      xi=hx*(ix-1)/lx
      yi=hy*(iy-1)/ly

      cx=dsin(2*pi*xi)
      cy=dsin(2*pi*yi)

      u(ix,iy)=1.0D0*cx*cy
      f(ix,iy)=((8.0D0*pi**2)+q)*u(ix,iy)+q
      u(ix,iy)=u(ix,iy)+1.0D0
   enddo
enddo

! Setting the type of the boundary conditions on each side of the rectangular domain:
! On the boundary laying on the line x=0(=ax) Dirichlet boundary condition
! will be used
! On the boundary laying on the line x=1(=bx) Dirichlet boundary condition
! will be used
! On the boundary laying on the line y=0(=ay) Neumann boundary condition will be used
! On the boundary laying on the line y=1(=by) Neumann boundary condition will be used
BCtype = 'DDNN'

! Setting the values of the boundary function G(x,y) that is equal to
! the TRUE solution in the mesh points laying on Dirichlet boundaries
do iy = 1,ny+1
   bd_ax(iy) = 1.0D0
   bd_bx(iy) = 1.0D0
enddo
! Setting the values of the boundary function g(x,y) that is equal to
! the normal derivative of the TRUE solution in the mesh points laying on
! Neumann boundaries
do ix = 1,nx+1
   bd_ay(ix) = -2.0*pi*dsin(2*pi*(ix-1)/nx)
   bd_by(ix) =  2.0*pi*dsin(2*pi*(ix-1)/nx)
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for 2D Helmholtz Solver
call d_init_Helmholtz_2D(ax, bx, ay, by, nx, ny, BCtype, q, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library for 2D Helmholtz Solver
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of 2D Helmholtz problem
! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by
! should not be changed between the Commit step and the subsequent call to
! the Solver routine! Otherwise the results may be wrong.
call d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by xhandle
call free_Helmholtz_2D(xhandle, ipar, stat)
if (stat.ne.0) goto 999
! Now we can use xhandle to solve another 2D Helmholtz problem

! Printing the results
write(*,10) nx
write(*,11) ny
print *, ''
! Watching the error along the line x=hx
ix=2
c1 = 0.0
do iy=1,ny+1
   write(*,12) (ix-1)*hx, (iy-1)*hy, f(ix,iy)-u(ix,iy)
   if (dabs(f(ix,iy)-u(ix,iy)).ge.c1) c1 = dabs(f(ix,iy)-u(ix,iy))
enddo
print *, ''

if (c1.ge.0.1D+0) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free Intel MKL memory if any was allocated
call mkl_free_buffers
! Success message to print if everything is OK
print *, ' Double precision 2D Helmholtz example has successfully PASSED'
print *, ' through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Double precision 2D Helmholtz example FAILED to compute the solution...'
    ! Free Intel MKL memory if any was allocated
    call mkl_free_buffers
    stop 1

10    format(1x,'The number of mesh intervals in x-direction is nx=',I1)
11    format(1x,'The number of mesh intervals in y-direction is ny=',I1)
12    format(1x,'In the mesh point (',F5.3,',',F5.3,') the error between ',&
                'the computed and the true solution is equal to ', E10.3)

! End of the example code
end
