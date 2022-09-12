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
!  Fortran-90 double precision example of solving 3D Laplace problem in a
!  parallelepiped domain using Intel(R) Math Kernel Library (Intel(R) MKL)
!  Poisson Library
!
!*******************************************************************************

program Laplace_3D_double_precision

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer nx, ny, nz
parameter(nx=4, ny=4, nz=4)

integer ix, iy, iz, i, stat
integer ipar(128)
double precision ax, bx, ay, by, az, bz, lx, ly, lz, hx, hy, hz, xi, yi, zi, c1
double precision dpar(13*(nx+ny)/2+9)
! Note that proper packing of data in right-hand side array f and boundary
! arrays bd_* is automatically provided by the following declarations of arrays
double precision f(nx+1,ny+1,nz+1), u(nx+1,ny+1,nz+1)
double precision bd_ax(ny+1,nz+1), bd_bx(ny+1,nz+1)
double precision bd_ay(nx+1,nz+1), bd_by(nx+1,nz+1)
double precision bd_az(nx+1,ny+1), bd_bz(nx+1,ny+1)
double precision q
type(DFTI_DESCRIPTOR), pointer :: xhandle, yhandle
character(6) BCtype

! Printing the header for the example
print *, ''
print *, ' Example of use of Intel MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of 3D Laplace equation'
print *, ' u_xx+u_yy+u_zz=0, 0<x<1, 0<y<1, 0<z<1,'
print *, ' with the following Dirichlet boundary conditions:'
print *, ' u(0,y,z)=  y^4+z^4-3* (y^2)*(z^2)'
print *, ' u(1,y,z)=1+y^4+z^4-3*((y^2)*(z^2)+y^2+z^2)'
print *, ' u(x,0,z)=x^4+  z^4-3* (x^2)*(z^2)'
print *, ' u(x,1,z)=x^4+1+z^4-3*((x^2)*(z^2)+x^2+z^2)'
print *, ' u(x,y,0)=x^4+y^4-  3* (x^2)*(y^2)'
print *, ' u(x,y,1)=x^4+y^4+1-3*((x^2)*(y^2)+x^2+y^2)'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/nx^2+1.0/ny^2+1.0/nz^2)'
print '(1x,a,I1)', ' For this example, the value of nx=ny=nz is ', nx
print *, ' The approximation error should be of order 0.2D-1, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, '                      DOUBLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the parallelepiped domain 0<x<1, 0<y<1, 0<z<1 for 3D Laplace Solver
ax=0.0D0
bx=1.0D0
ay=0.0D0
by=1.0D0
az=0.0D0
bz=1.0D0

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the first of two steps on the way to use Helmholtz Solver
! to solve Laplace problem!
!*******************************************************************************
q=0.0D0

! Computing the mesh size hx in x-direction
lx=bx-ax
hx=lx/nx
! Computing the mesh size hy in y-direction
ly=by-ay
hy=ly/ny
! Computing the mesh size hx in z-direction
lz=bz-az
hz=lz/nz

! Filling in the values of the TRUE solution
! u(x,y,z)=x^4+y^4+z^4-3*[(x^2)*(y^2)+(x^2)*(z^2)+(y^2)*(z^2)]
! in the mesh points into the array u
!*******************************************************************************
! Filling in the right-hand side f(x,y,z)=0.0 in the mesh points into the array f.
! Note that this the second of two steps on the way to use Helmholtz Solver
! to solve Laplace problem! Current implementation of Poisson Library requires
! the array f for the solution of Laplace problem
!*******************************************************************************
! Here we are using the mesh sizes hx and hy computed before to compute
! the coordinates (xi,yi,zi) of the mesh points
do iz=1,nz+1
   do iy=1,ny+1
      do ix=1,nx+1
         xi=hx*(ix-1)/lx
         yi=hy*(iy-1)/ly
         zi=hz*(iz-1)/lz

         u(ix,iy,iz)=xi**4+yi**4+zi**4-3.0D0*((xi**2)*(yi**2)+(xi**2)*(zi**2)+&
                                                                     (yi**2)*(zi**2))
         f(ix,iy,iz)=0.0D0
      enddo
   enddo
enddo

! Setting the type of the boundary conditions on each surface
! of the parallelepiped domain:
! On the boundary laying on the plane x=0(=ax) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane x=1(=bx) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane y=0(=ay) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane y=1(=by) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane z=0(=az) Dirichlet boundary condition
! will be used
! On the boundary laying on the plane z=1(=bz) Dirichlet boundary condition
! will be used
BCtype = 'DDDDDD'

! Setting the values of the boundary function G(x,y,z) that is equal to
! the TRUE solution in the mesh points laying on Dirichlet boundaries
do iy = 1,ny+1
   do iz = 1,nz+1
      yi=hy*(iy-1)/ly
      zi=hz*(iz-1)/lz
      bd_ax(iy,iz) = yi**4+zi**4-3.0D0*(yi**2)*(zi**2)
      bd_bx(iy,iz) = 1.0D0+yi**4+zi**4-3.0D0*(yi**2+zi**2+(yi**2)*(zi**2))
   enddo
enddo
do ix = 1,nx+1
   do iz = 1,nz+1
      xi=hx*(ix-1)/lx
      zi=hz*(iz-1)/lz
      bd_ay(ix,iz) = xi**4+zi**4-3.0D0*(xi**2)*(zi**2)
      bd_by(ix,iz) = xi**4+1.0D0+zi**4-3.0D0*(xi**2+(xi**2)*(zi**2)+zi**2)
   enddo
enddo
do ix = 1,nx+1
   do iy = 1,ny+1
      xi=hx*(ix-1)/lx
      yi=hy*(iy-1)/ly
      bd_az(ix,iy) = xi**4+yi**4-3.0D0*(xi**2)*(yi**2)
      bd_bz(ix,iy) = xi**4+yi**4+1.0D0-3.0D0*((xi**2)*(yi**2)+xi**2+yi**2)
   enddo
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for 3D Laplace Solver
call d_init_Helmholtz_3D(ax, bx, ay, by, az, bz, nx, ny, nz, BCtype, q, ipar, dpar,&
                                                                                 stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library for 3D Laplace Solver
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call d_commit_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                            yhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of 3D Laplace problem
! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz
! should not be changed between the Commit step and the subsequent call to
! the Solver routine! Otherwise the results may be wrong.
call d_Helmholtz_3D(f, bd_ax, bd_bx, bd_ay, bd_by, bd_az, bd_bz, xhandle,&
                                                            yhandle, ipar, dpar, stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by xhandle and yhandle
call free_Helmholtz_3D(xhandle, yhandle, ipar, stat)
if (stat.ne.0) goto 999
! Now we can use xhandle and yhandle to solve another 3D Laplace problem

! Printing the results
write(*,9)  nx
write(*,10) ny
write(*,11) nz
! Watching the error in the plane x=hx
ix=2
c1 = 0.0
do iy=1,ny+1
   do iz=1,nz+1
      write(*,12) (ix-1)*hx, (iy-1)*hy, (iz-1)*hz, f(ix,iy,iz)-u(ix,iy,iz)
      if (dabs(f(ix,iy,iz)-u(ix,iy,iz)).ge.c1) c1 = dabs(f(ix,iy,iz)-u(ix,iy,iz))
   enddo
enddo
print *, ''

if (c1.ge.0.2D-1) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free Intel MKL memory if any was allocated
call mkl_free_buffers
! Success message to print if everything is OK
print *, ' Double precision 3D Laplace example has successfully PASSED'
print *, ' through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Double precision 3D Laplace example FAILED to compute the solution...'
    ! Free Intel MKL memory if any was allocated
    call mkl_free_buffers
    stop 1

9     format(1x,'The number of mesh intervals in x-direction is nx=',I1)
10    format(1x,'The number of mesh intervals in y-direction is ny=',I1)
11    format(1x,'The number of mesh intervals in z-direction is nz=',I1)
12    format(1x,'In the mesh point (',F4.2,',',F4.2,',',F4.2,') the error between ',&
                'the computed and the true solution is equal to ', E10.3)

! End of the example code
end
