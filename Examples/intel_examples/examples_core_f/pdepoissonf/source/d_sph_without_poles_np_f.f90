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
!  Fortran-90 precision example of solving Helmholtz problem in a spherical
!  rectangle using Intel(R) Math Kernel Library (Intel(R) MKL) Poisson Library
!
!*******************************************************************************

program d_sph_without_poles_np_f

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer np,nt
parameter(np=8,nt=8)
double precision pi
parameter(pi=3.14159265358979324D0)

double precision ap,bp,hp,at,bt,ht,q,lp,lt,theta_i,phi_i,ct,cp,c1
double precision  u(np+1,nt+1),f(np+1,nt+1)
type(DFTI_DESCRIPTOR), pointer :: handle
integer stat
double precision dpar(5*np/2+nt+10)
integer ip,it,i
integer ipar(128)

! Printing the header for the example
print *, ''
print *, ' Example of use of Intel MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of Helmholtz problem in a sphercal'
print *, ' rectangle 0.1<p<2*pi-0.1, 0.1<t<pi-0.1, with Helmholtz coefficient q=1'
print *, ' and right-hand side f(p,t)=(p-0.1)(2*pi-0.1-p)(t-0.1)(pi-0.1-t)+'
print *, '                        (2-(pi-2t)*cos(t)/sin(t))(p-0.1)(2*pi-0.1-p)+'
print *, '                         2(t-0.1)(pi-0.1-t)/sin(t)/sin(t)'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/np^2+1.0/nt^2)'
print '(1x,a,I1)', ' For this example, the value of np=nt is ', np
print *, ' The approximation error should be of order 1.2D-1, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, '                      DOUBLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the rectangular domain on a sphere 0.1<p<2*pi-0.1, 0.1<t<pi-0.1
! for Helmholtz Solver on a sphere
! Poisson Library will automatically detect that this problem is not on
! a whole sphere and will use the homogeneous Dirichlet boundary conditions
! on the whole boundary
ap=0.1D0
bp=2*pi-0.1D0
at=0.1D0
bt=pi-0.1D0

! Setting the coefficient q to 1.0D0 for Helmholtz problem
! If you like to solve Poisson problem, please set q to 0.0D0
q=1.0D0

! Computing the mesh size hp in phi-direction
lp=bp-ap
hp=lp/np
! Computing the mesh size ht in theta-direction
lt=bt-at
ht=lt/nt

! Filling in the values of the TRUE solution
! u(p,t)=(p-0.1)(2*pi-0.1-p)(t-0.1)(pi-0.1-t)
! in the mesh points into the array u
! Filling in the right-hand side f(p,t)=(p-0.1)(2*pi-0.1-p)(t-0.1)(pi-0.1-t)+
! ( 2-(pi-2t)*cos(t)/sin(t) )(p-0.1)(2*pi-0.1-p)+2(t-0.1)(pi-0.1-t)/sin(t)/sin(t)
! in the mesh points into the array f.
! We choose the right-hand side to correspond to the TRUE solution of
! Helmholtz equation on a sphere.
! Here we are using the mesh sizes hp and ht computed before to compute
! the coordinates (phi_i,theta_i) of the mesh points
do it=1,nt+1
   do ip=1,np+1
      theta_i=ht*(it-1)+at
      phi_i=hp*(ip-1)+ap
      ct=(theta_i-at)*(bt-theta_i)
      cp=(  phi_i-ap)*(bp-  phi_i)
      u(ip,it)=ct*cp
      f(ip,it)=q*ct*cp+(2.0D0-(dcos(theta_i)/dsin(theta_i))*(at+bt-2*theta_i))*cp+ &
                                                2.0D0*ct/(dsin(theta_i)*dsin(theta_i))
   enddo
enddo

do ip=1,np+1
  f(ip,   1)=0.0D0
  f(ip,nt+1)=0.0D0
enddo

do it=1,nt+1
  f(1   ,it)=0.0D0
  f(np+1,it)=0.0D0
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for Helmholtz Solver on
! a sphere
! As we are looking for the solution in a spherical rectangle,
! this is a NON-PERIODIC problem
! Therefore, the routines ending with "_NP" are used to find the solution
call D_INIT_SPH_NP(ap,bp,at,bt,np,nt,q,ipar,dpar,stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library
! for Helmholtz Solver on a sphere
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call D_COMMIT_SPH_NP(f,handle,ipar,dpar,stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of Helmholtz problem in a spherical rectangle
call D_SPH_NP(f,handle,ipar,dpar,stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by handle
call FREE_SPH_NP(handle,ipar,stat)
if (stat.ne.0) goto 999
! Now we can use handle to solve another Helmholtz problem after
! proper initialization

! Printing the results
write(*,10) np
write(*,11) nt
print *, ''
! Watching the error along the line phi=hp
ip=2
c1 = 0.0
do it=1,nt+1
   write(*,12) (ip-1)*hp, (it-1)*ht, f(ip,it)-u(ip,it)
   if (dabs(f(ip,it)-u(ip,it)).ge.c1) c1 = dabs(f(ip,it)-u(ip,it))
enddo
print *, ''

if (c1.ge.1.2D-1) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free Intel MKL memory if any was allocated
call mkl_free_buffers
! Success message to print if everything is OK
print *, ' Double precision Helmholtz example in a spherical rectangle has '
print *, ' successfully PASSED through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Double precision Helmholtz example in a spherical rectangle FAILED ',&
             'to compute the solution...'
    ! Free Intel MKL memory if any was allocated
    call mkl_free_buffers
    stop 1

10    format(1x,'The number of mesh intervals in phi-direction is np=',I1)
11    format(1x,'The number of mesh intervals in theta-direction is nt=',I1)
12    format(1x,'In the mesh point (',F5.3,',',F5.3,') the error between ', &
                'the computed and the true solution is equal to ', E10.3)

! End of the example code
end