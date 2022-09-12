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
!  Fortran-90 single precision example of solving Helmholtz problem on a whole
!  sphere excluding poles using Intel(R) Math Kernel Library (Intel(R) MKL)
!  Poisson Library
!
!*******************************************************************************

program s_sph_without_poles_p_f

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none

integer np,nt
! Note that the size of the transform np must be even !!!
parameter(np=8,nt=8)
real pi
parameter(pi=3.14159265358979324E0)

real ap,bp,hp,at,bt,ht,q,lp,lt,theta_i,ct,c1
real  u(np+1,nt+1),f(np+1,nt+1)
type(DFTI_DESCRIPTOR), pointer :: handle_s, handle_c
integer stat
real spar(5*np/2+nt+10)
integer ip,it,i
integer ipar(128)

! Printing the header for the example
print *, ''
print *, ' Example of use of Intel MKL Poisson Library'
print *, ' **********************************************'
print *, ''
print *, ' This example gives the solution of Helmholtz problem on a sphere'
print *, ' excluding poles 0<p<2*pi, 0.1<t<pi-0.1, with Helmholtz coefficient q=1'
print *, ' and right-hand side f(p,t)=(t-0.1)*(pi-0.1-t)+2-(pi-2*t)*cos(t)/sin(t)'
print *, ' -----------------------------------------------------------------------'
print *, ' In general, the error should be of order O(1.0/np^2+1.0/nt^2)'
print '(1x,a,I1)', ' For this example, the value of np=nt is ', np
print *, ' The approximation error should be of order 0.35E-1, if everything is OK'
print *, ' -----------------------------------------------------------------------'
print *, ' Note that np should be even to solve the PERIODIC problem!'
print *, ' -----------------------------------------------------------------------'
print *, '                      SINGLE PRECISION COMPUTATIONS                     '
print *, ' ======================================================================='
print *, ''

! Defining the rectangular domain on a sphere 0<p<2*pi, 0.1<t<pi-0.1
! for Helmholtz Solver on a sphere
! Poisson Library will automatically detect that this problem is
! PERIODIC along phi-direction only
ap=0.0E0
bp=2*pi
at=0.1E0
bt=pi-0.1E0

! Setting the coefficient q to 1.0E0 for Helmholtz problem
! If you like to solve Poisson problem, please set q to 0.0E0
q=1.0E0

! Computing the mesh size hp in phi-direction
lp=bp-ap
hp=lp/np
! Computing the mesh size ht in theta-direction
lt=bt-at
ht=lt/nt

! Filling in the values of the TRUE solution u(p,t)=(t-0.1)*(pi-0.1-t)
! in the mesh points into the array u
! Filling in the right-hand side f(p,t)=(t-0.1)*(pi-0.1-t)+2-(pi-2*t)*cos(t)/sin(t)
! in the mesh points into the array f.
! We choose the right-hand side to correspond to the TRUE solution
! of Helmholtz equation.
! Here we are using the mesh sizes hp and ht computed before to compute
! the coordinates (phi_i,theta_i) of the mesh points
do it=1,nt+1
   do ip=1,np+1
      theta_i=ht*(it-1)+at
      ct=(theta_i-at)*(bt-theta_i)
      u(ip,it)=ct
      f(ip,it)=q*ct+2.0E0-(cos(theta_i)/sin(theta_i))*(at+bt-2*theta_i)
   enddo
enddo

do ip=1,np+1
   f(ip,   1)=0.0E0
   f(ip,nt+1)=0.0E0
enddo

! Initializing ipar array to make it free from garbage
do i=1,128
   ipar(i)=0
enddo

! Initializing simple data structures of Poisson Library for Helmholtz Solver
! on a sphere
! As we are looking for the solution on a whole interval over phi,
! this is a PERIDOC problem
! Therefore, the routines ending with "_P" are used to find the solution
call S_INIT_SPH_P(ap,bp,at,bt,np,nt,q,ipar,spar,stat)
if (stat.ne.0) goto 999

! Initializing complex data structures of Poisson Library for Helmholtz Solver
! on a sphere
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call S_COMMIT_SPH_P(f,handle_s,handle_c,ipar,spar,stat)
if (stat.ne.0) goto 999

! Computing the approximate solution of Helmholtz problem on a sphere without poles
call S_SPH_P(f,handle_s,handle_c,ipar,spar,stat)
if (stat.ne.0) goto 999

! Cleaning the memory used by handle_s and handle_c
call FREE_SPH_P(HANDLE_S,HANDLE_C,IPAR,STAT)
if (stat.ne.0) goto 999
! Now we can use handle_s and handle_c to solve another Helmholtz problem after
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
   if (abs(f(ip,it)-u(ip,it)).ge.c1) c1 = abs(f(ip,it)-u(ip,it))
enddo
print *, ''

if (c1.ge.0.35E-1) then
   print *, 'The computed solution seems to be inaccurate.'
   goto 999
endif

! Free Intel MKL memory if any was allocated
call mkl_free_buffers
! Success message to print if everything is OK
print *, ' Single precision Helmholtz example on a sphere without poles has'
print *, ' successfully PASSED through all steps of computation!'
stop 0
! Failure message to print if something went wrong
999 print *, 'Single precision Helmholtz example on a sphere without poles has',&
             'FAILED to compute the solution...'
    ! Free Intel MKL memory if any was allocated
    call mkl_free_buffers
    stop 1

10    format(1x,'The number of mesh intervals in phi-direction is np=',I1)
11    format(1x,'The number of mesh intervals in theta-direction is nt=',I1)
12    format(1x,'In the mesh point (',F5.3,',',F5.3,') the error between the ',&
                'computed and the true solution is equal to ', E10.3)

! End of the example code
end
