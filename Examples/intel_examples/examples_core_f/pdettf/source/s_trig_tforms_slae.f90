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
!  Single precision Fortran90 test example for trigonometric transforms
!*******************************************************************************
!
! This example gives the solution of the linear algebraic system,
! corresponding to the three-point approximation of the 1D differential
! problems with the equation  -u"+u=f(x), 0<x<1, and with 3 types of
! boundary conditions:
! u(0)=u(1)=0 (DD case), or u'(0)=u'(1)=0 (NN case), or u'(0)=u(1)=0 (ND case)

program s_tt_example_slae

  use mkl_dfti
  use mkl_trig_transforms

  implicit none

  include 'mkl_service.fi'

  integer n, i, k, tt_type
  integer ir, ipar(128)
  parameter (n=8)
  real pi, xi, c
  real c1, c2, c3, c4, c5, c6
! NOTE: This example uses shorter spar array of size 3n/2+2 instead of 5n/2+2
! as only sine, cosine, and staggered cosine transforms are used. More details
! can be found in Chapter 13 of Intel(R) Math Kernel Library (Intel(R) MKL) Manual.
  real u(n+1), f(n+1), spar(3*n/2+2), lambda(n+1)
  parameter (pi=3.14159265358979324E0)
  type(dfti_descriptor), pointer :: handle

! Printing the header for the example
  print *, ''
  print *, ' Example of use of Intel MKL Trigonometric Transforms'
  print *, ' **********************************************'
  print *, ''
  print *, ' This example gives the the exact solutions of the linear algebraic '
  print *, ' systems, corresponding to the three-point approximations of the 1D '
  print *, ' differential problems with the equation -u"+u=f(x), 0<x<1, '
  print *, ' and with 3 types of boundary conditions:'
  print *, ' DD case: u(0)=u(1)=0,'
  print *, ' NN case: u''(0)=u''(1)=0,'
  print *, ' ND case: u''(0)=u(1)=0.'
  print *, ' -----------------------------------------------------------------'
  print *, ' The error expected to be close to the machine zero in all 3 cases'
  print *, ' if everything is OK'
  print *, ' -----------------------------------------------------------------'
  print *, '                   SINGLE PRECISION COMPUTATIONS                  '
  print *, ' ================================================================='
  print *, ''

  do i=0,2
! Varying the type of the transform
    tt_type=i

! Computing test solutions u(x)
    do k=1,n+1
      xi=1.0E0*real(k-1,4)/real(n,4)
! The exact solution with Dirichlet-Dirichlet boundary conditions
      if (tt_type.eq.0) c=sin(pi*xi)*exp(xi)
! The exact solution with Neumann-Neumann boundary conditions
      if (tt_type.eq.1) c=(pi*cos(pi*xi)-sin(pi*xi))*exp(xi)
! The exact solution with Neumann-Dirichlet boundary conditions
      if (tt_type.eq.2) c=cos(0.5E0*pi*xi)*(exp(xi)-xi)
        u(k)=c
    end do
! Computing the right-hand side f(x) that corresponds to the finite element
! approximation and provides the exact solution computed above of the resulting
! algebraic system
    do k=2,n
      f(k)=(2.0E0+1.0E0/(real(n,4)**2))*u(k)-u(k-1)-u(k+1)
    end do
    if (tt_type.eq.0) then
! The Dirichlet boundary conditions
      f(1)=0.0E0
      f(n+1)=0.0E0
    end if
    if (tt_type.eq.1) then
! The Neumann boundary conditions
      f(1)=(2.0E0+1.0E0/(real(n,4)**2))*u(1)-2.0E0*u(2)
      f(n+1)=(2.0E0+1.0E0/(real(n,4)**2))*u(n+1)-2.0E0*u(n)
    end if
    if (tt_type.eq.2) then
! The mixed Neumann-Dirichlet boundary conditions
      f(1)=(2.0E0+1.0E0/(real(n,4)**2))*u(1)-2.0E0*u(2)
      f(n+1)=0.0E0
    end if

! Computing the eigenvalues for the three-point finite-difference problem
    if (tt_type.eq.0.or.tt_type.eq.1) then
      do k=1,n+1
        lambda(k)=(2.0E0*sin(0.5E0*pi*real(k-1,4)/real(n,4)))**2+1.0E0/(real(n,4)**2)
      end do
    end if
    if (tt_type.eq.2) then
      do k=1,n+1
        lambda(k)=(2.0E0*sin(0.25E0*pi*real(2*k-1,4)/real(n,4)))**2+1.0E0/(real(n,4)**2)
      end do
    end if

! Computing the solution of 1D problem using trigonometric transforms
! First we initialize the transform
    CALL S_INIT_TRIG_TRANSFORM(n,tt_type,ipar,spar,ir)
    if (ir.ne.0) goto 99
! Then we commit the transform. Note that the data in f will be changed at
! this stage !
! If you want to keep them, save them in some other array before the call to
! the routine
    CALL S_COMMIT_TRIG_TRANSFORM(f,handle,ipar,spar,ir)
    if (ir.ne.0) goto 99
! Now we can apply trigonometric transform
    CALL S_FORWARD_TRIG_TRANSFORM(f,handle,ipar,spar,ir)
    if (ir.ne.0) goto 99

! Scaling the solution by the eigenvalues
    do k=1,n+1
      f(k)=f(k)/lambda(k)
    end do

! Now we can apply trigonometric transform once again as ONLY input vector f
! has changed
    CALL S_BACKWARD_TRIG_TRANSFORM(f,handle,ipar,spar,ir)
    if (ir.ne.0) goto 99
! Cleaning the memory used by handle
! Now we can use handle for other kind of trigonometric transform
    CALL FREE_TRIG_TRANSFORM(handle,ipar,ir)
    if (ir.ne.0) goto 99

! Performing the error analysis
    c1=0.0E0
    c2=0.0E0
    c3=0.0E0
    do k=1,n+1
! Computing the absolute value of the exact solution
      c4=abs(u(k))
! Computing the absolute value of the computed solution
! Note that the solution is now in place of the former right-hand side !
       c5=abs(f(k))
! Computing the absolute error
      c6=abs(f(k)-u(k))
! Computing the maximum among the above 3 values c4-c6
      if (c4.gt.c1) c1=c4
      if (c5.gt.c2) c2=c5
      if (c6.gt.c3) c3=c6
    end do

! Printing the results
    if (tt_type.eq.0) then
      print *, 'The computed solution of DD problem is'
    end if
    if (tt_type.eq.1) then
      print *, 'The computed solution of NN problem is'
    end if
    if (tt_type.eq.2) then
      print *, 'The computed solution of ND problem is'
    endif
    print *, ''
    do k=1,n+1
      write(*,11) k,f(k)
    end do
    print *, ''
    write(*,12) c3/c1
    print *, ''
    if (c3/c1.ge.1.0E-5) then
      print *, 'The computed solution seems to be inaccurate.'
      goto 99
    endif
! End of the loop over the different kind of transforms and problems
  end do

! Free Intel MKL memory if any was allocated
  call mkl_free_buffers
! Success message to print if everything is OK
  print *, 'This example has successfully PASSED through all steps of computation!'
  stop 0

! Failure message to print if something went wrong
99    print *, 'FAILED to compute the solution(s)...'
      ! Free Intel MKL memory if any was allocated
      call mkl_free_buffers
      stop 1

! Print formats
11    format(1x,'u(',I1,')=',F6.3)
12    format(1x,'Error=',E10.3)

! End of the example code
end
