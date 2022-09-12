!===============================================================================
! Copyright 2011-2020 Intel Corporation.
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

! Content:
!       Example of using dfftw_plan_dft_r2c_2d function.
!
!*****************************************************************************

program dp_plan_dft_r2c_2d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 2D transform
  integer, parameter :: N1 = 32
  integer, parameter :: N2 = 10

  ! Arbitrary harmonic used to verify computation
  integer, parameter :: H1 = 1, H2 = 2

  ! Working precision is double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0

  ! Data arrays
  complex(WP), allocatable :: x_cmplx(:,:)
  real(WP), allocatable :: x_real(:,:)

  ! FFTW plan
  integer*8 :: plan_r2c = 0


  print *,"Example dp_plan_dft_r2c_2d"
  print *,"2D real-to-complex out-of-place transform"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0"]")', N1, N2
  print '("  H = ["I0","I0"]")', H1, H2

  print *,"Allocate data arrays"
  allocate ( x_real (N1    , N2), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx(INT(N1/2.0)+1, N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW real-to-complex plan"
  call dfftw_plan_dft_r2c_2d(plan_r2c, N1, N2, x_real, x_cmplx, FFTW_ESTIMATE)
  if (0 == plan_r2c) goto 999

  print *,"Initialize data for r2c transform"
  call init(x_real, N1, N2, H1, H2)

  print *,"Compute real-to-complex transform"
  call dfftw_execute(plan_r2c)

  print *,"Verify real-to-complex transform"
  status = verificate(x_cmplx, N1, N2, H1, H2)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plan"
  call dfftw_destroy_plan(plan_r2c)

  print *,"Free data arrays"
  deallocate(x_real, x_cmplx)

  if (status == 0) then
    print *, "TEST PASSED"
    call exit(0)
  else
    print *, "TEST FAILED"
    call exit(1)
  end if

999 print '("  Error, status = ",I0)', status
  status = 1
  goto 100

contains

  ! Compute mod(K*L,M) accurately
  pure integer*8 function moda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    moda = mod(k8*l,m)
  end function moda

  ! Initialize x(:,:) to harmonic H
  subroutine init(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)

    integer k1, k2
    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP

    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = 2 * cos( TWOPI * ( &
        real  (moda(H1,k1-1,N1),WP) / N1 &
        + real(moda(H2,k2-1,N2),WP) / N2)) / (N1*N2)
    end forall
    if (mod(H1,N1)==0 .and. mod(H2,N2)==0) then
      x(1:N1,1:N2) = x(1:N1,1:N2) / 2
    end if
  end subroutine init

  ! Verify that x(:,:) contains harmonic (H1,H2)
  integer function verificate(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0
    do k2 = 1, N2
      do k1 = 1, N1/2+1
        if (mod(k1-1-H1,N1)==0 .and. mod(k2-1-H2,N2)==0) then
          res_exp = 1.0
        else if (mod(-k1+1-H1,N1)==0 .and. mod(-k2+1-H2,N2)==0) then
          res_exp = 1.0
        else
          res_exp = 0.0
        end if
        res_got = x(k1, k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1,k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', x(k1,k2)
          print '(" err "G10.3)', err
          print *,"  Verification FAILED"
          verificate = 1
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program dp_plan_dft_r2c_2d
