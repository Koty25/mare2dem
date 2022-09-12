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
!       Example of using dfftw_plan_dft_r2c function.
!
!*****************************************************************************

program dp_plan_dft_r2c

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 4D transform
  integer :: N(4) = [32,16,8,4]

  ! Arbitrary harmonic used to test FFT
  integer :: H(4) = [1,-1,2,-2]

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0

  ! Data arrays
  complex(WP), allocatable :: x_cmplx(:,:,:,:)
  real(WP), allocatable :: x_real(:,:,:,:)

  ! FFTW plan handle
  integer*8 :: plan_r2c = 0

  print *,"Example dp_plan_dft_r2c"
  print *,"4D real-to-complex FFT"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0","I0"]")', N(1:4)
  print '("  H = ["I0","I0","I0","I0"]")', H(1:4)

  print *,"Allocate data arrays"
  allocate ( x_real (N(1),     N(2), N(3), N(4)), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx(N(1)/2+1, N(2), N(3), N(4)), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW real-to-complex plan for 4D FFT"
  call dfftw_plan_dft_r2c(plan_r2c, 4, N(1:4), x_real, x_cmplx, FFTW_ESTIMATE)
  if (0 == plan_r2c) goto 999

  print *,"Initialize data for real-to-complex FFT"
  call init(x_real, N, H)

  print *,"Compute real-to-complex transform"
  call dfftw_execute(plan_r2c)

  print *,"Verify real-to-complex transform"
  status = verificate(x_cmplx, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plan"
  call dfftw_destroy_plan(plan_r2c)

  print *,"Deallocate data arrays"
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

  ! Initialize x(N(:)) to produce unit peaks at x(H) and x(N-H)
  subroutine init(x, N, H)
    integer N(4), H(4)
    real(WP) :: x(:, :, :, :)

    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP
    integer k1,k2,k3,k4

    forall (k1=1:N(1), k2=1:N(2), k3=1:N(3), k4=1:N(4))
      x(k1,k2,k3,k4) = 2 * cos( TWOPI * ( &
        real  (moda(H(1),k1-1,N(1)),WP) / N(1) &
        + real(moda(H(2),k2-1,N(2)),WP) / N(2) &
        + real(moda(H(3),k3-1,N(3)),WP) / N(3) &
        + real(moda(H(4),k4-1,N(4)),WP) / N(4))) / product(N)
    end forall
    if (all(mod(H,N)==0)) then
      x(:,:,:,:) =  x(:,:,:,:) / 2
    end if
  end subroutine init

  ! Verify that x(N) is unit peak at H and N-H
  integer function verificate(x, N, H)
    integer :: N(4), H(4)
    complex(WP) :: x(:,:,:,:)

    integer k1, k2, k3, k4
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(product(N),WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0
    do k4 = 1, N(4)
      do k3 = 1, N(3)
        do k2 = 1, N(2)
          do k1 = 1, N(1)/2+1
            if (mod(k1-1-H(1),N(1))==0 .and. &
              mod  (k2-1-H(2),N(2))==0 .and. &
              mod  (k3-1-H(3),N(3))==0 .and. &
              mod  (k4-1-H(4),N(4))==0) then
              res_exp = 1.0
            else if (mod(-k1+1-H(1),N(1))==0 .and. &
              mod       (-k2+1-H(2),N(2))==0 .and. &
              mod       (-k3+1-H(3),N(3))==0 .and. &
              mod       (-k4+1-H(4),N(4))==0) then
              res_exp = 1.0
            else
              res_exp = 0.0
            end if
            res_got = x(k1,k2,k3,k4)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0","I0"):"$)', k1, k2, k3, k4
              print '(" expected ("G24.17","G24.17"),"$)', res_exp
              print '(" got ("G24.17","G24.17"),"$)', res_got
              print '(" err "G10.3)', err
              print *,"  Verification FAILED"
              verificate = 1
              return
            end if
          end do
        end do
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program dp_plan_dft_r2c
