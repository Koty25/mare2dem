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
!       Example of using dfftw_plan_dft_1d function.
!
!*****************************************************************************

program dp_plan_dft_1d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Size of 1D transform
  integer, parameter :: N = 7

  ! Arbitrary harmonic used to verify FFT
  integer :: H = 1

  ! Working precision is double precision (using dfftw_* functions)
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0

  ! FFTW plan variables
  integer*8 :: fwd = 0, bwd = 0

  ! The data array
  complex(WP), allocatable :: x(:)

  print *,"Example dp_plan_dft_1d"
  print *,"Forward and backward complex 1D in-place FFT"
  print *,"Configuration parameters:"
  print '("  N = "I0)', N
  print '("  H = "I0)', H

  print *,"Allocate data array"
  allocate ( x(N), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW plan for forward in-place transform"
  call dfftw_plan_dft_1d(fwd, N, x, x, FFTW_FORWARD, FFTW_ESTIMATE)
  if (0 == fwd) goto 999

  print *,"Create FFTW plan for backward in-place transform"
  call dfftw_plan_dft_1d(bwd, N, x, x, FFTW_BACKWARD, FFTW_ESTIMATE)
  if (0 == bwd) goto 999

  print *,"Initialize data for forward FFT"
  call init(x, N, H)

  print *,"Compute forward transform"
  call dfftw_execute(fwd)

  print *,"Verify the result of the forward transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Initialize data for backward FFT"
  call init(x, N, -H)

  print *,"Compute backward transform"
  call dfftw_execute(bwd)

  print *,"Verify the result of the backward transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plans"
  call dfftw_destroy_plan(fwd)
  call dfftw_destroy_plan(bwd)

  print *,"Deallocate data array"
  deallocate(x)

  if (status == 0) then
    print *,"TEST PASSED"
    call exit(0)
  else
    print *,"TEST FAILED"
    call exit(1)
  endif

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

  ! Initialize array x(N) with harmonic H
  subroutine init(x, N, H)
    integer N, H
    complex(WP) :: x(:)

    integer k
    complex(WP), parameter :: I_TWOPI = (0,6.2831853071795864769_WP)

    do k = 1, N
      x(k) = exp( I_TWOPI * real(moda(k-1, H, N), WP)/N ) / N
    end do
  end subroutine init

  ! Verify that x(N) is unit peak at x(H)
  integer function verificate(x, N, H)
    integer N, H
    complex(WP) :: x(:)

    integer k
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N, WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0
    do k = 1, N
      if (mod(k-1-H,N)==0) then
        res_exp = 1
      else
        res_exp = 0
      end if
      res_got = x(k)
      err = abs(res_got - res_exp)
      maxerr = max(err,maxerr)
      if (.not.(err < errthr)) then
        print '("  x("I0"): "$)', k
        print '(" expected "G24.17","$)', res_exp
        print '(" got "G24.17","$)', res_got
        print '(" err "G10.3)', err
        print *," Verification FAILED"
        verificate = 1
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program dp_plan_dft_1d
