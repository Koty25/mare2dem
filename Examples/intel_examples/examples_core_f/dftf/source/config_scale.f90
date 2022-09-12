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
! An example of using configuration parameters DFTI_FORWARD_SCALE and
! DFTI_BACKWARD_SCALE. Default value for each of these parameters is 1.0.
!
!*****************************************************************************

program config_scale

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Size of 1D transform
  integer, parameter :: N = 32

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H = 1

  ! Working precision is double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! The data array
  complex(WP), allocatable :: x (:)

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  ! Forward and backward scale factors
  real(WP) :: fscale, bscale

  hand => null()

  print *,"Example config_scale"
  print *,"Forward and backward FFT with scaling"
  print *,"Configuration parameters:"
  print '("  DFTI_PRECISION      = DFTI_DOUBLE")'
  print '("  DFTI_FORWARD_DOMAIN = DFTI_COMPLEX")'
  print '("  DFTI_DIMENSION      = 1")'
  print '("  DFTI_LENGTHS        = /"I0"/" )', N

  print *,"Create DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  if (0 /= status) goto 999

  ! Arbitrary scaling factors may be used
  fscale = 1.0_WP / real(N,WP)
  bscale = 1.0_WP / sqrt( real(N,WP) )

  print '(" Set forward  scale = ",G24.17)', fscale
  status = DftiSetValue(hand, DFTI_FORWARD_SCALE, fscale)
  if (0 /= status) goto 999

  print '(" Set backward scale = ",G24.17)', bscale
  status = DftiSetValue(hand, DFTI_BACKWARD_SCALE, bscale)
  if (0 /= status) goto 999

  print *,"Commit the descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Allocate array for input data"
  allocate ( x(N), STAT = status)
  if (0 /= status) goto 999

  print *,"Initialize input for forward transform"
  call init(x, N, H, fscale)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x)
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Initialize input for backward transform"
  call init(x, N, -H, bscale)

  print *,"Compute backward transform"
  status = DftiComputeBackward(hand, x)
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Release the DFTI descriptor"
  ignored_status = DftiFreeDescriptor(hand)

  if (allocated(x)) then
      print *,"Deallocate data array"
      deallocate(x)
  endif

  if (status == 0) then
    print *,"TEST PASSED"
    call exit(0)
  else
    print *,"TEST FAILED"
    call exit(1)
  endif

999 print '("  Error, status = ",I0)', status
  goto 100

contains

  ! Compute mod(K*L,M) accurately
  pure real(WP) function moda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    moda = real(mod(k8*l,m),WP)
  end function moda

  ! Initialize input so that FFT scaled with'scale' produces unit peak at H
  subroutine init(x, N, H, scale)
    integer N, H
    complex(WP) :: x(:)
    real(WP) :: scale

    integer k
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)
    real(WP) :: factor

    factor = 1.0 / scale
    do k = 1, N
      x(k) = factor * exp( I_TWOPI * moda(k-1, H, N) / N ) / N
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

    maxerr = 0.0_WP
    do k = 1, N
      if (mod(k-1-H,N)==0) then
        res_exp = 1.0_WP
      else
        res_exp = 0.0_WP
      end if
      res_got = x(k)
      err = abs(res_got - res_exp)
      maxerr = max(err,maxerr)
      if (.not.(err < errthr)) then
        print '("  x("I0"): "$)', k
        print '(" expected ("G24.17","G24.17"),"$)', res_exp
        print '(" got ("G24.17","G24.17"),"$)', res_got
        print '(" err "G10.3)', err
        print *," Verification FAILED"
        verificate = 100
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program config_scale
