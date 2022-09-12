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
! A simple example of double-precision real-to-complex out-of-place 1D
! FFT using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
!
!*****************************************************************************
program basic_dp_real_dft_1d

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Size of 1D transform
  integer, parameter :: N = 7

  ! Arbitrary harmonic used to test FFT
  integer, parameter :: H = 1

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! Data arrays
  real(WP), allocatable :: x_real (:)
  complex(WP), allocatable :: x_cmplx (:)

  type(DFTI_DESCRIPTOR), POINTER :: hand

  hand => null()

  print *,"Example basic_dp_real_dft_1d"
  print *,"Forward-Backward double-precision real-to-complex",               &
    &      " out-of-place 1D transform"
  print *,"Configuration parameters:"
  print *,"DFTI_PRECISION              = DFTI_DOUBLE"
  print *,"DFTI_FORWARD_DOMAIN         = DFTI_REAL"
  print *,"DFTI_DIMENSION              = 1"
  print '(" DFTI_LENGTHS                = /"I0"/" )', N
  print *,"DFTI_PLACEMENT              =  DFTI_NOT_INPLACE"
  print *,"DFTI_CONJUGATE_EVEN_STORAGE = DFTI_COMPLEX_COMPLEX"

  print *,"Create DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 1, N)
  if (0 /= status) goto 999

  print *,"Set out-of-place"
  status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  if (0 /= status) goto 999

  print *,"Set CCE storage"
  status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,                   &
    &                    DFTI_COMPLEX_COMPLEX)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Allocate data arrays"
  allocate ( x_real(N), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx(INT(N/2.0) + 1), STAT = status)
  if (0 /= status) goto 999

  print *,"Initialize data for real-to-complex FFT"
  call init_r(x_real, N, H)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x_real, x_cmplx)
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verify_c(x_cmplx, N, H)
  if (0 /= status) goto 999

  print *,"Initialize data for complex-to-real FFT"
  call init_c(x_cmplx, N, H)

  print *,"Compute backward transform"
  status = DftiComputeBackward(hand, x_cmplx, x_real)
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verify_r(x_real, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Release the DFTI descriptor"
  ignored_status = DftiFreeDescriptor(hand)

  if (allocated(x_real) .or. allocated(x_cmplx)) then
    print *,"Deallocate data arrays"
  endif
  if (allocated(x_real))     deallocate(x_real)
  if (allocated(x_cmplx))    deallocate(x_cmplx)

  if (status == 0) then
    print *, "TEST PASSED"
    call exit(0)
  else
    print *, "TEST FAILED"
    call exit(1)
  end if

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

  ! Initialize x to be inverse of unit peak at H
  subroutine init_r(x, N, H)
    integer N, H
    real(WP) :: x(:)

    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP
    integer k

    forall (k=1:N)
      x(k) = 2 * cos( TWOPI * moda(H,k-1,N) / N) / N
    end forall
    if (mod(2*(N-H),N)==0) x(1:N) =  x(1:N) / 2
  end subroutine init_r

  ! Initialize x to be inverse of unit peak at H
  subroutine init_c(x, N, H)
    integer N, H
    complex(WP) :: x(:)

    complex(WP), parameter :: ITWOPI = (0.0_WP,6.2831853071795864769_WP)
    integer k

    forall (k=1:N/2+1)
      x(k) = exp( -ITWOPI * moda(H,k-1,N) / N) / N
    end forall
  end subroutine init_c

  ! Verify that x(N) is unit peak at H and N-H
  integer function verify_c(x, N, H)
    integer N, H
    complex(WP) :: x(:)

    integer k
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k = 1, N/2+1
      if (mod(k-1-H,N)==0 .or. mod(-k+1-H,N)==0) then
        res_exp = 1.0_WP
      else
        res_exp = 0.0_WP
      end if
      res_got = x(k)
      err = abs(res_got - res_exp)
      maxerr = max(err,maxerr)
      if (.not.(err < errthr)) then
        print '("  x("I0"):"$)', k
        print '(" expected ("G24.17","G24.17"),"$)', res_exp
        print '(" got ("G24.17","G24.17"),"$)', res_got
        print '(" err "G10.3)', err
        print *," Verification FAILED"
        verify_c = 1
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verify_c = 0
  end function verify_c

  ! Verify that x is unit peak at H
  integer function verify_r(x, N, H)
    integer N, H
    real(WP) :: x(:)

    integer k
    real(WP) err, errthr, maxerr, res_exp

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k = 1, N
      if (mod(k-1-H, N)==0) then
        res_exp = 1.0_WP
      else
        res_exp = 0.0_WP
      end if
      err = abs(x(k) - res_exp)
      maxerr = max(err,maxerr)
      if (.not.(err < errthr)) then
        print '("  x("I0"):"$)', k
        print '(" expected "G24.17","$)', res_exp
        print '(" got "G24.17","$)', x(k)
        print '(" err "G10.3)', err
        print *," Verification FAILED"
        verify_r = 100
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verify_r = 0
  end function verify_r

end program basic_dp_real_dft_1d
