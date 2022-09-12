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
! A simple example of double-precision real-to-complex out-of-place 2D
! FFT using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
!
!*****************************************************************************
program basic_dp_real_dft_2d

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Sizes of 2D transform
  integer, parameter :: N1 = 7
  integer, parameter :: N2 = 13

  ! Arbitrary harmonic used to test FFT
  integer, parameter :: H1 = 3
  integer, parameter :: H2 = 1

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! Strides define data layout for real and complex domains
  integer cstrides(3), rstrides(3)

  ! Data arrays
  real(WP), allocatable :: x_real (:,:)
  complex(WP), allocatable :: x_cmplx (:,:)

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  hand => null()

  print *,"Example basic_dp_real_dft_2d"
  print *,"Forward-Backward double-precision real out-of-place 2D transform"
  print *,"Configuration parameters:"
  print *,"DFTI_PRECISION              = DFTI_DOUBLE"
  print *,"DFTI_FORWARD_DOMAIN         = DFTI_REAL"
  print *,"DFTI_DIMENSION              = 2"
  print '(" DFTI_LENGTHS                = /"I0","I0"/" )', N1, N2
  print *,"DFTI_PLACEMENT              =  DFTI_NOT_INPLACE"
  print *,"DFTI_CONJUGATE_EVEN_STORAGE = DFTI_COMPLEX_COMPLEX"

  print *,"Create DFTI descriptor for real transform"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 2, (/N1, N2/))
  if (0 /= status) goto 999

  print *,"Set out-of-place"
  status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  if (0 /= status) goto 999

  print *,"Set CCE storage"
  status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE,                   &
    &                    DFTI_COMPLEX_COMPLEX)
  if (0 /= status) goto 999

  rstrides = [0, 1, N1]
  cstrides = [0, 1, INT(N1/2.0) + 1]

  print '(" Set input  strides = "3(I0:", "))', rstrides
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides)
  if (0 /= status) goto 999

  print '(" Set output strides = "3(I0:", "))', cstrides
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Allocate data arrays"
  allocate ( x_real(N1, N2), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx(INT(N1/2.0) + 1, N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Initialize data for real-to-complex FFT"
  call init_r(x_real, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x_real(:,1), x_cmplx(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verify_c(x_cmplx, N1, N2, H1, H2)
  if (0 /= status) goto 999


  print *,"Reconfigure DFTI descriptor for backward transform"

  print '(" Set input  strides = "3(I0:", "))', cstrides
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, cstrides)
  if (0 /= status) goto 999

  print '(" Set output strides = "3(I0:", "))', rstrides
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, rstrides)
  if (0 /= status) goto 999

  print *,"Recommit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize data for complex-to-real FFT"
  call init_c(x_cmplx, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeBackward(hand, x_cmplx(:,1), x_real(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verify_r(x_real, N1, N2, H1, H2)
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

  ! Initialize x(:,:) to harmonic H
  subroutine init_r(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)

    integer k1, k2
    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP

    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = 2 * cos( TWOPI * ( moda(H1,k1-1,N1) / N1 &
    &    +                           moda(H2,k2-1,N2) / N2)) / (N1*N2)
    end forall
    if (mod(2*(N1-H1),N1)==0 .and. mod(2*(N2-H2),N2)==0) then
      x(1:N1,1:N2) = x(1:N1,1:N2) / 2
    end if
  end subroutine init_r

  ! Initialize x(:,:) to produce unit peak at x(H1,H2)
  subroutine init_c(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)
    integer k1,k2

    forall (k1=1:N1/2+1, k2=1:N2)
      x(k1,k2) = exp( -I_TWOPI * ( moda(H1,k1-1,N1) / N1 &
    &    +                          moda(H2,k2-1,N2) / N2)) / (N1*N2)
    end forall
  end subroutine init_c


  ! Verify that x(:,:) has unit peak at (H1,H2)
  integer function verify_c(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k2 = 1, N2
      do k1 = 1, N1/2+1
        if (mod(k1-1-H1,N1)==0 .and. mod(k2-1-H2,N2)==0) then
          res_exp = 1.0_WP
        else if (mod(-k1+1-H1,N1)==0 .and. mod(-k2+1-H2,N2)==0) then
          res_exp = 1.0_WP
        else
          res_exp = 0.0_WP
        end if
        res_got = x(k1, k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1,k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', res_got
          print '(" err "G10.3)', err
          print *,"  Verification FAILED"
          verify_c = 100
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verify_c = 0
  end function verify_c

  ! Verify that x(:,:) is unit peak at x(H1,H2)
  integer function verify_r(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)

    integer k1, k2
    real(WP) err, errthr, maxerr, res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k2 = 1, N2
      do k1 = 1, N1
        if (mod(k1-1-H1, N1)==0 .AND. mod(k2-1-H2,N2)==0) then
          res_exp = 1.0_WP
        else
          res_exp = 0.0_WP
        end if
        res_got = x(k1, k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1, k2
          print '(" expected "G24.17","$)', res_exp
          print '(" got "G24.17","$)', res_got
          print '(" err "G10.3)', err
          print *,"  Verification FAILED"
          verify_r = 100
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verify_r = 0
  end function verify_r

end program basic_dp_real_dft_2d
