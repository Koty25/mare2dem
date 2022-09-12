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
! A simple example of double-precision complex-to-complex in-place 2D
! FFT using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
!
!*****************************************************************************

program basic_dp_complex_dft_2d

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Sizes of 2D transform
  integer, parameter :: N1 = 7
  integer, parameter :: N2 = 13

  ! Arbitrary harmonic to test the FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2

  ! need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  complex(WP), allocatable :: x (:,:)

  type(DFTI_DESCRIPTOR), POINTER :: hand

  hand => null()

  print *,"Example basic_dp_complex_dft_2d"
  print *,"Forward and backward double-precision complex-to-complex ",       &
    &    " in-place 2D transform"
  print *,"Configuration parameters:"
  print *,"DFTI_PRECISION      = DFTI_DOUBLE"
  print *,"DFTI_FORWARD_DOMAIN = DFTI_COMPLEX"
  print *,"DFTI_DIMENSION      = 2"
  print '(" DFTI_LENGTHS        = /"I0","I0"/" )', N1, N2

  print *,"Create DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,[N1,N2])
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Allocate array for input data"
  allocate ( x(N1, N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Initialize input for forward transform"
  call init(x, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verificate(x, N1, N2, H1, H2)
  if (0 /= status) goto 999

  print *,"Initialize input for backward transform"
  call init(x, N1, N2, -H1, -H2)

  print *,"Compute backward transform"
  status = DftiComputeBackward(hand, x(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verificate(x, N1, N2, H1, H2)
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

  ! Initialize array with harmonic /H1, H2/
  subroutine init(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = exp( I_TWOPI * ( &
    &    moda(  k1-1,H1, N1) / N1  &
    &    + moda(k2-1,H2, N2) / N2 )) / (N1*N2)
    end forall
  end subroutine init

  ! Verify that x(N1,N2) is unit peak at x(H1,H2)
  integer function verificate(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    real(WP) err, errthr, maxerr
    complex(WP) :: res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k2 = 1, N2
      do k1 = 1, N1
        if (mod(k1-1-H1, N1)==0 .AND. mod(k2-1-H2, N2)==0) then
          res_exp = 1.0_WP
        else
          res_exp = 0.0_WP
        end if
        res_got = x(k1,k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1, k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', res_got
          print '(" err "G10.3)', err
          print *," Verification FAILED"
          verificate = 100
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program basic_dp_complex_dft_2d
