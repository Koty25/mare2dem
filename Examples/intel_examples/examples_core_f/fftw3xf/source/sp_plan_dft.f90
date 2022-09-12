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
!       Example of using sfftw_plan_dft function.
!
!*****************************************************************************

program sp_plan_dft

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 4D transform
  integer, parameter :: N(4) = [7,13,5,4]

  ! Arbitrary harmonic used to verify FFT
  integer :: H(4) = [1,2,3,4]

  ! Working precision is single precision (sfftw_* functions)
  integer, parameter :: WP = selected_real_kind(6,37)

  ! FFTW plan handles
  integer*8 :: fwd = 0, bwd = 0

  ! Data array
  complex(WP), allocatable :: x(:,:,:,:)

  ! Execution status
  integer :: status = 0

  print *,"Example sp_plan_dft"
  print *,"Compute forward and backward 4D in-place FFT"
  print *,"Configuration parameters:"
  print '("  N = ["3(I0",")I0"]")', N
  print '("  H = ["3(I0",")I0"]")', H

  print *,"Allocate data array"
  allocate ( x(N(1),N(2),N(3),N(4)), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW plan for forward transform"
  call sfftw_plan_dft(fwd,4,N(1:4),x,x,FFTW_FORWARD,FFTW_ESTIMATE)
  if (0 == fwd) goto 999

  print *,"Create FFTW plan for backward transform"
  call sfftw_plan_dft(bwd,4,N(1:4),x,x,FFTW_BACKWARD,FFTW_ESTIMATE)
  if (0 == bwd) goto 999

  print *,"Initialize data for forward FFT"
  call init(x, N(1:4), H(1:4))

  print *,"Compute forward transform"
  call sfftw_execute(fwd)

  print *,"Verify the result of forward transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Initialize data for backward FFT"
  H(1:4) = -H(1:4)
  call init(x, N(1:4), H(1:4))

  print *,"Compute backward transform"
  call sfftw_execute_dft(bwd, x, x)

  print *,"Verify the result of backward transform"
  H(1:4) = -H(1:4)
  status = verificate(x, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Deallocate data array"
  deallocate(x)

  print *,"Destroy FFTW plans"
  call sfftw_destroy_plan(fwd)
  call sfftw_destroy_plan(bwd)

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

  ! Initialize array x with harmonic H
  subroutine init(x, N, H)
    integer :: N(4), H(4)
    complex(WP) :: x(:,:,:,:)

    integer :: k1,k2,k3,k4
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

    forall (k1=1:N(1), k2=1:N(2), k3=1:N(3), k4=1:N(4))
      x(k1,k2,k3,k4) = exp( I_TWOPI * (       &
    &    real(  moda(k1-1,H(1),N(1)),WP)/real(N(1),WP)  &
    &    + real(moda(k2-1,H(2),N(2)),WP)/real(N(2),WP)  &
    &    + real(moda(k3-1,H(3),N(3)),WP)/real(N(3),WP)  &
    &    + real(moda(k4-1,H(4),N(4)),WP)/real(N(4),WP)))&
    &    / cmplx(product(N))
    end forall
  end subroutine init

  ! Verify that x(N) is unit peak at x(H)
  integer function verificate(x, N, H)
    integer :: N(4), H(4)
    complex(WP) :: x(:,:,:,:)

    integer k1,k2,k3,k4
    real(WP) err, errthr, maxerr
    complex(WP) :: res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * sum(log(real(N,WP))) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k4 = 1, N(4)
      do k3 = 1, N(3)
        do k2 = 1, N(2)
          do k1 = 1, N(1)
            if (mod(k1-1-H(1),N(1))==0 .AND.                            &
    &          mod  (k2-1-H(2),N(2))==0 .AND.                           &
    &          mod  (k3-1-H(3),N(3))==0 .AND.                           &
    &          mod  (k4-1-H(4),N(4))==0) then
              res_exp = 1.0_WP
            else
              res_exp = 0.0_WP
            end if
            res_got = x(k1,k2,k3,k4)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("3(I0",")I0"): "$)', k1,k2,k3,k4
              print '(" expected "G14.7","$)', res_exp
              print '(" got "G14.7","$)', res_got
              print '(" err "G10.3)', err
              print *," Verification FAILED"
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

end program sp_plan_dft
