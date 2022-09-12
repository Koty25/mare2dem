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
!       Example of using dfftw_plan_dft_2d function.
!
!*****************************************************************************

program dp_plan_dft_2d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 2D transform
  integer, parameter :: N1 = 8
  integer, parameter :: N2 = 16

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = -N2/2

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0

  ! Data array
  complex(WP), allocatable :: x(:,:)

  ! FFTW plan
  integer*8 :: fwd = 0, bwd = 0

  print *,"Example dp_plan_dft_2d"
  print *,"Forward and backward 2D complex in-place FFT"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0"]")', N1, N2
  print '("  H = ["I0","I0"]")', H1, H2

  print *,"Allocate array"
  allocate ( x(N1,N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW forward transform plan"
  call dfftw_plan_dft_2d(fwd, N1, N2, x, x, FFTW_FORWARD, FFTW_ESTIMATE)
  if (0 == fwd) goto 999

  print *,"Create FFTW backward transform plan"
  call dfftw_plan_dft_2d(bwd, N1, N2, x, x, FFTW_BACKWARD, FFTW_ESTIMATE)
  if (0 == bwd) goto 999

  print *,"Initialize data for forward transform"
  call init(x, N1, N2, H1, H2)

  print *,"Compute forward transform"
  call dfftw_execute(fwd)

  print *,"Verify the result of the forward transform"
  status = verificate(x, N1, N2, H1, H2)
  if (0 /= status) goto 999

  print *,"Initialize data for backward transform"
  call init(x, N1, N2, -H1, -H2)

  print *,"Compute backward transform"
  call dfftw_execute_dft(bwd, x, x)

  print *,"Verify the result of the backward transform"
  status = verificate(x, N1, N2, H1, H2)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plans"
  call dfftw_destroy_plan(fwd)
  call dfftw_destroy_plan(bwd)

  print *,"Deallocate array"
  deallocate(x)

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

  ! Initialize array with harmonic /H1, H2/
  subroutine init(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    complex(WP), parameter :: I_TWOPI = (0,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = exp( I_TWOPI * ( &
        real(  moda(k1-1,H1, N1), WP)/N1 &
        + real(moda(k2-1,H2, N2), WP)/N2 )) / (N1*N2)
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

    maxerr = 0
    do k2 = 1, N2
      do k1 = 1, N1
        if (mod(k1-1-H1, N1)==0 .AND. mod(k2-1-H2, N2)==0) then
          res_exp = 1
        else
          res_exp = 0
        end if
        res_got = x(k1,k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1, k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', x(k1,k2)
          print '(" err "G10.3)', err
          print *," Verification FAILED"
          verificate = 1
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program dp_plan_dft_2d
