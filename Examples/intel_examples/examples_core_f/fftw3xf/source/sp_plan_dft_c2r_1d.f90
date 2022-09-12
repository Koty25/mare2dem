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
!       Example of using sfftw_plan_dft_c2r_1d function.
!
!*****************************************************************************

program sp_plan_dft_c2r_1d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Size of 1D transform
  integer, parameter :: N = 7

  ! Arbitrary harmonic used to verify FFT
  integer :: H = 1

  ! Need single precision
  integer, parameter :: WP = selected_real_kind(6,37)

  ! Execution status
  integer :: status = 0

  ! Data arrays
  complex(WP), allocatable :: x_cmplx(:)
  real(WP), allocatable :: x_real(:)

  ! FFTW plan handle
  integer*8 :: plan_c2r = 0

  print *,"Example sp_plan_dft_c2r_1d"
  print *,"1D complex-to-real out-of-place FFT"
  print *,"Configuration parameters:"
  print '("  N = "I0)', N
  print '("  H = "I0)', H

  print *,"Allocate data arrays"
  allocate ( x_real (N)    , STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx(INT(N/2.0)+1), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW complex-to-real plan for 1D FFT"
  call sfftw_plan_dft_c2r_1d(plan_c2r, N, x_cmplx, x_real, FFTW_ESTIMATE)
  if (0 == plan_c2r) goto 999

  print *,"Initialize data"
  call init(x_cmplx, N, H)

  print *,"Compute complex-to-real transform"
  call sfftw_execute(plan_c2r)

  print *,"Verify complex-to-real transform"
  status = verificate(x_real, N, H)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(plan_c2r)

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

  ! Initialize x(N) to produce unit peak at H
  subroutine init(x, N, H)
    integer :: N, H
    complex(WP) :: x(:)

    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)
    integer k

    forall (k=1:N/2+1)
      x(k) = exp( -I_TWOPI * (real(moda(H,k-1,N),WP)/real(N,WP)) ) / cmplx(N)
    end forall
  end subroutine init

  ! Verify that x is unit peak at H
  integer function verificate(x, N, H)
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
        print '(" expected "G14.7","$)', res_exp
        print '(" got "G14.7","$)', x(k)
        print '(" err "G10.3)', err
        print *," Verification FAILED"
        verificate = 1
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program sp_plan_dft_c2r_1d
