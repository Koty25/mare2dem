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
!       Example of using sfftw_plan_many_dft function.
!
!*****************************************************************************
program sp_plan_many_dft

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 3D transform and the number of transforms
  integer, parameter :: N1 = 10
  integer, parameter :: N2 = 12
  integer, parameter :: N3 = 18
  integer, parameter :: M  = 4

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2
  integer, parameter :: H3 = 3

  ! Leading stride, sizes of embedding array, and distance
  ! for multiple in-place transform
  integer :: stride, nembed(3), dist

  ! FFTW plans
  integer*8 :: fwd = 0, bwd = 0

  ! Need single precision
  integer, parameter :: WP = selected_real_kind(6,37)
  complex(WP), allocatable :: x(:,:,:,:)

  ! Execution status
  integer :: status = 0

  print *,"Example sp_plan_many_dft"
  print *,"Forward and backward multiple 3D FFT"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0"]")', N1, N2, N3
  print '("  M = "I0)', M
  print '("  H = ["I0","I0","I0"]")', H1, H2, H3

  print *, "Allocate data array"
  ! Arbitrary padding may be added here
  allocate ( x(N1+1, N2+2, N3+3, M), STAT = status)
  if (0 /= status) goto 999

  print *,"Setup layout"
  stride = 1
  nembed(1) = size(x,dim=1)
  nembed(2) = size(x,dim=2)
  nembed(3) = size(x,dim=3)
  dist = product(nembed)

  print *,"Create plan for forward FFT"
  call sfftw_plan_many_dft(fwd, 3, (/N1, N2, N3/), M, &
    &                       x, nembed, stride, dist,   &
    &                       x, nembed, stride, dist,   &
    &                       FFTW_FORWARD, FFTW_ESTIMATE)
  if (0 == fwd) goto 999

  call sfftw_plan_many_dft(bwd, 3, (/N1, N2, N3/), M, &
    &                       x, nembed, stride, dist,   &
    &                       x, nembed, stride, dist,   &
    &                       FFTW_BACKWARD, FFTW_ESTIMATE)
  if (0 == bwd) goto 999

  print *,"Initialize data for forward transform"
  call init(x, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute Forward transform"
  call sfftw_execute(fwd)

  print *,"Verify forward transform"
  status = verificate(x, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

  print *,"Initialize data for backward transform"
  call init(x, N1, N2, N3, M, -H1, -H2, -H3)

  print *,"Compute backward transform"
  call sfftw_execute_dft(bwd, x, x)

  print *,"Verify backward transform"
  status = verificate(x, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

100 continue

  ! Destroy FFTW plans
  call sfftw_destroy_plan(fwd)
  call sfftw_destroy_plan(bwd)

  ! Deallocate array
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

  ! Initialize array with harmonic /H1, H2, H3/
  subroutine init(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3, M
    complex(WP) :: x(:,:,:,:)

    integer k1, k2, k3, mm
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2, k3=1:N3, mm=1:M)
      x(k1,k2,k3,mm) = exp( I_TWOPI * ( &
    &    real(  moda(k1-1,H1, N1), WP)/real(N1,WP) &
    &    + real(moda(k2-1,H2, N2), WP)/real(N2,WP) &
    &    + real(moda(k3-1,H3, N3), WP)/real(N3,WP) )) / cmplx(N1*N2*N3)
    end forall
  end subroutine init

  ! Verify that x is unit peak at x(H1,H2,H3,:)
  integer function verificate(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3, M
    complex(WP) :: x(:,:,:,:)

    integer k1, k2, k3, mm
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do mm = 1, M
      do k3 = 1, N3
        do k2 = 1, N2
          do k1 = 1, N1
            if (mod(k1-1-H1, N1)==0 .AND.                                &
    &          mod  (k2-1-H2, N2)==0 .AND.                                &
    &          mod  (k3-1-H3, N3)==0) then
              res_exp = 1.0_WP
            else
              res_exp = 0.0_WP
            end if
            res_got = x(k1,k2,k3,mm)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0"):"$)', k1, k2, k3
              print '(" expected ("G14.7","G14.7"),"$)', res_exp
              print '(" got ("G14.7","G14.7"),"$)', res_got
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

end program sp_plan_many_dft
