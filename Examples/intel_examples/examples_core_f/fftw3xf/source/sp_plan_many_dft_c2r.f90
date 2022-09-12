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
!       Example of using sfftw_plan_many_dft_c2r function.
!
!*****************************************************************************
program sp_plan_many_dft_c2r

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 3D transform and number of transforms
  integer, parameter :: N1 = 7
  integer, parameter :: N2 = 13
  integer, parameter :: N3 = 5
  integer, parameter :: M = 4 ! howmany

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2
  integer, parameter :: H3 = 3

  ! These parameters describe data layout for real and complex domains
  integer :: rembed(3), rstride, rdist
  integer :: cembed(3), cstride, cdist

  ! FFTW plans
  integer*8 :: plan_r2c = 0, plan_c2r = 0

  ! Need single precision
  integer, parameter :: WP = selected_real_kind(6,37)

  ! Data arrays
  complex(WP), allocatable :: x_cmplx(:,:,:,:)
  real(WP), allocatable :: x_real(:,:,:,:)

  ! Execution status
  integer :: status = 0

  print *,"Example sp_plan_many_dft_c2r"
  print *,"Multiple 3D real-to-complex and complex-to-real FFT"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0"]")', N1, N2, N3
  print '("  M = "I0)', M
  print '("  H = ["I0","I0","I0"]")', H1, H2, H3

  print *,"Allocate arrays"
  ! Arbitrary padding may be added here
  allocate ( x_real (  N1     +1, N2+2, N3+3, M), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_cmplx( (INT(N1/2.0)+1)+4, N2+5, N3+6, M), STAT = status)
  if (0 /= status) goto 999

  print *,"Setup embedding array sizes, strides, and distances"
  rembed(1) = size(x_real, dim=1)
  rembed(2) = size(x_real, dim=2)
  rembed(3) = size(x_real, dim=3)
  rstride   = 1
  rdist     = product(rembed)

  cembed(1) = size(x_cmplx, dim=1)
  cembed(2) = size(x_cmplx, dim=2)
  cembed(3) = size(x_cmplx, dim=3)
  cstride   = 1
  cdist     = product(cembed)

  print '("  rembed=["I0","I0","I0"], rstride="I0", rdist="I0)', rembed, rstride, rdist
  print '("  cembed=["I0","I0","I0"], cstride="I0", cdist="I0)', cembed, cstride, cdist

  print *,"Create FFTW plan for r2c transform"
  call sfftw_plan_many_dft_r2c(plan_r2c, 3, (/N1, N2, N3/), M, &
    x_real,  rembed, rstride, rdist,                           &
    x_cmplx, cembed, cstride, cdist, FFTW_ESTIMATE)
  if (0 == plan_r2c) goto 999

  print *,"Create FFTW plan for c2r transform"
  call sfftw_plan_many_dft_c2r(plan_c2r, 3, (/N1, N2, N3/), M, &
    x_cmplx, cembed, cstride, cdist,                           &
    x_real,  rembed, rstride, rdist, FFTW_ESTIMATE)
  if (0 == plan_c2r) goto 999

  print *,"Initialize input data for r2c transform"
  call init_r(x_real, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute real-to-complex transform"
  call sfftw_execute(plan_r2c)

  print *,"Verify the result of r2c transform"
  status = verificate_complex(x_cmplx, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

  print *,"Initialize input data for c2r transform"
  call init_c(x_cmplx, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute complex-to-real transform"
  call sfftw_execute(plan_c2r)

  print *,"Verify the result of c2r transform"
  status = verificate_real(x_real, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

100 continue

  ! Destroy FFTW plans
  call sfftw_destroy_plan(plan_r2c)
  call sfftw_destroy_plan(plan_c2r)

  ! Deallocate arrays
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

  ! Initialize x to be inverse of unit peak at H1,H2,H3
  subroutine init_r(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
    real(WP) :: x(:, :, :, :)

    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP
    integer k1,k2,k3, mm

    forall (k1=1:N1, k2=1:N2, k3=1:N3, mm=1:M)
      x(k1,k2,k3,mm) = 2.0_WP * cos( TWOPI * ( &
    &    real  (moda(H1,k1-1,N1),WP) / real(N1,WP) &
    &    + real(moda(H2,k2-1,N2),WP) / real(N2,WP) &
    &    + real(moda(H3,k3-1,N3),WP) / real(N3,WP))) / real((N1*N2*N3),WP)
    end forall
    if (mod(H1,N1)==0 .and. mod(H2,N2)==0 .and. mod(H3,N3)==0) then
      x(:,:,:,:) =  x(:,:,:,:) / 2
    end if
  end subroutine init_r

  ! Initialize x to be inverse of unit peak at H1,H2,H3
  subroutine init_c(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
    complex(WP) :: x(:, :, :, :)

    complex(WP), parameter :: ITWOPI = (0.0_WP,6.2831853071795864769_WP)
    integer k1,k2,k3, mm

    forall (k1=1:N1/2+1, k2=1:N2, k3=1:N3, mm=1:M)
      x(k1,k2,k3,mm) = exp( -ITWOPI * ( &
    &    real  (moda(H1,k1-1,N1),WP) / real(N1,WP) &
    &    + real(moda(H2,k2-1,N2),WP) / real(N2,WP) &
    &    + real(moda(H3,k3-1,N3),WP) / real(N3,WP))) / cmplx(N1*N2*N3)
    end forall
  end subroutine init_c

  ! Verify that x is unit peak at (H1,H2,H3,mm)
  integer function verificate_complex(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
    complex(WP) :: x(:, :, :, :)

    integer k1,k2,k3, mm
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do mm = 1, M
      do k3 = 1, N3
        do k2 = 1, N2
          do k1 = 1, N1/2+1
            if (mod(k1-1-H1,N1)==0 .and. &
    &          mod  (k2-1-H2,N2)==0 .and. &
    &          mod  (k3-1-H3,N3)==0) then
              res_exp = 1.0_WP
            else if (mod(-k1+1-H1,N1)==0 .and. &
    &          mod       (-k2+1-H2,N2)==0 .and. &
    &          mod       (-k3+1-H3,N3)==0) then
              res_exp = 1.0_WP
            else
              res_exp = 0.0_WP
            end if
            res_got = x(k1,k2,k3,mm)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0", "I0"):"$)', k1,k2,k3,mm
              print '(" expected ("G14.7","G14.7"),"$)', res_exp
              print '(" got ("G14.7","G14.7"),"$)', res_got
              print '(" err "G10.3)', err
              print *," Verification FAILED"
              verificate_complex = 1
              return
            end if
          end do
        end do
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate_complex = 0
  end function verificate_complex

  ! Verify that x is unit peak at (H1,H2,H3,mm)
  integer function verificate_real(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
    real(WP) :: x(:, :, :, :)

    integer k1,k2,k3, mm
    real(WP) err, errthr, maxerr
    real(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do mm = 1, M
      do k3 = 1, N3
        do k2 = 1, N2
          do k1 = 1, N1/2+1
            if (mod(k1-1-H1,N1)==0 .and. &
    &          mod  (k2-1-H2,N2)==0 .and. &
    &          mod  (k3-1-H3,N3)==0) then
              res_exp = 1.0_WP
            else
              res_exp = 0.0_WP
            end if
            res_got = x(k1,k2,k3,mm)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0", "I0"):"$)', k1,k2,k3,mm
              print '(" expected "G14.7","$)', res_exp
              print '(" got "G14.7","$)', res_got
              print '(" err "G10.3)', err
              print *," Verification FAILED"
              verificate_real = 1
              return
            end if
          end do
        end do
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate_real = 0
  end function verificate_real

end program sp_plan_many_dft_c2r
