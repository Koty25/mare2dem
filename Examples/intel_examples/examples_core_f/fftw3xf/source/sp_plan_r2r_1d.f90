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
!       Example of using sfftw_plan_r2r_1d function.
!
!*****************************************************************************

program sp_plan_r2r_1d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Physical size of 1D r2r transform, N>1
  integer, parameter :: N = 7

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H = 1

  ! Working precision is single precision
  integer, parameter :: WP = selected_real_kind(6,37)

  ! Execution status
  integer :: status = 0

  ! Data array
  real(WP), allocatable :: x(:)

  ! FFTW plan
  integer*8 :: r2r = 0

  integer i

  print *,"Example sp_plan_r2r_1d"
  print *,"1D real-to-real transforms"
  print '("  N = "I0)', N
  print '("  H = "I0)', H

  print *,"Allocate data array"
  ! Extra element is needed for MKL_RODFT00 kind
  allocate ( x(N+1), STAT = status)
  if (0 /= status) goto 999

  print *,"========== REDFT00 (DCT1) ========="

  print '(" Create FFTW_REDFT00 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_REDFT00, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for REDFT00, H="I0)', H
  call init_redft00(x, N, H)

  print *,"Compute FFTW_REDFT00 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_REDFT00 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== REDFT10 (DCT2) ========="

  print '(" Create FFTW_REDFT10 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_REDFT10, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for REDFT10, H="I0)', H
  call init_redft10(x, N, H)

  print *,"Compute FFTW_REDFT10 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_REDFT10 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== REDFT01 (DCT3) ========="

  print '(" Create FFTW_REDFT01 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_REDFT01, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for REDFT01, H="I0)', H
  call init_redft01(x, N, H)

  print *,"Compute FFTW_REDFT01 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_REDFT01 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== REDFT11 (DCT4) ========="

  print '(" Create FFTW_REDFT11 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_REDFT11, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for REDFT11, H="I0)', H
  call init_redft11(x, N, H)

  print *,"Compute FFTW_REDFT11 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_REDFT11 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== RODFT00 (DST1) ========="

  print '(" Create FFTW_RODFT00 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_RODFT00, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for RODFT00, H="I0)', H
  call init_rodft00(x, N, H)

  print *,"Compute FFTW_RODFT00 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_RODFT00 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== RODFT10 (DST2) ========="

  print '(" Create FFTW_RODFT10 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_RODFT10, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for RODFT10, H="I0)', H
  call init_rodft10(x, N, H)

  print *,"Compute FFTW_RODFT10 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_RODFT10 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== RODFT01 (DST3) ========="

  print '(" Create FFTW_RODFT01 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_RODFT01, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for RODFT01, H="I0)', H
  call init_rodft01(x, N, H)

  print *,"Compute FFTW_RODFT01 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_RODFT01 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== RODFT11 (DST4) ========="

  print '(" Create FFTW_RODFT11 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, FFTW_RODFT11, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for RODFT11, H="I0)', H
  call init_rodft11(x, N, H)

  print *,"Compute FFTW_RODFT11 transform"
  !print '("X="10(F10.4:","))',x(1:N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:N)

  print *,"Verify the result of FFTW_RODFT11 transform"
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"========== MKL_RODFT00 (DST1) ========="

  print '(" Create MKL_RODFT00 plan for N="I0)', N
  call sfftw_plan_r2r_1d(r2r, N, x, x, MKL_RODFT00, FFTW_ESTIMATE)
  if (0 == r2r) goto 999

  print '(" Initialize data for MKL_RODFT00, H="I0)', H
  call init_rodft00(x, N, H)
  ! Unlike FFTW_RODFT00, MKL_RODFT00 requires extra zero element at x(1)
  ! Some compilers cannot do this well: x(2:1+N) = x(1:N)
  do i = N, 1, -1
    x(i+1) = x(i)
  end do
  x(1) = 0.0_WP

  print *,"Compute MKL_RODFT00 transform"
  !print '("X="10(F10.4:","))',x(1:1+N)
  call sfftw_execute(r2r)
  !print '("Y="10(F10.4:","))',x(1:1+N)

  print *,"Verify the result of MKL_RODFT00 transform"
  ! Ignore the extra zero element at x(1)
  x(1:N) = x(2:1+N)
  status = verificate(x, N, H)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(r2r)

  print *,"Deallocate data array"
  deallocate(x)

100 continue

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
  real(WP) function amoda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    amoda = real( mod(k8*l,m), WP )
  end function amoda

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_redft00(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*(N-1) ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * cos( M_PI * amoda(j,k-1,NL) / real((N-1),WP) )
    end do
    if ( j==0 .or. j==N-1 ) x(:) = x(:) / 2.0
    x(:) = x(:) / real(NL,WP)
  end subroutine init_redft00

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_redft10(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * cos( M_PI * amoda(j,2*k-1,2*NL) / real((2*N),WP) )
    end do
    if (j==0) x(:) = x(:) / 2.0
    x(:) = x(:) / real(NL,WP)
  end subroutine init_redft10

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_redft01(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * cos( M_PI * amoda(2*j+1,k-1,2*NL) / real((2*N),WP) )
    end do
    x(:) = x(:) / real(NL,WP)
  end subroutine init_redft01

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_redft11(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * cos( M_PI * amoda(2*j+1,2*k-1,4*NL) / real((4*N),WP) )
    end do
    x(:) = x(:) / real(NL,WP)
  end subroutine init_redft11

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_rodft00(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*(N+1) ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * sin( M_PI * amoda(j+1,k,NL) / real((N+1),WP) )
    end do
    x(:) = x(:) / real(NL,WP)
  end subroutine init_rodft00

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_rodft10(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * sin( M_PI * amoda(j+1,2*k-1,2*NL) / real((2*N),WP) )
    end do
    if (j==N-1) x(:) = x(:) / 2.0
    x(:) = x(:) / real(NL,WP)
  end subroutine init_rodft10

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_rodft01(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * sin( M_PI * amoda(2*j+1,k,2*NL) / real((2*N),WP) )
    end do
    x(:) = x(:) / real(NL,WP)
  end subroutine init_rodft01

  ! Initialize x(:) to be inverse of (n==H)
  subroutine init_rodft11(x, N, H)
    integer N, H
    real(WP) :: x(:)
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP
    integer j, k, NL

    j = modulo(H,N)
    NL = 2*N ! logical size of transform
    do k = 1, N
      x(k) = 2.0 * sin( M_PI * amoda(2*j+1,2*k-1,4*NL) / real((4*N),WP) )
    end do
    x(:) = x(:) / real(NL,WP)
  end subroutine init_rodft11

  ! Verify that x(:) is unit peak at x(modulo(H,N))
  integer function verificate(x, N, H)
    integer N, H
    real(WP) :: x(:)

    integer j,k
    real(WP) err, errthr, maxerr, res_exp
    real(WP), parameter:: M_PI = 3.14159265358979323846_WP

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N, WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr
    j = modulo(H,N) + 1
    maxerr = 0.0_WP
    do k = 1,N
      if (k==j) then
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
        print *,"  Verification FAILED"
        verificate = 1
        return
      end if
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program sp_plan_r2r_1d
