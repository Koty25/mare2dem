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
!       Example of using dfftw_plan_guru_split_dft function.
!
!*****************************************************************************
program dp_plan_guru_split_dft

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 3D transform, and how many of them to execute per call
  integer, parameter :: N1 = 4
  integer, parameter :: N2 = 8
  integer, parameter :: N3 = 16
  integer, parameter :: M  = 32

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2
  integer, parameter :: H3 = 3

  ! Strides and distance describe data layout
  integer :: stride(3), dist

  ! FFTW plan
  integer*8 :: plan = 0

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Data arrays
  real(WP), allocatable :: x_re(:,:,:,:), x_im(:,:,:,:)

  ! Execution status
  integer :: status = 0

  print *,"Example dp_plan_guru_split_dft"
  print *,"Multiple 3D split-complex transform using FFTW guru interface"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0"]")', N1, N2, N3
  print '("  M = "I0)', M
  print '("  H = ["I0","I0","I0"]")', H1, H2, H3

  print *,"Allocate arrays"
  ! Arbitrary but identical for x_re and x_im padding may be added here
  allocate ( x_re(N1, N2, N3, M), STAT = status)
  if (0 /= status) goto 999
  allocate ( x_im(N1, N2, N3, M), STAT = status)
  if (0 /= status) goto 999

  print *,"Setup strides and distance"
  stride(1) = 1
  stride(2) = size(x_re, dim=1)
  stride(3) = stride(2)*size(x_re, dim=2)
  dist      = stride(3)*size(x_re, dim=3)

  print *,"Create FFTW plan for in-place multiple 3D FFT"
  call dfftw_plan_guru_split_dft(plan, &
    3, (/N1, N2, N3/), stride, stride, &
    1, M, dist, dist,                  &
    x_re, x_im, x_re, x_im, FFTW_ESTIMATE)
  if (0 == plan) goto 999

  print *,"Initialize input data for forward transform"
  call init(x_re, x_im, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute forward transform"
  call dfftw_execute(plan)

  print *,"Verify the result of the forward transform"
  status = verificate(x_re, x_im, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

  print *,"Initialize input for backward transform"
  call init(x_re, x_im, N1, N2, N3, M, -H1, -H2, -H3)

  print *,"Compute backward transform (x_re <-> x_im)"
  call dfftw_execute_split_dft(plan, x_im, x_re, x_im, x_re)

  print *,"Verify the result of backward transform"
  status = verificate(x_re, x_im, N1, N2, N3, M, H1, H2, H3)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plan"
  call dfftw_destroy_plan(plan)

  print *,"Deallocate arrays"
  deallocate(x_re, x_im)

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
  subroutine init(x_re, x_im, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3, M
    real(WP) :: x_re(:,:,:,:), x_im(:,:,:,:)

    integer k1,k2,k3,mm
    real(WP), parameter:: TWOPI = 2*3.14159265358979323846_WP

    forall (k1=1:N1, k2=1:N2, k3=1:N3, mm = 1:M)
      x_re(k1,k2,k3,mm) = cos( TWOPI * (   &
        real(  moda(H1,k1-1,N1),WP) / N1   &
        + real(moda(H2,k2-1,N2),WP) / N2   &
        + real(moda(H3,k3-1,N3),WP) / N3 )) / (N1*N2*N3)
      x_im(k1,k2,k3,mm) = sin( TWOPI * (   &
        real(  moda(H1,k1-1,N1),WP) / N1   &
        + real(moda(H2,k2-1,N2),WP) / N2   &
        + real(moda(H3,k3-1,N3),WP) / N3 )) / (N1*N2*N3)
    end forall
  end subroutine init

  ! Verify that x is unit peak at (H1,H2,H3,mm)
  integer function verificate(x_re, x_im, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3, M
    real(WP) :: x_re(:,:,:,:), x_im(:,:,:,:)

    integer k1,k2,k3,mm
    real(WP) err, errthr, maxerr, re_exp, im_exp, re_got, im_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0
    do mm = 1, M
      do k3 = 1, N3
        do k2 = 1, N2
          do k1 = 1, N1
            if (mod(k1-1-H1,N1)==0 .AND. &
              mod  (k2-1-H2,N2)==0 .AND. &
              mod  (k3-1-H3,N3)==0) then
              re_exp = 1
            else
              re_exp = 0
            end if
            im_exp = 0
            re_got = x_re(k1,k2,k3,mm)
            im_got = x_im(k1,k2,k3,mm)
            err = (abs(re_got-re_exp) + abs(im_got - im_exp))
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0", "I0"):"$)', k1,k2,k3,mm
              print '(" expected ("G24.17","G24.17"),"$)', re_exp, im_exp
              print '(" got ("G24.17","G24.17"),"$)', re_got, im_got
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

end program dp_plan_guru_split_dft
