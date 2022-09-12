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
!       Example of using sfftw_plan_guru_dft function.
!
!*****************************************************************************
program sp_plan_guru_dft

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 3D transform
  integer, parameter :: N1 = 7
  integer, parameter :: N2 = 13
  integer, parameter :: N3 = 5
  integer, parameter :: M = 4 ! howmany

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2
  integer, parameter :: H3 = 3

  ! Input/output strides for 3D FFT
  integer :: istride(3)
  integer :: ostride(3)

  ! Input/output distance for multiple FFT
  integer :: idist
  integer :: odist

  ! FFTW plan
  integer*8 :: fwd = 0

  ! Need single precision
  integer, parameter :: WP = selected_real_kind(6,37)

  ! Data array
  complex(WP), allocatable :: x(:,:,:,:)
  complex(WP), allocatable :: y(:,:,:,:)

  ! Execution status
  integer :: status = 0

  print *,"Example sp_plan_guru_dft"
  print *,"Multiple out-of-place 3D complex FFT"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0"]")', N1, N2, N3
  print '("  M = "I0)', M
  print '("  H = ["I0","I0","I0"]")', H1, H2, H3

  print *,"Allocate data arrays"
  ! Arbitrary padding may be added here
  allocate ( x(N1, N2, N3, M), STAT = status)
  if (0 /= status) goto 999
  allocate ( y(M, N3, N2, N1), STAT = status)
  if (0 /= status) goto 999

  print *,"Set strides and distances"
  istride(1) = 1
  istride(2) = size(x,dim=1)
  istride(3) = istride(2)*size(x,dim=2)
  idist      = istride(3)*size(x,dim=3)

  odist = 1
  ostride(3) = size(y,dim=1)
  ostride(2) = ostride(3)*size(y,dim=2)
  ostride(1) = ostride(2)*size(y,dim=3)
  print '("  istride=["I0","I0","I0"], idist="I0)',istride,idist
  print '("  ostride=["I0","I0","I0"], odist="I0)',ostride,odist

  print *,"Create FFTW plan for the forward transform"
  call sfftw_plan_guru_dft(fwd, 3, (/N1,N2,N3/), istride, ostride,  &
                           1, M, idist, odist , x, y,               &
                           FFTW_FORWARD, FFTW_ESTIMATE)
  if (0 == fwd) goto 999

  print *,"Initialize input for forward transform"
  call init(x, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute forward transform"
  call sfftw_execute_dft(fwd, x, y)

  print *,"Verify the result of forward transform"
  status = verificate(y, M, N3, N2, N1, H3, H2, H1)
  if (0 /= status) goto 999

100 continue

  print *,"Destroy FFTW plan"
  call sfftw_destroy_plan(fwd)

  print *,"Deallocate arrays"
  deallocate(x, y)

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

  ! Verify that x(mm,k3,k2,k1) is unit peak at x(mm,H3,H2,H1)
  integer function verificate(x, M, N3, N2, N1, H3, H2, H1)
    integer N1, N2, N3, H1, H2, H3, M
    complex(WP) :: x(:,:,:,:)

    integer k1, k2, k3, mm
    real(WP) err, errthr, maxerr
    complex(WP) x_exp, x_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k1 = 1, N1
      do k2 = 1, N2
        do k3 = 1, N3
          do mm = 1, M
            if (mod(k1-1-H1,N1)==0 .AND.                                    &
    &            mod(k2-1-H2,N2)==0 .AND.                                    &
    &            mod(k3-1-H3,N3)==0) then
              x_exp = 1.0_WP
            else
              x_exp = 0.0_WP
            end if
            x_got = x(mm, k3, k2, k1)
            err = abs(x_got - x_exp)
            maxerr =  max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0", "I0","I0","I0"):"$)', mm,k3,k2,k1
              print '(" expected ("G14.7","G14.7")"$)', x_exp
              print '(" got ("G14.7","G14.7")"$)', x_got
              print '(" err "G10.3)', err
              print *,"  Verification FAILED"
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

end program sp_plan_guru_dft
