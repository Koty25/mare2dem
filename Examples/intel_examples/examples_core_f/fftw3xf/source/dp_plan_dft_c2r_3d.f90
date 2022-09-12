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
!       Example of using dfftw_plan_dft_c2r_3d function.
!
!*****************************************************************************

program dp_plan_dft_c2r_3d

  include 'fftw3.f'
  include 'fftw3_mkl.f'

  ! Sizes of 3D transform
  integer, parameter :: N1 = 7
  integer, parameter :: N2 = 13
  integer, parameter :: N3 = 5

  ! Arbitrary harmonic used to verify FFT
  integer :: H1 = 1, H2 = 6, H3 = 1

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0

  ! Data array for in-place transform
  complex(WP), allocatable :: x_cmplx(:,:,:)

  ! FFTW plan
  integer*8 :: c2r = 0

  print *,"Example dp_plan_dft_c2r_3d"
  print *,"3D complex-to-real in-place transform"
  print *,"Configuration parameters:"
  print '("  N = ["I0","I0","I0"]")', N1, N2, N3
  print '("  H = ["I0","I0","I0"]")', H1, H2, H3

  print *,"Allocate data array for in-place transform"
  allocate ( x_cmplx(INT(N1/2.0)+1, N2, N3), STAT = status)
  if (0 /= status) goto 999

  print *,"Create FFTW complex-to-real plan"
  call dfftw_plan_dft_c2r_3d(c2r, N1, N2, N3, x_cmplx, x_cmplx, FFTW_ESTIMATE)
  if (0 == c2r) goto 999

  print *,"Initialize data for complex-to-real FFT"
  call init(x_cmplx, N1, N2, N3, H1, H2, H3)

  print *,"Compute complex-to-real in-place transform"
  call dfftw_execute(c2r)

  print *,"Verify complex-to-real in-place transform"
  status = verificate(x_cmplx, N1, N2, N3, H1, H2, H3)
  if (0 /= status) goto 999

  print *,"Destroy FFTW plan"
  call dfftw_destroy_plan(c2r)

  print *,"Free data array"
  deallocate(x_cmplx)

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
  pure integer*8 function moda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    moda = mod(k8*l,m)
  end function moda

  ! Initialize x to produce unit peak at (H1,H2,H3)
  subroutine init(x, N1, N2, N3, H1, H2, H3)
    integer ::  N1, N2, N3, H1, H2, H3
    complex(WP) :: x(:,:,:)

    complex(WP), parameter :: I_TWOPI = (0,6.2831853071795864769_WP)
    integer k1,k2,k3

    forall (k1=1:N1/2+1, k2=1:N2, k3=1:N3)
      x(k1,k2,k3) = exp( -I_TWOPI * (    &
        real  (moda(H1,k1-1,N1),WP) / N1 &
        + real(moda(H2,k2-1,N2),WP) / N2 &
        + real(moda(H3,k3-1,N3),WP) / N3)) / (N1*N2*N3)
    end forall
  end subroutine init

  ! Verify that x(:,:,:) is unit peak at x(H1,H2,H3)
  integer function verificate(xc, N1, N2, N3, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3
    complex(WP) :: xc(:,:,:)

    integer k1, k2, k3
    real(WP) err, errthr, maxerr, res_exp, res_got

    ! Use Cray pointer to associate real view with complex array
    real(WP) :: x( 2*size(xc,dim=1), size(xc,dim=2), size(xc,dim=3) )
    pointer (p,x)
    p = loc(xc)

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0
    do k3 = 1, N3
      do k2 = 1, N2
        do k1 = 1, N1
          if (mod(k1-1-H1, N1)==0 .AND.                                     &
              mod(k2-1-H2, N2)==0 .AND.                                     &
              mod(k3-1-H3, N3)==0) then
            res_exp = 1
          else
            res_exp = 0
          end if
          res_got = x(k1, k2, k3)
          err = abs(res_got - res_exp)
          maxerr = max(err,maxerr)
          if (.not.(err < errthr)) then
            print '("  x("I0","I0","I0"):"$)', k1, k2, k3
            print '(" expected "G24.17","$)', res_exp
            print '(" got "G24.17","$)', res_got
            print '(" err "G10.3)', err
            print *,"  Verification FAILED"
            verificate = 1
            return
          end if
        end do
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program dp_plan_dft_c2r_3d
