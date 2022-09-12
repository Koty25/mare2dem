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
! An example of using DFTI_NUMBER_OF_TRANSFORMS configuration parameter.
! The parameter defines how many identical transforms are computed by one call
! of DftiComputeForward or DftiComputeBackward function.
!
! Values:
! Any positive integer (default 1)
!
!*****************************************************************************
program config_number_of_transforms

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Sizes of 3D transform
  integer, parameter :: N1 = 8
  integer, parameter :: N2 = 12
  integer, parameter :: N3 = 16

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2
  integer, parameter :: H3 = 3

  ! Number of transforms
  integer, parameter :: M = 100

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! Data array
  complex(WP), allocatable :: x (:,:,:,:)

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  ! Distance between first elements for multiple transforms
  integer :: dist

  hand => null()

  print *,"Example config_number_of_transforms"
  print *,"Forward and backward multiple 3D FFT"
  print *,"Configuration parameters:"
  print '(" DFTI_PRECISION            = DFTI_DOUBLE")'
  print '(" DFTI_FORWARD_DOMAIN       = DFTI_COMPLEX")'
  print '(" DFTI_DIMENSION            = 3")'
  print '(" DFTI_LENGTHS              = /"I0","I0","I0"/" )', N1, N2, N3
  print '(" DFTI_NUMBER_OF_TRANSFORMS = " I0)', M


  print *,"Allocate data array"
  allocate ( x(N1, N2, N3, M), STAT = status)
  if (0 /= status) goto 999

  print *,"Create DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, [N1,N2,N3])
  if (0 /= status) goto 999

  print *,"Set number of transforms"
  status = DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, M)
  if (0 /= status) goto 999

  dist = size(x,dim=1) * size(x,dim=2) * size(x,dim=3)

  print '(" Set input/output distance = "I0)', dist
  status = DftiSetValue(hand, DFTI_INPUT_DISTANCE, dist)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize input for forward transform"
  call init(x, N1, N2, N3, M, H1, H2, H3)

  print *,"Compute forward transform:"
  status = DftiComputeForward(hand, x(:,1,1,1))
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verificate(x, N1, N2, N3, M, H1, H2, H3)
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

  ! Initialize arrays with harmonic /H1, H2, H3/
  subroutine init(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
    complex(WP) :: x(:,:,:,:)

    integer k1, k2, k3, mm
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2, k3=1:N3, mm=1:M)
      x(k1,k2,k3,mm) = exp( I_TWOPI * ( &
        moda(  k1-1,H1, N1)/N1 &
        + moda(k2-1,H2, N2)/N2 &
        + moda(k3-1,H3, N3)/N3 )) / (N1*N2*N3)
    end forall
  end subroutine init

  ! Verify that x(:,:,:,mm) are unit peaks at /H1, H2, H3/
  integer function verificate(x, N1, N2, N3, M, H1, H2, H3)
    integer N1, N2, N3, M, H1, H2, H3
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
              mod  (k2-1-H2, N2)==0 .AND.                                &
              mod  (k3-1-H3, N3)==0) then
              res_exp = 1.0_WP
            else
              res_exp = 0.0_WP
            end if
            res_got = x(k1,k2,k3,mm)
            err = abs(res_got - res_exp)
            maxerr = max(err,maxerr)
            if (.not.(err < errthr)) then
              print '("  x("I0","I0","I0", "I0"):"$)', k1, k2, k3, mm
              print '(" expected ("G24.17","G24.17"),"$)', res_exp
              print '(" got ("G24.17","G24.17"),"$)', res_got
              print '(" err "G10.3)', err
              print *," Verification FAILED"
              verificate = 100
              return
            end if
          end do
        end do
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program config_number_of_transforms