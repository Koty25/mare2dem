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
! An example of using Intel(R) Math Kernel Library (Intel(R) MKL) DFTI
! configuration parameter DFTI_PLACEMENT.
! The parameter defines if the result overwrites the input data or not.
!
! Values:
!   DFTI_INPLACE (default) - result overwrites input data
!   DFTI_NOT_INPLACE       - result is placed in a separate array
!
! Note 1: In many cases in-place transforms need storage association of input
!         and output. When input and output have different type (e.g. real
!         transforms), the data can be storage associated by using Cray
!         pointers extension of Fortran.
!
! Note 2: When storage data types of forward and backward domains are
!         the same, the configuration parameters for the layout of input are
!         also used for the layout of output (e.g. output strides are
!         ignored). Otherwise, both input and output layout shall be defined
!         (for example real transform with conjugate even storage
!         set to DFTI_COMPLEX_COMPLEX).
!
! This example computes an in-place and an out-of-place 2D real transform.
!
!*****************************************************************************

program config_placement

  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Sizes of 2D transform
  integer, parameter :: N1 = 100
  integer, parameter :: N2 = 200

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = 2

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! Data arrays
  real(WP), allocatable :: xr(:,:)
  complex(WP), allocatable :: xc(:,:)

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  ! Strides define data layout for real and complex domains
  integer cstrides(3), rstrides(3)

  hand => null()

  print *,"Example config_placement"
  print *,"In-place and out-of-place 2D real FFT"
  print *,"Configuration parameters:"
  print '(" DFTI_FORWARD_DOMAIN = DFTI_REAL")'
  print '(" DFTI_PRECISION      = DFTI_DOUBLE")'
  print '(" DFTI_DIMENSION      = 2")'
  print '(" DFTI_LENGTHS        = /"I0","I0"/" )', N1, N2

  print *,"======= In-place 2D real FFT ======="

  print *,"Allocate data array"
  allocate ( xr(2*(N1/2+1), N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Create DFTI descriptor"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 2, [N1,N2])
  if (0 /= status) goto 999

  print *,"Set CCE storage"
  status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, &
    &                    DFTI_COMPLEX_COMPLEX)
  if (0 /= status) goto 999

  rstrides = [0, 1, size(xr,dim=1)   ]
  cstrides = [0, 1, size(xr,dim=1)/2 ]

  print '(" Set input  strides = "3(I0:", "))', rstrides
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides)
  if (0 /= status) goto 999

  ! This must be done because type of elements in forward
  ! and backward domain is different
  print '(" Set output strides = "3(I0:", "))', cstrides
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize data for real-to-complex FFT"
  call init_r(xr, N1, N2, H1, H2)

  print *,"Compute forward in-place transform"
  status = DftiComputeForward(hand, xr(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result by viewing xr as a complex array"
  status = verify_r_as_c(xr, N1, N2, H1, H2)
  if (0 /= status) goto 999

  print *,"======= Reconfigure for out-on-place 2D real FFT ======="

  print *,"Allocate data array for output"
  allocate ( xc(N1/2+1, N2), STAT = status)
  if (0 /= status) goto 999

  print *,"Set out-of-place"
  status = DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  if (0 /= status) goto 999

  rstrides = [0, 1, size(xr,dim=1) ]
  cstrides = [0, 1, size(xc,dim=1) ]

  print '(" Set input  strides = "3(I0:", "))', rstrides
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, rstrides)
  if (0 /= status) goto 999

  print '(" Set output strides = "3(I0:", "))', cstrides
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, cstrides)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize data for real-to-complex FFT"
  call init_r(xr, N1, N2, H1, H2)

  print *,"Compute forward out-of-place transform"
  status = DftiComputeForward(hand, xr(:,1), xc(:,1) )
  if (0 /= status) goto 999

  print *,"Verify the result"
  status = verify_c(xc, N1, N2, H1, H2)
  if (0 /= status) goto 999

100 continue

  print *,"Release the DFTI descriptor"
  ignored_status = DftiFreeDescriptor(hand)

  if (allocated(xr) .or. allocated(xc)) then
    print *,"Deallocate data arrays"
  endif
  if (allocated(xr))  deallocate(xr)
  if (allocated(xc))  deallocate(xc)

  if (status == 0) then
    print *, "TEST PASSED"
    call exit(0)
  else
    print *, "TEST FAILED"
    call exit(1)
  end if

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

  ! Initialize x(:,:) to harmonic H
  subroutine init_r(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)

    integer k1, k2
    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP
    real(WP) :: factor

    factor = 2
    if (mod(2*(N1-H1),N1)==0 .and. mod(2*(N2-H2),N2)==0) factor = 1
    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = factor * cos(TWOPI*( moda(H1,k1-1,N1) / N1   &
    &    +                             moda(H2,k2-1,N2) / N2)) / (N1*N2)
    end forall
  end subroutine init_r

  ! Verify that r(:,:) has unit peak at (H1,H2)
  integer function verify_r_as_c(r, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: r(:,:)

    ! Use Cray pointer to associate complex view with the real array r.
    ! The first dimension of r must be divisible by 2.
    complex(WP) :: c( 1:size(r,dim=1)/2, 1:size(r,dim=2) )
    pointer (p,c)

    ! Associate complex view with r
    p = loc(r)

    verify_r_as_c = verify_c(c, N1, N2, H1, H2)
  end function verify_r_as_c

  ! Verify that x(:,:) has unit peak at (H1,H2)
  integer function verify_c(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    complex(WP) :: x(:,:)

    integer k1, k2
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k2 = 1, N2
      do k1 = 1, N1/2+1
        if (mod(k1-1-H1,N1)==0 .and. mod(k2-1-H2,N2)==0) then
          res_exp = 1.0_WP
        else if (mod(-k1+1-H1,N1)==0 .and. mod(-k2+1-H2,N2)==0) then
          res_exp = 1.0_WP
        else
          res_exp = 0.0_WP
        end if
        res_got = x(k1, k2)
        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1,k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', res_got
          print '(" err "G10.3)', err
          print *,"  Verification FAILED"
          verify_c = 100
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verify_c = 0
  end function verify_c

end program config_placement
