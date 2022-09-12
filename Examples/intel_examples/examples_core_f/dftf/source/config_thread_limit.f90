!===============================================================================
! Copyright 2012-2020 Intel Corporation.
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
! An example of using DFTI_THREAD_LIMIT configuration parameter.
! The parameter specifies maximum number of OpenMP threads FFT can use.
!
! Values:
!   0 (default) = use number of threads specified by
!                 mkl_[domain_]set_num_threads()
!   Any positive integer N = use not more than N threads
!
!*****************************************************************************
program config_thread_limit

#if defined(_OPENMP)
  use omp_lib
#endif

  include 'mkl_service.fi'

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Number of parallel user threads
#if defined(_OPENMP)
  integer, parameter :: NUT = 2
#endif

  integer :: failed = 0
  integer :: thr_failed, me, team
  integer*4 :: flag = 0

  print *,"Example config_thread_limit"

  ! Enable nested parallel OpenMP sections (maybe oversubscribed)
  ! Use kind(4) literals because some compilers fail on this in i8 mode
#if defined(_OPENMP)
  call omp_set_nested(.true._4)
  call omp_set_dynamic(.false._4)
#endif
 
  ! Enable threading of Intel(R) Math Kernel Library (Intel(R) MKL) called from OpenMP parallel sections
  call mkl_set_dynamic(flag)

#if defined(_OPENMP)  
  print '(" Run parallel FFTs on "I0" parallel threads")', NUT
!$omp parallel num_threads(NUT) default(shared) private(thr_failed, me)
#else
  print '(" Run parallel FFT on a single thread")'
#endif

#if defined(_OPENMP)
  me = omp_get_thread_num()
  team = omp_get_num_threads()
#else
  me = 0
  team = 1
#endif

  if (me == 0) then
    print '(" Thread "I0": parallel team is "I0" threads")', me, team
  endif

  if (me == 0) then
    thr_failed = run_dft(me, 2, 100,200,300, -1,-2,-3)
  else
    thr_failed = run_dft(me, 3, 200,300,100, -1,-2,-3)
  endif
  if (0 /= thr_failed) failed = thr_failed

#if defined(_OPENMP)
!$omp end parallel
#endif

  if (failed == 0) then
    print *,"TEST PASSED"
    call exit(0)
  else
    print *,"TEST FAILED"
    call exit(1)
  endif

contains

  integer function run_dft(tid,tlimit,n1,n2,n3,h1,h2,h3)

    use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

    integer :: tid              ! Id of this thread
    integer :: tlimit           ! Thread limit
    integer :: N1,N2,N3         ! Sizes of 3D transform
    integer :: H1,H2,H3         ! Arbitrary harmonic used to test the FFT


    ! Execution status
    integer :: status = 0, ignored_status

    ! DFTI descriptor handle
    type(DFTI_DESCRIPTOR), POINTER :: hand

    ! Data array
    complex(WP), allocatable :: x (:,:,:)

    ! Temporary variable
    integer :: tl

    print '(" Thread "I0": Create DFTI descriptor for "I0"x"I0"x"I0" FFT")', &
      & tid, N1, N2, N3
    status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, [N1,N2,N3])
    if (0 /= status) goto 999

    print '(" Thread "I0": Set number of user threads "I0)', tid, tlimit
    status = DftiSetValue(hand, DFTI_THREAD_LIMIT, tlimit)
    if (0 /= status) goto 999

    ! If tlimit > 1 check if we linked with sequential MKL
    if (tlimit > 1) then
      ! Get thread limit of uncommitted descriptor
      status = DftiGetValue(hand, DFTI_THREAD_LIMIT, tl)
      if (0 /= status) goto 999
      print '(" Thread "I0": uncommitted descriptor thread limit "I0)', tid, tl
    endif

    print '(" Thread "I0": Commit DFTI descriptor")', tid
    status = DftiCommitDescriptor(hand)
    if (0 /= status) goto 999

    ! Get thread limit of committed descriptor
    status = DftiGetValue(hand, DFTI_THREAD_LIMIT, tl)
    if (0 /= status) goto 999
    print '(" Thread "I0": committed descriptor thread limit "I0)', tid, tl

    print '(" Thread "I0": Allocate data array")', tid
    allocate ( x(N1, N2, N3), STAT = status)
    if (0 /= status) goto 999

    print '(" Thread "I0": Initialize input for forward transform")', tid
    call init(x, N1, N2, N3, H1, H2, H3)

    print '(" Thread "I0": Compute forward transform")', tid
    status = DftiComputeForward(hand, x(:, 1, 1))
    if (0 /= status) goto 999

    print '(" Thread "I0": Verify the result")', tid
    status = verificate(x, N1, N2, N3, H1, H2, H3)
    if (0 /= status) goto 999

100 continue

    print '(" Thread "I0": Release the DFTI descriptor")', tid
    ignored_status = DftiFreeDescriptor(hand)

    if (allocated(x)) then
      print '(" Thread "I0": Deallocate data array")', tid
      deallocate(x)
    endif

    if (status == 0) then
      print '(" Thread "I0": Subtest Passed")', tid
    else
      print '(" Thread "I0": Subtest Failed")', tid
    endif

    run_dft = status
    return

999 print '(" Thread "I0":  Error, status = ",I0)', tid, status
    goto 100

  end function run_dft

  ! Compute mod(K*L,M) accurately
  pure real(WP) function moda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    moda = real(mod(k8*l,m),WP)
  end function moda

  ! Initialize arrays with harmonic /H1, H2, H3/
  subroutine init(x, N1, N2, N3, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3
    complex(WP) :: x(:,:,:)

    integer k1, k2, k3
    complex(WP), parameter :: I_TWOPI = (0.0_WP,6.2831853071795864769_WP)

    forall (k1=1:N1, k2=1:N2, k3=1:N3)
      x(k1,k2,k3) = exp( I_TWOPI * ( &
    &    moda(  k1-1,H1, N1)/N1 &
    &    + moda(k2-1,H2, N2)/N2 &
    &    + moda(k3-1,H3, N3)/N3 )) / (N1*N2*N3)
    end forall
  end subroutine init

  ! Verify that x(:,:,:) are unit peaks at /H1, H2, H3/
  integer function verificate(x, N1, N2, N3, H1, H2, H3)
    integer N1, N2, N3, H1, H2, H3
    complex(WP) :: x(:,:,:)

    integer k1, k2, k3
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 5.0 * log(real(N1*N2*N3,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k3 = 1, N3
      do k2 = 1, N2
        do k1 = 1, N1
          if (mod(k1-1-H1, N1)==0 .AND.                                &
    &        mod  (k2-1-H2, N2)==0 .AND.                                &
    &        mod  (k3-1-H3, N3)==0) then
            res_exp = 1.0_WP
          else
            res_exp = 0.0_WP
          end if
          res_got = x(k1,k2,k3)
          err = abs(res_got - res_exp)
          maxerr = max(err,maxerr)
          if (.not.(err <= errthr)) then
            print '("  x("I0","I0","I0"):"$)', k1, k2, k3
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
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end program config_thread_limit
