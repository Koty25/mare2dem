!===============================================================================
! Copyright 2020 Intel Corporation.
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
!  Content:
!      Intel(R) oneMKL FORTRAN vsrnggaussian example for VSL openMP offload
!*******************************************************************************

include "mkl_omp_offload.f90"
include "errcheck.inc"

program vsrnggaussian_example
    use onemkl_vsl_omp_offload
    use MKL_VSL_TYPE

    type (VSL_STREAM_STATE) :: stream
    integer :: i, j, n, nn
    integer :: brng,method,seed
    integer :: errcode
    real a,sigma
    real,allocatable :: r(:)
    real :: tM,tD,tQ,tD2
    real :: sM,sD
    real :: sum, sum2
    real :: s
    real :: DeltaM,DeltaD

    n=1000
    nn=10

    brng=VSL_BRNG_PHILOX4X32X10
    method=VSL_RNG_METHOD_GAUSSIAN_ICDF
    seed=1

    a=0.0
    sigma=1.0

    allocate(r(n))
    if (.not. allocated(r)) goto 998

    ! Initialize stream
    errcode=vslnewstream(stream, brng, seed)
    call CheckVslError(errcode)

    ! Calling vsrnggaussian on the GPU
    !$omp target data map(r)
    !$omp target variant dispatch use_device_ptr(r)
    errcode=vsrnggaussian( method, stream, n, r, a, sigma)
    !$omp end target variant dispatch
    !$omp end target data
    call CheckVslError(errcode)

 !     ***** Theoretical moments *****
    tM=a
    tD=sigma*sigma
    tQ=720.0*sigma*sigma*sigma*sigma

!     ***** Sample moments *****
    sum=0.0
    sum2=0.0
    do i=1,n
    sum=sum+r(i)
    sum2=sum2+r(i)*r(i)
    end do
    sM=sum/n
    sD=sum2/n-sM*sM

!     ***** Comparison of theoretical and sample moments *****
    tD2=tD*tD
    s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n))
    DeltaM=(tM-sM)/sqrt(tD/n)
    DeltaD=(tD-sD)/sqrt(s)

!     ***** Printing results *****
    print *,"Sample of vsRngGaussian."
    print *,"------------------------"
    print *,""
    print *,"Parameters:"
    print 11,"    a=",a
    print 11,"    sigma=",sigma

    print *,""
    print *,"Results (first 10 of 1000):"
    print *,"---------------------------"
    do i=1,nn
    print 10,r(i)
    end do

    print *,""
    if (abs(DeltaM)>3.0 .OR. abs(DeltaD)>3.0) then
        print 12,"Error: sample moments (mean=",                        &
    &    sM,", variance=",sD,                                          &
    &    ") disagree with theory (mean=",                              &
    &    tM,", variance=",tD,")."
        deallocate(r)
        stop 1
    else
        print 12,"Sample moments (mean=",sM,                            &
    &    ", variance=",sD,") agree with theory (mean=",                &
    &    tM,", variance=",tD,")."
    end if

    ! Deinitialize stream
    errcode=vsldeletestream(stream)


    deallocate(r)

    stop

    998 print *, 'Error: cannot allocate memory' 
    999 stop 1

    10    format(F7.3)
    11    format(A,F5.3)
    12    format(A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)

end program
