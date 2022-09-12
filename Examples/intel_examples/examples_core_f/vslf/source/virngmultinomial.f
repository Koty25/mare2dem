!===============================================================================
! Copyright 2018-2020 Intel Corporation.
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
!    viRngMultinomial  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=4) i,j,nn
      integer n
      integer(kind=4) errcode

      integer(kind=4) ntrial
      integer(kind=4) k
      real(kind=8) p(5)
      integer(kind=4) r(5000)
      integer brng,method,seed

      real(kind=8) tM,tD,tQ,tD2
      real(kind=8) sM,sD
      real(kind=8) sum, sum2
      real(kind=8) s
      real(kind=8) DeltaM,DeltaD

      TYPE (VSL_STREAM_STATE) :: stream

      n=1000
      nn=10

      brng=VSL_BRNG_MCG31
      method=VSL_RNG_METHOD_MULTINOMIAL_MULTPOISSON
      seed=1

      ntrial=10
      k=5
      do i=1,k
        p(i)=0.2
      end do

!     ***** Initialize *****
      errcode=vslnewstream( stream, brng,  seed )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode=virngmultinomial( method, stream, n, r, ntrial, k, p )
      call CheckVslError(errcode)

!     ***** Printing results *****
      print *,"Sample of viRngMultinomial."
      print *,"------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    ntrial=",ntrial
      print 11,"    k=",k
      print *,""
      print *,"    probability vector:"

      do i=1,k
        write(*, 13, advance="no") "   ",p(i)
      end do
      print *,""

      print *,""
      print *,"Results (first 10 of 1000):"
      print *,"---------------------------"
      do i=1,nn
        do j=1,k
          write(*, 10, advance="no") r(i*k+j)
        end do
        print *,""
      end do

      do i=1,k

!       !***** Theoretical moments *****
        tM=ntrial*p(i)
        tD=ntrial*p(i)*(1.0-p(i))
        tQ=ntrial*(1.0-p(i))*(4.0*ntrial*p(i)-
     &  4.0*ntrial*p(i)*p(i)+4.0*p(i)-3.0)

!       ***** Sample moments *****
        sum=0.0
        sum2=0.0
        do j=1,n-1
          sum=sum+r(j*k+i)
          sum2=sum2+r(j*k+i)*r(j*k+i)
        end do
        sM=sum/n
        sD=sum2/n-sM*sM

!       ***** Comparison of theoretical and sample moments *****
        tD2=tD*tD
        s=((tQ-tD2)/n)-(2*(tQ-2*tD2)/(n*n))+((tQ-3*tD2)/(n*n*n))
        DeltaM=(tM-sM)/sqrt(tD/n)
        DeltaD=(tD-sD)/sqrt(s)

        print *,""
        if (abs(DeltaM)>3.0 .OR. abs(DeltaD)>3.0) then
          print 12,"Error: Category=",i," sample moments (mean=",         &
     &    sM,", variance=",sD,                                            &
     &    ") disagree with theory (mean=",                                &
     &    tM,", variance=",tD,")."
          stop 1
        else
          print 12,"Category=",i," sample moments (mean=",sM,              &
     &    ", variance=",sD,") agree with theory (mean=",                  &
     &    tM,", variance=",tD,")."
        end if
      end do

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)

10    format(I5)
11    format(A,I5)
12    format(A,I5,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)
13    format(A,F5.2)

      end
