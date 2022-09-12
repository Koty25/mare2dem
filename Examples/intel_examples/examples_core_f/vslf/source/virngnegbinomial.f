!===============================================================================
! Copyright 2003-2020 Intel Corporation.
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
!    viRngNegbinomial  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=4) i,nn
      integer n
      integer(kind=4) errcode

      real(kind=8) a,p
      integer(kind=4) r(1000)
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
      method=VSL_RNG_METHOD_NEGBINOMIAL_NBAR
      seed=1

      a=1.0
      p=0.75

!     ***** Initialize *****
      errcode=vslnewstream( stream, brng,  seed )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode=virngnegbinomial( method, stream, n, r, a, p )
      call CheckVslError(errcode)

!     ***** Theoretical moments *****
      tM=a*(1-p)/p
      tD=a*(1-p)/(p*p)
      tQ=(a*(1-p)*(p*p-3*a*p-6*p+3*a+6))/(p*p*p*p)

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
      print *,"Sample of viRngNegbinomial."
      print *,"---------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    a=",a
      print 11,"    p=",p

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
        stop 1
      else
        print 12,"Sample moments (mean=",sM,                            &
     &    ", variance=",sD,") agree with theory (mean=",                &
     &    tM,", variance=",tD,")."
      end if

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)

10    format(I5)
11    format(A,F5.2)
12    format(A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)

      end
