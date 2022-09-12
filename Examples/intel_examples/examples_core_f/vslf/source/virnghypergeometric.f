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
!    viRngHypergeometric  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=4) i,nn
      integer n
      integer(kind=4) errcode

      integer(kind=4) l,ss,m
      integer(kind=4) r(1000)
      integer brng,method,seed

      real(kind=8) tM,tD,tQ,tD2
      real(kind=8) sM,sD
      real(kind=8) sum, sum2
      real(kind=8) s
      real(kind=8) DeltaM,DeltaD
      real(kind=8) K,L2,L3,L4,L5,L6,KL,KL4,S2,S3,S4,M2,M3,M4
      real(kind=8) r_l,r_ss,r_m

      TYPE (VSL_STREAM_STATE) :: stream

      n=1000
      nn=10

      brng=VSL_BRNG_MCG31
      method=VSL_RNG_METHOD_HYPERGEOMETRIC_H2PE
      seed=1

      l=100
      ss=10
      m=30

      r_l=l
      r_ss=ss
      r_m=m

!     ***** Initialize *****
      errcode=vslnewstream( stream, brng,  seed )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode=virnghypergeometric( method, stream, n, r, l, ss, m )
      call CheckVslError(errcode)

!     ***** Theoretical moments *****
      K = (r_l-1)*(r_l-2)*(r_l-3)
      L2 = r_l*r_l
      L3 = L2*r_l
      L4 = L2*L2
      L5 = L3*L2
      L6 = L3*L3
      KL = K*r_l
      KL4 = K*L4
      S2 = r_ss*r_ss
      S3 = S2*r_ss
      S4 = S2*S2
      M2 = r_m*r_m
      M3 = M2*r_m
      M4 = M2*M2

      tM=r_m*r_ss/r_l
      tD=(r_m*r_ss*(r_l-r_m)*(r_l-r_ss))/(r_l*r_l*(r_l-1))
      tQ=( (3*r_l+18)    *S4/KL4 - (6*L2+36*r_l)  *S3/KL4               &
     &   + (3*L3+24*L2)   *S2/KL4 - 6        *r_ss/KL  ) * M4           &
     &   + ( (-6*L2-36*r_l)*S4/KL4 + (12*L3+72*L2)  *S3/KL4             &
     &   - (6*L4+38*L3)   *S2/KL4 + 12       *r_ss/K   ) * M3           &
     &   + ( (3*L3+24*L2)  *S4/KL4 - (6*L4+48*L3)   *S3/KL4             &
     &   + (31*L4+3*L5+L3)*S2/KL4 - (L4+7*L5)*r_ss/KL4 ) * M2           &
     &   + ( -6            *S4/KL  + 12             *S3/K               &
     &   - (4*L4+7*L5)    *S2/KL4 + (L6+L5)  *r_ss/KL4 ) * r_m

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
      print *,"Sample of viRngHypergeometric."
      print *,"------------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    l=",l
      print 11,"    s=",ss
      print 11,"    m=",m

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
11    format(A,I5)
12    format(A,F5.2,A,F5.2,A,F5.2,A,F5.2,A)

      end
