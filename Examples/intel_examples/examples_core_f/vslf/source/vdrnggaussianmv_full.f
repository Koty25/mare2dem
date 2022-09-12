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
!    vdRngGaussianMV  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "statcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=4) i
      integer(kind=4) errcode

      integer(kind=4) nn
      integer ndim,info
      integer n

      parameter(n=1000,nn=5,ndim=3)

      integer brng,method,seed
      integer me

      real(kind=8) c(ndim,ndim),t(ndim,ndim),a(ndim)
      real(kind=8) r(ndim,n)
      real(kind=8) dbS(ndim),dbS2(ndim),dbMean(ndim),dbVar(ndim)
      real(kind=8) dbCovXY,dbCovXZ,dbCovYZ

      real(kind=8) S(ndim),D2(ndim),Q(ndim)
      real(kind=8) DeltaM(ndim),DeltaD(ndim)

      TYPE (VSL_STREAM_STATE) :: stream

      brng=VSL_BRNG_MCG31
      seed=7777777
      method=VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
      me=VSL_MATRIX_STORAGE_FULL

!     Variance-covariance matrix for test
!     (should be symmetric,positive-definite)

!     This is full storage for dpotrf subroutine
      c(1,1)=16.0D0
      c(1,2)=8.0D0
      c(1,3)=4.0D0

      c(2,1)=8.0D0
      c(2,2)=13.0D0
      c(2,3)=17.0D0

      c(3,1)=4.0D0
      c(3,2)=17.0D0
      c(3,3)=62.0D0

      a(1)=3.0D0
      a(2)=5.0D0
      a(3)=2.0D0

      t  = c

      print *,'Variance-covariance matrix C'
      write (*,'(3F7.3)') c
      print *,''

      print *,'Mean vector a:'
      write (*,'(3F7.3)') a
      print *,''

      print *,'VSL_MATRIX_STORAGE_FULL'
      print *,'-----------------------'
      print *,''

      call dpotrf('U',ndim,t,ndim,info)

!     Stream initialization
      errcode=vslnewstream(stream,brng,seed)
      call CheckVslError(errcode)

!     Generating random numbers
!     from multivariate normal distribution
      errcode=vdrnggaussianmv(method,stream,n,r,ndim,me,a,t)
      call CheckVslError(errcode)

!     Printing random numbers
      print 11,' Results (first ',nn,' of ',n,')'
      print *,'--------------------------'
11    format(A,I1,A,I5,A)

      do i=1,nn
        print 12,' r(',i,')=(',r(:,i),')'
      end do
12    format(A,I1,A,3F8.3,A)
      print *,''

      call dCalculateGaussianMVSampleCharacteristics(ndim, n, r,        &
     &      dbS, dbS2, dbMean, dbVar, dbCovXY, dbCovXZ, dbCovYZ)

!     Printing
      print *,'Sample characteristics:'
      print *,'-----------------------'
      print *,'      Sample             Theory'
      print 13,' Mean :(',dbMean(1),dbMean(2),dbMean(3),                &
     &         ')  (',a(1),a(2),a(3),')'
      print 13,' Var. :(',dbVar(1),dbVar(2),dbVar(3),                   &
     &         ')  (',c(1,1),c(2,2),c(3,3),')'
      print 14,' CovXY: ',dbCovXY,'          ',c(1,2)
      print 14,' CovXZ: ',dbCovXZ,'          ',c(1,3)
      print 14,' CovYZ: ',dbCovYZ,'          ',c(2,3)
13    format(A,F5.1,F5.1,F5.1,A,F5.1,F5.1,F5.1,A)
14    format(A,F6.1,A,F6.1)
      print *,''

      errcode=dGaussianMVCheckResults(ndim, n, a, c, dbMean, dbVar, S,  &
     &      D2, Q, DeltaM, DeltaD)

      if (errcode /= 0) then
        print *,"Error: sample moments"
        print *,"disagree with theory"
        print 15, "    DeltaM: ", DeltaM(1), DeltaM(2), DeltaM(3)
        print 15, "    DeltaD: ", DeltaD(1), DeltaD(2), DeltaD(3)
        print *,  "   ( at least one of the Deltas > 3.0) "
        stop 1
      else
        print *,"Sample moments"
        print *,"agree with theory"
        print 15, "    DeltaM: ", DeltaM(1), DeltaM(2), DeltaM(3)
        print 15, "    DeltaD: ", DeltaD(1), DeltaD(2), DeltaD(3)
        print *,  "   ( all Deltas < 3.0) "
      end if
15    format(A,F7.3,F7.3,F7.3)
      print *,''

!     Stream finalization
      errcode=vslDeleteStream(stream)
      call CheckVslError(errcode)

      end
