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
!    Calculation of cross-product matrix Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 4        ! Task dimension
      integer,parameter :: NN  = 1000     ! Number of observations

      real(kind=4),parameter :: P_THRESHOLD = 0.01

      real(kind=4) C(DIM*DIM)
      data C / 1.0, 0.0, 0.0, 0.0,                                      &
     &         0.0, 1.0, 0.0, 0.0,                                      &
     &         0.0, 0.0, 1.0, 0.0,                                      &
     &         0.0, 0.0, 0.0, 1.0 /

      real(kind=4) a(DIM)
      data a / 5.0, 5.0, 5.0, 5.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer cov_storage
      integer cp_storage
      real(kind=4) x(DIM,NN)
      real(kind=4) mean(DIM)
      real(kind=4) cov(DIM,DIM), cp(DIM,DIM)
      integer i, j, order
      integer(kind=4) errcode
      integer errnums, method
      integer(kind=8) estimate

      real(kind=4) pval_cov(DIM,DIM)

!     ***** Initializing parameters for Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_COLS
      cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
      cp_storage      = VSL_SS_MATRIX_STORAGE_FULL

!     ***** Generate data set using VSL GaussianMV RNG *****
      errcode = sGenerateGaussianMVData( p, n, x, a, C )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )


!    ***** Initialization of the task parameters related
!          to cross-product matrix computation *****/
      errcode = vslssseditcp(task, mean, cp=cp, cp_storage=cp_storage)
      call CheckVslError( errcode )

!    ***** Cross-product matrix is included in the list of estimates
!          to compute *****
      estimate = VSL_SS_CP
      method   = VSL_SS_METHOD_FAST

!     ***** Compute the estimates using FAST method *****
      errcode = vslssscompute( task, estimate, method )
      call CheckVslError( errcode )

!     ***** Edit task parameters for computation of covariance *****
      errcode = vslsssedittask( task,VSL_SS_ED_COV,cov )
      call CheckVslError( errcode )
      errcode = vslissedittask( task,VSL_SS_ED_COV_STORAGE,cov_storage)
      call CheckVslError( errcode )

!    ***** Convert cross-product matrix into correlation matrix *****
      estimate = VSL_SS_COV
      method = VSL_SS_METHOD_CP_TO_COVCOR
      errcode = vslssscompute( task, estimate, method );
      call CheckVslError( errcode )


!     ***** Testing stat characteristics of the computed estimates *****
!     ***** Compute p-values for covariance computed from cross-product *
      call sComputePvalsVariance( p, n, cov, C, pval_cov )
      call sComputePvalsCovariance( p, n, cov, C, pval_cov )

!     ***** Checking the validity of p-values for all estimates *****
      errnums = 0
      do i = 1, p
        do j = 1, i
          if (pval_cov(i, j) < P_THRESHOLD) then
            errnums = errnums + 1
          end if
        end do
      end do

!     ***** Printing results *****
      print 9, 'Task dimension :         ', p
      print 9, 'Number of observations : ', n
      print *, ''

!     ***** Printing computed cross-produt matrix *****
      print 10, '  Computed cross-product matrix               '
      print *, ''
      do i = 1, p
        write (*, 14) cp(i, 1:p)
        print *, ''
      end do
      print *, ''
      print *, ''


!     ***** Printing p-values for covariance matrix estimate *****
      print *, 'P-values for covariance obtained from cross-product'
      do i = 1, p
        write (*, 14) pval_cov(i, 1:p)
        print *, ''
      end do
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' All the computed estimates agree with theory'
      else
        print *, ' Error: At least one of the computed estimates',      &
     &           ' disagrees with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

9     format(A,I4)
10    format(A,$)
11    format(A,I1,A,$)
12    format(F11.6,A,$)
13    format(F11.6,A,F11.6,A,F11.6,A,$)
14    format(4F14.6,$)

      end
