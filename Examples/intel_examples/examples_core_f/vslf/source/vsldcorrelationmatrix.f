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
!    Calculation of correlation matrix  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3        ! Task dimension
      integer,parameter :: NN  = 10000    ! Number of observations

      real(kind=8),parameter :: P_THRESHOLD = 0.01

      real(kind=8) C(DIM*DIM)             ! Exact covariance matrix
      data C  / 16.0,  8.0,  4.0,                                       &
     &           8.0, 13.0, 17.0,                                       &
     &           4.0, 17.0, 62.0 /

      real(kind=8) a(DIM)                 ! Exact mean
      data a  / 3.0, 5.0, 2.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer cov_storage
      integer cor_storage
      real(kind=8) x(DIM,NN)
      real(kind=8) cov(DIM,DIM), cor(DIM,DIM)
      real(kind=8) mean(DIM), temp
      integer i, j
      integer(kind=4) errcode
      integer errnums
      integer(kind=8) task_params
      integer task_method

      real(kind=8) pval_mean(DIM)
      real(kind=8) pval_cov(DIM,DIM)

!     ***** Initialize parameters of Summary Statistics task *****
      p = DIM
      n = NN
      x_storage   = VSL_SS_MATRIX_STORAGE_COLS
      cov_storage = VSL_SS_MATRIX_STORAGE_FULL
      cor_storage = VSL_SS_MATRIX_STORAGE_FULL
      task_params = ior( VSL_SS_COV, VSL_SS_COR )
      task_method = VSL_SS_METHOD_FAST
      errcode     = 0

!     ***** Generate data set using VSL GaussianMV RNG *****
      errcode = dGenerateGaussianMVData( p, n, x, a, C )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters using FULL_STORAGE
!           for covariance/correlation matrices *****
      errcode = vsldsseditcovcor( task, mean, cov, cov_storage,         &
     &                            cor, cor_storage )
      call CheckVslError( errcode )

!     ***** Compute covariance/correlation matrices using FAST method  *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Testing stat characteristics of mean and covariance matrix *****
      call dComputePvalsMean( p, n, mean, a, C, pval_mean )
      call dComputePvalsVariance( p, n, cov, C, pval_cov )
      call dComputePvalsCovariance( p, n, cov, C, pval_cov )

      errnums = 0
      do i = 1, p
        if ( pval_mean(i) < P_THRESHOLD ) then
          errnums = errnums + 1
        end if
        do j = 1, i
          if ( pval_cov(i,j) < P_THRESHOLD ) then
            errnums = errnums + 1
          end if
        end do
      end do

!     ***** Printing results *****
      print *, 'Task dimension :         ', p
      print *, 'Number of observations : ', n
      print *, ''

!     ***** Print the exact mean, covariance and correlation matrices *****
      print *, 'Exact means'
      do i = 1, p
         print 5, a(i), ' '
      end do

      print *, ''
      print *, ''

      print *, ' Exact covariance matrix         ',                     &
     &         ' Exact correlation matrix'
      do i = 1, p
        do j = 1, p
            print 5, C((j - 1)*p + i), ' '
        end do
        print 6, '   '
        do j = 1, p
          temp = 1.0
          if ( j /= i ) then
            temp = sqrt( C((i - 1)*p + i) * C((j - 1)*p + j) )
          end if
          print 5, C((j - 1)*p + i) / temp, ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Print the computed mean, covariance and correlation matrices *****
      print *, ' Computed means'
      do i = 1, p
         print 5, mean(i), ' '
      end do

      print *, ''
      print *, ''

      print *, ' Computed covariance matrix      ',                     &
     &         ' Computed correlation matrix'
      do i = 1, p
        do j = 1, p
          print 5, cov(i, j), ' '
        end do
        print 6, '   '
        do j = 1, p
          print 5, cor(i, j), ' '
        end do
        print *, ''
      end do

      print *,''

      print *, 'P-values of the computed means'
      do i = 1, p
        print 5, pval_mean(i), ' '
      end do

      print *, ''
      print *, ''

      print *, 'P-values of the computed covariance matrix'
      do i = 1, p
        do j = 1, i
          print 5, pval_cov(i, j), ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Mean and covariance estimates agree with theory'
      else
        print *, ' Error: Mean and/or covariance estimates',            &
     &           ' disagree with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

5     format (F9.6,A,$)
6     format (A,$)

      end
