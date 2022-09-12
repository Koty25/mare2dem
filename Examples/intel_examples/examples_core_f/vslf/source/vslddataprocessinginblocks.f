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
!    Data processing in blocks  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM      = 5        ! Task dimension
      integer,parameter :: NN       = 100000   ! Number of observations
      integer,parameter :: NNBLOCKS = 100      ! Number of data portions
      integer,parameter :: BLOCK_LENGTH = NN / NNBLOCKS

      real(kind=8),parameter :: P_THRESHOLD = 0.05

!     Exact covariance matrix
      real(kind=8) C(DIM*DIM)
      data C  /  1.0, 0.05, 0.05, 0.05, 0.05,                           &
     &          0.05,  1.0, 0.05, 0.05, 0.05,                           &
     &          0.05, 0.05,  1.0, 0.05, 0.05,                           &
     &          0.05, 0.05, 0.05,  1.0, 0.05,                           &
     &          0.05, 0.05, 0.05, 0.05,  1.0 /

!     Exact means vector
      real(kind=8) a(DIM)
      data a  / -1.0,  1.0,  2.0, -3.0,  4.0 /

!     Array of accumulated weights
      real(kind=8) W(2)
      data W  / 0.0, 0.0 /

      integer indices(DIM)
      data indices / 0, 1, 1, 0, 1 /

      TYPE(VSL_SS_TASK) task
      TYPE(VSL_STREAM_STATE) stream
      integer p
      integer n
      integer nblocks, block_size
      integer x_storage
      integer cov_storage
      real(kind=8) x(DIM,BLOCK_LENGTH)
      real(kind=8) weights(NN)
      real(kind=8) cov(DIM,DIM)
      real(kind=8) mean(DIM)
      real(kind=8) r2(DIM), c2(DIM)
      real(kind=8) min_est(DIM), max_est(DIM)
      real(kind=8) T(DIM,DIM)
      integer(kind=8) estimate
      integer task_method
      integer i, j, k, order
      integer(kind=4) errcode
      integer errnums

      real(kind=8) pval_mean(DIM)
      real(kind=8) pval_cov(DIM,DIM)
      real(kind=8) pval_r2(DIM), pval_c2(DIM)

!     ***** Initializing parameters for Summary Statistics task *****
      p           = DIM
      n           = NN
      nblocks     = NNBLOCKS
      block_size  = BLOCK_LENGTH
      x_storage   = VSL_SS_MATRIX_STORAGE_COLS
      cov_storage = VSL_SS_MATRIX_STORAGE_FULL
      task_method = VSL_SS_METHOD_FAST
      errcode     = 0
      errnums     = 0

      do i = 1, n
        weights(i) = 1.0
      end do

!     ***** Generate data set using VSL GaussianMV RNG *****
      errcode = dInitGaussianMVDataGenerator( p, C, T, stream )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, block_size, x_storage, x,       &
     &                         weights, indices )
      call CheckVslError( errcode )

!     ***** Register array of weights in the task ******
      errcode = vsldssedittask( task, VSL_SS_ED_ACCUM_WEIGHT, W )
      call CheckVslError( errcode )

!     ***** Edit task parameters for computing of mean estimate and
!           2nd raw and central moments estimates *****
      errcode = vsldsseditmoments( task, mean, r2, c2m = c2 )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters using FULL_STORAGE
!           for covariance matrix computation *****
      errcode = vsldsseditcovcor( task, mean, cov, cov_storage )
      call CheckVslError( errcode )

!     ***** Edit task parameters for min and max computation *****
      errcode = vsldssedittask( task, VSL_SS_ED_MAX, max_est )
      call CheckVslError( errcode )

      errcode = vsldssedittask( task, VSL_SS_ED_MIN, min_est )
      call CheckVslError( errcode )

!     ***** Minimum and maximum are included in the list of estimates
!           to compute *****
      estimate = ior( VSL_SS_MIN, VSL_SS_MAX )

!     ***** Mean and 2nd raw and central moments are included in the list
!           of estimates to compute *****
      estimate = ior( estimate, ior( VSL_SS_MEAN,                       &
     &                ior( VSL_SS_2R_MOM, VSL_SS_2C_MOM ) ) )

!     ***** Covariance matrix is included in the list of estimates
!           to compute *****
      estimate = ior( estimate, VSL_SS_COV )

      do i = 1, nblocks
!       ***** Generate new portion of data using VSL GaussianMV RNG *****
        errcode = dGenerateGaussianMVDataBlock( p, block_size, x,       &
     &                                          stream, a, T )
        call CheckVslError( errcode )

        if ( i == 0 ) then
          do j = 1, p
            min_est(j) = x(j,1)
            max_est(j) = x(j,1)
          end do
        end if

!       ***** Compute the estimates using FAST method *****
        errcode = vsldsscompute( task, estimate, task_method )
        call CheckVslError( errcode )

!       Comparison of observations with min and max estimates
        do k = 1, p
          if( indices(k) /= 0 ) then
            do j = 1, block_size
              if (x(k, j) < min_est(k)) then
                errnums = errnums + 1
              end if
              if (x(k, j) > max_est(k)) then
                errnums = errnums + 1
              end if
            end do
          end if
        end do
      end do

!     ***** Testing stat characteristics of the computed estimates *****
!     Compute p-values for mean estimates
      call dComputePvalsMean( p, n, mean, a, C, pval_mean )
!     Compute p-values for variance estimates
      call dComputePvalsVariance( p, n, cov, C, pval_cov )
!     Compute p-values for covariance estimates
      call dComputePvalsCovariance( p, n, cov, C, pval_cov )
!     Compute p-values for raw moments estimates
      order = 2
      call dComputePvalsRawMoments( p, n, r2, order, a, C, pval_r2 )
!     Compute p-values for central moments estimates
      call dComputePvalsCentralMoments( p, n, c2, order, a, C,          &
     &                                  pval_c2 )

!     ***** Checking the validity of p-values for all estimates *****
      do i = 1, p
        if (indices(i) /= 0) then
          if (pval_mean(i) < P_THRESHOLD) then
            errnums = errnums + 1
          end if
          if (pval_r2(i) < P_THRESHOLD) then
            errnums = errnums + 1
          end if
          if (pval_c2(i) < P_THRESHOLD) then
            errnums = errnums + 1
          end if

          do j = 1, i
            if (indices(j) /= 0) then
              if (pval_cov(i, j) < P_THRESHOLD) then
                errnums = errnums + 1
              end if
            end if
          end do
        end if
      end do

!     ***** Printing results *****
      print 10, 'Task dimension :         ', p
      print 10, 'Number of observations : ', n
      print 10, 'Number of blocks :       ', nblocks
      print *, ''

!     ***** Printing computed minimum, maximum, mean
!           and moments estimates *****
      print 11, '             Min        Max        Mean       '
      print 11, '2nd_raw    2nd_cen'
      print *, ''

      do i = 1, p
        if (indices(i) /= 0) then
          print 12, 'Variable #', i
          print 13, min_est(i), ' '
          print 13, max_est(i), ' '
          print 13, mean(i), ' '
          print 13, r2(i), ' '
          print 13, c2(i), ' '
          print *, ''
        end if
      end do
      print *, ''

!     ***** Printing computed covariance matrix *****
      print *, 'Computed covariance matrix'
      do i = 1, p
        if (indices(i) /= 0) then
          print 12, 'Variable #', i

          do j = 1, p
            if (indices(j) /= 0) then
              print 13, cov(i,j), ''
            end if
          end do
          print *, ''
        end if
      end do
      print *, ''
      print *, ''

!     ***** Printing p-values for mean and moments estimates *****
      print *, 'P-values of the computed estimates'
      print *, ''

      print *, '            Mean       2nd_raw    2nd_cen'
      do i = 1, p
        if (indices(i) /= 0) then
          print 12, 'Variable #', i
          print 13, pval_mean(i), ' '
          print 13, pval_r2(i), ' '
          print 13, pval_c2(i), ' '
          print *, ''
        end if
      end do
      print *, ''

!     ***** Printing p-values for covariance matrix estimate *****
      print *, 'Covariance matrix'
      do i = 1, p
        if (indices(i) /= 0) then
          print 12, 'Variable #', i

          do j = 1, p
            if (indices(j) /= 0) then
              print 13, pval_cov(i,j), ''
            end if
          end do
          print *, ''
        end if
      end do
      print *, ''
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

10    format(A,I6)
11    format(A,$)
12    format(A,I1,$)
13    format(F10.6,A,$)
14    format(5F8.5)

      end
