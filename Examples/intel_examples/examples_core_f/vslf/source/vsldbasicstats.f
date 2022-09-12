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
!    Basic functionality Example Program Text
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

      real(kind=8),parameter :: P_THRESHOLD = 0.01

      real(kind=8) C(DIM*DIM)
      data C / 1.0, 0.0, 0.0, 0.0,                                      &
     &         0.0, 1.0, 0.0, 0.0,                                      &
     &         0.0, 0.0, 1.0, 0.0,                                      &
     &         0.0, 0.0, 0.0, 1.0 /

      real(kind=8) a(DIM)
      data a / 5.0, 5.0, 5.0, 5.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer cov_storage
      integer cor_storage
      real(kind=8) x(DIM,NN)
      real(kind=8) cov(DIM,DIM), cor(DIM,DIM)
      real(kind=8) mean(DIM)
      real(kind=8) min_estimate(DIM), max_estimate(DIM)
      real(kind=8) raw2(DIM), raw3(DIM), raw4(DIM)
      real(kind=8) cen2(DIM), cen3(DIM), cen4(DIM)
      real(kind=8) skewness(DIM), kurtosis(DIM), variation(DIM)
      integer i, j, order
      integer(kind=4) errcode
      integer errnums
      integer(kind=8) estimate
      integer task_method

      real(kind=8) pval_mean(DIM)
      real(kind=8) pval_cov(DIM,DIM)
      real(kind=8) pval_raw2(DIM), pval_raw3(DIM), pval_raw4(DIM)
      real(kind=8) pval_cen2(DIM), pval_cen3(DIM), pval_cen4(DIM)
      real(kind=8) pval_kurt(DIM), pval_skew(DIM), pval_var(DIM)

!     ***** Initializing parameters for Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_COLS
      cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
      cor_storage     = VSL_SS_MATRIX_STORAGE_FULL
      task_method     = VSL_SS_METHOD_FAST

!     ***** Generate data set using VSL GaussianMV RNG *****
      errcode = dGenerateGaussianMVData( p, n, x, a, C )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Edit task parameters for min and max computation *****
      errcode = vsldssedittask( task, VSL_SS_ED_MIN, min_estimate )
      call CheckVslError( errcode )

      errcode = vsldssedittask( task, VSL_SS_ED_MAX, max_estimate )
      call CheckVslError( errcode )

!     ***** Edit task parameters for computating of mean estimate and 2nd, 3rd
!           and 4th raw and central moments estimates *****
      errcode = vsldsseditmoments( task, mean, raw2, raw3, raw4,        &
     &                             cen2, cen3, cen4 )
      call CheckVslError( errcode )

!     ***** Edit task parameters for kurtosis, skewness and variation
!           computation *****
      errcode = vsldssedittask( task, VSL_SS_ED_KURTOSIS, kurtosis )
      call CheckVslError( errcode )

      errcode = vsldssedittask( task, VSL_SS_ED_SKEWNESS, skewness )
      call CheckVslError( errcode )

      errcode = vsldssedittask( task, VSL_SS_ED_VARIATION,              &
     &                          variation )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters using FULL_STORAGE
!           for covariance/correlation matrices computation *****
      errcode = vsldsseditcovcor( task, mean, cov, cov_storage,         &
     &                            cor, cor_storage )
      call CheckVslError( errcode )

!     ***** Minimum and maximum are included in the list of estimates
!           to compute *****
      estimate = ior( VSL_SS_MIN, VSL_SS_MAX )

!     ***** Mean and 2nd, 3rd and 4th raw and central moments are included
!            in the list of estimates to compute *****

      estimate = ior( estimate, ior( VSL_SS_MEAN,                       &
     &           ior( VSL_SS_2R_MOM, ior( VSL_SS_3R_MOM,                &
     &           ior( VSL_SS_4R_MOM, ior( VSL_SS_2C_MOM,                &
     &           ior( VSL_SS_3C_MOM, VSL_SS_4C_MOM ) ) ) ) ) ) )

!     ***** Kurtosis, skewness and variation are included in the list
!           of estimates to compute *****
      estimate = ior( estimate, ior( VSL_SS_KURTOSIS,                   &
     &           ior( VSL_SS_SKEWNESS, VSL_SS_VARIATION ) ) )

!     ***** Covariance and correlation matrices are included in the list
!           of estimates to compute *****
      estimate = ior( estimate, ior( VSL_SS_COV, VSL_SS_COR ) )

!     ***** Compute the estimates using FAST method *****
      errcode = vsldsscompute( task, estimate, task_method )
      call CheckVslError( errcode )

!     ***** Testing stat characteristics of computed estimates *****
      errnums = 0

!     Comparison of observations with min and max estimates
      do i = 1, p
        do j = 1, n
          if (x(i, j) < min_estimate(i)) then
            errnums = errnums + 1
          end if
          if (x(i, j) > max_estimate(i)) then
            errnums = errnums + 1
          end if
        end do
      end do

!     Compute p-values for mean estimates
      call dComputePvalsMean( p, n, mean, a, C, pval_mean )
!     Compute p-values for variance estimates
      call dComputePvalsVariance( p, n, cov, C, pval_cov )
!     Compute p-values for covariance estimates
      call dComputePvalsCovariance( p, n, cov, C, pval_cov )
!     Compute p-values for raw moments estimates
      order = 2
      call dComputePvalsRawMoments( p, n, raw2, order, a, C, pval_raw2 )
      order = 3
      call dComputePvalsRawMoments( p, n, raw3, order, a, C, pval_raw3 )
      order = 4
      call dComputePvalsRawMoments( p, n, raw4, order, a, C, pval_raw4 )
!     Compute p-values for central moments estimates
      order = 2
      call dComputePvalsCentralMoments( p, n, cen2, order, a, C,        &
     &                                  pval_cen2 )
      order = 3
      call dComputePvalsCentralMoments( p, n, cen3, order, a, C,        &
     &                                  pval_cen3 )
      order = 4
      call dComputePvalsCentralMoments( p, n, cen4, order, a, C,        &
     &                                  pval_cen4 )
!     Compute p-values for kurtosis, skewness and variation estimates
      call dComputePvalsKurtosis( p, n, kurtosis, C, pval_kurt )
      call dComputePvalsSkewness( p, n, skewness, C, pval_skew )
      call dComputePvalsVariation( p, n, variation, a, C, pval_var )

!     ***** Checking the validity of p-values for all estimates *****
      errnums = 0
      do i = 1, p
        if (pval_mean(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_raw2(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_raw3(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_raw4(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cen2(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cen3(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cen4(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_kurt(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_skew(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_var(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if

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

!     ***** Printing computed minimum, maximum, mean and moments estimates *****
      print 10, '               '
      print 10, 'Min         Max         Mean       '
      print 10, '2nd_raw     3rd_raw     4th_raw      '
      print 10, '2nd_cen     3rd_cen     4th_cen'
      print *, ''

      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 12, min_estimate(i), ' '
        print 12, max_estimate(i), ' '
        print 12, mean(i), ' '
        print 13, raw2(i), ' ', raw3(i), ' ', raw4(i), ' '
        print 13, cen2(i), ' ', cen3(i), ' ', cen4(i), ''
        print *, ''
      end do
      print *, ''

!     ***** Printing computed kurtosis, skewness and variation estimates *****
      print *, '              Kurtosis    Skewness    Variation'
      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 13, kurtosis(i), ' ', skewness(i), ' ', variation(i), ' '
        print *, ''
      end do
      print *, ''

!     ***** Printing computed covariance and correlation matrices *****
      print 10, '  Computed covariance matrix               '
      print 10, '  Computed correlation matrix'
      print *, ''
      do i = 1, p
        write (*, 14) cov(i, 1:p)
        print 10, '   '
        write (*, 14) cor(i, 1:p)
        print *, ''
      end do
      print *, ''
      print *, ''

!     ***** Printing p-values for mean and moments estimates *****
      print *, 'P-values of the computed estimates'
      print *, ''
      print *, ''
      print 10, '               Mean       '
      print 10, '2nd_raw    3rd_raw    4th_raw    '
      print 10, '2nd_cen    3rd_cen    4th_cen'
      print *, ''

      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 12, pval_mean(i), ''
        print 13, pval_raw2(i), '', pval_raw3(i), '', pval_raw4(i), ''
        print 13, pval_cen2(i), '', pval_cen3(i), '', pval_cen4(i), ''
        print *, ''
      end do
      print *, ''

!     ***** Printing p-values for kurtosis, skewness and variation
!           estimates *****
      print *, '              Kurtosis    Skewness    Variation'
      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 13, pval_kurt(i), ' ', pval_skew(i), ' ', pval_var(i), ' '
        print *, ''
      end do
      print *, ''

!     ***** Printing p-values for covariance matrix estimate *****
      print *, 'Covariance matrix'
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
14    format(4F10.6,$)

      end
