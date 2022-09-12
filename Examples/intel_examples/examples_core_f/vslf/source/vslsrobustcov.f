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
!    Computation of robust covariance  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 5        ! Task dimension
      integer,parameter :: NN  = 5000     ! Number of observation

      real(kind=4),parameter :: P_THRESHOLD = 0.001

      integer,parameter      :: RATIO = 2     ! Ratio of outliers in the dataset
      real(kind=4),parameter :: M     = 100.0 ! Mean of the outliers
      real(kind=4),parameter :: COEFF = 1.0   ! Coefficient to compute covarince
                                              ! of outliers

!     ***** Robust method parameters *****
      real(kind=4),parameter :: BD_POINT = 0.4
      real(kind=4),parameter :: ARP      = 0.001
      real(kind=4),parameter :: ACCURACY = 0.001
      integer,parameter      :: ITER_NUM = 5

!     ***** Parameters for major distribution *****
      real(kind=4) C(DIM*DIM)
      data C  / 1.0, 0.1, 0.1, 0.1, 0.1,                                &
     &          0.1, 2.0, 0.1, 0.1, 0.1,                                &
     &          0.1, 0.1, 1.0, 0.1, 0.1,                                &
     &          0.1, 0.1, 0.1, 2.0, 0.1,                                &
     &          0.1, 0.1, 0.1, 0.1, 1.0 /

      real(kind=4) a(DIM)
      data a  / 0.0, 0.0, 0.0, 0.0, 0.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer cov_storage
      integer rcov_storage
      real(kind=4) x(NN,DIM)
      real(kind=4) mean(DIM), rmean(DIM)
      real(kind=4) cov(DIM*(DIM+1)/2), rcov(DIM*(DIM+1)/2)
      real(kind=4) cov_buf(DIM,DIM)
      real(kind=4) pval_c(DIM,DIM), pval_r(DIM,DIM)
      integer i, j, k, l, k1, k2
      integer(kind=4) errcode
      integer errnums
      integer task_method
      integer(kind=8) task_params

      integer robust_params_n
      real(kind=4) robust_method_params(VSL_SS_TBS_PARAMS_N)

!     ***** Initializing parameters for Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_ROWS
      cov_storage     = VSL_SS_MATRIX_STORAGE_L_PACKED
      rcov_storage    = VSL_SS_MATRIX_STORAGE_L_PACKED
      robust_params_n = VSL_SS_TBS_PARAMS_N
      task_params     = ior( VSL_SS_COV, VSL_SS_ROBUST_COV )
      task_method     = ior( VSL_SS_METHOD_FAST, VSL_SS_METHOD_TBS )
      errcode         = 0
      errnums         = 0

!     ***** Generate dataset *****
      errcode = sGenerateContaminatedDataset( p, n, x, a, C,            &
     &                                        RATIO, M, COEFF )

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters using L_PACKED_STORAGE
!           for covariance matrix*****
      errcode = vslssseditcovcor( task, mean, cov, cov_storage )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters
!           for robust covariance estimator *****
      robust_method_params(1) = BD_POINT
      robust_method_params(2) = ARP
      robust_method_params(3) = ACCURACY
      robust_method_params(4) = ITER_NUM

      errcode = vslssseditrobustcovariance( task, rcov_storage,         &
     &      robust_params_n, robust_method_params, rmean, rcov )
      call CheckVslError( errcode )

!     ***** Compute covariance and robust covariance matrices using FAST
!           method  *****
      errcode = vslssscompute( task, task_params , task_method )
      call CheckVslError( errcode )


!     ***** Printing results *****
      print 10, 'Task dimension :         ', p
      print 10, 'Number of observations : ', n
      print *, ''

      print *, 'Exact covariance matrix:'
      do i = 1, p
        do j = 1, p
          print 12, C((i - 1)*p + j)
        end do
        print *, ''
      end do
      print *, ''

      print *, 'Exact means vector:'
      write (*, 11) a
      print *, ''

      print *, 'Classical covariance estimate:'
      k2 = 0
      do i = 1, p
        k1 = k2 + 1
        k2 = k1 + i - 1
        write (*, 11) cov(k1:k2)
      end do
      print *, ''

      print *, 'Classical mean estimate:'
      write (*, 11) mean
      print *, ''

      print *, 'Robust covariance estimate:'
      k2 = 0
      do i = 1, p
          k1 = k2 + 1
          k2 = k1 + i - 1
          write (*, 11) rcov(k1:k2)
      end do
      print *, ''

      print *, 'Robust mean estimate:'
      write (*, 11) rmean
      print *, ''

!     ***** Testing stat characteristics of classic and robust
!           covariance matrices *****
      k = 1
      do i = 1, p
        do j = 1, i - 1
          cov_buf(i, j) = cov(k)
          cov_buf(j, i) = cov(k)
          k = k + 1
        end do

        cov_buf(i, i) = cov(k)
        k = k + 1
      end do
      call sComputePvalsVariance( p, n, cov_buf, C, pval_c )
      call sComputePvalsCovariance( p, n, cov_buf, C, pval_c )

      k = 1
      do i = 1, p
        do j = 1, i - 1
          cov_buf(i, j) = rcov(k)
          cov_buf(j, i) = rcov(k)
          k = k + 1
        end do

        cov_buf(i, i) = rcov(k)
        k = k + 1
      end do
      call sComputePvalsVariance( p, n, cov_buf, C, pval_r )
      call sComputePvalsCovariance( p, n, cov_buf, C, pval_r )

      errnums = 0
      do i = 1, p
        do j = 1, i
          if (pval_r(i, j) < P_THRESHOLD) then
            errnums = errnums + 1
          end if
        end do
      end do

      print *, 'P-values of the computed classic covariance matrix'
      do i = 1, p
        write (*, 11) pval_c(i, 1:i)
      end do
      print *, ''

      print *, 'P-values of the computed robust covariance matrix'
      do i = 1, p
        write (*, 11) pval_r(i, 1:i)
      end do
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Robust covariance estimate agrees with theory'
      else
        print *, ' Error: Robust covariance estimate',                  &
     &           ' disagrees with theory'
        stop 1
      end if

      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

10    format(A,I5)
11    format(5F12.6)
12    format(F12.6,$)

      end
