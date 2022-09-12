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
!    Calculation of group/poopled covariances  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3       ! Task dimension
      integer,parameter :: NN  = 10000   ! Number of observations
      integer,parameter :: GG  = 2       ! Number of groups
      integer,parameter :: GGN = 2       ! Number of group covariance matrices

      real(kind=8),parameter :: P_THRESHOLD = 0.005

      real(kind=8) C(DIM*DIM)            ! Exact covariance matrix
      data C  / 1.0, 0.0, 0.0,                                          &
     &          0.0, 1.0, 0.0,                                          &
     &          0.0, 0.0, 1.0 /

      real(kind=8) m(DIM)                ! Exact means vector
      data m  / 0.0, 0.0, 0.0 /

      integer group_matrix_indices(GG)
      data group_matrix_indices / 1, 1 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer cov_storage
      integer pld_cov_storage
      integer grp_cov_storage
      real(kind=8) x(NN,DIM)
      real(kind=8) mean(DIM), pld_mean(DIM)
      real(kind=8) cov(DIM*DIM), pld_cov(DIM*DIM)
      real(kind=8),target :: grp_mean(DIM*GGN)
      real(kind=8),target :: grp_cov(DIM*DIM*GGN)
      real(kind=8),pointer :: grp_mean_ptr(:), grp_cov_ptr(:)
      real(kind=8) a, sigma
      integer group_indices(NN)
      integer i, j
      integer(kind=4) errcode
      integer errnums
      integer x_len
      integer task_method
      integer(kind=8) task_params

      real(kind=8) pval_pld_mean(DIM)
      real(kind=8) pval_pld_cov(DIM*DIM)
      real(kind=8),target :: pval_grp_mean(DIM*GGN)
      real(kind=8),target :: pval_grp_cov(DIM*DIM*GGN)
      real(kind=8),pointer :: pval_grp_mean_ptr(:), pval_grp_cov_ptr(:)

!     ***** Initialize parameters of Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_ROWS
      cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
      pld_cov_storage = VSL_SS_MATRIX_STORAGE_FULL
      grp_cov_storage = VSL_SS_MATRIX_STORAGE_FULL
      a               = 0.0
      sigma           = 1.0
      task_params     = ior( ior( VSL_SS_COV, VSL_SS_POOLED_COV ),      &
     &                       VSL_SS_GROUP_COV )
      task_method     = VSL_SS_METHOD_1PASS
      errcode         = 0

!     ***** Generate data set using VSL Gaussian RNG
!           with a = 0 and sigma = 1 *****
      errcode = dGenerateGaussianData(  p, n, x, a, sigma )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode  = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters for pooled covariance
!           estimator *****
      errcode = vslissedittask( task, VSL_SS_ED_POOLED_COV_STORAGE,     &
     &                          pld_cov_storage )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters for group covariance
!           estimator *****
      status = vslissedittask( task, VSL_SS_ED_GROUP_COV_STORAGE,       &
     &                         grp_cov_storage )
      call CheckVslError( errcode )

      errcode = vsldsseditcovcor( task, mean, cov, cov_storage )
      call CheckVslError( errcode )

!     Dividing elements into odd and even
      do i = 1, n
        group_indices(i) = mod(i - 1, 2)
      end do

      errcode = vsldsseditpooledcovariance( task, group_indices,        &
     &      pld_mean, pld_cov, group_matrix_indices, grp_mean, grp_cov )
      call CheckVslError( errcode )

!     ***** Compute covariance matrices using FAST method  *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Testing stat characteristics of mean and covariance matrices *****
!     Compute p-values for group mean estimates
      call dComputePvalsMean( p, n, grp_mean, m, C, pval_grp_mean )
      grp_mean_ptr      => grp_mean( p+1 : 2*p )
      pval_grp_mean_ptr => pval_grp_mean( p+1 : 2*p )
      call dComputePvalsMean( p, n, grp_mean_ptr, m, C,                 &
     &                        pval_grp_mean_ptr )
!     Compute p-values for group variance estimates
      call dComputePvalsVariance( p, n, grp_cov, C, pval_grp_cov )
      grp_cov_ptr       => grp_cov( p*p + 1 : GGN*p*p )
      pval_grp_cov_ptr  => pval_grp_cov( p*p + 1 : GGN*p*p )
      call dComputePvalsVariance( p, n, grp_cov_ptr, C,                 &
     &                            pval_grp_cov_ptr )
!     Compute p-values for group variance estimates
      call dComputePvalsCovariance( p, n, grp_cov, C, pval_grp_cov )
      call dComputePvalsCovariance( p, n, grp_cov_ptr, C,               &
     &                              pval_grp_cov_ptr )
!     Compute p-values for pooled mean estimates
      call dComputePvalsMean( p, n, pld_mean, m, C, pval_pld_mean )
!     Compute p-values for pooled variance estimates
      call dComputePvalsVariance( p, n, pld_cov, C, pval_pld_cov )
!     Compute p-values for pooled covariance estimates
      call dComputePvalsCovariance( p, n, pld_cov, C, pval_pld_cov )

      errnums = 0
      do i = 1, p
        if( pval_grp_mean(i) < P_THRESHOLD ) then
          errnums = errnums + 1
        end if
        if( pval_grp_mean(i + p) < P_THRESHOLD ) then
          errnums = errnums + 1
        end if
        if( pval_pld_mean(i) < P_THRESHOLD ) then
          errnums = errnums + 1
        end if
        do j = 1, p
          if( pval_grp_cov((i - 1)*p + j) < P_THRESHOLD ) then
            errnums = errnums + 1
          end if
          if( pval_grp_cov(p*p + (i - 1)*p + j) < P_THRESHOLD ) then
            errnums = errnums + 1
          end if
          if( pval_pld_cov((i - 1)*p + j) < P_THRESHOLD ) then
            errnums = errnums + 1
          end if
        end do
      end do

!     ***** Printing results *****
      print *, 'Task dimension :         ', p
      print *, 'Number of observations : ', n
      print *, ''

!     ***** Print exact covariance matrix and mean *****
      print *, 'Exact covariance matrix:'
      do i = 1, p
        do j = 1, p
          print 5, C((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'Exact means:'
      do i = 1, p
        print 5, m(i), ' '
      end do
      print *, ''
      print *, ''

!     ***** Print computed covariance matrix and mean estimates *****
      print *, 'Computed covariance matrix:'
      do i = 1, p
        do j = 1, p
          print 5, cov((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'Computed means:'
      do i = 1, p
        print 5, mean(i), ' '
      end do
      print *, ''
      print *, ''

!     ***** Print group covariance matrices and mean estimates *****
      print *, 'Group covariance matrices:'
      do i = 1, p
        do j = 1, p
          print 5, grp_cov((i - 1)*p + j), ' '
        end do

        print 6, '     '

        do j = 1, p
          print 5, grp_cov_ptr((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'Group means:'
      do i = 1, p
        print 5, grp_mean(i), ' '
      end do

      print 6, '     '

      do i = 1, p
        print 5, grp_mean_ptr(i), ' '
      end do
      print *, ''
      print *, ''

!     ***** Print pooled covariance matrix and mean estimates *****
      print *,'Pooled covariance matrix:'
      do i = 1, p
        do j = 1, p
            print 5, pld_cov((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'Pooled means:'
      do i = 1, p
        print 5, pld_mean(i), ' '
      end do
      print *, ''
      print *, ''
      print *, ''

!     ***** Print P-values of the group covariance matrices *****
      print *, 'P-values of the computed group covariance matrices:'
      do i = 1, p
        do j = 1, p
          print 5, pval_grp_cov((i - 1)*p + j), ' '
        end do

        print 6, '     '

        do j = 1, p
          print 5, pval_grp_cov_ptr((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'P-values of the computed group mean:'
      do i = 1, p
        print 5, pval_grp_mean(i), ' '
      end do

      print 6, '     '

      do i = 1, p
        print 5, pval_grp_mean_ptr(i), ' '
      end do
      print *, ''
      print *, ''

!     ***** Printing P-values of the pooled covariance matrix *****
      print *,'P-values of the computed pooled covariance matrix:'
      do i = 1, p
        do j = 1, p
            print 5, pval_pld_cov((i - 1)*p + j), ' '
        end do
        print *, ''
      end do
      print *, ''

      print *, 'P-values of the computed pooled mean:'
      do i = 1, p
        print 5, pval_pld_mean(i), ' '
      end do
      print *, ''
      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Pooled and group covariance matrices',               &
     &           ' and mean estimates agree with theory'
      else
        print *, ' Error: Pooled and group covariance matrices',        &
     &           ' and mean estimates disagree with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

5     format (F9.6,A,$)
6     format (A,$)

      end
