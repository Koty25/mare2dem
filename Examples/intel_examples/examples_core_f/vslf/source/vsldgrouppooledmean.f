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
!    Calculation of group/poopled means  Example Program Text
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
      integer,parameter :: GGN = 2       ! Number of group means

      real(kind=8),parameter :: P_THRESHOLD = 0.005

      real(kind=8) C(DIM*DIM)            ! Exact covariance matrix
      data C  / 1.0, 0.0, 0.0,                                          &
     &          0.0, 1.0, 0.0,                                          &
     &          0.0, 0.0, 1.0 /

      real(kind=8) m(DIM)                ! Exact means vector
      data m  / 0.0, 0.0, 0.0 /

      integer group_mean_indices(GG)
      data group_mean_indices / 1, 1 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer pld_cov_storage
      integer grp_cov_storage
      real(kind=8) x(NN,DIM)
      real(kind=8) pld_mean(DIM)
      real(kind=8),target :: grp_mean(DIM*GGN)
      real(kind=8),pointer :: grp_mean_ptr(:)
      real(kind=8) a, sigma
      integer group_indices(NN)
      integer i, j
      integer(kind=4) errcode
      integer errnums
      integer x_len
      integer task_method
      integer(kind=8) task_params

      real(kind=8) pval_pld_mean(DIM)
      real(kind=8),target :: pval_grp_mean(DIM*GGN)
      real(kind=8),pointer :: pval_grp_mean_ptr(:)

!     ***** Generate data set using VSL Gaussian RNG
!           with a = 0 and sigma = 1 *****

      p               = DIM
      n               = NN
      a               = 0.0
      sigma           = 1.0
      errcode = dGenerateGaussianData(  p, n, x, a, sigma )
      call CheckVslError( errcode )

!     Dividing elements into odd and even
      do i = 1, n
        group_indices(i) = mod(i - 1, 2)
      end do

!     ***** Initialize parameters of Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_ROWS
      task_params     = ior( VSL_SS_GROUP_MEAN, VSL_SS_POOLED_MEAN )
      task_method     = VSL_SS_METHOD_1PASS
      errcode         = 0

!     ***** Create Summary Statistics task *****
      errcode  = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

      errcode = vsldsseditpooledcovariance( task, group_indices,        &
     &      pld_mean, grp_cov_indices=group_mean_indices,               &
     &      grp_means=grp_mean )
      call CheckVslError( errcode )

!     ***** Compute means using 1PASS method  *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Testing stat characteristics of means *****
!     Compute p-values for group mean estimates
      call dComputePvalsMean( p, n, grp_mean, m, C, pval_grp_mean )
      grp_mean_ptr      => grp_mean( p+1 : 2*p )
      pval_grp_mean_ptr => pval_grp_mean( p+1 : 2*p )
      call dComputePvalsMean( p, n, grp_mean_ptr, m, C,                 &
     &                        pval_grp_mean_ptr )

!     Compute p-values for pooled mean estimates
      call dComputePvalsMean( p, n, pld_mean, m, C, pval_pld_mean )

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
      end do

!     ***** Printing results *****
      print *, 'Task dimension :         ', p
      print *, 'Number of observations : ', n
      print *, ''

!     ***** Print exact mean *****
      print *, 'Exact means:'
      do i = 1, p
        print 5, m(i), ' '
      end do
      print *, ''
      print *, ''

!     ***** Print group mean estimates *****
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

!     ***** Print pooled mean estimates *****
      print *, 'Pooled means:'
      do i = 1, p
        print 5, pld_mean(i), ' '
      end do
      print *, ''
      print *, ''
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

      print *, 'P-values of the computed pooled mean:'
      do i = 1, p
        print 5, pval_pld_mean(i), ' '
      end do
      print *, ''
      print *, ''
      print *, ''

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Pooled and group mean estimates',                    &
     &           ' agree with theory'
      else
        print *, ' Pooled and group mean estimates',                    &
     &           ' disagree with theory'
        stop 1
      end if

5     format (F9.6,A,$)
6     format (A,$)

      end
