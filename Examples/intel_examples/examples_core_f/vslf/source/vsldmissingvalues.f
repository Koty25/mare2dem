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
!  Support of missing values, Multiple Imputation algorithm Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM     = 4     ! Task dimension
      integer,parameter :: NN      = 1000  ! Number of observations
      integer,parameter :: NNPATT  = 9     ! Number of different patterns
      integer,parameter :: MM      = 5     ! Number of sets of imputed values
      real(kind=8),parameter :: EPSILON = 0.2   ! Missing values percentage

      real(kind=8),parameter :: dNAN    = VSL_SS_DNAN

      integer,parameter :: EM_ITER_NUM = 100
      integer,parameter :: DA_ITER_NUM = 10
      integer,parameter :: EM_ACCURACY = 0.001

!     Number of simulated missing values
      integer,parameter :: MI_SIMVALS_N = 20000
!     Initial number of estimates
      integer,parameter :: MI_INIT_EST_N = DIM + DIM * (DIM + 1) / 2
!     Number of estimates
      integer,parameter :: MI_EST_N =                                   &
     &                     MM * DA_ITER_NUM * DIM * (DIM + 3) / 2

!     Exact covariance matrix
      real(kind=8) C(DIM*DIM)
      data C  /   1.0, 0.125, 0.125, 0.125,                             &
     &          0.125,   1.0, 0.125, 0.125,                             &
     &          0.125, 0.125,   1.0, 0.125,                             &
     &          0.125, 0.125, 0.125,   1.0  /

!     Exact mean
      real(kind=8) a(DIM)
      data a  /   0.0,   0.0,   0.0,   0.0  /

!     Matrix of patterns
      integer patt(DIM*NNPATT)
      data patt  / 1, 0, 0, 0,                                          &
     &             0, 1, 0, 0,                                          &
     &             0, 0, 1, 0,                                          &
     &             0, 0, 0, 1,                                          &
     &             1, 0, 1, 0,                                          &
     &             0, 0, 1, 1,                                          &
     &             1, 1, 0, 1,                                          &
     &             1, 0, 1, 1,                                          &
     &             1, 1, 1, 1  /

      real(kind=8) W(2)
      data W  / 0.0, 0.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer m
      integer n_miss_vals
      integer index
      integer x_storage
      integer q_storage
      integer npatt
      real(kind=8) eps
      real(kind=8) x(NN,DIM)
      real(kind=8) mean(DIM,MM), raw2(DIM,MM), variance(DIM,MM)
      real(kind=8) mean_buf(DIM), raw2_buf(DIM), variance_buf(DIM)
      real(kind=8) qmean(DIM), qvariance(DIM)
      real(kind=8) B(DIM), U(DIM), T(DIM), sqrtT(DIM)
      real(kind=8) cint_left(DIM), cint_right(DIM)
      integer miss_vals(NN)
      integer task_method
      integer(kind=8)task_params
      integer i, j, k, k1, k2, ipatt
      integer(kind=4) errcode
      integer errnums

      integer nparams
      integer n_simvals, n_init_est, n_est, n_prior
      real(kind=8) params(VSL_SS_MI_PARAMS_SIZE)
      real(kind=8) simvals(MI_SIMVALS_N), init_est(MI_INIT_EST_N)
      real(kind=8) est(MI_EST_N), prior(1)

!     ***** Initializing parameters for Summary Statistics task *****
      p           = DIM
      n           = NN
      m           = MM
      x_storage   = VSL_SS_MATRIX_STORAGE_ROWS
      q_storage   = VSL_SS_MATRIX_STORAGE_COLS
      npatt       = NNPATT
      eps         = EPSILON
      nparams     = VSL_SS_MI_PARAMS_SIZE
      n_simvals   = MI_SIMVALS_N
      n_est       = MI_EST_N
      n_init_est  = MI_INIT_EST_N
      n_prior     = 0
      task_params = VSL_SS_MISSING_VALS
      task_method = VSL_SS_METHOD_MI
      errcode     = 0
      errnums     = 0

!     ***** Define covariance matrix and mean of input data *****
      k = 1
      do i = 1, p
        init_est(i) = a(i)
        do j = i, p
          init_est(p + k) = C((i - 1)*p + j)
          k = k + 1
        end do
      end do

!     ***** Generate dataset with missing values *****
      do i = 1, n
        miss_vals(i) = -1
      end do

      errcode = dGenerateMissingValuesInput( p, n, x, eps, npatt, patt, &
     &                                       a, C, miss_vals,           &
     &                                       n_miss_vals )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Initialization of the task parameters for missing values
!           estimator *****
      n_simvals = n_miss_vals * m

      params(1) = EM_ITER_NUM
      params(2) = DA_ITER_NUM
      params(3) = EM_ACCURACY
      params(4) = m
      params(5) = n_miss_vals

      errcode = vsldsseditmissingvalues( task, nparams, params,         &
     &                                   n_init_est, init_est,          &
     &                                   n_prior, prior,                &
     &                                   n_simvals, simvals,            &
     &                                   n_est, est )
      call CheckVslError( errcode )

!     ***** Compute missing values estimates using MULTIPLE_IMPUTATION
!           method *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Compute mean and variance for all imputations *****
      index = 1
      do k = 1, m
!       Recover of the input matrix
        do i = 1, n
          ipatt = miss_vals(i)
          if (ipatt > 0) then
            do j = 1, p
              if (patt((ipatt - 1) * p + j) == 1) then
                x(i, j) = simvals(index)
                index = index + 1
              end if
            end do
          end if
        end do

        W(1) = 0.0
        W(2) = 0.0

!       Initializing of the task parameters for mean, 2nd raw moment and
!       variance estimators computation
        errcode = vsldsseditmoments( task, mean_buf, raw2_buf,          &
     &                               c2m = variance_buf )
        call CheckVslError( errcode )

!       Compute mean and variance using fast method
        task_params = ior( VSL_SS_MEAN , ior( VSL_SS_2C_MOM,            &
     &                     VSL_SS_2R_MOM ) )
        task_method = VSL_SS_METHOD_FAST
        errcode = vsldsscompute( task, task_params, task_method )
        call CheckVslError( errcode )

        do i = 1, p
          mean(i, k) = mean_buf(i)
          variance(i, k) = variance_buf(i)
        end do
      end do

!     ***** Testing stat characteristics of computed estimates *****

!     Register matrix of means estimates as observations matrix in
!     Summary Statistics task
      errcode = vsldssedittask( task, VSL_SS_ED_OBSERV, mean )
      call CheckVslError( errcode )

!     Set proper number of observations
      errcode = vslissedittask( task, VSL_SS_ED_OBSERV_N, m )
      call CheckVslError( errcode )

!     Register proper observation storage format in the task
      errcode = vslissedittask( task, VSL_SS_ED_OBSERV_STORAGE,         &
     &                          q_storage )
      call CheckVslError( errcode )

!     Initializing of the task parameters for mean, 2nd raw moment and
!     variance estimators computation
      errcode = vsldsseditmoments( task, qmean, raw2, c2m = B )
      call CheckVslError( errcode )

      W(1) = 0.0
      W(2) = 0.0

!     ***** Compute mean, 2nd raw moments and variance est for
!           matrix of means est for all datasets *****
      task_params = ior( VSL_SS_MEAN, VSL_SS_2C_MOM )
      task_method = VSL_SS_METHOD_FAST
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     Calculate average of variance esimates
      do j = 1, p
        qvariance(j) = 0.0
        do i = 1, m
          qvariance(j) = qvariance(j) + variance(j, i)
        end do
        qvariance(j) = qvariance(j) / m
      end do

!     Calculate "within imputation" variance
      do j = 1, p
        U(j) = qvariance(j) / n
      end do

!     Calculate "total" variance
      do j = 1, p
        T(j) = U(j) + (1.0 + 1.0 / m) * B(j)
      end do

!     Compute borders of 95% confidence interval for mean estimates
      call vdsqrt( p, T, sqrtT )
      errnums = 0
      do j = 1, p
        cint_right(j) = qmean(j) + 2*sqrtT(j)
        cint_left(j)  = qmean(j) - 2*sqrtT(j)

        if ( cint_right(j) < a(j) .or. cint_left(j) > a(j) ) then
          errnums = errnums + 1
        end if
      end do

!     ***** Printing results *****
      print 10, 'Task dimension :         ', p
      print 10, 'Number of observations : ', n
      print *, ''
      print 10, 'Number of initial estimates :        ', n_init_est
      print 10, 'Number of simulated missing values : ', n_simvals
      print 10, 'Number of estimates :                ', n_est
      print 10, 'Number of missing values :           ', n_miss_vals
      print *, ''

!     ***** Print exact mean and variance *****
      print *, 'Exact mean'
      write (*, 11) a
      print *, ''

      print *, 'Exact variance'
      do i = 1, p
        print 12, C((i - 1)*p + i), ''
      end do
      print *, ''
      print *, ''

!     ***** Print computes means and variances for all sets
!           of imputed values *****
      print *,  ' Set:     Mean:'
      do i = 1, m
        print 13, '   ', i, '    '
        write (*, 11) mean(1:p, i)
      end do
      print *, ''

      print 14, 'Average '
      write (*, 11) qmean
      print *, ''
      print *, ''

      print *,  ' Set:     Variance:'
      do i = 1, m
        print 13, '   ', i, '    '
        write (*, 11) variance(1:p, i)
      end do
      print *, ''

      print 14, 'Average '
      write (*, 11) qvariance
      print *, ''

!     ***** Print between-imputation, within-imputation and total
!           variances *****
      print *, 'Between-imputation variance:'
      print *, ''
      write (*, 11) B
      print *, ''

      print *, 'Within-imputation variance:'
      print *, ''
      write (*, 11) U
      print *, ''

      print *, 'Total variance:'
      print *, ''
      write (*, 11) T
      print *, ''

!     ***** Print borders of 95% confidence intervals for mean estimates *****
      print *, '95%% confidence interval:'
      print *, ''
      print 14, ' right    '
      write (*, 11) cint_right
      print 14, ' left     '
      write (*, 11) cint_left
      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Computed missing values estimates',                  &
     &           ' agree with theory'
      else
        print *, ' Error: Computed missing values estimates',           &
     &           ' disagree with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

10    format(A,I4)
11    format(4F10.6)
12    format(F10.6,A,$)
13    format(A,I1,A,$)
14    format(A,$)

      end
