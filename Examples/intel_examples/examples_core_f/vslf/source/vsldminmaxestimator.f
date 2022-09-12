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
!    Calculation of min/max estimators  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3       ! Task dimension
      integer,parameter :: NN  = 1000    ! Number of observations

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      real(kind=8) x(DIM,NN)
      real(kind=8) min_est(DIM), max_est(DIM)
      real(kind=8) a, sigma
      integer i, j
      integer(kind=4) errcode
      integer task_method
      integer(kind=8) task_params
      integer errnums

!     ***** Initialize parameters of Summary Statistics task *****
      p           = DIM
      n           = NN
      x_storage   = VSL_SS_MATRIX_STORAGE_COLS
      a           = 0.0
      sigma       = 1.0
      task_params = IOR( VSL_SS_MAX, VSL_SS_MIN )
      task_method = VSL_SS_METHOD_FAST
      errcode     = 0
      errnums     = 0

!     ***** Generate data set using VSL Gaussian RNG
!           with a = 0 and sigma = 1 *****
      errcode = dGenerateGaussianData( p, n, x, a, sigma )
      call CheckVslError( errcode )

!     ***** Set initial values of the estimates *****
      do i = 1, p
        min_est(i) = x(i, 1)
        max_est(i) = x(i, 1)
      end do

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Edit task parameters for min and max computation *****
      errcode = vsldssedittask( task, VSL_SS_ED_MIN, min_est )
      call CheckVslError( errcode )

      errcode = vsldssedittask( task, VSL_SS_ED_MAX, max_est )
      call CheckVslError( errcode )

!     ***** Compute min and max estimates using FAST method *****
      errcode = vsldsscompute( task, task_params , task_method )
      call CheckVslError( errcode )

!     ***** Comparison of observations with min and max estimates *****
      do i = 1, p
         do j = 1, n
            if ( x(i, j) < min_est(i) ) then
               errnums = errnums + 1
            end if
            if ( x(i, j) > max_est(i) ) then
               errnums = errnums + 1
            end if
         end do
      end do

!     ***** Printing results *****
      print 9, 'Task dimension :         ', p
      print 9, 'Number of observations : ', n
      print *, ''
      print *,'                 Min        Max'
      do i = 1, p
         print 10, 'Variable #', i ,':  ', min_est(i), max_est(i)
      end do

      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
         print *, ' All observations are',                              &
     & ' within ranges for all dimensions'
      else
         print 11, ' Error: There are ', errnums,                       &
     & ' observations beyond the ranges'
         stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

9     format (A,I4)
10    format (A,I2,A,F10.6,F10.6)
11    format (A,I4,A)

      end
