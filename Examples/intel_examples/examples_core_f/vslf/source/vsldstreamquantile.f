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
!    Calculation of quantiles, streaming data case  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3          ! Task dimension
      integer,parameter :: NN  = 1000       ! Number of observations
      integer,parameter :: MM  = 9          ! Number of deciles

      real(kind=8),parameter :: EPS = 0.01  ! Accuracy of quantile computation

      integer indices(DIM)
      data indices / 1, 0, 0 /

      TYPE(VSL_SS_TASK) task

      integer p
      integer n
      integer m
      integer x_storage
      integer o_storage
      integer params_n
      real(kind=8) x(NN,DIM)
      real(kind=8) xT(DIM,NN)
      real(kind=8) o_quant(MM)
      real(kind=8) o_stat(NN,DIM)
      real(kind=8) eq(MM,DIM), sq(MM,DIM)
      real(kind=8) min_abs_eq, min_abs_sq
      real(kind=8) params(1)
      real(kind=8) a, sigma
      integer l1, l2
      integer i, j, k
      integer(kind=4) errcode
      integer numRight, numLeft
      integer errnums
      integer task_method
      integer(kind=8) task_params

!     ***** Initialize parameters of Summary Statistics task *****
      p           = DIM
      n           = NN
      m           = MM
      x_storage   = VSL_SS_MATRIX_STORAGE_ROWS
      o_storage   = VSL_SS_MATRIX_STORAGE_ROWS
      params_n    = 1
      params(1)   = EPS
      a           = 0.0
      sigma       = 1.0
      task_params = VSL_SS_STREAM_QUANTS
      task_method = VSL_SS_METHOD_SQUANTS_ZW
      errcode     = 0
      errnums     = 0

!     ***** Generate transposed data set using VSL Gaussian RNG
!           with mean a = 0 and stdev = 1 *****
      errcode = dGenerateGaussianData( p, n, xT, a, sigma )
      call CheckVslError( errcode )

      do j = 1, p
        do i = 1, n
          x(i, j) = xT(j, i)
        end do
      end do

      do i = 1, m
        o_quant(i) = i / ( m + 1.0 )
      end do

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x,                &
     &                         indices = indices )
      call CheckVslError( errcode )

      errcode = vsldsseditstreamquantiles( task, m, o_quant, sq,        &
     &                                     params_n, params )
      call CheckVslError( errcode )

!     ***** Compute streaming quantiles with accuracy EPS *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Edit task parameters for standard deciles and order statistics
!           computation *****
      errcode = vsldsseditquantiles( task, m, o_quant, eq, o_stat,      &
     &                               o_storage )
      call CheckVslError( errcode )

!     ***** Compute standard quantiles and order statistics using
!           FAST method *****
      task_params = ior( VSL_SS_QUANTS, VSL_SS_ORDER_STATS )
      task_method = VSL_SS_METHOD_FAST
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Check the correctness of computed streaming quantiles *****
      errnums = 0

      do i = 1, p
        if ( indices(i) == 1 ) then
          do j = 1, m
            min_abs_eq = abs( o_stat(1, i) - eq(m, i) )
            min_abs_sq = abs( o_stat(1, i) - sq(m, i) )
            l1 = 1
            l2 = 1
            do k = 2, n
              if ( abs( o_stat(k, i) - eq(m, i) ) < min_abs_eq ) then
                min_abs_eq = abs( o_stat(k, i) - eq(m, i) )
                l1 = k
              end if

              if ( abs( o_stat(k, i) - sq(m, i) ) < min_abs_sq ) then
                min_abs_sq = abs( o_stat(k, i) - sq(m, i) )
                l2 = k
              end if
            end do

            if ( abs( l1 - l2 ) >= params(1) * n ) then
              errnums = errnums + 1
            end if
          end do
        end if
      end do

!     ***** Printing results *****
      print *, 'Task dimension :         ', p
      print *, 'Number of observations : ', n
      print *, ''

!     ***** Printing part of the initial matrix of observations *****
      print *,'1st 4 and last 4 observations in source matrix'
      do j = 1, p
        do i = 1, 4
          write(*, 5) x(i, j), ' '
        end do
        print 6, '     ...      '
        do i = n - 4, n
          write(*, 5) x(i, j), ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Printing computed standard quantiles *****
      print *, 'Standard deciles of the observations:'
      print 6, '             '
      do i = 1, m
        print 7, '  D', i, '   '
      end do

      print *, ''

      k = 1
      do j = 1, p
        if ( indices(j) == 1 ) then
          print 7, 'Variable #', j, ' '
          do i = 1, m
            write(*, 5) eq(i, k), ' '
          end do
          print *, ''
          k = k + 1
        end if
      end do
      print *, ''
      print *, ''

!     ***** Printing computed streaming quantiles *****
      print *, 'Streaming deciles of the observations:'
      print 6, '             '
      do i = 1, m
        print 7, '  D', i, '   '
      end do

      print *, ''

      k = 1
      do j = 1, p
        if ( indices(j) == 1 ) then
          print 7, 'Variable #', j, ' '
          do i = 1, m
            write(*, 5) sq(i, k), ' '
          end do
          print *, ''
          k = k + 1
        end if
      end do
      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if( errnums == 0 ) then
        print *, ' Computed streaming quantiles agree with theory'
      else
        print *, ' Error: Computed streaming quantiles',                &
     &           ' disagree with theory'
        stop 1
      end if


!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

5     format (F6.3,A,$)
6     format (A,$)
7     format (A,I1,A,$)

      end
