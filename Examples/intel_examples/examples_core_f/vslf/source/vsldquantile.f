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
!    Calculation of quantiles and odered statistic  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3        ! Task dimension
      integer,parameter :: NN  = 1000     ! Number of observations
      integer,parameter :: MM  = 9        ! Number of deciles

      TYPE(VSL_SS_TASK) task
      integer p
      integer n, n2
      integer m
      integer x_storage
      integer o_storage
      real(kind=8) x(NN,DIM)
      real(kind=8) xT(DIM,NN)
      real(kind=8) o_stat(NN,DIM)
      real(kind=8) o_quant(MM)
      real(kind=8) quantiles(MM,DIM)
      real(kind=8) a, sigma
      integer i, j
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
      a           = 0.0
      sigma       = 1.0
      task_params = ior( VSL_SS_QUANTS, VSL_SS_ORDER_STATS )
      task_method = VSL_SS_METHOD_FAST
      errcode     = 0
      errnums     = 0

      n2 = n / 2
      if ( n > n2 * 2 ) then
        n2 = n2 + 1
      end if

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
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

      errcode = vsldsseditquantiles( task, m, o_quant, quantiles,       &
     &                               o_stat, o_storage )
      call CheckVslError( errcode )

!     ***** Compute quantiles and order statistics using FAST method *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Check the correctness of computed quantiles and order
!           statistics *****

      errnums = 0

      do j = 1, p
        do i = 1, n-1
          if( o_stat(i, j) >= o_stat(i + 1, j) ) then
            errnums = errnums + 1
          end if
        end do

        numRight = 0
        numLeft  = 0

        do i = 1, n
          if( x(i, j) >= quantiles((m / 2) + 1, j) ) then
            numRight = numRight + 1
          end if

          if( x(i, j) <= quantiles((m / 2) + 1, j) ) then
            numLeft = numLeft + 1
          end if
        end do

        if ( numRight < n2 .or. numLeft < n2 ) then
          errnums = errnums + 1
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
          write(*, 5) x(i,j), ' '
        end do
        print 6, '     ...      '
        do i = n-4, n
          write(*, 5) x(i,j), ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Printing computed quantiles *****
      print *,'Deciles of the observations for all variables:'
      do i = 1, m
        print 7, '  D', i, '   '
      end do

      print *, ''

      do j = 1, p
        do i = 1, m
          write(*, 5) quantiles(i, j), ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Printing part of the order statistics matrix *****
      print *,'1st 4 and last 4 observations in order statistics matrix'
      do j = 1, p
        do i = 1, 4
          write(*, 5) o_stat(i,j),' '
         end do
         print 6, '     ...      '
        do i = n-4, n
          write(*, 5) o_stat(i,j),' '
         end do

         print *, ''
      end do

      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if( errnums == 0 ) then
        print *, ' Computed quantiles and order statistics',            &
     &           ' agree with theory'
      else
        print *, ' Error: Computed quantiles and/or order statistics',  &
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
