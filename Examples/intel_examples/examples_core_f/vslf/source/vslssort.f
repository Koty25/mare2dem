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
!    Sorting data array  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3      ! Task dimension
      integer,parameter :: NN  = 10     ! Number of observations

      TYPE(VSL_STREAM_STATE) :: stream
      integer,parameter :: SEED   = 1
      integer,parameter :: BRNG   = VSL_BRNG_MCG31
      integer,parameter :: METHOD = VSL_RNG_METHOD_UNIFORM_STD

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      integer y_storage
      real(kind=4) x(DIM, NN)
      real(kind=4) y(DIM, NN)
      real(kind=4) lBound, rBound
      integer i, j
      integer(kind=4) errcode
      integer errnums
      integer task_method
      integer edit_value
      integer(kind=8) task_params

!     ***** Initialize parameters *****
      p           = DIM
      n           = NN
      x_storage   = VSL_SS_MATRIX_STORAGE_COLS
      y_storage   = VSL_SS_MATRIX_STORAGE_COLS
      lBound      = 0.0
      rBound      = 10.0
      task_params = VSL_SS_SORTED_OBSERV
      task_method = VSL_SS_METHOD_RADIX
      edit_value  = VSL_SS_ED_SORTED_OBSERV_STORAGE
      errcode     = 0
      errnums     = 0

!     ***** Generate data set using Uniforn RNG *****
!     ***** Initialize *****
      errcode = vslnewstream( stream, BRNG,  SEED )
      call CheckVslError(errcode)

!     ***** Call RNG *****
      errcode = vsrnguniform(METHOD, stream, n * p, x, lBound, rBound)
      call CheckVslError(errcode)

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

      errcode = vslsssedittask( task, VSL_SS_ED_SORTED_OBSERV, y )
      call CheckVslError( errcode )

      errcode = vslissedittask( task, edit_value, y_storage )
      call CheckVslError( errcode )

!     ***** Sort data using RADIX method  *****
      errcode = vslssscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Check the correctness of sorting *****
      do j = 1, n - 1
        do i = 1, p-1
          if( y(i, j) >= y(i, j + 1) ) then
            errnums = errnums + 1
          end if
        end do
      end do

!     ***** Printing results *****
      print *, ''
      print *, 'Task dimension :        ', p
      print *, 'Number of observations :', n
      print *, ''

!     ***** Printing of the initial matrix *****
      print *,'Initial matrix'
      do j = 1, n
        do i = 1, p
           write(*, 5) x(i,j), ' '
        end do
        print *, ''
      end do


!     ***** Printing of the sorted matrix *****
      print *,'Sorted matrix:'
      do j = 1, n
        do i = 1, p
          write(*, 5) y(i,j), ' '
        end do
        print *, ''
      end do

      print *, ''

!     ***** Printing summary of the test *****
      if( errnums == 0 ) then
        print *, 'Sorting is correct'
      else
        print *, 'Error: Sorting is incorrect'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

5     format (F6.3,A,$)

      end
