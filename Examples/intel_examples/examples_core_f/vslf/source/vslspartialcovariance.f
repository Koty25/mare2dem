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
!    Calculation of partial covariance matrix  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM      = 4        ! Task dimension
      integer,parameter :: PART_DIM = 2        ! Partial covariance dimension

      real(kind=4),parameter :: EPSILON = 1e-6

      real(kind=4) cov(DIM,DIM)
      data cov  /  1.0,  0.1, 0.15,  0.1,                               &
     &             0.1,  2.0,  0.1,  0.1,                               &
     &            0.15,  0.1,  1.0,  0.1,                               &
     &             0.1,  0.1,  0.1,  1.0 /

      integer pcov_index(DIM)
      data pcov_index  /  1,  1, -1, -1 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer dummy_n
      integer part_p
      integer dummy_x_storage
      integer cov_storage
      integer pcov_storage
      real(kind=4) dummy_x(DIM,1)
      real(kind=4) cp_cov(DIM,DIM)
      real(kind=4) pcov(PART_DIM,PART_DIM)
      real(kind=4) th_pcov(PART_DIM,PART_DIM)
      integer task_method
      integer(kind=8) task_params
      integer i, i1, j, j1
      integer(kind=4) errcode
      integer errnums

!     ***** Initializing parameters for Summary Statistics task *****
      p               = DIM
      part_p          = PART_DIM
      dummy_n         = 1
      dummy_x_storage = VSL_SS_MATRIX_STORAGE_COLS
      cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
      pcov_storage    = VSL_SS_MATRIX_STORAGE_FULL
      task_params     = VSL_SS_PARTIAL_COV
      task_method     = VSL_SS_METHOD_FAST
      errcode         = 0
      errnums         = 0

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, dummy_n, dummy_x_storage,       &
     &                         dummy_x )
      call CheckVslError( errcode )

!     ***** Edit task parameters for partial covariance matrix computation *****
      errcode = vslssseditpartialcovcor( task, pcov_index, cov,         &
     &                                   cov_storage, p_cov = pcov,     &
     &                                   p_cov_storage = pcov_storage,  &
     &                                   p_cor_storage = pcov_storage )
      call CheckVslError( errcode )

!     ***** Compute partial covariance matrix using FAST method *****
      errcode = vslssscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Printing results *****
      print 10, 'Task dimension : ', p
      print *, ''

!     Print input covariance matrix
      print *, 'Covariance matrix'
      do i = 1, p
        write (*, 11) cov(i, 1:p)
      end do
      print *, ''

!     Print computed partial covariance matrix estimate
      print *, 'Computed partial covariance matrix'
      do i = 1, part_p
        write (*, 12) pcov(i, 1:part_p)
      end do
      print *, ''

!     ***** Testing stat characteristics of partial covariance matrix *****
!     Compute theoretical partial covariance estimate using sweep operator
      do i = 1, p
        do j = 1, p
          cp_cov(i, j) = cov(i, j)
        end do
      end do

      do i = 1, p
        if (pcov_index(i) == -1) then
          call sSweep( i, p, cp_cov )
        end if
      end do

      i1 = 1
      j1 = 1
      do i = 1, p
        if (pcov_index(i) == 1) then
          j1 = 1
          do j = 1, p
            if (pcov_index(j) == 1) then
              th_pcov(i1, j1) = cp_cov(i, j)
              j1 = j1 + 1
            end if
          end do
          i1 = i1 + 1
        end if
      end do

!     Print theoretical partial covariance estimate
      print *, 'Theoretical partial covariance matrix'
      do i = 1, part_p
        write (*, 12) th_pcov(i, 1:part_p)
      end do
      print *, ''
      print *, ''

!     Check the correctness of computed partial covariance matrix
      errnums = 0
      do i = 1, part_p
        do j = 1, part_p
          if (abs(pcov(i, j) - th_pcov(i, j)) > EPSILON) then
            errnums = errnums + 1
          end if
        end do
      end do

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' Computed partial covariance matrix estimate',        &
     &           ' agrees with theory'
      else
        print *, ' Error: Computed partial covariance matrix estimate', &
     &           ' disagrees with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

10    format (A,I1)
11    format (4F9.6)
12    format (2F9.6)

      end
