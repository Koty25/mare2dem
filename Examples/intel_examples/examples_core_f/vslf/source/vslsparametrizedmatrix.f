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
!  Parametrization of correlation matrix  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=8),parameter :: DIM    = 3      ! Task dimension
      integer(kind=8),parameter :: LLWORK = 2*25*DIM

      real(kind=4),parameter :: TEST_THRESHOLD = -1e-6

!     Distorted correlation matrix
      real(kind=4) cor(DIM*DIM)
      data cor /  1.0, 0.95,  0.7,                                      &
     &           0.95,  1.0, 0.29,                                      &
     &            0.7, 0.29,  1.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer dummy_n
      integer dummy_x_storage
      integer cor_storage
      integer pcor_storage
      real(kind=4) dummy_x(1)
      real(kind=4) pcor(DIM*DIM)
      real(kind=4) copy_cor(DIM*DIM)
      integer i, k1, k2
      integer(kind=4) errcode
      integer errnums
      integer(kind=8) task_params
      integer task_method

!     ***** Following variables are used in SSYEV routine which finds
!           eigenvalues of simmetric matrix *****
      real(kind=4) eigenvals(DIM), work(LLWORK)
      integer(kind=8) lwork
      integer info
      character jobz, uplo



!     ***** Initialize parameters of Summary Statistics task *****
      p               = DIM
      dummy_n         = 1
      dummy_x_storage = VSL_SS_MATRIX_STORAGE_COLS
      cor_storage     = VSL_SS_MATRIX_STORAGE_FULL
      pcor_storage    = VSL_SS_MATRIX_STORAGE_FULL
      task_params     = VSL_SS_PARAMTR_COR
      task_method     = VSL_SS_METHOD_SD
      errcode         = 0
      errnums         = 0

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, dummy_n, dummy_x_storage,       &
     &                         dummy_x )
      call CheckVslError( errcode )

!     ***** Edit task parameters for parameterization of correlation *****
      errcode = vslssseditcorparameterization( task, cor,               &
     &                           cor_storage, pcor, pcor_storage )
      call CheckVslError( errcode )

!     ***** Parametrize correlation *****
      errcode = vslssscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Compute eigenvalues of distorted correlation matrix *****
      copy_cor = cor

      lwork = LLWORK
      uplo = 'U'
      jobz = 'N'
      call ssyev( jobz, uplo, p, copy_cor, p, eigenvals, work, lwork,   &
     &            info )
      call CheckVslError(int(info,kind=4))

!     ***** Printing results *****
      print 9, 'Task dimension :     ', p
      print *, ''

!     ***** Print distorted correlation matrix and it's eigenvalues *****
      print *, 'Distorted correlation matrix:'
      k1 = 1
      k2 = p
      do i = 1, p
        write (*, 10) cor(k1:k2)
        k1 = k1 + p
        k2 = k2 + p
      end do
      print *, ''

      print *, 'Eigenvalues of the distorted correlation matrix'
      write (*, 10) eigenvals
      print *, ''
      print *, ''

!     ***** Compute eigenvalues of parameterized correlation matrix *****
      copy_cor = pcor

      uplo = 'U'
      jobz = 'N'
      call ssyev( jobz, uplo, p, copy_cor, p, eigenvals, work, lwork,   &
     &            info )
      call CheckVslError(int(info,kind=4))

!     ***** Print parameterized correlation matrix and it's eigenvalues *****
      print *, 'Parameterized correlation matrix:'
      k1 = 1
      k2 = p
      do i = 1, p
        write (*, 10) pcor(k1:k2)
        k1 = k1 + p
        k2 = k2 + p
      end do
      print *, ''

      print *, 'Eigenvalues of the parameterized correlation matrix'
      write (*, 10) eigenvals
      print *, ''
      print *, ''

      errnums = 0
      do i = 1, p
        if (eigenvals(i) < TEST_THRESHOLD) then
          errnums = errnums + 1
        end if
      end do

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' All eigenvalues of parametrized correlation are',    &
     &           ' non negative'
      else
        print 11, ' Error: Parameterized correlation matrix has',       &
     &            errnums, ' negative eigenvalue(s)'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

!     ***** Print parametrized correlation matrix and it's eigenvalues *****

9     format(A,I2)
10    format(3F8.5)
11    format(A,I2,A)

      end
