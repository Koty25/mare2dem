!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) LAPACK95 TPPACK Example
!*******************************************************************************
      PROGRAM LA_TPPACK_EXAMPLE

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: MKL_TPPACK
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), AP(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'TPPACK Example Program Results.'
      READ (*,*) N
      ALLOCATE( A(N,N), AP(N,N) )
      A = 0
      AP =0

      DO J=1,N
      DO I=1,J
         READ(*,'(F2.0)') A(I,J)
      ENDDO
      ENDDO

      WRITE(*,*) 'The matrix A:'
      DO I=1,N; WRITE(*,"(4(I3,1X),I3,1X)") INT(A(I,:)); ENDDO

      WRITE(*,*) 'CALL MKL_TPPACK( AP, A )'

      CALL MKL_TPPACK( AP, 1, 1, N, N, A, 'U', 'N', INFO )

      WRITE(*,*) 'AP on exit:'
      DO I=1,N; WRITE(*,"(4(I3,1X),I3,1X)") INT(AP(I,:)); ENDDO

      WRITE(*,*)'INFO on exit:', INFO

      END PROGRAM LA_TPPACK_EXAMPLE