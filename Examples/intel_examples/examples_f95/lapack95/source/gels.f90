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
!      Example shows how to migrate from NETLIB LAPACK95 to
!      Intel(R) Math Kernel Library (Intel(R) MKL) LAPACK95
!*******************************************************************************
      PROGRAM SGELS_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GELS
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, M, N, NRHS
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), BB(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GELS Example Program Results'
      M=6; N = 4; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE( A(M,N), AA(M,N), B(M,NRHS), BB(M,NRHS) )

      DO J=1,N
         DO I=1,M
            READ(*,'(F2.0)') A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(24(F9.5))") A(I,:);
      ENDDO

      AA=A

      DO J=1,NRHS
         DO I=1,M
            READ(*,'(F3.0)') B(I,J)
         ENDDO
      ENDDO

      BB=B
      WRITE(*,*)'Matrix B :'
      DO I=1,M;
         WRITE(*,"(3(F9.5))") B(I,:);
      ENDDO

      WRITE(*,*) 'CALL GELS( A, B )'
      CALL GELS(  A, B )

      WRITE(*,*) 'The matrix B on exit :'
      DO I=1,M;
         WRITE(*,"(3(E13.6))") B(I,:);
      ENDDO
      WRITE(*,*)
      WRITE(*,*)' * EXAMPLE 2 * '

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(24(F9.5))") AA(I,:);
      ENDDO

      WRITE(*,*)'Matrix B(:,1) :'
      WRITE(*,"(F9.5)") BB(:,1);

      WRITE(*,*) "CALL GELS( A, B(:,1), 'T', INFO )"
      CALL GELS( AA, BB(:,1), 'T', INFO )

      WRITE(*,*) 'The matrix B on exit :'
      WRITE(*,"(E13.6)") BB(:,1);

      WRITE(*,*) 'INFO = ', INFO

      END PROGRAM SGELS_MAIN
