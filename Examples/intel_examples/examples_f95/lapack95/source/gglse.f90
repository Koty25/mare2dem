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
!     S G G L S E  Example Program Text
!*******************************************************************************

      PROGRAM SGGLSE_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
!      USE lapack95, ONLY: GGLSE, LANGE
      USE lapack95, ONLY: GGLSE
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" .
      REAL(WP) ::  R1
      INTEGER :: I, J, INFO, M, N, P
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), C(:), D(:), X(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GGLSE Example Program Results'
      M=5; N = 3; P = 2
      ALLOCATE( A(M,N), B(P,N), C(M), D(P), X(N) )

      A = 0
      DO I=1,M
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      B = 0
      DO I=1,P
         DO J=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      C = 0
      DO I=1,M
         READ(*,*) C(I)
      ENDDO

      D = 0
      DO I=1,P
         READ(*,*) D(I)
      ENDDO

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(3(F9.5))") A(I,:);
      ENDDO

      WRITE(*,*)'Matrix B : '
      DO I=1,P;
         WRITE(*,"(3(F9.5))") B(I,:);
      ENDDO

      WRITE(*,*)'Vector C : '
      DO I=1,M;
         WRITE(*,"(F9.5)") C(I);
      ENDDO

      WRITE(*,*) 'Vector D : '
      DO I=1,P;
         WRITE(*,"(F9.5)") D(I);
      ENDDO
      WRITE(*,*)
      WRITE(*,*) 'CALL GGLSE( A, B, C, D, X, INFO )'
      CALL GGLSE( A, B, C, D, X, INFO )

      WRITE(*,*) 'C on exit : '
      DO I=1,M;
         WRITE(*,"(E14.6)") C(I);
      ENDDO

      WRITE(*,*) 'X on exit : '
      DO I=1,N;
         WRITE(*,"(E14.5)") X(I);
      ENDDO

      WRITE(*,*)'INFO = ',INFO

      DEALLOCATE(A, B, C, D, X)

      END PROGRAM SGGLSE_MAIN
