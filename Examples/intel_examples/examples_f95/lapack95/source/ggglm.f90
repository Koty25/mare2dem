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
!     S G G G L M  Example Program Text
!*******************************************************************************

      PROGRAM SGGGLM_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GGGLM
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" .
      INTEGER :: I, J, INFO, M, N, P
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), D(:), X(:), Y(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GGGLM Example Program Results'
      M=5; N = 4; P = 2
      ALLOCATE( A(M,N), B(M,P), D(M), X(N), Y(P) )

      A = 0
      DO I=1,M
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      B = 0
      DO I=1,M
         DO J=1,P
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      DO I=1,M
         READ(*,*) D(I)
      ENDDO

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(4(F9.5))") A(I,:);
      ENDDO
      WRITE(*,*)'Matrix B :'
      DO I=1,M;
         WRITE(*,"(2(F9.5))") B(I,:);
      ENDDO
      WRITE(*,*)'Vector D :'
      DO I=1,M;
         WRITE(*,"(F9.5)") D(I);
      ENDDO

      WRITE(*,*) 'CALL GGGLM( A, B, D, X, Y, INFO )'
      CALL GGGLM( A, B, D, X, Y, INFO )
      WRITE(*,*) 'X on exit:'
      DO I=1,N; WRITE(*,"(E14.6)") X(I); ENDDO

      WRITE(*,*) 'Y on exit:'
      DO I=1,P; WRITE(*,"(E14.6)") Y(I); ENDDO

      WRITE(*,*)'INFO = ',INFO

      DEALLOCATE(A, B, D, X, Y)

      END PROGRAM SGGGLM_MAIN
