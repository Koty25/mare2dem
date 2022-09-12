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
!     C G G S V D  Example Program Text
!*******************************************************************************

      PROGRAM CGGSVD_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GGSVD
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, M, N, P, K, L
      REAL(WP), ALLOCATABLE :: ALPHA(:), BETA(:)
!  .. "Local Arrays" ..
      COMPLEX(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), BB(:,:), U(:,:), V(:,:), Q(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC REAL, AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'GGSVD Example Program Results'
      N = 3; M=5; P=2
      ALLOCATE( A(M,N), AA(M,N), B(P,N), BB(P,N), ALPHA(N), BETA(N), U(M,M), V(P,P), Q(N,N) )

      A = 0
      DO J=1,N
         DO I=1,M
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      AA=A

      WRITE(*,*)'Matrix A : '
      DO I=1,M
         WRITE(*,"(3('('(I3,1X,',',I3)')',1X,1X))") INT(A(I,1)),INT(AIMAG(A(I,1))), &
              INT(A(I,2)),INT(AIMAG(A(I,2))), &
              INT(A(I,3)),INT(AIMAG(A(I,3)))
      ENDDO

      B = 0
      DO J=1,N
         DO I=1,P
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      BB=B
      WRITE(*,*)
      WRITE(*,*)'Matrix B : '
      DO I=1,P
         WRITE(*,"(3('('(I3,1X,',',I3)')',1X,1X))") INT(B(I,1)), INT(AIMAG(B(I,1))), &
              INT(B(I,2)),INT(AIMAG(B(I,2))), &
              INT(B(I,3)),INT(AIMAG(B(I,3)))
      ENDDO

      WRITE(*,*)
      WRITE(*,*) "CALL GGSVD( A, B, ALPHA, BETA, K, L ) "
      CALL GGSVD( A, B, ALPHA, BETA, K, L )

      WRITE(*,*)'A on exit : '
      DO I=1,M;
         WRITE(*,"(3('('(E13.6,1X,',',E13.6)')',1X))") REAL(A(I,1)),AIMAG(A(I,1)), &
              REAL(A(I,2)),AIMAG(A(I,2)), &
              REAL(A(I,3)),AIMAG(A(I,3))
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'ALPHA on exit : '
      DO I=1,N
         WRITE(*,"((E14.6,1X))") ALPHA(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'BETA on exit : '
      DO I=1,N
         WRITE(*,"((E14.6,1X))") BETA(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)' K = ', K
      WRITE(*,*)' L = ', L

      WRITE(*,*)
      WRITE(*,*)' * EXAMPLE 2 * '
      WRITE(*,*)
      WRITE(*,*) "CALL GGSVD( A, B, ALPHA, BETA, K, L, U, V, Q, INFO=INFO )"
      CALL GGSVD(  AA, BB, ALPHA, BETA, K, L, U, V, Q , INFO=INFO )

      WRITE(*,*)'U on exit : '
      DO I=1,M;
         WRITE(*,"(5('('(E13.6,1X,',',1X,E13.6)')',1X))") REAL(U(I,1)), AIMAG(U(I,1)), &
              REAL(U(I,2)), AIMAG(U(I,2)), &
              REAL(U(I,3)), AIMAG(U(I,3)), &
              REAL(U(I,4)), AIMAG(U(I,4)), &
              REAL(U(I,5)), AIMAG(U(I,5))
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'V on exit : '
      DO I=1,P;
         WRITE(*,"(2('('(E13.6,1X,',',E13.6)')',1X,1X))") REAL(V(I,1)), AIMAG(V(I,1)), &
              REAL(V(I,2)), AIMAG(V(I,2))
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'Q on exit : '
      DO I=1,N;
         WRITE(*,"(3('('(E13.6,1X,',',E13.6)')',1X))") REAL(Q(I,1)), AIMAG(Q(I,1)), &
              REAL(Q(I,2)), AIMAG(Q(I,2)), &
              REAL(Q(I,3)), AIMAG(Q(I,3))
      ENDDO
      WRITE(*,*)
      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, AA, B, BB, ALPHA, BETA, U, V, Q)

      END PROGRAM CGGSVD_MAIN
