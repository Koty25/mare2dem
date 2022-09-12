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
!     S P O S V  Example Program Text
!*******************************************************************************

      PROGRAM SPOSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: POSV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NRHS
      CHARACTER(LEN=1) :: UPLO
!  .. "Local Arrays" ..
 !    INTEGER, DIMENSION(NP) :: IPUT
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), AA(:,:),BB(:,:), S(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPOSV Example Program Results'
      N = 5; NRHS = 3

      ALLOCATE ( A(N,N), B(N,NRHS), AA(N,N), BB(N,NRHS), S(N) )

      A = 0
      DO I=1,N
         DO J=I,N
         READ(*,'(F3.0)') A(I,J);
         ENDDO
      ENDDO;

      WRITE(*,*) 'The array A :'
      DO I=1,N
      WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(A(I,1:N));
      ENDDO

      DO I = 1, N
        DO J = 1, NRHS
        B(I,J) = (SUM(A(I,I:N)) + SUM(A(1:I-1,I)))*J
        ENDDO
      ENDDO

      AA = TRANSPOSE(A)
      BB = B
      WRITE(*,*) 'The array B :'
      DO I=1,N
      WRITE(*,'(3(I3,1X,1X))') INT(B(I,1:NRHS));
      ENDDO

      WRITE(*,*) 'CALL POSV( A, B )'
      CALL POSV(  A, B )

      WRITE(*,*)'A on exit: '
      DO I=1,N; WRITE(*,"(5(E15.6))") A(I,1:N);
      ENDDO

      WRITE(*,*)'B on exit: '
      DO I=1,N; WRITE(*,"(3(F9.5))") B(I,1:NRHS);
      ENDDO

      WRITE(*,*)' * Example 2 * '

      WRITE(*,*) 'The array A :'
      DO J=1,N
      WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(AA(J,:));
      ENDDO

      WRITE(*,*) 'The array B(:,1) :'
      Do I=1,N
      WRITE(*,'(3(I3,1X,1X))') INT(BB(I,1));
      ENDDO

      WRITE(*,*) "CALL POSV( A, B(:,1), 'L' )"
      UPLO='L'
      CALL POSV(  AA, BB(:,1), 'L' )

      WRITE(*,*)'A on exit: '
      DO J=1,N; WRITE(*,"(6(E13.6))") AA(J,:); ENDDO

      WRITE(*,*)'B on exit: '
      DO I=1,N; WRITE(*,"(3(F8.5))") BB(I,1); ENDDO

      DEALLOCATE(A, B, AA, BB, S)

      END PROGRAM SPOSV_MAIN
