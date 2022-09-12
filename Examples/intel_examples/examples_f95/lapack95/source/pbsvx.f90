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
!     S P B S V X  Example Program Text
!*******************************************************************************

      PROGRAM SPBSVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: PBSVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Parameters" ..
      CHARACTER(LEN=1) :: EQUED
!  .. "Local Scalars" ..
      INTEGER :: KD, I, J, INFO, N, NRHS
      REAL(WP) :: RCOND
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), X(:,:), S(:)
      REAL(WP), ALLOCATABLE :: FERR(:), BERR(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPBSVX Example Program Results.'
      N = 7; KD = 3; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE ( A(KD+1,N), B(N,NRHS) )
      ALLOCATE ( X(N,NRHS), S(N), FERR(NRHS), BERR(NRHS) )

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,'(F3.0)') A(I,J);
         ENDDO
      ENDDO


      B(:,1) = 0.0
      DO I = 1, N
         DO J = MAX(1,-N+I+KD+1), KD
            B(I,1) = A(J,I-J+KD+1) + B(I,1)
         ENDDO
         DO J = MAX(1,KD+2-I), KD+1
            B(I,1) = A(J,I) + B(I,1)
         ENDDO
      ENDDO
      DO J = 2, NRHS; B(:,J) = B(:,1)*J; ENDDO

      WRITE(*,*) 'AB on entry:'
      DO I = 1, KD+1; WRITE (*,'(7(F9.5))') A(I,:); ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO J = 1, N; WRITE (*,'(3(F10.5))') B(J,:); ENDDO

print*, "***"
      CALL PBSVX(A, B, X, FACT='E', EQUED=EQUED, S=S, FERR=FERR, BERR=BERR, RCOND=RCOND, INFO=INFO)
print*, "***"
      WRITE(*,*)'EQUED = ',EQUED
      WRITE(*,*)'S on exit :'
      DO I=1,N; WRITE(*,"(5(F8.5))") S(I);
      ENDDO
      WRITE(*,*)'X on exit :'
      DO I=1,N; WRITE(*,"(5(F8.5))") X(I,:);
      ENDDO
      WRITE(*,*)'FERR on exit :'
      DO I=1,NRHS; WRITE(*,"(5(E13.6))") FERR(I);
      ENDDO
      WRITE(*,*)'BERR = '
      DO I=1,NRHS; WRITE(*,"(5(E13.6))") BERR(I);
      ENDDO
      WRITE(*,*)'RCOND = ',RCOND

      WRITE(*,*)'INFO = ' ,INFO

      WRITE(*,*) '\noindent'
      WRITE(*,*) 'The solution of the system $ A\,X = B $ is:'
      WRITE(*,*) '$$ X = \left( \begin{array}{rrr}'
      DO I=1,N; WRITE(*,"(2(F9.5,' & '),F9.5,' \\')") X(I,:); ENDDO
      WRITE(*,*) '\end{array} \right). $$'

      DEALLOCATE(A, B, X, S, FERR, BERR)

      END PROGRAM SPBSVX_MAIN
