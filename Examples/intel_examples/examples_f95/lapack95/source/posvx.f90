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
!     S P O S V X  Example Program Text
!*******************************************************************************

      PROGRAM SPOSVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: POSVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NRHS
      CHARACTER(LEN=1) :: EQUED
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), X(:,:), S(:), &
                               FERR(:), BERR(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPOSVX Example Program Results'
      N = 5; NRHS = 3
      ALLOCATE ( A(N,N), B(N,NRHS), X(N,NRHS), S(N), FERR(NRHS), &
                 BERR(NRHS) )

      A = 0
      DO I=1,N
         DO J=I,N
         READ(*,'(F3.0)') A(I,J);
         ENDDO
      ENDDO;

      A(:,1)=1E-6*A(:,1);A(1,2:N)=1E-6*A(1,2:N)

      DO J = 1, NRHS; DO I = 1, N
      B(I,J) = (SUM(A(I,I:N)) + SUM(A(1:I-1,I)))*J
      ENDDO; ENDDO

      WRITE(*,*) 'The array B :'
      DO J=1,NRHS; DO I = 1, N
       WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(B(I,J));
      ENDDO; ENDDO

      WRITE(*,*) 'CALL POSVX( A, B, X, FACT="E", EQUED=EQUED, S=S )'
      CALL POSVX(  A, B, X, FACT='E', EQUED=EQUED, S=S )

      WRITE(*,*)'EQUED = ', EQUED
      WRITE(*,*)'S = ', S
      WRITE(*,*)'FERR = ', FERR
      WRITE(*,*)'BERR = ', BERR

      WRITE(*,*) 'The solution of the system $ A\,X = B $ is:'
      WRITE(*,*) '$$ X = \left( \begin{array}{rrr}'
      DO I=1,N; WRITE(*,"(2(F9.5,' & '),F9.5,' \\')") X(I,:); ENDDO
      WRITE(*,*) '\end{array} \right). $$'

      DEALLOCATE(A, B, X, S, FERR)

      END PROGRAM SPOSVX_MAIN
