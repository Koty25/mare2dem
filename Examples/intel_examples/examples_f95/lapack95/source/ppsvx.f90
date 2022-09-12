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
!     S P P S V X  Example Program Text
!*******************************************************************************

      PROGRAM SPPSVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: PPSVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NN, NRHS
      REAL(WP) :: RCOND
!  .. Local Arrays ..
      REAL(WP), ALLOCATABLE :: A(:), B(:,:),X(:,:), FERR(:), BERR(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPPSVX Example Program Results.'
      N = 5; NRHS = 1
      NN = N*(N+1)/2
      ALLOCATE ( A(NN), B(N,NRHS), X(N,NRHS),FERR(NRHS), BERR(NRHS) )

      A = 0
      DO I=1,NN
         READ(*,'(F3.0)') A(I);
      ENDDO;

      WRITE(*,*) 'The array A :'
      DO I=1,NN
      WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(A(I));
      ENDDO

      B = 0
      DO I=1,N
      DO J=1,NRHS
         READ(*,'(F3.0)') B(I,J);
      ENDDO;
      ENDDO;

      WRITE(*,*) 'The array B :'
      DO I=1,N
      WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(B(I,:));
      ENDDO;

      WRITE(*,*) 'CALL PPSVX( A, B, X, "L", FERR=FERR, BERR=BERR, RCOND=RCOND )'
      CALL PPSVX(  A, B, X, 'L', FERR=FERR, BERR=BERR, RCOND=RCOND )

      WRITE(*,*)'X on exit :'
      DO I=1,N; WRITE(*,"(5(E13.6))") X(I,:);
      ENDDO
      WRITE(*,*)'FERR on exit :'
      DO I=1,NRHS; WRITE(*,"(5(E13.6))") FERR(I);
      ENDDO
      WRITE(*,*)'BERR = '
      DO I=1,NRHS; WRITE(*,"(5(E13.6))") BERR(I);
      ENDDO
      WRITE(*,*)'RCOND = ', RCOND

      WRITE(*,*) '\noindent'
      WRITE(*,*) 'The solution of the system $ A\,X = B $ is:'
      WRITE(*,*) '$$ X = \left( \begin{array}{rrr}'
      DO I=1,N; WRITE(*,"(2(F9.5,' & '),F9.5,' \\')") X(I,:); ENDDO
      WRITE(*,*) '\end{array} \right). $$'

      DEALLOCATE(A, B, X, FERR, BERR)

      END PROGRAM SPPSVX_MAIN
