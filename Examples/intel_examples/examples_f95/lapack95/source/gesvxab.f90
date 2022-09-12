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
!     S G E S V X  Example Program Text
!*******************************************************************************

      PROGRAM SGESVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GESVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NRHS, INFO
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP) :: RCOND, RPVGRW
      REAL(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), X(:,:),BB(:,:), FERR(:), BERR(:)
      REAL(WP), ALLOCATABLE :: RR(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SGESVX Example Program Results.'
      N = 2
      NRHS = 1
      ALLOCATE( A(N,N), AA(N,N), B(N,NRHS), X(N,NRHS),BB(N,NRHS), IPIV(N), RR(N,N), FERR(NRHS), BERR(NRHS) )

      AA(1,1) = 0; AA(1,2) = 1; AA(2,1) = 0; AA(2,2) = 1

      DO J = 1, NRHS; BB(:,J) = SUM( AA, DIM=2)*J; ENDDO

      WRITE(*,*) 'The matrix A:'
      DO I=1,N; WRITE(*,"(4(I3,1X),I3,1X)") INT(AA(I,:)); ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO I=1,N; WRITE(*,"(2(I3,1X),I3,1X)") INT(BB(I,:)); ENDDO

      WRITE(*,*) 'CALL GESVX( A, B, X, FERR=FERR, BERR=BERR, ', &
                 'RCOND=RCOND, RPVGRW=RPVGRW, INFO =INFO )'
      A=AA; B=BB
      CALL GESVX( A, B, X, FERR=FERR, BERR=BERR, RCOND=RCOND, &
                     RPVGRW=RPVGRW, INFO=INFO )

      WRITE(*,*) 'INFO = ', INFO
      WRITE(*,*) 'FERR = ', FERR
      WRITE(*,*) 'BERR = ', BERR
      WRITE(*,*) 'RCOND = ', RCOND
      WRITE(*,*) 'RPVGRW = ', RPVGRW

      WRITE(*,*) '\noindent'
      WRITE(*,*) 'The solution of the system $ A\,X = B $ is:'
      WRITE(*,*) '$$ X = \left( \begin{array}{rrr}'
      DO I=1,N; WRITE(*,"(2(F9.5,' & '),F9.5,' \\')") X(I,:); ENDDO
      WRITE(*,*) '\end{array} \right). $$'

      WRITE(*,*) 'The matrix A on exit:'
      DO I=1,N; WRITE(*,*) AA(I,:); ENDDO

      DEALLOCATE(A, AA, B, X, BB, IPIV, RR, FERR, BERR)

      END PROGRAM SGESVX_MAIN
