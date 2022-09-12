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
!     S G B S V X  Example Program Text
!*******************************************************************************

      PROGRAM SGBSVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GBSVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: K, KL, KU, I, J, N, NRHS
      CHARACTER(LEN=1) :: EQUED
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: AB(:,:), B(:,:), X(:,:), R(:), C(:), &
                               FERR(:), BERR(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GBSVX Example Program Results.'
      N = 6; KL = 2; KU = 1; NRHS = 2
      ALLOCATE ( AB(  KL+KU+1,N), B(N,NRHS), X(N,NRHS), R(N), C(N), &
                 FERR(NRHS),BERR(NRHS) )

      AB = 0
!     DO I=KL+1,2*KL+KU+1
      DO I=   1,  KL+KU+1
      DO J=1,N
         READ(*,'(F2.0)') AB(I,J)
      ENDDO
      ENDDO

      AB(:,1:1)=1e-6*AB(:,1:1)
      DO J = 1, N
        DO K = 1, 1+KU-MIN(J-1,KU)-1; AB(K,J) = HUGE(1.0_WP); END DO
        DO K = 1+ KU+MIN(N-J,KL)+KL+1, 1+KU+KL; AB(K,J) = HUGE(1.0_WP); END DO
      END DO
      WRITE(*,*) 'The array AB:'
      DO I=1,N; WRITE(*,*) AB(I,:); ENDDO

      B = 0.0_WP
      DO I = 1, NRHS;
        DO J = 1, N
          DO K = 1+KU-MIN(J-1,KU), 1+ KU+MIN(N-J,KL)
            B(J,I) = AB(K,J)+B(J,I);
          ENDDO
        ENDDO;
        B(:,I) = B(:,I)*I;
      ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO I=1,N; WRITE(*,"(1(I3,1X),I3,1X)") INT(B(I,:)); ENDDO

      WRITE(*,*) "CALL GBSVX( AB, B, X, 2, TRANS='T', EQUED=EQUED,", &
                 " R=R, C=C, FERR=FERR, BERR=BERR )"

      CALL GBSVX( AB, B, X, 2, FACT='E',TRANS='T',EQUED=EQUED, &
                     R=R, C=C, FERR=FERR, BERR=BERR )

      WRITE(*,*)'X on exit: '
      DO I=1,N; WRITE(*,"(1(E12.6,1X),E12.6,1X)") X(I,:); ENDDO
      WRITE(*,*)'EQUED = ', EQUED
      WRITE(*,*)'R = ', R
      WRITE(*,*)'C = ', C
      WRITE(*,*)'FERR = ', FERR
      WRITE(*,*)'BERR = ', BERR

      DEALLOCATE(AB, B, X, R, C, FERR, BERR)

      END PROGRAM SGBSVX_MAIN
