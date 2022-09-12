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
!     S P B S V  Example Program Text
!*******************************************************************************

      PROGRAM SPBSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: PBSV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Parameters" ..
      REAL(WP), PARAMETER :: ZERO = 0.0_WP, ONE = 1.0_WP
!  .. "Local Scalars" ..
      INTEGER :: KD, I, J, INFO, N, NRHS
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: AB(:,:), B(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPBSV Example Program Results.'
      N = 7; KD = 3; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE ( AB(KD+1,N), B(N,NRHS) )

      AB = 0
      DO I=1,KD+1
         DO J=1,N
         READ(*,'(F3.0)') AB(I,J);
         ENDDO
      ENDDO

      B(:,1) = ZERO
      DO I = 1, N
         DO J = MAX(1,-N+I+KD+1), KD
         B(I,1) = AB(J,I-J+KD+1) + B(I,1)
         ENDDO
         DO J = MAX(1,KD+2-I), KD+1
         B(I,1) = AB(J,I) + B(I,1)
         ENDDO
      ENDDO

      DO J = 2, NRHS; B(:,J) = B(:,1)*J; ENDDO

      WRITE(*,*) 'AB on entry:'
      DO I = 1,KD+1;
         WRITE (*,'(7(F8.5))') AB(I,:);
      ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO J = 1, N; WRITE (*,*) B(J,:);
      ENDDO

      CALL PBSV(AB, B, INFO=INFO)

      WRITE(*,*)'AB on exit :'
      DO I=1,KD+1; WRITE(*,"(10(E12.6))") AB(I,:);
      ENDDO

      WRITE(*,*)'B on exit :'
      DO J=1,N; WRITE(*,"(10(F8.5))") B(J,:);
      ENDDO

      WRITE(*,*)'INFO = ' ,INFO

      DEALLOCATE(AB, B)

      END PROGRAM SPBSV_MAIN
