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
!      Example shows how to migrate from NETLIB LAPACK95 to
!      Intel(R) Math Kernel Library (Intel(R) MKL) LAPACK95
!*******************************************************************************
      PROGRAM SGBSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GBSV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: K, KL, KU, I, J, N, NRHS, INFO
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: AB(:,:), B(:,:)
!  .. "Executable Statements" ..
       WRITE (*,*) 'SGBSV Example Program Results.'
       N = 6; KL = 2; KU = 1; NRHS = 2
       ALLOCATE ( AB(2*KL+KU+1,N), B(N,NRHS), IPIV(N))
       AB=0.0; B=0.0; IPIV=0

       DO I=KL+1,2*KL+KU+1
       DO J=1,N
          READ(*,'(F2.0)') AB(I,J)
       ENDDO
       ENDDO

       WRITE(*,*) 'The matrix AB:'
       DO I=1,N; WRITE(*,"(5(I3,1X,1X),I3,1X)") INT(AB(I,:));
       ENDDO

       DO I = 1, NRHS
       DO J = 1, N
       DO K = MAX(1,J-KL),MIN(J+KU,N); B(J,I) = AB(KL+KU+1+J-K,K)+B(J,I);
       ENDDO
       ENDDO;
       B(:,I) = B(:,I)*I;
       ENDDO

       WRITE(*,*) 'The RHS matrix B:'
       DO I=1,N; WRITE(*,"(1(I3,1X),I3,1X)") INT(B(I,:)); ENDDO
       WRITE(*,*) 'CALL GBSV( AB, B, 2, IPIV, INFO)'

       CALL GBSV( AB, B, 2, IPIV, INFO )

       WRITE(*,*)'AB on exit: '
       DO I=1,2*KL+KU+1; WRITE(*,"(5(E12.6,1X),E12.6,1X)") AB(I,:); ENDDO
       WRITE(*,*)'B on exit: '
       DO I=1,N; WRITE(*,"(1(E12.6,1X),E12.6,1X)") B(I,:); ENDDO
       WRITE(*,*)'IPIV on exit: ', IPIV
       WRITE(*,*)'INFO on exit: ', INFO

       END PROGRAM SGBSV_MAIN





