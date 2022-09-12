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
!     S G E L S S  Example Program Text
!*******************************************************************************

      PROGRAM SGELSS_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GELSS
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" .
      REAL(WP) ::  R1, R2, R3
      INTEGER :: RANK, I, J, INFO, M, N, NRHS
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), S(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GELSS Example Program Results'
      M=6; N = 4; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE( A(M,N), B(M,NRHS), S(N) )

      A = 0
      B = 0
      DO I=1,M
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      DO I=1,M
         DO J=1,NRHS
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(24(F9.5))") A(I,:);
      ENDDO

      WRITE(*,*)'Matrix B :'
      DO I=1,M;
         WRITE(*,"(3(F9.5))") B(I,:);
      ENDDO

      WRITE(*,*) 'CALL GELSS( A, B, RANK, S, 0.00001_WP, INFO )'
      CALL GELSS( A, B, RANK, S, 0.00001_WP, INFO=INFO )

      WRITE(*,*) ' A on exit : '
      DO I=1,M;
         WRITE(*,"(4(E14.6))") A(I,:);
      ENDDO

      WRITE(*,*) ' B on exit : '
      DO I=1,M;
         WRITE(*,"(3(E14.6))") B(I,:);
      ENDDO

      WRITE(*,*) ' S on exit : '
      DO I=1,N;
         WRITE(*,"(E14.6)") S(I);
      ENDDO

      WRITE(*,*) 'RANK = ', RANK
      WRITE(*,*)'INFO = ',INFO

      DEALLOCATE(A, B, S)

      END PROGRAM SGELSS_MAIN
