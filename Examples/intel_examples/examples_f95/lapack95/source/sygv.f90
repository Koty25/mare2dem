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
!     S S Y G V  Example Program Text
!*******************************************************************************

      PROGRAM SSYGV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYGV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), BB(:,:), W(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SYEV Example Program Results'
      N = 5
      ALLOCATE( A(N,N), AA(N,N), B(N,N), BB(N,N), W(N) )

      A = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      AA=TRANSPOSE(A)

      WRITE(*,*)'Matrix A : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      B = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      BB=TRANSPOSE(B)

      WRITE(*,*)'Matrix B : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(B(I,:));
      ENDDO

      WRITE(*,*) 'CALL SYGV( A, B, W )'
      CALL SYGV( A, B, W )

      WRITE(*,*)'Matrix B on exit:'
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") B(I,:);
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO

      WRITE(*,*)
      WRITE(*,*) ' * EXAMPLE 2 * '

      WRITE(*,*)'Matrix A : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(AA(I,:));
      ENDDO

      WRITE(*,*)'Matrix B : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(BB(I,:));
      ENDDO

      WRITE(*,*) "CALL SYGV( A, B, W, '2', 'V', 'L', INFO )"
      CALL SYGV( AA, BB, W, 2, 'V', 'L', INFO )

      WRITE(*,*)'Matrix A on exit:'
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") AA(I,:);
      ENDDO

      WRITE(*,*)'Matrix B on exit:'
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") BB(I,:);
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(F11.5,1X))") W(I)
      ENDDO

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, AA, B, BB, W)

      END PROGRAM SSYGV_MAIN
