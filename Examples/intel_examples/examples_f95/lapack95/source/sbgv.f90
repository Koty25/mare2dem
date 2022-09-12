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
!     S S B V G  Example Program Text
!*******************************************************************************

      PROGRAM SSBGV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SBGV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N, KD
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), W(:), Z(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SBGV Example Program Results'
      N = 5; KD = 2;
      ALLOCATE( A(KD+1,N), B(KD+1,N), W(N), Z(N,N) )

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix AB : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      B = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix BB : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(B(I,:));
      ENDDO

      WRITE(*,*) 'CALL SBGV( AB, BB, W )'
      CALL SBGV(  A, B, W)

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO

      WRITE(*,*) 'BB on exit : '
      DO I=1,KD+1;
         WRITE(*,"(5(E14.6,1X))") B(I,:);
      ENDDO

      WRITE(*,*)
      WRITE(*,*) ' * EXAMPLE 2 * '

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix AB : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO


      B = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix BB : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(B(I,:));
      ENDDO

      WRITE(*,*) "CALL SBGV( AB, BB, W, 'L', Z, INFO )"
      CALL SBGV( A, B, W, 'L', Z, INFO )

      WRITE(*,*)'W on exit : '
      DO I=1,N
         WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO

      WRITE(*,*) 'Z on exit : '
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") Z(I,:);
      ENDDO

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, B, W, Z)

      END PROGRAM SSBGV_MAIN
