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
!     S S B E V  Example Program Text
!*******************************************************************************

      PROGRAM SSBEV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SBEV, SBEVD
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N, KD
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), W(:), Z(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SYEV Example Program Results'
      N = 5; KD = 2;
      ALLOCATE( A(KD+1,N), W(N), Z(N,N) )

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      WRITE(*,*) 'CALL SBEV( A, W )'

      CALL SBEV(  A, W)

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(F9.5))") W(I)
      ENDDO

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)
      WRITE(*,*) ' * EXAMPLE 2 * '

      WRITE(*,*)'Matrix A : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      WRITE(*,*) "CALL SBEVD( A, W, 'L', Z, INFO )"
      CALL SBEVD( A, W, 'L', Z, INFO )

      WRITE(*,*) 'Z on exit : '
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") Z(I,:);
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(F9.5))") W(I)
      ENDDO

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, W, Z)

      END PROGRAM SSBEV_MAIN
