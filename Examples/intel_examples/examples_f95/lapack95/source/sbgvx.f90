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
!     S S B V G X  Example Program Text
!*******************************************************************************

      PROGRAM SSBGVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SBGVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N, KD
      REAL(WP) :: VL, VU
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), W(:), Z(:,:), Q(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SBGVX Example Program Results'
      N = 5; KD = 2;
      ALLOCATE( A(KD+1,N), B(KD+1,N), W(N), Z(N,N), Q(N,N) )

      A = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      B = 0
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix B : '
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(B(I,:));
      ENDDO

      WRITE(*,*)
      WRITE(*,*) 'CALL SBGVX( A, B, W, Z=Z, VL=0.0_WP, VU=100.0_WP, M=M, Q=Q )'
      VL=0; VU=100
      CALL SBGVX(  A, B, W, Z=Z, VL=0.0_WP, VU=100.0_WP, M=M, Q=Q )

      WRITE(*,*) ' W on exit : '
      DO I=1,N
      WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO

      WRITE(*,*) ' M = ',M

      WRITE(*,*) ' Matrix Q on exit : '
      DO I=1,N
         WRITE(*,'(5(E14.6,1X))') Q(I,:)
      ENDDO

      DEALLOCATE(A, B, W, Z, Q)

      END PROGRAM SSBGVX_MAIN
