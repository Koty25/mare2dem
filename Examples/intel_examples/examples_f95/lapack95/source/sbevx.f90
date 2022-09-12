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
      PROGRAM SSBEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SBEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N, KD
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), W(:), Z(:,:), Q(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SBEVX Example Program Results'
      N = 5; KD = 2;
      ALLOCATE( A(KD+1,N), W(N), Z(N,N), Q(N,N) )

 !     OPEN(UNIT=21,FILE='sbevu.ma',STATUS='UNKNOWN')
      DO I=1,KD+1
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,KD+1;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      WRITE(*,*) 'CALL SBEVX( A, W, Z=Z, VL=-4.0_WP, VU=100.0_WP, M=M, Q=Q )'
      CALL SBEVX( A, W, Z=Z, VL=-4.0_WP, VU=100.0_WP, M=M, Q=Q)

      WRITE(*,*) 'W on exit : '
      DO I=1,M
      WRITE(*,"(5(F9.5))") W(I)
      ENDDO
      WRITE(*,*) 'M = ',M

      WRITE(*,*) 'The matrix Q'
      WRITE(*,'(5(F9.5))') ((Q(I,J),J=1,N),I=1,N)

      END PROGRAM SSBEVX_MAIN


