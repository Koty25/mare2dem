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
      PROGRAM SSYGVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYGVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IFAIL(:)
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:), W(:), Z(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SYEVX Example Program Results'
      N = 5
      ALLOCATE( A(N,N), B(N,N), W(N), IFAIL(N), Z(N,N) )

      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(A(I,:));
      ENDDO

      DO J=1,N
         DO I=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix B:'
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(B(I,:));
      ENDDO

      WRITE(*,*) 'CALL SYGVX( A, B, W, 3, Z=Z, VL=-10.0_WP, VU=10.0_WP, M=M, IFAIL=IFAIL )'
      CALL SYGVX(  A, B, W, 3, Z=Z, VL=-10.0_WP, VU=10.0_WP, M=M, IFAIL=IFAIL )

      WRITE(*,*)'Matrix B on exit:'
      DO I=1,N;
         WRITE(*,"(5(F12.7,1X))") B(I,:);
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,M
      WRITE(*,"(5(F12.7,1X))") W(I)
      ENDDO

      WRITE(*,*) 'M=',M

      WRITE(*,*)'IFAIL = ', IFAIL

      WRITE(*,*) 'Z on exit : '
      DO J=1,M
      WRITE(*,"(5(F12.6,1X))") (Z(I,J), I = 1,N)
      ENDDO

      DEALLOCATE( A, B, W, IFAIL, Z )

      END PROGRAM SSYGVX_MAIN


