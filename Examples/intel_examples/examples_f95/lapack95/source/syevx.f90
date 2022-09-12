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
!     S S Y E V X  Example Program Text
!*******************************************************************************

      PROGRAM SSYEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), W(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SYEVX Example Program Results'
      N = 5
      ALLOCATE( A(N,N), W(N) )

      A = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N;
         WRITE(*,"(5(F9.5))") A(I,:);
      ENDDO

      WRITE(*,*) 'CALL SYEVX( A, W, VL=-7.0_WP, VU=7.0_WP, M=M, ABSTOL=1.0E-4_WP )'
      W = 0
      M = 0
      CALL SYEVX( A, W, VL=-7.0_WP, VU=7.0_WP, M=M, ABSTOL=1.0E-4_WP )
      WRITE(*,*) 'W on exit'
      DO I=1,N
      WRITE(*,"(5(F9.5))") W(I)
      ENDDO
      WRITE(*,*)' M =',M

      DEALLOCATE(A, W)

      END PROGRAM SSYEVX_MAIN
