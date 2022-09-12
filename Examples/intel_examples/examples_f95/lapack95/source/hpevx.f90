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
!     C H P E V X  Example Program Text
!*******************************************************************************

      PROGRAM CHPEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HPEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N
      REAL(WP), ALLOCATABLE :: W(:)
!  .. "Local Arrays" ..
      COMPLEX(WP), ALLOCATABLE :: A(:), Z(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'HPEVX Example Program Results'
      N = 5
      ALLOCATE( A(N*(N+1)/2), W(N), Z(N,N) )

      A = 0
      DO J=1,N*(N+1)/2
         READ(*,*) A(J)
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N*(N+1)/2
         WRITE(*,"((I3,1X,'+',1X,I3,'i',1X))") INT(A(I)), INT(AIMAG(A(I)))
      ENDDO


      WRITE(*,*) "CALL HPEVX( A, W, IL=2, IU=5, ABSTOL=1.0E-5_WP ) "
      W = 0
      CALL HPEVX( A, W, IL=2, IU=5, ABSTOL=1.0E-5_WP )

      WRITE(*,*) 'W on exit : '
      DO I=1,N
      WRITE(*,"(5(F9.5))") W(I)
      ENDDO

      DEALLOCATE(A, W, Z)

      END PROGRAM CHPEVX_MAIN
