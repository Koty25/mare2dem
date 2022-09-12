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
      PROGRAM CHPEV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HPEV, HPEVD
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, INFO, N
      REAL(WP), ALLOCATABLE :: W(:)
!  .. "Local Arrays" ..
      COMPLEX(WP), ALLOCATABLE :: A(:), Z(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC  AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'HPEV Example Program Results'
      N = 5
      ALLOCATE( A(N*(N+1)/2), W(N), Z(N,N) )

      DO I=1,N*(N+1)/2
         READ(*,*) A(I)
      ENDDO

      WRITE(*,*)'Matrix A : '
      DO I=1,N*(N+1)/2
         WRITE(*,"((I3,1X,'+',1X,I3,'i',1X))") INT(A(I)),INT(AIMAG(A(I)))
      ENDDO

      WRITE(*,*) "CALL HPEV( A, W) "
      CALL HPEV( A, W )

      WRITE(*,*) 'W on exit : '
      DO I=1,N
      WRITE(*,"(5(E14.6))") W(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) ' * EXAMPLE 2 * '

      DO I=1,N*(N+1)/2
         READ(*,*) A(I)
      ENDDO

      WRITE(*,*)'Matrix A : '
      DO I=1,N*(N+1)/2
         WRITE(*,"((I3,1X,'+',1X,I3,'i',1X))") INT(A(I)),INT(AIMAG(A(I)))
      ENDDO

      WRITE(*,*) "CALL HPEVD( A, W, 'L', Z, INFO )"
      CALL HPEVD( A, W, 'L', Z, INFO )

      WRITE(*,*)'Z on exit: '
      DO I=1,N; WRITE(*,*) Z(I,:);
      ENDDO

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, W, Z)

      END PROGRAM CHPEV_MAIN
