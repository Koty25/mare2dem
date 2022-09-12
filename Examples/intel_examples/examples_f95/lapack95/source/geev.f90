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
!     C G E E V  Example Program Text
!*******************************************************************************

      PROGRAM CGEEV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GEEV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N
      REAL(WP), ALLOCATABLE ::  WR(:), WI(:)
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:)
      COMPLEX(WP), ALLOCATABLE :: AA(:,:), W(:), VL(:,:), VR(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC  REAL, AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'GEEV Example Program Results'
      N = 5
      ALLOCATE( A(N,N), AA(N,N), W(N), WR(N), WI(N), VL(N,N), VR(N,N) )

      A = 0
      AA = 0
      DO I=1,N
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N
         WRITE(*,"(5(I3,1X))") INT(A(I,:))
      ENDDO

      WRITE(*,*) "CALL GEEV( A, WR, WI ) "
      CALL GEEV( A, WR, WI )

      WRITE(*,*) 'WR on exit : '
      DO I=1,N
         WRITE(*,"((F9.6,1X))") REAL(WR(I))
      ENDDO

      WRITE(*,*) 'WI on exit : '
      DO I=1,N
         WRITE(*,"((F9.5,1X))") REAL(WI(I))
      ENDDO

      WRITE(*,*)
      WRITE(*,*) ' * EXAMPLE 2 * '

      DO I=1,N
         DO J=1,N
            READ(*,*) AA(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A : '
      DO J=1,N
         WRITE(*,"(5(I3,1X,'+',I3,'i',1X,1X,1X))") INT(AA(J,1)),INT(AIMAG(AA(J,1))), &
              INT(AA(J,2)),INT(AIMAG(AA(J,2))), &
              INT(AA(J,3)),INT(AIMAG(AA(J,3))), &
              INT(AA(J,4)),INT(AIMAG(AA(J,4))), &
              INT(AA(J,5)),INT(AIMAG(AA(J,5)))
      ENDDO

      WRITE(*,*) "CALL GEEV( A, W, VL, VR, INFO )"
      CALL GEEV( AA, W, VL, VR, INFO )

      WRITE(*,*) ' W on exit : '
      DO I=1,N;
         WRITE(*,"('('E14.7,1X,',',1X,E14.7,')',1X)")  W(I)
      ENDDO

      WRITE(*,*)' VL on exit : '
      DO I=1,N;
         WRITE(*,"(5('('E14.7,1X,',',1X,E14.7,')',1X))")  VL(I,:)
      ENDDO;

      WRITE(*,*) ' VR on exit : '
      DO I=1,N;
         WRITE(*,"(5('('E14.7,1X,',',1X,E14.7,')',1X))")  VR(I,:)
      ENDDO;

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE( A, AA, W, WR, WI, VL, VR )

      END PROGRAM CGEEV_MAIN
