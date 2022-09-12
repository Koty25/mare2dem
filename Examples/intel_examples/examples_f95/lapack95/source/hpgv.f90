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
!     C H P G V  Example Program Text
!*******************************************************************************

      PROGRAM CHPGV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HPGV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N
      REAL(WP), ALLOCATABLE :: W(:)
!  .. "Local Arrays" ..
      COMPLEX(WP), ALLOCATABLE :: A(:), AA(:), B(:), BB(:), Z(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC REAL, AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'HPGV Example Program Results'
      N = 5
      ALLOCATE( A(N*(N+1)/2), AA(N*(N+1)/2), B(N*(N+1)/2), BB(N*(N+1)/2), W(N), Z(N,N) )

      A = 0
      DO J=1,N*(N+1)/2
         READ(*,*) A(J)
      ENDDO

      WRITE(*,*)'Matrix AP : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,',',I3)')')") INT(A(I)), INT(AIMAG(A(I)))
      ENDDO

      B = 0
      DO J=1,N*(N+1)/2
         READ(*,*) B(J)
      ENDDO

      WRITE(*,*)'Matrix B : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,',',I3)')')") INT(B(I)), INT(AIMAG(B(I)))
      ENDDO

      WRITE(*,*) "CALL HPGV( AP, BP, W) "
      CALL HPGV( A, B, W )

      WRITE(*,*)'BP on exit : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(E14.6,1X,','E14.6)')')") REAL(B(I)), AIMAG(B(I))
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,N
         WRITE(*,"(5(E14.6))") W(I)
      ENDDO

      WRITE(*,*)
      WRITE(*,*)' * EXAMPLE 2 * '

      AA = 0
      DO J=1,N*(N+1)/2
         READ(*,*) AA(J)
      ENDDO

      WRITE(*,*)'Matrix AP : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,',',I3)')')") INT(AA(I)), INT(AIMAG(AA(I)))
      ENDDO

      BB = 0
      DO J=1,N*(N+1)/2
         READ(*,*) BB(J)
      ENDDO

      WRITE(*,*)'Matrix BP : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,',',I3)')')") INT(BB(I)), INT(AIMAG(BB(I)))
      ENDDO

      WRITE(*,*)
       WRITE(*,*) "CALL HPGV( A, B, W, 3, 'L', Z, INFO)"
       CALL HPGV( AA, BB, W, 3, 'L', Z, INFO )

      WRITE(*,*)'Z on exit : '
      DO I=1,N
         DO J=1,N-1
            WRITE(*,"('('(E14.6,1X,',',1X,E14.6)')')") REAL(Z(I,J)), AIMAG(Z(I,J))
         ENDDO
         WRITE(*,"('('(E14.6,1X,',',1X,E14.6)')')") REAL(Z(I,N)), AIMAG(Z(I,N))
      ENDDO
      WRITE(*,*)'BP on exit : '
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(E14.6,1X,',',1X,E14.6)')')") REAL(BB(I)), AIMAG(BB(I))
      ENDDO

      WRITE(*,*) 'W on exit : '
      DO I=1,N
      WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO
      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, AA, B, BB, W, Z)

      END PROGRAM CHPGV_MAIN
