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
!     S G E S V D  Example Program Text
!*******************************************************************************

      PROGRAM SGESVD_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GESVD
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, M, N
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), AA(:,:), S(:), U(:,:), V(:,:), VT(:,:), WW(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GESVD Example Program Results'
      N = 5; M=3
      ALLOCATE( A(M,N), AA(M,N), S(MIN(M,N)), U(M,M), V(N,N), VT(N,N), WW(MIN(M,N)-1))

      A = 0
      DO J=1,N
         DO I=1,M
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      AA=A

      WRITE(*,*)'Matrix A:'
      DO I=1,M; WRITE(*,"(5(I3,1X))") INT(A(I,:)); ENDDO

      WRITE(*,*)
      WRITE(*,*) 'CALL GESVD( A, S )'
      CALL GESVD( A, S )

      WRITE(*,*) 'S on exit : '
      WRITE(*,"(3(F10.5,1X))") S(:)

      WRITE(*,*)
      WRITE(*,*)' * EXAMPLE 2 * '
      WRITE(*,*) "CALL GESVD( A, S, VT=VT, WW=WW, JOB='U', INFO=INFO )"
      CALL GESVD( AA, S, VT=VT, WW=WW, JOB='U', INFO=INFO )

      WRITE(*,*) 'A on exit : '
      DO I=1,M;
         WRITE(*,"(5(E14.6,1X))") AA(I,:);
      ENDDO

      WRITE(*,*) 'VT on exit : '
      DO I=1,N;
         WRITE(*,"(5(E14.6,1X))") VT(I,:);
      ENDDO

      WRITE(*,*)'WW on exit : '
      WRITE(*,"(5(E14.6,1X))") WW(:)

      WRITE(*,*) ' INFO = ', INFO

      DEALLOCATE(A, AA, S, U, V, VT, WW)

      END PROGRAM SGESVD_MAIN
