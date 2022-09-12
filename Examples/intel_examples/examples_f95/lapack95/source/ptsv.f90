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
!     S P T S V  Example Program Text
!*******************************************************************************

      PROGRAM SPTSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: PTSV
!  .. "IMPLICIT STATEMENT" ..
      IMPLICIT NONE
!  .. "Parameters" ..
      REAL(WP), PARAMETER :: ONE = 1.0_WP, THREE = 3.0_WP
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N, NRHS
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:), D(:), E(:), B(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPTSV Example Program Results.'
      N = 5; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE (A(N), D(N), E(N-1), B(N,NRHS) )

      A = 0
      DO I=1,N
         READ(*,'(F3.0)') A(I);
      ENDDO;


      E = THREE; D(1) = 2*E(1)+ONE; D(2:N) = 2*E+ONE
      B(1,1) = D(1)+E(1); B(N,1) = E(N-1)+D(N)
      B(2:N-1,1) = 2*E + D(2:N-1)
      DO I = 2, NRHS; B(:,I) = B(:,1)*I;
      ENDDO
      WRITE(*,*) '  on entry:'
      WRITE (*,'(2HD:,8(1X,F9.2))') D
      WRITE (*,'(2HE:,8(1X,F9.2))') E
      WRITE(*,*) 'The RHS matrix B:'
      DO J = 1, N; WRITE (*,*) B(J,:);
      ENDDO

      WRITE(*,*)' CALL PTSV( D, E, B, INFO )'
      CALL PTSV( D, E, B, INFO )

      WRITE(*,*) 'Vector D on exit:'
      DO I = 1, N; WRITE (*,'(F8.5)') D(I);
      ENDDO

      WRITE(*,*) 'Vector E on exit:'
      DO I = 1, N-1; WRITE (*,'(E15.5)') E(I);
      ENDDO


      WRITE(*,*) 'The matrix B on exit:'
      DO J = 1, N; WRITE (*,*) B(J,:);
      ENDDO

      WRITE(*,*)'INFO = ',INFO

      DEALLOCATE(A, D, E, B)

      END PROGRAM SPTSV_MAIN
