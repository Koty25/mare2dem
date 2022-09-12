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
      PROGRAM SGESV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999
!
!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GESV, GETRS, MKL_GETRFNPI
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N, NRHS
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), BB(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SGESV Example Program Results.'
      READ (*, * ) N
      READ (*, * ) NRHS
      ALLOCATE( A(N,N), AA(N,N), B(N,NRHS), BB(N,NRHS), IPIV(N) )
      DO J=1,N
      DO I=1,N
         READ(*,'(F2.0)') AA(I,J)
      ENDDO
      IPIV(J)=J
      ENDDO

      DO J = 1, NRHS; BB(:,J) = SUM( AA, DIM=2)*J; ENDDO

      WRITE(*,*) 'The matrix A:'
      DO I=1,N; WRITE(*,"(4(I3,1X),I3,1X)") INT(AA(I,:)); ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO I=1,N; WRITE(*,"(2(I3,1X),I3,1X)") INT(BB(I,:)); ENDDO

      WRITE(*,*) 'CALL GESV( A, B )'
      A=AA; B=BB

      CALL MKL_GETRFNPI(A)
      CALL GETRS(A, IPIV, B)

      WRITE(*,*) 'B - the solution vectors computed by GESV'
      DO I=1,N; WRITE(*,"(2(E12.6,1X),E12.6,1X)") B(I,:); ENDDO

      WRITE(*,*) 'CALL GESV( A, B(:,1), IPIV, INFO )'

      CALL GESV(  AA, BB(:,1), IPIV, INFO )

      WRITE(*,*) ' A on exit:'
      DO I=1,N; WRITE(*,"(4(E12.6,1X),E12.6,1X)") AA(I,:); ENDDO

      WRITE(*,*) 'B on exit:'
      DO I=1,N; WRITE(*,"(4(E12.6,1X),E12.6,1X)") BB(I,1); ENDDO

      WRITE(*,*)'IPIV on exit:', IPIV

      WRITE(*,*)'INFO on exit:', INFO

      END PROGRAM SGESV_MAIN









