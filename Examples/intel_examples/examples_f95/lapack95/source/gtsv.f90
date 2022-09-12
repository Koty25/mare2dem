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
!     S G T S V  Example Program Text
!*******************************************************************************

      PROGRAM SGTSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GTSV
!  .."Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, INFO, N, NRHS
      INTEGER, ALLOCATABLE :: IPIV(:)
!  .. Local Arrays ..
      REAL(WP), ALLOCATABLE :: DL(:), D(:), DU(:), B(:,:),X(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GTSV Example Program Results.'
      N = 6; NRHS = 3

      ALLOCATE(DL(N-1),D(N),DU(N-1),B(N,NRHS), X(N,NRHS), IPIV(N) )

      DL = 0
      D = 0
      DU = 0
          DO I=1,N-1
          READ(*,'(F2.0)') DL(I)
          ENDDO
          DO I=1,N
          READ(*,'(F2.0)') D(I)
          ENDDO
          DO I=1,N-1
          READ(*,'(F2.0)') DU(I)
          ENDDO

      WRITE(*,*)'DU :'
      WRITE(*,"(I3,1X)") INT(DU(:));

      WRITE(*,*)'D :'
      WRITE(*,"(I3,1X)") INT(D(:));

      WRITE(*,*)'DL :'
      WRITE(*,"(I3,1X)") INT(DL(:));

      DO I = 2, N-1; B(I,:) = DL(I-1) + D(I) + DU(I); ENDDO
      B(1,:) = D(1) + DU(1);B(N,:) = DL(N-1) + D(N)
      DO I = 1, NRHS; B(:,I) = B(:,I)*I; ENDDO

      WRITE(*,*) 'The RHS matrix B:'
      DO I=1,N; WRITE(*,"(3(I3,1X))") INT(B(I,:)); ENDDO

      WRITE(*,*) ' CALL GTSV(DL, D, DU, B, INFO)'
      CALL GTSV(DL, D, DU, B, INFO )

      WRITE(*,*)'DL on exit: '
      WRITE(*,"(F8.5)") DL(:);
      WRITE(*,*)'D on exit: '
      WRITE(*,"(F8.5)") D(:);
      WRITE(*,*)'DU on exit: '
      WRITE(*,"(F8.5)") DU(:);
      WRITE(*,*)'B on exit: '
      DO I=1,N; WRITE(*,"(6(F8.5))") B(I,:); ENDDO

      WRITE(*,*)'INFO = ',INFO

      DEALLOCATE(DL, D, DU, B, X, IPIV)

      END PROGRAM SGTSV_MAIN
