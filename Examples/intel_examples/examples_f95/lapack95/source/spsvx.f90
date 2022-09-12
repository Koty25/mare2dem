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
!     S S P S V X  Example Program Text
!*******************************************************************************

      PROGRAM SSPSVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SPSVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, INFO, N, NN, NRHS
!  .. Local Arrays ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: B(:,:), AP(:), X(:,:), AFP(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SSPSVX Example Program Results.'
      N = 5; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      NN = N*(N+1)/2
      ALLOCATE ( AP(NN), B(N,NRHS), X(N,NRHS), AFP(NN), IPIV(N) )

      AP = 0
      DO I=1,NN
            READ(*,'(F3.0)') AP(I)
         ENDDO

      WRITE(*,*)'Matrix AP :'
      DO I=1,NN; WRITE(*,"(15(I3,1X,1X),I3,1X)") INT(AP(I));
      ENDDO

      B = 0
      DO I=1,N
         READ(*,'(F3.0)') B(I,1)
      ENDDO

      WRITE(*,*)"CALL SPSVX( AP, B, X, 'L', AFP, IPIV, INFO=INFO)"
      CALL SPSVX( AP, B, X, 'L', AFP, IPIV, INFO=INFO)

      WRITE(*,*)'AFP = '
      DO I=1,NN; WRITE(*,"(F8.5)") AFP(I);
      ENDDO
      WRITE(*,*)'IPIV = ', IPIV

      DEALLOCATE(AP, B, X, AFP, IPIV)

      END PROGRAM SSPSVX_MAIN
