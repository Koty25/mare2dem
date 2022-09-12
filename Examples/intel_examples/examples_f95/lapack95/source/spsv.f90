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
!     S S P S V  Example Program Text
!*******************************************************************************

      PROGRAM SSPSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SPSV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, N, NN, NRHS
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: B(:,:), AP(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SSPSV Example Program Results.'
      N = 5; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      NN = N*(N+1)/2
      ALLOCATE ( AP(NN), B(N,NRHS), IPIV(N) )

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

      WRITE(*,*)'Matrix B :'
      DO I=1,N; WRITE(*,"(10(I3,1X,1X),I3,1X)") INT(B(I,1));
      ENDDO

      WRITE(*,*)" CALL SPSV( AP, B, 'L', IPIV )"

      CALL SPSV( AP, B, 'L', IPIV )

      WRITE(*,*)'AP on exit: '
      DO I=1,NN; WRITE(*,"(15(E13.5))") AP(I);
      ENDDO

      WRITE(*,*)'Matrix B on exit :'
      DO I=1,N; WRITE(*,"(F9.5)") B(I,1);
      ENDDO
      WRITE(*,*)'IPIV = ', IPIV

      DEALLOCATE(AP, B, IPIV)

      END PROGRAM SSPSV_MAIN
