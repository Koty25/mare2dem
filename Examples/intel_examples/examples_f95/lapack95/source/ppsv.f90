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
!     S P P S V  Example Program Text
!*******************************************************************************

      PROGRAM SPPSV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: PPSV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NN, NRHS
!  .. Local Arrays ..
      REAL(WP), ALLOCATABLE :: AP(:) ,B(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SPPSV Example Program Results.'
      N = 5; NRHS = 1
      NN = N*(N+1)/2
      ALLOCATE ( AP(NN), B(N,NRHS) )

      AP = 0
      DO I=1,NN
         READ(*,'(F3.0)') AP(I);
      ENDDO;

      WRITE(*,*) 'The array AP :'
      DO I=1,NN
      WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(AP(I));
      ENDDO

      B = 0
      DO I=1,N
      DO J=1,NRHS
         READ(*,'(F3.0)') B(I,J);
      ENDDO;
      ENDDO;

      WRITE(*,*) 'The array B :'
      DO I=1,N
       WRITE(*,'(10(I3,1X,1X),I3,1X)') INT(B(I,:));
      ENDDO;

      WRITE(*,*) "CALL PPSV( AP, B, 'L' )"
      CALL PPSV(  AP, B, 'L' )

      WRITE(*,*)'AP on exit :'
      DO I=1,NN; WRITE(*,"(5(E13.6))") AP(I);
      ENDDO

      WRITE(*,*)'B on exit :'
      DO I=1,N; WRITE(*,"(5(E13.6))") B(I,:);
      ENDDO

      DEALLOCATE(AP, B)

      END PROGRAM SPPSV_MAIN
