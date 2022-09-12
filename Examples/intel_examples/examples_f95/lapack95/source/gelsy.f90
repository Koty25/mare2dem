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
!     S G E L S Y  Example Program Text
!*******************************************************************************

      PROGRAM SGELSY_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GELSY
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N, NRHS, RANK
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: JPVT(:)
      REAL(WP), ALLOCATABLE :: A(:,:), B(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GELSY Example Program Results'
      M=6; N = 4; NRHS = 3
      WRITE(*,'(5H M = , I4, 5H N = , I4, 9H; NRHS = , I4)') M, N, NRHS
      ALLOCATE( A(M,N), B(MAX(M,N),NRHS), JPVT(N) )

      A = 0
      B = 0
      DO I=1,M
         DO J=1,N
            READ(*,'(F2.0)') A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A :'
      DO I=1,M;
         WRITE(*,"(24(F9.5))") A(I,:);
      ENDDO

      DO I=1,MAX(M,N)
         DO J=1,NRHS
            READ(*,'(F3.0)') B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix B :'
      DO I=1,M;
         WRITE(*,"(3(F9.5))") B(I,:);
      ENDDO

      JPVT(1)=0; JPVT(2)=0; JPVT(3)=1; JPVT(4)=0;

      WRITE(*,*) 'CALL GELSY( A, B, RANK, JPVT, 1.0E-5_WP )'

      CALL GELSY( A, B, RANK, JPVT, 1.0E-5_WP )

      WRITE(*,*) 'The matrix B on exit :'
      DO I=1,M;
         WRITE(*,"(3(E15.7),1X)") B(I,:);
      ENDDO


      WRITE(*,*)'RANK = ', RANK
      WRITE(*,*)'JPVT = ', JPVT

      DEALLOCATE(A, B, JPVT)

      END PROGRAM SGELSY_MAIN
