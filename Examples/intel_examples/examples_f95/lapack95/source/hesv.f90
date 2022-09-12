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
!     S H E S V  Example Program Text
!*******************************************************************************

      PROGRAM SHESV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HESV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, NRHS
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      COMPLEX(WP), ALLOCATABLE :: A(:,:), B(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SHESV Example Program Results'
      N = 5; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      ALLOCATE ( A(N,N), B(N,NRHS), IPIV(N) )

      A = 0
      DO I=1,N
         DO J=I,N
         READ(*,*) A(I,J);
         ENDDO
      ENDDO;

      WRITE(*,*) 'The matrix A :'
      DO I=1,N
      WRITE(*,"(5(I3,1X,'+',1X,I3,'i',1X))") INT(A(I,1)),INT(AIMAG(A(I,1))), &
                                             INT(A(I,2)),INT(AIMAG(A(I,2))), &
                                             INT(A(I,3)),INT(AIMAG(A(I,3))), &
                                             INT(A(I,4)),INT(AIMAG(A(I,4))), &
                                             INT(A(I,5)),INT(AIMAG(A(I,5)))
      ENDDO

      B = 0
      DO I=1,N
         READ(*,*) B(I,1);
      ENDDO;

      WRITE(*,*) 'The array B :'
      DO I=1,N
      WRITE(*,"((I3,1X,'+',1X,I3,'i',1X))") INT(B(I,1)),INT(AIMAG(B(I,1)))
      ENDDO

      WRITE(*,*)' CALL HESV( A, B, IPIV=IPIV )'
      CALL HESV( A, B, IPIV=IPIV )

      WRITE(*,*)'A on exit: '
      DO I=1,N; WRITE(*,*) A(I,:);
      ENDDO

      WRITE(*,*)'B on exit: '
      DO I=1,N; WRITE(*,*) (B(I,1));
      ENDDO

      WRITE(*,*)'IPIV on exit : ', IPIV

      DEALLOCATE(A, B, IPIV)

      END PROGRAM SHESV_MAIN
