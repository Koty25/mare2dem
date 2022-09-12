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
      PROGRAM HETRS2_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HETRS2, HETRS, HETRF

      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, LDA, LDB, NRHS, INFO, LWORK
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      COMPLEX(WP), ALLOCATABLE :: A1(:,:), A2(:,:), B1(:,:), B2(:,:), C(:,:), WORK(:)

!  .. "Executable Statements" ..
       WRITE (*,*) 'HETRS2 Example Program Results.'
      N = 5; NRHS = 3
      WRITE(*,'(5H N = , I4, 9H; NRHS = , I4)') N, NRHS
      LDA = N; LWORK =5000; LDB=N
      ALLOCATE ( A1(LDA,N), A2(LDA,N), B1(LDB,NRHS), B2(LDB,NRHS), C(LDB,NRHS), IPIV(N), WORK(LWORK))
      A1=0.0; A2=0.0; B1=0.0; B2=0.0; C=0.0; IPIV=0;

       A1 = 0
       DO I=1,N
          DO J=I,N
          READ(*,*) A1(I,J);
          ENDDO
       ENDDO;

       B1 = 0
       DO I=1,N
         READ(*,*) B1(I,1);
       ENDDO;
       B2=B1

       CALL HETRF( A1, 'U', IPIV, INFO )
       WRITE(*,*)'INFO on exit from HETRF: ', INFO
       A2=A1

       WRITE(*,*) 'CALL HETRS( AB, B, IPIV, L, INFO)'
       CALL HETRS( A1, B1, IPIV, 'U', INFO )

       WRITE(*,*) 'CALL HETRS2( AB, B, IPIV, L, INFO)'
       CALL HETRS2( A2, B2, IPIV, 'U', INFO )

       WRITE(*,*)'B on exit from HETRS: '
       DO I=1,LDB; WRITE(*,"(5(E12.6,1X),E12.6,1X)") B1(I,:); ENDDO
       WRITE(*,*)'B on exit from HETRS2: '
       DO I=1,LDB; WRITE(*,"(1(E12.6,1X),E12.6,1X)") B2(I,:); ENDDO

       DO I=1,LDB
          DO J=1,NRHS
             C(I,J)= B1(I,J) - B2(I,J)
          ENDDO
       ENDDO

       WRITE(*,*)'B1-B2 on exit: '
       DO I=1,LDB; WRITE(*,"(1(E12.6,1X),E12.6,1X)") C(I,:); ENDDO

       WRITE(*,*)'IPIV on exit: ', IPIV
       WRITE(*,*)'INFO on exit: ', INFO

       END PROGRAM HETRS2_MAIN
