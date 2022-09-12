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

      PROGRAM SSYTRI2_EXAMPLE

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYTRI2, SYTRI, SYTRF
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, K, N, LDA, INFO
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: A1(:,:) , A2(:,:), B(:,:)
!  .. "Executable Statements" ..
       WRITE (*,*) 'SYTRI2 Example Program Results.'
       N = 5; LDA = N
       ALLOCATE ( A1(LDA,N), A2(LDA,N), B(N,N), IPIV(N))
       A1=0.0; B=0.0

      A1 = 0
      DO I=1,N
         DO J=I,N
         READ(*,'(F3.0)') A1(I,J)
         ENDDO
      ENDDO
       
       CALL SYTRF(A1,'U',IPIV,INFO )
       WRITE(*,*)'INFO on exit from SYTRF: ', INFO

       A2 = A1

       WRITE(*,*) 'The matrix A1:'
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") A1(I,:)
       ENDDO

       WRITE(*,*) 'CALL SYTRI2( A1, IPIV, U, INFO )'

       CALL SYTRI2(A1,IPIV,'U',INFO)
       CALL SYTRI(A2,IPIV,'U',INFO)

       WRITE(*,*)'A1 on exit from SYTRI2: '
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") A1(I,:); ENDDO
       
       WRITE(*,*)'A2 on exit from SYTRI: '
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") A2(I,:); ENDDO


       do I=1,N
          do j=1,N
             B(I,J) = A2(I,J) - A1(I,J)
          enddo
       enddo

       WRITE(*,*)'A1- A2 is: '
       DO I=1,N; WRITE(*,"(5(E12.6,1X),E12.6,1X)") B(I,:); ENDDO

       WRITE(*,*)'IPIV on exit: ', IPIV
       WRITE(*,*)'INFO on exit: ', INFO

       DEALLOCATE(A1,A2,B,IPIV)

       END PROGRAM SSYTRI2_EXAMPLE
