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
      PROGRAM SSYCONV_EXAMPLE

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYCONV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, LDA, INFO
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IPIV(:)
      REAL(WP), ALLOCATABLE :: AB(:,:)
      REAL(WP), ALLOCATABLE :: E(:)

!  .. "Executable Statements" ..
       WRITE (*,*) 'SYCONV Example Program Results.'
       N = 5; LDA = N
       ALLOCATE ( AB(LDA,N), IPIV(N), E(N) )
       AB=0.0; IPIV=0; E=0.0

       DO J=1,N
          IPIV(J) = - J
       ENDDO

      AB = 0
      DO I=1,N
         DO J=I,N
         READ(*,'(F3.0)') AB(I,J)
         ENDDO
      ENDDO

       WRITE(*,*) 'The matrix AB:'
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") AB(I,:)
       ENDDO

       WRITE(*,*) 'CALL SYCONV( AB, IPIV, E, L, C, INFO )'

       CALL SYCONV( AB, IPIV, E, 'U', 'C', INFO )
!       CALL SSYCONV( 'L','C',N,AB,LDA, IPIV, WORK, INFO )

       WRITE(*,*)'AB on exit: '
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") AB(I,:); ENDDO
       WRITE(*,*)'IPIV on exit: ', IPIV
       WRITE(*,*)'E on exit:'; WRITE(*,"(8(F8.5,1X))") E(:);
       WRITE(*,*)'INFO on exit: ', INFO

       DEALLOCATE(AB, IPIV, E)

       END PROGRAM SSYCONV_EXAMPLE
