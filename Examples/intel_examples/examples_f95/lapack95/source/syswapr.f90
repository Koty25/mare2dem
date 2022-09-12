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
      PROGRAM SYSWAPR_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements"
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYSWAPR
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, LDA, I1, I2
!  .. "Local Arrays" ..

      REAL(WP), ALLOCATABLE :: AB(:,:)
!  .. "Executable Statements" ..
       WRITE (*,*) 'SYSWAPR Example Program Results.'
       N = 5; LDA = N
       I1 = 2; I2 = 5

       ALLOCATE ( AB(LDA,N))
       AB=0.0
       DO J=1,N
          DO I=1,N
             READ(*,*) AB(I,J)
          ENDDO
       ENDDO


       WRITE(*,*) 'The matrix AB:'
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") AB(I,:)
       ENDDO
       WRITE (*,*) "I1 = ",I1,"I2 = ",I2

       WRITE(*,*) 'CALL SYSWAPR( AB, I1, I2 )'

       CALL SYSWAPR( AB, I1, I2 )

       WRITE(*,*)'AB on exit: '
       DO I=1,LDA; WRITE(*,"(5(E12.6,1X),E12.6,1X)") AB(I,:); ENDDO
       DEALLOCATE( AB )

       END PROGRAM SYSWAPR_MAIN






