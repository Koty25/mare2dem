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
!     S S T E V R  Example Program Text
!*******************************************************************************

      PROGRAM SSTEVR_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: STEVR
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: D(:), E(:), Z(:,:), W(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'STEVR Program Results'
      N = 5;
      ALLOCATE( D(N), E(N), Z(N,N), W(N) )

      D = 0
      DO J=1,N
         READ(*,*) D(J)
      ENDDO

      E = 0
      DO J=1,N
         READ(*,*) E(J)
      ENDDO

      WRITE(*,*)' Vector D : '
      DO I=1,N
      WRITE(*,"(5(I5,1X))") INT(D(I))
      ENDDO

      WRITE(*,*)' Vector E : '
      DO I=1,N
      WRITE(*,"(5(I5,1X))") INT(E(I))
      ENDDO

      WRITE(*,*) 'CALL STEVR( D, E, W, Z, -5.0_WP, 5.0_WP, M=M )'
      W = 0
      Z = 0
      M = 0
      CALL STEVR( D, E, W, Z, -5.0_WP, 5.0_WP, M=M )

      WRITE(*,*) ' W on exit :'
      DO I=1,N
      WRITE(*,"(5(F12.5,1X))") W(I)
      ENDDO

      WRITE(*,*)'Z on exit:'
      DO I=1,N
         WRITE(*,'(5(F9.5))') Z(I,:)
      ENDDO

      WRITE(*,*)'M = ', M

      DEALLOCATE(D, E, Z, W)

      END PROGRAM SSTEVR_MAIN
