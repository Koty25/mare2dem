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
!     S S T E V X  Example Program Text
!*******************************************************************************

      PROGRAM SSTEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: STEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, IL, IU
      REAL(WP) :: ABSTOL
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: D(:), E(:), Z(:,:), W(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'STEV Example Program Results'
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

      WRITE(*,*) 'CALL STEVX( D, E, W, IL=1, IU=3, ABSTOL=1E-2 )'
      IL=1; IU=3; ABSTOL=1E-2
      W = 0
      CALL STEVX( D, E, W, IL=IL, IU=IU, ABSTOL=ABSTOL )

      WRITE(*,*) ' W on exit : '
      DO I=1,N
      WRITE(*,"(5(F12.5,1X))") W(I)
      ENDDO

      DEALLOCATE(D, E, Z, W)

      END PROGRAM SSTEVX_MAIN
