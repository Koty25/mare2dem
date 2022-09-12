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
!     S G E E V X  Example Program Text
!*******************************************************************************

      PROGRAM SGEEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GEEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, ILO, IHI
      REAL(WP) :: ABNRM
      REAL(WP), ALLOCATABLE :: WR(:), WI(:)
!  .. "Local Arrays" ..
      REAL(WP), ALLOCATABLE :: A(:,:), VL(:,:), VR(:,:), SCALE(:)
      REAL(WP), ALLOCATABLE :: RCONDE(:), RCONDV(:)
!  .. "Intrinsic Functions" ..
      INTRINSIC REAL, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'GEEVX Example Program Results'
      N = 5
      ALLOCATE( A(N,N), WR(N), WI(N), VL(N,N), VR(N,N), SCALE(N) )
      ALLOCATE( RCONDE(N), RCONDV(N) )

      A = 0
      DO I=1,N
         DO J=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N
         WRITE(*,"(5(I3,1X))") INT(A(I,:))
      ENDDO

      WRITE(*,*) "CALL GEEVX( A, WR, WI, 'B', ILO, IHI, SCALE,", &
                               " ABNRM, RCONDE, RCONDV )"
      CALL GEEVX( A, WR, WI, BALANC='B', ILO=ILO, IHI=IHI, &
           SCALE=SCALE, ABNRM=ABNRM, RCONDE=RCONDE, RCONDV=RCONDV )

      WRITE(*,*) "WR on exit : "
      DO I=1,N
         WRITE(*,"((E14.6,1X))") REAL(WR(I))
      ENDDO

      WRITE(*,*) "WI on exit : "
      DO I=1,N
         WRITE(*,"((E14.6,1X))") REAL(WI(I))
      ENDDO

      WRITE(*,*)'ILO : ',ILO
      WRITE(*,*)'IHI : ',IHI
      WRITE(*,*)'SCALE : ',SCALE
      WRITE(*,*)'ABNRM on exit :', ABNRM
      WRITE(*,*)'RCONDE on exit :', RCONDE
      WRITE(*,*)'RCONDV on exit :', RCONDV

      DEALLOCATE(A, WR, WI, VL, VR, SCALE, RCONDE, RCONDV)

      END PROGRAM SGEEVX_MAIN
