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
!     C H P G V X  Example Program Text
!*******************************************************************************

      PROGRAM CHPGVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: HPGVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N, ITYPE, IL, IU, M, INFO
      REAL(WP), ALLOCATABLE :: W(:)
      REAL(WP) :: ABSTOL
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE :: IFAIL(:)
      COMPLEX(WP), ALLOCATABLE :: A(:), B(:), Z(:,:)
!  .. "Intrinsic Functions" ..
      INTRINSIC AIMAG, INT
!  .. "Executable Statements" ..
      WRITE (*,*) 'HPGVX Example Program Results'
      N = 5
      ALLOCATE( A(N*(N+1)/2), B(N*(N+1)/2), W(N), Z(N,N), IFAIL(N) )

      A = 0
      DO J=1,N*(N+1)/2
         READ(*,*) A(J)
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,1X,',',I3)')')") INT(A(I)), INT(AIMAG(A(I)))
      ENDDO


      B = 0
      DO J=1,N*(N+1)/2
         READ(*,*) B(J)
      ENDDO

      WRITE(*,*)'Matrix B:'
      DO I=1,N*(N+1)/2
         WRITE(*,"('('(I3,1X,',',I3)')')") INT(B(I)), INT(AIMAG(B(I)))
      ENDDO

      WRITE(*,*) "CALL HPGVX( A, B, W, 2, Z=Z, IL=4, IU=5, M=M, ", &
                 "IFAIL=IFAIL, ABSTOL=1.0E-3_WP, INFO =INFO ) "
      ITYPE=2; IL=4; IU=5; ABSTOL=1e-3
      CALL HPGVX( A, B, W, 2, Z=Z, IL=4, IU=5, M=M, &
                     IFAIL=IFAIL, ABSTOL=1.0E-3_WP, INFO = INFO )

      WRITE(*,*) 'W on exit:'
      DO I=1,N
         WRITE(*,"(5(E14.6,1X))") W(I)
      ENDDO

      WRITE(*,*) 'IFAIL on exit : ',IFAIL

      WRITE(*,*) 'M and INFO on exit:', M, INFO
      WRITE(*,*) 'Z on exit : '
      DO I=1,N
      WRITE(*,"(5(1H(,2(F12.10,1X),1H),1X))") (Z(I,J), J = 1,M)
      ENDDO

      DEALLOCATE(A, B, W, Z, IFAIL)

      END PROGRAM CHPGVX_MAIN
