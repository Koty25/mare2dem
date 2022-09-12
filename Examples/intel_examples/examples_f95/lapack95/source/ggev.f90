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
!     C G E E V  Example Program Text
!*******************************************************************************

      PROGRAM CGGEV_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GGEV
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, INFO, N
!  .. "Local Arrays" ..
      COMPLEX(WP), ALLOCATABLE :: A(:,:), AA(:,:), B(:,:), BB(:,:), ALPHA(:)
      COMPLEX(WP), ALLOCATABLE :: VL(:,:), VR(:,:), BETA(:)
      REAL(WP), ALLOCATABLE :: AR(:,:), BR(:,:), ALPHAR(:), ALPHAI(:), BETAR(:), &
                               VLR(:,:), VRR(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GGEV Example Program Results'
      N = 5
      ALLOCATE( A(N,N), AA(N,N), B(N,N), BB(N,N), ALPHA(N), BETA(N), VL(N,N), VR(N,N) )
      ALLOCATE( AR(N,N), BR(N,N), ALPHAR(N), ALPHAI(N), BETAR(N), VLR(N,N), VRR(N,N) )

      WRITE (*,*) 'Example 1, Real Example'
      AR = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) AR(I,J)
         ENDDO
      ENDDO
      WRITE(*,*)'Matrix AR : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(AR(I,:));
      ENDDO

      BR = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) BR(I,J)
         ENDDO
      ENDDO
      WRITE(*,*)'Matrix BR : '
      DO I=1,N;
         WRITE(*,"(5(I3,1X))") INT(BR(I,:));
      ENDDO
      WRITE(*,*)
      WRITE(*,*) 'CALL GGEV( A, B, ALPHAR, ALPHAI, BETAR, VLR, VRR )'
      CALL GGEV( AR, BR, ALPHAR, ALPHAI, BETAR, VLR, VRR )

      WRITE(*,*); WRITE(*,*)'ALPHAR on exit : '; WRITE(*,*) ALPHAR(1:N)
      WRITE(*,*); WRITE(*,*)'ALPHAI on exit : '; WRITE(*,*) ALPHAI(1:N)
      WRITE(*,*); WRITE(*,*)'BETA on exit : '; WRITE(*,*) BETAR(1:N)
      WRITE(*,*)'Array VL:'; DO I =1,N; WRITE(*,*)I, VLR(I,1:N); ENDDO
      WRITE(*,*)'Array VR:'; DO I =1,N; WRITE(*,*)I, VRR(I,1:N); ENDDO

      WRITE(*,*)
      WRITE(*,*)' Generalized eigenvalues : '
      DO I=1,N
         WRITE(*,*) '(',ALPHAR(I)/BETAR(I),',',ALPHAI(I)/BETAR(I),')'
      ENDDO

      WRITE (*,*) 'Example 2, Complex Example'
      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      AA=A

      WRITE(*,*)'Matrix A : '
      DO I=1,N
         WRITE(*,"(5('('(I3,1X,',',I3)')',1X,1X))") INT(A(I,:)), INT(AIMAG(A(I,:)))
      ENDDO

      DO J=1,N
         DO I=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      BB=B

      WRITE(*,*)'Matrix B : '
      DO I=1,N
         WRITE(*,"(5('('(I3,1X,',',I3)')',1X,1X))") INT(B(I,:)), INT(AIMAG(B(I,:)))
      ENDDO
      WRITE(*,*)
      WRITE(*,*) 'CALL GGEV( A, B, ALPHA, BETA, INFO=INFO )'
      CALL GGEV( A, B, ALPHA, BETA,  INFO=INFO )

      WRITE(*,*)
      WRITE(*,*)'ALPHA on exit : '
      DO I=1,N
         WRITE(*,*) ALPHA(I)
      ENDDO

      WRITE(*,*)
      WRITE(*,*)'BETA on exit : '
      DO I=1,N
         WRITE(*,*) BETA(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'INFO = ', INFO

      WRITE(*,*)
      WRITE(*,*)' Generalized eigenvalues : '
      DO I=1,N
         WRITE(*,*)  ALPHA(I)/BETA(I)
      ENDDO

      WRITE(*,*)
      WRITE(*,*) 'CALL GGEV( A, B, ALPHAR, ALPHAI, BETA, VL, VR )'
      CALL GGEV( AA, BB, ALPHA, BETA, VL, VR )

      WRITE(*,*)
      WRITE(*,*)'Matrix VL on exit : '
      DO I=1,N;
         WRITE(*,"(5('('(F8.5,1X,',',F8.5)')'))") REAL(VL(I,:)), AIMAG(VL(I,:))
      ENDDO

      WRITE(*,*)
      WRITE(*,*)'Matrix VR on exit : '
      DO I=1,N;
         WRITE(*,"(5('('(F8.5,1X,',',F8.5)')'))") REAL(VR(I,:)), AIMAG(VR(I,:))
      ENDDO

      DEALLOCATE(A, AA, B, BB, ALPHA, BETA, VL, VR)
      DEALLOCATE(AR, BR, ALPHAR, ALPHAI, BETAR, VLR, VRR)

      END PROGRAM CGGEV_MAIN
