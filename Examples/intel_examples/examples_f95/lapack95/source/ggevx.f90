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
!     C G G E V X  Example Program Text
!*******************************************************************************

      PROGRAM CGGEVX_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: GGEVX
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, N
      REAL :: ABNRM, BBNRM
!  .. "Local Arrays" ..
      CHARACTER(LEN=*), PARAMETER :: FMT = '(4(1X,1H(,F7.3,1H,,F7.3,1H):))'
      CHARACTER(LEN=1) :: BALANC
      COMPLEX(WP), ALLOCATABLE :: A(:,:), B(:,:), ALPHA(:)
      COMPLEX(WP), ALLOCATABLE :: BETA(:), LAMBDA(:)
      REAL(WP), ALLOCATABLE :: LSCALE(:), RSCALE(:), RCONDE(:), RCONDV(:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'GGEVX Example Program Results'
      N = 5
      ALLOCATE( A(N,N), B(N,N), ALPHA(N), BETA(N), LSCALE(N), &
&       RSCALE(N), RCONDE(N), RCONDV(N), LAMBDA(N) )

      A = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A : '
      DO I=1,N
         WRITE(*,"(5('('(I3,1X,',',I3)')',1X,1X))") INT(A(I,:)), INT(AIMAG(A(I,:)))
      ENDDO

      B = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) B(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix B : '
      DO I=1,N
         WRITE(*,"(5('('(I3,1X,',',I3)')',1X,1X))") INT(B(I,:)), INT(AIMAG(B(I,:)))
      ENDDO

      WRITE(*,*)
      WRITE(*,*) "CALL GGEVX( A, B, BALANC='B', LSCALE=LSCALE, RSCALE=RSCALE,"
      WRITE(*,*)"ABNRM=ABNRM, BBNRM=BBNRM, RCONDE=RCONDE, RCONDV=RCONDV )"
      BALANC = 'B'
      CALL GGEVX( A, B, ALPHA, BETA, BALANC=BALANC, LSCALE=LSCALE, RSCALE=RSCALE, ABNRM=ABNRM, &
           BBNRM=BBNRM, RCONDE=RCONDE, RCONDV=RCONDV )

      WRITE(*,*)
      WRITE(*,*)'LSCALE : '
      DO I=1,N
         WRITE(*,'(F9.5)') LSCALE(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'RSCALE : '
      DO I=1,N
         WRITE(*,'(F9.5)') RSCALE(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'ABNRM = ', ABNRM
      WRITE(*,*)
      WRITE(*,*)'BBNRM = ', BBNRM
      WRITE(*,*)
      WRITE(*,*)'RCONDE : '
      DO I=1,N
         WRITE(*,'(F9.5)') RCONDE(I)
      ENDDO
      WRITE(*,*)
      WRITE(*,*)'RCONDV : '
      DO I=1,N
         WRITE(*,'(F9.5)') RCONDV(I)
      ENDDO

      WRITE(*,*)'Matrix ALPHA : '
      WRITE(*,FMT) ALPHA

      WRITE(*,*)'Matrix BETA : '
      WRITE(*,FMT) BETA
      lambda = alpha/beta
      print *,'lambda = ', LAMBDA
      WRITE(*,FMT) LAMBDA

      DEALLOCATE(A, ALPHA, BETA, LSCALE)
      DEALLOCATE(RSCALE, RCONDE, RCONDV, LAMBDA)

      END PROGRAM CGGEVX_MAIN
