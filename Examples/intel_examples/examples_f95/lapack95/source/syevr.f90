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
!     S S Y E V R  Example Program Text
!*******************************************************************************

      PROGRAM SSYEVR_MAIN

!  -- LAPACK95 EXAMPLE DRIVER ROUTINE (VERSION 1.0) --
!     UNI-C, DENMARK
!     DECEMBER, 1999

!  .. "Use Statements" ..
      USE f95_precision, ONLY: WP => SP
      USE lapack95, ONLY: SYEVR
!  .. "Implicit Statement" ..
      IMPLICIT NONE
!  .. "Local Scalars" ..
      INTEGER :: I, J, M, N, IL, IU
!  .. "Local Arrays" ..
      INTEGER, ALLOCATABLE ::  ISUPPZ(:)
      REAL(WP), ALLOCATABLE :: A(:,:),W(:), Z(:,:)
!  .. "Executable Statements" ..
      WRITE (*,*) 'SYEVR Example Program Results'
      N = 5
      ALLOCATE( A(N,N), W(N), Z(N,N), ISUPPZ(2*N)  )

      A = 0
      DO J=1,N
         DO I=1,N
            READ(*,*) A(I,J)
         ENDDO
      ENDDO

      WRITE(*,*)'Matrix A:'
      DO I=1,N;
         WRITE(*,"(5(F9.5))") A(I,:);
      ENDDO

      WRITE(*,*) 'CALL SYEVR( A, W, Z=Z, ISUPPZ=ISUPPZ, IL=1, IU=2, M=M )'
      IL=1; IU=2

      CALL SYEVR( A, W, IL=1, IU=2, M=M )

      WRITE(*,*)'Matrix A on exit :'
      DO I=1,N;
         WRITE(*,"(5(F9.5))") A(I,:);
      ENDDO


      WRITE(*,*)'M = ', M
      WRITE(*,*) 'W on exit :'
      DO I=1,M
         WRITE(*,"(5(F9.5))") W(I)
      ENDDO

      WRITE(*,*)'Z on exit :'
      DO I=1,N
         WRITE(*,"(5(E14.6))") Z(I,1:M)
      ENDDO
      WRITE(*,*)'ISUPPZ on exit:'
      WRITE(*,*) ISUPPZ(1:2*max(1,m))

      DEALLOCATE(A, W, Z, ISUPPZ)

      END PROGRAM SSYEVR_MAIN
