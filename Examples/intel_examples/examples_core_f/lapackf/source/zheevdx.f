!===============================================================================
! Copyright 2009-2020 Intel Corporation.
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

*
*  ZHEEVD Example.
*  ==============
*
*  Program computes all eigenvalues and eigenvectors of a complex Hermitian
*  matrix A using divide and conquer algorithm, where A is:
*
*  (  3.40,  0.00) ( -2.36, -1.93) ( -4.68,  9.55) (  5.37, -1.23)
*  ( -2.36,  1.93) (  6.94,  0.00) (  8.13, -1.47) (  2.07, -5.78)
*  ( -4.68, -9.55) (  8.13,  1.47) ( -2.14,  0.00) (  4.68,  7.44)
*  (  5.37,  1.23) (  2.07,  5.78) (  4.68, -7.44) ( -7.42,  0.00)
*
*  Description.
*  ============
*
*  The routine computes all eigenvalues and, optionally, eigenvectors of an
*  n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies
*
*  A*v(j) = lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The computed eigenvectors are
*  orthonormal.
*  If the eigenvectors are requested, then this routine uses a divide and
*  conquer algorithm to compute eigenvalues and eigenvectors.
*
*  Example Program Results.
*  ========================
*
* ZHEEVD Example Program Results
* 
* Eigenvalues
* -21.97  -0.05   6.46  16.34
* 
* Eigenvectors (stored columnwise)
* (  0.41,  0.00) ( -0.34,  0.00) ( -0.69,  0.00) (  0.49,  0.00)
* (  0.02, -0.30) (  0.32, -0.21) ( -0.57, -0.22) ( -0.59, -0.21)
* (  0.18,  0.57) ( -0.42, -0.32) (  0.06,  0.16) ( -0.35, -0.47)
* ( -0.62, -0.09) ( -0.58,  0.35) ( -0.15, -0.31) ( -0.10, -0.12)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, LIWORK, LRWORK
*
*     .. Local Arrays ..
      INTEGER          IWORK( LWMAX )
      DOUBLE PRECISION W( N ), RWORK( LWMAX )
      COMPLEX*16       A( LDA, N ), WORK( LWMAX )
      DATA             A/
     $ ( 3.40, 0.00),(-2.36, 1.93),(-4.68,-9.55),( 5.37, 1.23),
     $ ( 0.00, 0.00),( 6.94, 0.00),( 8.13, 1.47),( 2.07, 5.78),
     $ ( 0.00, 0.00),( 0.00, 0.00),(-2.14, 0.00),( 4.68,-7.44),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),(-7.42, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZHEEVD
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'ZHEEVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      LIWORK = -1
      LRWORK = -1
      CALL ZHEEVD( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL ZHEEVD( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print eigenvalues.
*
      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A, 
     $                   LDA )
      STOP
      END
*
*     End of ZHEEVD Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing a real matrix.
*
      SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END
