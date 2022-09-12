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
*  CHEEV Example.
*  ==============
*
*  Program computes all eigenvalues and eigenvectors of a complex Hermitian
*  matrix A:
*
*  (  9.14,  0.00) ( -4.37, -9.22) ( -1.98, -1.72) ( -8.96, -9.50)
*  ( -4.37,  9.22) ( -3.35,  0.00) (  2.25, -9.51) (  2.57,  2.40)
*  ( -1.98,  1.72) (  2.25,  9.51) ( -4.82,  0.00) ( -3.24,  2.04)
*  ( -8.96,  9.50) (  2.57, -2.40) ( -3.24, -2.04) (  8.44,  0.00)
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
*
*  Example Program Results.
*  ========================
*
* CHEEV Example Program Results
* 
* Eigenvalues
* -16.00  -6.76   6.67  25.51
* 
* Eigenvectors (stored columnwise)
* (  0.34,  0.00) ( -0.55,  0.00) (  0.31,  0.00) ( -0.70,  0.00)
* (  0.44, -0.54) (  0.26,  0.18) (  0.45,  0.29) (  0.22, -0.28)
* ( -0.48, -0.37) ( -0.52, -0.02) ( -0.05,  0.57) (  0.15,  0.08)
* (  0.10, -0.12) ( -0.50,  0.28) ( -0.23, -0.48) (  0.34, -0.49)
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
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
*     RWORK dimension should be at least MAX(1,3*N-2)
      REAL             W( N ), RWORK( 3*N-2 )
      COMPLEX          A( LDA, N ), WORK( LWMAX )
      DATA             A/
     $ ( 9.14, 0.00),(-4.37, 9.22),(-1.98, 1.72),(-8.96, 9.50),
     $ ( 0.00, 0.00),(-3.35, 0.00),( 2.25, 9.51),( 2.57,-2.40),
     $ ( 0.00, 0.00),( 0.00, 0.00),(-4.82, 0.00),(-3.24,-2.04),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 8.44, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CHEEV
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CHEEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL CHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,
     $            INFO )
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
*     End of CHEEV Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX          A( LDA, * )
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
      REAL             A( LDA, * )
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
