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
*  CHEEVX Example.
*  ==============
*
*  Program computes eigenvalues specified by a selected range of values
*  and corresponding eigenvectors of a complex Hermitian matrix A:
*
*  (  6.51,  0.00) ( -5.92,  9.53) ( -2.46,  2.91) (  8.84,  3.21)
*  ( -5.92, -9.53) ( -1.73,  0.00) (  6.50,  2.09) (  1.32,  8.81)
*  ( -2.46, -2.91) (  6.50, -2.09) (  6.90,  0.00) ( -0.59,  2.47)
*  (  8.84, -3.21) (  1.32, -8.81) ( -0.59, -2.47) ( -2.85,  0.00)
*
*  Description.
*  ============
*
*  The routine computes selected eigenvalues and, optionally, eigenvectors of
*  an n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies
*
*  A*v(j) = lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The computed eigenvectors are
*  orthonormal.
*  Eigenvalues and eigenvectors can be selected by specifying either a range
*  of values or a range of indices for the desired eigenvalues.
*
*  Example Program Results.
*  ========================
*
* CHEEVX Example Program Results
*
* The total number of eigenvalues found: 3
* 
* Selected eigenvalues
*   0.09   9.53  18.75
* 
* Selected eigenvectors (stored columnwise)
* (  0.18,  0.00) ( -0.54,  0.00) (  0.67,  0.00)
* ( -0.40, -0.31) ( -0.21, -0.17) ( -0.30, -0.43)
* (  0.60,  0.40) ( -0.35, -0.28) ( -0.39, -0.34)
* ( -0.34,  0.26) ( -0.57,  0.35) (  0.05,  0.05)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA, LDZ
      PARAMETER        ( LDA = N, LDZ = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, IL, IU, M
      REAL             ABSTOL, VL, VU
*
*     .. Local Arrays ..
*     IWORK dimension should be at least 5*N
      INTEGER          IWORK( 5*N ), IFAIL( N )
*     RWORK dimension should be at least 7*N
      REAL             W( N ), RWORK( 7*N )
      COMPLEX          A( LDA, N ), Z( LDZ, N ), WORK( LWMAX )
      DATA             A/
     $ ( 6.51, 0.00),(-5.92,-9.53),(-2.46,-2.91),( 8.84,-3.21),
     $ ( 0.00, 0.00),(-1.73, 0.00),( 6.50,-2.09),( 1.32,-8.81),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 6.90, 0.00),(-0.59,-2.47),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),(-2.85, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CHEEVX
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CHEEVX Example Program Results'
*     Negative ABSTOL means using the default value
      ABSTOL = -1.0
*     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
      VL = 0.0
      VU = 100.0
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CHEEVX( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK,
     $             IFAIL, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL CHEEVX( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK,
     $             IFAIL, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print the number of eigenvalues found.
*
      WRITE(*,'(/A,I2)')' The total number of eigenvalues found:', M
*
*     Print eigenvalues.
*
      CALL PRINT_RMATRIX( 'Selected eigenvalues', 1, M, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Selected eigenvectors (stored columnwise)',
     $                   N, M, Z, LDZ )
      STOP
      END
*
*     End of CHEEVX Example.
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
