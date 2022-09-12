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
*  ZGEEV Example.
*  ==============
*
*  Program computes the eigenvalues and left and right eigenvectors of a general
*  rectangular matrix A:
*
*  ( -3.84,  2.25) ( -8.94, -4.75) (  8.95, -6.53) ( -9.87,  4.82)
*  ( -0.66,  0.83) ( -4.40, -3.82) ( -3.50, -4.26) ( -3.15,  7.36)
*  ( -3.99, -4.73) ( -5.88, -6.60) ( -3.36, -0.40) ( -0.75,  5.23)
*  (  7.74,  4.18) (  3.66, -7.53) (  2.58,  3.60) (  4.59,  5.41)
*
*  Description.
*  ============
*
*  The routine computes for an n-by-n complex nonsymmetric matrix A, the
*  eigenvalues and, optionally, the left and/or right eigenvectors. The right
*  eigenvector v(j) of A satisfies
*
*  A*v(j)= lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
*
*  u(j)H*A = lambda(j)*u(j)H
*
*  where u(j)H denotes the conjugate transpose of u(j). The computed
*  eigenvectors are normalized to have Euclidean norm equal to 1 and
*  largest component real.
*
*  Example Program Results.
*  ========================
*
* ZGEEV Example Program Results
* 
* Eigenvalues
* ( -9.43,-12.98) ( -3.44, 12.69) (  0.11, -3.40) (  5.76,  7.13)
* 
* Left eigenvectors
* (  0.24, -0.18) (  0.61,  0.00) ( -0.18, -0.33) (  0.28,  0.09)
* (  0.79,  0.00) ( -0.05, -0.27) (  0.82,  0.00) ( -0.55,  0.16)
* (  0.22, -0.27) ( -0.21,  0.53) ( -0.37,  0.15) (  0.45,  0.09)
* ( -0.02,  0.41) (  0.40, -0.24) (  0.06,  0.12) (  0.62,  0.00)
* 
* Right eigenvectors
* (  0.43,  0.33) (  0.83,  0.00) (  0.60,  0.00) ( -0.31,  0.03)
* (  0.51, -0.03) (  0.08, -0.25) ( -0.40, -0.20) (  0.04,  0.34)
* (  0.62,  0.00) ( -0.25,  0.28) ( -0.09, -0.48) (  0.36,  0.06)
* ( -0.23,  0.11) ( -0.10, -0.32) ( -0.43,  0.13) (  0.81,  0.00)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA, LDVL, LDVR
      PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
*     RWORK dimension should be at least 2*N
      DOUBLE PRECISION RWORK( 2*N )
      COMPLEX*16       A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ), 
     $                 W( N ), WORK( LWMAX )
      DATA             A/
     $ (-3.84, 2.25),(-0.66, 0.83),(-3.99,-4.73),( 7.74, 4.18),
     $ (-8.94,-4.75),(-4.40,-3.82),(-5.88,-6.60),( 3.66,-7.53),
     $ ( 8.95,-6.53),(-3.50,-4.26),(-3.36,-0.40),( 2.58, 3.60),
     $ (-9.87, 4.82),(-3.15, 7.36),(-0.75, 5.23),( 4.59, 5.41)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZGEEV
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'ZGEEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL ZGEEV( 'Vectors', 'Vectors', N, A, LDA, W, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, RWORK, INFO )
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
      CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print left eigenvectors.
*
      CALL PRINT_MATRIX( 'Left eigenvectors', N, N, VL, LDVL )
*
*     Print right eigenvectors.
*
      CALL PRINT_MATRIX( 'Right eigenvectors', N, N, VR, LDVR )
      STOP
      END
*
*     End of ZGEEV Example.
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
