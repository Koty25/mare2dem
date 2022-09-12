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
*  CSYSV Example.
*  ==============
*
*  The program computes the solution to the system of linear equations
*  with a complex symmetric matrix A and multiple right-hand sides B,
*  where A is the coefficient matrix:
*
*  (  9.99, -4.73) ( -5.68, -0.80) ( -8.94,  1.32) ( -9.42,  2.05)
*  ( -5.68, -0.80) ( -8.01,  4.61) (  1.64, -6.29) (  6.79, -2.17)
*  ( -8.94,  1.32) (  1.64, -6.29) (  9.04,  3.96) ( -4.51, -7.54)
*  ( -9.42,  2.05) (  6.79, -2.17) ( -4.51, -7.54) (  0.40,  4.06)
*
*  and B is the right-hand side matrix:
*
*  (  5.71, -1.20) (  2.84, -0.18)
*  ( -7.70,  6.47) ( -8.29, -1.72)
*  (  3.77, -7.40) ( -4.28, -8.25)
*  ( -3.78,  0.33) ( -2.70, -0.39)
*
*  Description.
*  ============
*
*  The routine solves for X the complex system of linear equations A*X = B,
*  where A is an n-by-n symmetric matrix, the columns of matrix B are
*  individual right-hand sides, and the columns of X are the corresponding
*  solutions.
*
*  The diagonal pivoting method is used to factor A as A = U*D*UT or
*  A = L*D*LT , where U (or L) is a product of permutation and unit upper
*  (lower) triangular matrices, and D is symmetric and block diagonal with
*  1-by-1 and 2-by-2 diagonal blocks.
*
*  The factored form of A is then used to solve the system of equations A*X = B.
*
*  Example Program Results.
*  ========================
*
* CSYSV Example Program Results
* 
* Solution
* (  0.13,  0.13) (  0.63,  0.34)
* (  0.32, -0.07) (  0.61,  0.21)
* ( -0.26, -0.44) ( -0.01, -0.10)
* ( -0.40,  0.51) (  0.21,  0.02)
* 
* Details of factorization
* (-16.42,  1.69) ( -0.53,  0.35) (  0.36,  0.41) ( -0.78,  0.49)
* (  0.00,  0.00) (  3.69,  0.64) (-16.58, -1.61) ( -0.10, -0.65)
* (  0.00,  0.00) (  0.00,  0.00) (  1.02, -3.74) ( -0.73, -0.52)
* (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  9.04,  3.96)
* 
* Pivot indices
*      1     -1     -1      3
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 4, NRHS = 2 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 100 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      INTEGER          IPIV( N )
      COMPLEX          A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $ ( 9.99,-4.73),( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),
     $ (-5.68,-0.80),(-8.01, 4.61),( 0.00, 0.00),( 0.00, 0.00),
     $ (-8.94, 1.32),( 1.64,-6.29),( 9.04, 3.96),( 0.00, 0.00),
     $ (-9.42, 2.05),( 6.79,-2.17),(-4.51,-7.54),( 0.40, 4.06)
     $                  /
      DATA             B/
     $ ( 5.71,-1.20),(-7.70, 6.47),( 3.77,-7.40),(-3.78, 0.33),
     $ ( 2.84,-0.18),(-8.29,-1.72),(-4.28,-8.25),(-2.70,-0.39)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CSYSV
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CSYSV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CSYSV( 'Upper', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL CSYSV( 'Upper', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
     $            INFO )
*
*     Check for the exact singularity.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The element of the diagonal factor '
         WRITE(*,*)'D(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'D is singular; the solution could not be computed.'
         STOP
      END IF
*
*     Print solution.
*
      CALL PRINT_MATRIX( 'Solution', N, NRHS, B, LDB )
*
*     Print details of factorization.
*
      CALL PRINT_MATRIX( 'Details of factorization', N, N, A, LDA )
*
*     Print pivot indices.
*
      CALL PRINT_INT_VECTOR( 'Pivot indices', N, IPIV )
      STOP
      END
*
*     End of CSYSV Example.
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
*     Auxiliary routine: printing a vector of integers.
*
      SUBROUTINE PRINT_INT_VECTOR( DESC, N, A )
      CHARACTER*(*)    DESC
      INTEGER          N
      INTEGER          A( N )
*
      INTEGER          I
*
      WRITE(*,*)
      WRITE(*,*) DESC
      WRITE(*,9999) ( A( I ), I = 1, N )
*
 9999 FORMAT( 11(:,1X,I6) )
      RETURN
      END
