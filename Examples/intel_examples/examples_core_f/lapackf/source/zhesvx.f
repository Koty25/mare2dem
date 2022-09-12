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
*  ZHESV Example.
*  ==============
*
*  The program computes the solution to the system of linear equations
*  with a Hermitian matrix A and multiple right-hand sides B,
*  where A is the coefficient matrix:
*
*  ( -2.90,  0.00) (  0.31,  4.46) (  9.66, -5.66) ( -2.28,  2.14)
*  (  0.31, -4.46) ( -7.93,  0.00) (  9.55, -4.62) ( -3.51,  3.11)
*  (  9.66,  5.66) (  9.55,  4.62) (  0.30,  0.00) (  9.33, -9.66)
*  ( -2.28, -2.14) ( -3.51, -3.11) (  9.33,  9.66) (  2.40,  0.00)
*
*  and B is the right-hand side matrix:
*
*  ( -5.69, -8.21) ( -2.83,  6.46)
*  ( -3.57,  1.99) ( -7.64,  1.10)
*  (  8.42, -9.83) ( -2.33, -4.23)
*  ( -5.00,  3.85) (  6.48, -3.81)
*
*  Description.
*  ============
*
*  The routine solves for X the complex system of linear equations A*X = B,
*  where A is an n-by-n Hermitian matrix, the columns of matrix B are
*  individual right-hand sides, and the columns of X are the corresponding
*  solutions.
*
*  The diagonal pivoting method is used to factor A as A = U*D*UH or
*  A = L*D*LH, where U (or L) is a product of permutation and unit upper
*  (lower) triangular matrices, and D is Hermitian and block diagonal with
*  1-by-1 and 2-by-2 diagonal blocks.
*
*  The factored form of A is then used to solve the system of equations A*X = B.
*
*  Example Program Results.
*  ========================
*
* ZHESV Example Program Results
* 
* Solution
* (  0.22, -0.95) ( -1.13,  0.18)
* ( -1.42, -1.30) (  0.70,  1.13)
* ( -0.65, -0.40) (  0.04,  0.07)
* ( -0.48,  1.35) (  1.15, -0.27)
* 
* Details of factorization
* (  3.17,  0.00) (  7.32,  3.28) ( -0.36,  0.06) (  0.20, -0.82)
* (  0.00,  0.00) (  0.03,  0.00) ( -0.48,  0.03) (  0.25, -0.76)
* (  0.00,  0.00) (  0.00,  0.00) (  0.30,  0.00) (  9.33, -9.66)
* (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  2.40,  0.00)
* 
* Pivot indices
*     -1     -1     -3     -3
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
      COMPLEX*16       A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $ (-2.90, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),
     $ ( 0.31, 4.46),(-7.93, 0.00),( 0.00, 0.00),( 0.00, 0.00),
     $ ( 9.66,-5.66),( 9.55,-4.62),( 0.30, 0.00),( 0.00, 0.00),
     $ (-2.28, 2.14),(-3.51, 3.11),( 9.33,-9.66),( 2.40, 0.00)
     $                  /
      DATA             B/
     $ (-5.69,-8.21),(-3.57, 1.99),( 8.42,-9.83),(-5.00, 3.85),
     $ (-2.83, 6.46),(-7.64, 1.10),(-2.33,-4.23),( 6.48,-3.81)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZHESV
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'ZHESV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL ZHESV( 'Upper', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL ZHESV( 'Upper', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
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
*     End of ZHESV Example.
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
