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
*  CPOSV Example.
*  ==============
*
*  The program computes the solution to the system of linear
*  equations with a Hermitian positive-definite matrix A and multiple
*  right-hand sides B, where A is the coefficient matrix:
*
*  (  5.96,  0.00) (  0.40, -1.19) ( -0.83, -0.48) ( -0.57,  0.40)
*  (  0.40,  1.19) (  7.95,  0.00) (  0.33,  0.09) (  0.22,  0.74)
*  ( -0.83,  0.48) (  0.33, -0.09) (  4.43,  0.00) ( -1.09,  0.32)
*  ( -0.57, -0.40) (  0.22, -0.74) ( -1.09, -0.32) (  3.46,  0.00)
*
*  and B is the right-hand side matrix:
*
*  ( -2.94,  5.79) (  8.44,  3.07)
*  (  8.12, -9.12) (  1.00, -4.62)
*  (  9.09, -5.03) (  3.64, -2.33)
*  (  7.36,  6.77) (  8.04,  2.87)
*
*  Description.
*  ============
*
*  The routine solves for X the complex system of linear equations 
*  A*X = B, where A is an n-by-n Hermitian positive-definite 
*  matrix, the columns of matrix B are individual right-hand sides, 
*  and the columns of X are the corresponding solutions.
*
*  The Cholesky decomposition is used to factor A as
*  A = UH*U, if uplo = 'U' or A = L*LH, if uplo = 'L',
*  where U is an upper triangular matrix and L is a lower triangular matrix.
*  The factored form of A is then used to solve the system of equations A*X = B.
*
*  Example Program Results.
*  ========================
*
* CPOSV Example Program Results
* 
* Solution
* (  0.80,  1.62) (  2.52,  0.61)
* (  1.26, -1.78) (  0.01, -1.38)
* (  3.38, -0.29) (  2.42, -0.52)
* (  3.46,  2.92) (  3.77,  1.37)
* 
* Details of Cholesky factorization
* (  2.44,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00)
* (  0.16,  0.49) (  2.77,  0.00) (  0.00,  0.00) (  0.00,  0.00)
* ( -0.34,  0.20) (  0.10, -0.10) (  2.06,  0.00) (  0.00,  0.00)
* ( -0.23, -0.16) (  0.12, -0.30) ( -0.57, -0.20) (  1.71,  0.00)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 4, NRHS = 2 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )
*
*     .. Local Scalars ..
      INTEGER          INFO
*
*     .. Local Arrays ..
      COMPLEX          A( LDA, N ), B( LDB, NRHS )
      DATA             A/
     $ ( 5.96, 0.00),( 0.40, 1.19),(-0.83, 0.48),(-0.57,-0.40),
     $ ( 0.00, 0.00),( 7.95, 0.00),( 0.33,-0.09),( 0.22,-0.74),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 4.43, 0.00),(-1.09,-0.32),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 3.46, 0.00)
     $                  /
      DATA             B/
     $ (-2.94, 5.79),( 8.12,-9.12),( 9.09,-5.03),( 7.36, 6.77),
     $ ( 8.44, 3.07),( 1.00,-4.62),( 3.64,-2.33),( 8.04, 2.87)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CPOSV
      EXTERNAL         PRINT_MATRIX
*
*     .. Executable Statements ..
      WRITE(*,*)'CPOSV Example Program Results'
*
*     Solve the equations A*X = B.
*
      CALL CPOSV( 'Lower', N, NRHS, A, LDA, B, LDB, INFO )
*
*     Check for the exact singularity.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The leading minor of order ',INFO,' is not positive'
         WRITE(*,*)'definite; the solution could not be computed.'
         STOP
      END IF
*
*     Print solution.
*
      CALL PRINT_MATRIX( 'Solution', N, NRHS, B, LDB )
*
*     Print details of Cholesky factorization.
*
      CALL PRINT_MATRIX( 'Details of Cholesky factorization', N, N, A, 
     $                   LDA )
      STOP
      END
*
*     End of CPOSV Example.
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
