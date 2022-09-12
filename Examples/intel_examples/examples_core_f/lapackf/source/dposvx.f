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
*  DPOSV Example.
*  ==============
*
*  The program computes the solution to the system of linear
*  equations with a symmetric positive-definite matrix A and multiple
*  right-hand sides B, where A is the coefficient matrix:
*
*    3.14   0.17  -0.90   1.65  -0.72
*    0.17   0.79   0.83  -0.65   0.28
*   -0.90   0.83   4.53  -3.70   1.60
*    1.65  -0.65  -3.70   5.32  -1.37
*   -0.72   0.28   1.60  -1.37   1.98
*
*  and B is the right-hand side matrix:
*
*   -7.29   6.11   0.59
*    9.25   2.90   8.88
*    5.99  -5.05   7.57
*   -1.94  -3.80   5.57
*   -8.30   9.66  -1.67
*
*  Description.
*  ============
*
*  The routine solves for X the real system of linear equations 
*  A*X = B, where A is an n-by-n symmetric positive-definite 
*  matrix, the columns of matrix B are individual right-hand sides, 
*  and the columns of X are the corresponding solutions.
*
*  The Cholesky decomposition is used to factor A as
*  A = UT*U, if uplo = 'U' or A = L*LT, if uplo = 'L',
*  where U is an upper triangular matrix and L is a lower triangular matrix.
*  The factored form of A is then used to solve the system of equations A*X = B.
*
*  Example Program Results.
*  ========================
*
* DPOSV Example Program Results
* 
* Solution
*  -6.02   3.95  -3.14
*  15.62   4.32  13.05
*   3.02  -8.25   4.91
*   3.25  -4.83   6.11
*  -8.78   9.04  -3.57
* 
* Details of Cholesky factorization
*   1.77   0.10  -0.51   0.93  -0.41
*   0.00   0.88   0.99  -0.84   0.36
*   0.00   0.00   1.81  -1.32   0.57
*   0.00   0.00   0.00   1.42   0.05
*   0.00   0.00   0.00   0.00   1.16
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 5, NRHS = 3 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )
*
*     .. Local Scalars ..
      INTEGER          INFO
*
*     .. Local Arrays ..
      DOUBLE PRECISION A( LDA, N ), B( LDB, NRHS )
      DATA             A/
     $  3.14, 0.00, 0.00, 0.00, 0.00,
     $  0.17, 0.79, 0.00, 0.00, 0.00,
     $ -0.90, 0.83, 4.53, 0.00, 0.00,
     $  1.65,-0.65,-3.70, 5.32, 0.00,
     $ -0.72, 0.28, 1.60,-1.37, 1.98
     $                  /
      DATA             B/
     $ -7.29, 9.25, 5.99,-1.94,-8.30,
     $  6.11, 2.90,-5.05,-3.80, 9.66,
     $  0.59, 8.88, 7.57, 5.57,-1.67
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DPOSV
      EXTERNAL         PRINT_MATRIX
*
*     .. Executable Statements ..
      WRITE(*,*)'DPOSV Example Program Results'
*
*     Solve the equations A*X = B.
*
      CALL DPOSV( 'Upper', N, NRHS, A, LDA, B, LDB, INFO )
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
*     End of DPOSV Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
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
