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
*  DSYSV Example.
*  ==============
*
*  The program computes the solution to the system of linear equations
*  with a real symmetric matrix A and multiple right-hand sides B,
*  where A is the coefficient matrix:
*
*   -5.86   3.99  -5.93  -2.82   7.69
*    3.99   4.46   2.58   4.42   4.61
*   -5.93   2.58  -8.52   8.57   7.69
*   -2.82   4.42   8.57   3.72   8.07
*    7.69   4.61   7.69   8.07   9.83
*
*  and B is the right-hand side matrix:
*
*    1.32  -6.33  -8.77
*    2.22   1.69  -8.33
*    0.12  -1.56   9.54
*   -6.41  -9.49   9.56
*    6.33  -3.67   7.48
*
*  Description.
*  ============
*
*  The routine solves for X the real system of linear equations A*X = B,
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
* DSYSV Example Program Results
* 
* Solution
*   1.17   0.52  -0.86
*  -0.71   1.05  -4.90
*  -0.63  -0.52   0.99
*  -0.33   0.43   1.22
*   0.83  -1.22   1.96
* 
* Details of factorization
*  -5.86   0.00   0.00   0.00   0.00
*  -0.68   7.18   0.00   0.00   0.00
*   1.01  -0.20  -2.82   0.00   0.00
*   0.48   0.35  11.93   4.21   0.00
*  -1.31   1.37   0.02   0.16   6.22
* 
* Pivot indices
*      1      2     -4     -4      5
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 5, NRHS = 3 )
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
      DOUBLE PRECISION A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $ -5.86, 3.99,-5.93,-2.82, 7.69,
     $  0.00, 4.46, 2.58, 4.42, 4.61,
     $  0.00, 0.00,-8.52, 8.57, 7.69,
     $  0.00, 0.00, 0.00, 3.72, 8.07,
     $  0.00, 0.00, 0.00, 0.00, 9.83
     $                  /
      DATA             B/
     $  1.32, 2.22, 0.12,-6.41, 6.33,
     $ -6.33, 1.69,-1.56,-9.49,-3.67,
     $ -8.77,-8.33, 9.54, 9.56, 7.48
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR
      EXTERNAL         DSYSV
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DSYSV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DSYSV( 'Lower', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
     $            INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL DSYSV( 'Lower', N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, 
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
*     End of DSYSV Example.
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
