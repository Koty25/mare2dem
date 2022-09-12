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
*  ZGESV Example.
*  ==============
*
*  The program computes the solution to the system of linear
*  equations with a square matrix A and multiple
*  right-hand sides B, where A is the coefficient matrix:
*
*  (  1.23, -5.50) (  7.91, -5.38) ( -9.80, -4.86) ( -7.32,  7.57) 
*  ( -2.14, -1.12) ( -9.92, -0.79) ( -9.18, -1.12) (  1.37,  0.43) 
*  ( -4.30, -7.10) ( -6.47,  2.52) ( -6.51, -2.67) ( -5.86,  7.38) 
*  (  1.27,  7.29) (  8.90,  6.92) ( -8.82,  1.25) (  5.41,  5.37) 
*
*  and B is the right-hand side matrix:
*
*  (  8.33, -7.32) ( -6.11, -3.81) 
*  ( -6.18, -4.80) (  0.14, -7.71) 
*  ( -5.71, -2.80) (  1.41,  3.40) 
*  ( -1.60,  3.08) (  8.54, -4.05) 
*
*  Description.
*  ============
*
*  The routine solves for X the system of linear equations A*X = B, 
*  where A is an n-by-n matrix, the columns of matrix B are individual 
*  right-hand sides, and the columns of X are the corresponding 
*  solutions.
*
*  The LU decomposition with partial pivoting and row interchanges is 
*  used to factor A as A = P*L*U, where P is a permutation matrix, L 
*  is unit lower triangular, and U is upper triangular. The factored 
*  form of A is then used to solve the system of equations A*X = B.
*
*  Example Program Results.
*  ========================
*
* ZGESV Example Program Results
* 
* Solution
* ( -1.09, -0.18) (  1.28,  1.21)
* (  0.97,  0.52) ( -0.22, -0.97)
* ( -0.20,  0.19) (  0.53,  1.36)
* ( -0.59,  0.92) (  2.22, -1.00)
* 
* Details of LU factorization
* ( -4.30, -7.10) ( -6.47,  2.52) ( -6.51, -2.67) ( -5.86,  7.38)
* (  0.49,  0.47) ( 12.26, -3.57) ( -7.87, -0.49) ( -0.98,  6.71)
* (  0.25, -0.15) ( -0.60, -0.37) (-11.70, -4.64) ( -1.35,  1.38)
* ( -0.83, -0.32) (  0.05,  0.58) (  0.93, -0.50) (  2.66,  7.86)
* 
* Pivot indices
*      3      3      3      4
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
      INTEGER          IPIV( N )
      COMPLEX*16       A( LDA, N ), B( LDB, NRHS )
      DATA             A/
     $ ( 1.23,-5.50),(-2.14,-1.12),(-4.30,-7.10),( 1.27, 7.29),
     $ ( 7.91,-5.38),(-9.92,-0.79),(-6.47, 2.52),( 8.90, 6.92),
     $ (-9.80,-4.86),(-9.18,-1.12),(-6.51,-2.67),(-8.82, 1.25),
     $ (-7.32, 7.57),( 1.37, 0.43),(-5.86, 7.38),( 5.41, 5.37)
     $                  /
      DATA             B/
     $ ( 8.33,-7.32),(-6.18,-4.80),(-5.71,-2.80),(-1.60, 3.08),
     $ (-6.11,-3.81),( 0.14,-7.71),( 1.41, 3.40),( 8.54,-4.05)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZGESV
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR
*
*     .. Executable Statements ..
      WRITE(*,*)'ZGESV Example Program Results'
*
*     Solve the equations A*X = B.
*
      CALL ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*     Check for the exact singularity.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element of the triangular factor of A,'
         WRITE(*,*)'U(',INFO,',',INFO,') is zero, so that'
         WRITE(*,*)'A is singular; the solution could not be computed.'
         STOP
      END IF
*
*     Print solution.
*
      CALL PRINT_MATRIX( 'Solution', N, NRHS, B, LDB )
*
*     Print details of LU factorization.
*
      CALL PRINT_MATRIX( 'Details of LU factorization', N, N, A, LDA )
*
*     Print pivot indices.
*
      CALL PRINT_INT_VECTOR( 'Pivot indices', N, IPIV )
      STOP
      END
*
*     End of ZGESV Example.
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
