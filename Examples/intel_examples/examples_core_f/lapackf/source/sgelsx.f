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
*  SGELS Example.
*  ==============
*
*  Program computes the least squares solution to the overdetermined linear
*  system A*X = B with full rank matrix A using QR factorization,
*  where A is the coefficient matrix:
*
*    1.44  -7.84  -4.39   4.53
*   -9.96  -0.28  -3.24   3.83
*   -7.55   3.24   6.27  -6.64
*    8.34   8.09   5.28   2.06
*    7.08   2.52   0.74  -2.47
*   -5.45  -5.70  -1.19   4.70
*
*  and B is the right-hand side matrix:
*
*    8.58   9.35
*    8.26  -4.43
*    8.48  -0.70
*   -5.28  -0.26
*    5.72  -7.36
*    8.93  -2.52
*
*  Description.
*  ============
*
*  The routine solves overdetermined or underdetermined real linear systems
*  involving an m-by-n matrix A, or its transpose, using a QR or LQ
*  factorization of A. It is assumed that A has full rank.
*
*  Several right hand side vectors b and solution vectors x can be handled
*  in a single call; they are stored as the columns of the m-by-nrhs right
*  hand side matrix B and the n-by-nrhs solution matrix X.
*
*  Example Program Results.
*  ========================
*
* SGELS Example Program Results
* 
* Solution
*  -0.45   0.25
*  -0.85  -0.90
*   0.71   0.63
*   0.13   0.14
* 
* Residual sum of squares for the solution
* 195.36 107.06
* 
* Details of QR factorization
* -17.54  -4.76  -1.96   0.42
*  -0.52  12.40   7.88  -5.84
*  -0.40  -0.14  -5.75   4.11
*   0.44  -0.66  -0.20  -7.78
*   0.37  -0.26  -0.17  -0.15
*  -0.29   0.46   0.41   0.24
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N, NRHS
      PARAMETER        ( M = 6, N = 4, NRHS = 2 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = M, LDB = M )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 100 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      REAL             A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $  1.44,-9.96,-7.55, 8.34, 7.08,-5.45,
     $ -7.84,-0.28, 3.24, 8.09, 2.52,-5.70,
     $ -4.39,-3.24, 6.27, 5.28, 0.74,-1.19,
     $  4.53, 3.83,-6.64, 2.06,-2.47, 4.70
     $                  /
      DATA             B/
     $  8.58, 8.26, 8.48,-5.28, 5.72, 8.93,
     $  9.35,-4.43,-0.70,-0.26,-7.36,-2.52
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         SGELS
      EXTERNAL         PRINT_MATRIX, PRINT_VECTOR_NORM
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'SGELS Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL SGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK, 
     $            LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL SGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK, 
     $            LWORK, INFO )
*
*     Check for the full rank.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
         WRITE(*,*)'factor of A is zero, so that A does not have full '
         WRITE(*,*)'rank; the least squares solution could not be '
         WRITE(*,*)'computed.'
         STOP
      END IF
*
*     Print least squares solution.
*
      CALL PRINT_MATRIX( 'Least squares solution', N, NRHS, B, LDB )
*
*     Print residual sum of squares for the solution
*
      CALL PRINT_VECTOR_NORM( 
     $     'Residual sum of squares for the solution', M-N, NRHS,
     $     B( N+1, 1 ), LDB )
*
*     Print details of QR factorization.
*
      CALL PRINT_MATRIX( 'Details of QR factorization', M, N, A, LDA )
      STOP
      END
*
*     End of SGELS Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
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
*
*     Auxiliary routine: printing norms of matrix columns.
*
      SUBROUTINE PRINT_VECTOR_NORM( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      REAL             A( LDA, * )
*
      REAL             TEMP
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO J = 1, N
         TEMP = 0.0
         DO I = 1, M
            TEMP = TEMP + A( I, J )*A( I, J )
         END DO
         WRITE(*,9998,ADVANCE='NO') TEMP
      END DO
      WRITE(*,*)
*
 9998 FORMAT( 11(:,1X,F6.2) )
      RETURN
      END
