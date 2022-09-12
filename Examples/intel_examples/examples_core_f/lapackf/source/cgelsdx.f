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
*  CGELSD Example.
*  ==============
*
*  Program computes the minimum norm-solution to a complex linear least squares
*  problem using the singular value decomposition of A,
*  where A is the coefficient matrix:
*
*  (  4.55, -0.32) ( -4.36, -4.76) (  3.99, -6.84) (  8.03, -6.47)
*  (  8.87, -3.11) (  0.02,  8.43) (  5.43, -9.30) (  2.28,  8.94)
*  ( -0.74,  1.16) (  3.80, -6.12) ( -7.24,  0.72) (  2.21,  9.52)
*
*  and B is the right-hand side matrix:
*
*  ( -8.25,  7.98) (  2.91, -8.81)
*  ( -5.04,  3.33) (  6.19,  0.19)
*  (  7.98, -4.38) ( -5.96,  7.18)
*
*  Description.
*  ============
*
*  The routine computes the minimum-norm solution to a complex linear least
*  squares problem: minimize ||b - A*x|| using the singular value
*  decomposition (SVD) of A. A is an m-by-n matrix which may be rank-deficient.
*
*  Several right hand side vectors b and solution vectors x can be handled
*  in a single call; they are stored as the columns of the m-by-nrhs right
*  hand side matrix B and the n-by-nrhs solution matrix X.
*
*  The effective rank of A is determined by treating as zero those singular
*  values which are less than rcond times the largest singular value.
*
*  Example Program Results.
*  ========================
*
* CGELSD Example Program Results
* 
* Minimum norm solution
* ( -0.08,  0.09) (  0.04,  0.16)
* ( -0.17,  0.10) (  0.17, -0.47)
* ( -0.92, -0.01) (  0.71, -0.41)
* ( -0.47, -0.26) (  0.69,  0.02)
*
* Effective rank =      3
* 
* Singular values
*  20.01  18.21   7.88
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N, NRHS
      PARAMETER        ( M = 3, N = 4, NRHS = 2 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = M, LDB = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, RANK
      REAL             RCOND
*
*     .. Local Arrays ..
*     IWORK dimension should be at least 3*MIN(M,N)*NLVL + 11*MIN(M,N),
*     RWORK dimension should be at least 10*MIN(M,N)+2*MIN(M,N)*SMLSIZ+
*     +8*MIN(M,N)*NLVL+3*SMLSIZ*NRHS+(SMLSIZ+1)**2,
*     where NLVL = MAX( 0, INT( LOG_2( MIN(M,N)/(SMLSIZ+1) ) )+1 )
*     and SMLSIZ = 25
      INTEGER          IWORK( 3*M*0+11*M )
      REAL             S( M ), RWORK( 10*M+2*M*25+8*M*0+3*25*NRHS+26*26
     $                 )
      COMPLEX          A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $ ( 4.55,-0.32),( 8.87,-3.11),(-0.74, 1.16),
     $ (-4.36,-4.76),( 0.02, 8.43),( 3.80,-6.12),
     $ ( 3.99,-6.84),( 5.43,-9.30),(-7.24, 0.72),
     $ ( 8.03,-6.47),( 2.28, 8.94),( 2.21, 9.52)
     $                  /
      DATA             B/
     $ (-8.25, 7.98),(-5.04, 3.33),( 7.98,-4.38),( 0.00, 0.00),
     $ ( 2.91,-8.81),( 6.19, 0.19),(-5.96, 7.18),( 0.00, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CGELSD
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CGELSD Example Program Results'
*     Negative RCOND means using default (machine precision) value
      RCOND = -1.0
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, 
     $             LWORK, RWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, 
     $             LWORK, RWORK, IWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge;'
         WRITE(*,*)'the least squares solution could not be computed.'
         STOP
      END IF
*
*     Print minimum norm solution.
*
      CALL PRINT_MATRIX( 'Minimum norm solution', N, NRHS, B, LDB )
*
*     Print effective rank.
*
      WRITE(*,'(/A,I6)')' Effective rank = ', RANK
*
*     Print singular values.
*
      CALL PRINT_RMATRIX( 'Singular values', 1, M, S, 1 )
      STOP
      END
*
*     End of CGELSD Example.
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
