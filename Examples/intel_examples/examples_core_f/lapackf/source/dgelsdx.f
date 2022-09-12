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
*  DGELSD Example.
*  ==============
*
*  Program computes the minimum norm-solution to a real linear least squares
*  problem using the singular value decomposition of A,
*  where A is the coefficient matrix:
*
*    0.12  -8.19   7.69  -2.26  -4.71
*   -6.91   2.22  -5.12  -9.08   9.96
*   -3.33  -8.94  -6.72  -4.40  -9.98
*    3.97   3.33  -2.74  -7.92  -3.20
*
*  and B is the right-hand side matrix:
*
*    7.30   0.47  -6.28
*    1.33   6.58  -3.42
*    2.68  -1.71   3.46
*   -9.62  -0.79   0.41
*
*  Description.
*  ============
*
*  The routine computes the minimum-norm solution to a real linear least
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
* DGELSD Example Program Results
*
* Minimum norm solution
*  -0.69  -0.24   0.06
*  -0.80  -0.08   0.21
*   0.38   0.12  -0.65
*   0.29  -0.24   0.42
*   0.29   0.35  -0.30
*
* Effective rank =      4
*
* Singular values
*  18.66  15.99  10.01   8.51
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N, NRHS
      PARAMETER        ( M = 4, N = 5, NRHS = 3 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = M, LDB = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, RANK
      DOUBLE PRECISION RCOND
*
*     .. Local Arrays ..
*     IWORK dimension should be at least 3*MIN(M,N)*NLVL + 11*MIN(M,N),
*     where NLVL = MAX( 0, INT( LOG_2( MIN(M,N)/(SMLSIZ+1) ) )+1 )
*     and SMLSIZ = 25
      INTEGER          IWORK( 3*M*0+11*M )
      DOUBLE PRECISION A( LDA, N ), B( LDB, NRHS ), S( M ), 
     $                 WORK( LWMAX )
      DATA             A/
     $  0.12,-6.91,-3.33, 3.97,
     $ -8.19, 2.22,-8.94, 3.33,
     $  7.69,-5.12,-6.72,-2.74,
     $ -2.26,-9.08,-4.40,-7.92,
     $ -4.71, 9.96,-9.98,-3.20
     $                  /
      DATA             B/
     $  7.30, 1.33, 2.68,-9.62, 0.00,
     $  0.47, 6.58,-1.71,-0.79, 0.00,
     $ -6.28,-3.42, 3.46, 0.41, 0.00
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DGELSD
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DGELSD Example Program Results'
*     Negative RCOND means using default (machine precision) value
      RCOND = -1.0
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, 
     $             LWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, 
     $             LWORK, IWORK, INFO )
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
      CALL PRINT_MATRIX( 'Singular values', 1, M, S, 1 )
      STOP
      END
*
*     End of DGELSD Example.
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
