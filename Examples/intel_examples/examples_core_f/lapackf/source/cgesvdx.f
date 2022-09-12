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
*  CGESVD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular complex matrix A:
*
*  (  5.91, -5.69) (  7.09,  2.72) (  7.78, -4.06) ( -0.79, -7.21)
*  ( -3.15, -4.08) ( -1.89,  3.27) (  4.57, -2.07) ( -3.88, -3.30)
*  ( -4.89,  4.20) (  4.10, -6.70) (  3.28, -3.84) (  3.84,  1.19)
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a complex
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written as
*
*  A = U*SIGMA*VH
*
*  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
*  diagonal elements, U is an m-by-m unitary matrix and VH (V conjugate
*  transposed) is an n-by-n unitary matrix. The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and are
*  returned in descending order. The first min(m, n) columns of U and V are
*  the left and right singular vectors of A.
*
*  Note that the routine returns VH, not V.
*
*  Example Program Results.
*  ========================
*
* CGESVD Example Program Results
* 
* Singular values
*  17.63  11.61   6.78
* 
* Left singular vectors (stored columnwise)
* ( -0.86,  0.00) (  0.40,  0.00) (  0.32,  0.00)
* ( -0.35,  0.13) ( -0.24, -0.21) ( -0.63,  0.60)
* (  0.15,  0.32) (  0.61,  0.61) ( -0.36,  0.10)
* 
* Right singular vectors (stored rowwise)
* ( -0.22,  0.51) ( -0.37, -0.32) ( -0.53,  0.11) (  0.15,  0.38)
* (  0.31,  0.31) (  0.09, -0.57) (  0.18, -0.39) (  0.38, -0.39)
* (  0.53,  0.24) (  0.49,  0.28) ( -0.47, -0.25) ( -0.15,  0.19)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N
      PARAMETER        ( M = 3, N = 4 )
      INTEGER          LDA, LDU, LDVT
      PARAMETER        ( LDA = M, LDU = M, LDVT = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
*     RWORK dimension should be at least MAX( 1, 5*MIN(M,N) )
      REAL             S( M ), RWORK( 5*M )
      COMPLEX          A( LDA, N ), U( LDU, M ), VT( LDVT, N ), 
     $                 WORK( LWMAX )
      DATA             A/
     $ ( 5.91,-5.69),(-3.15,-4.08),(-4.89, 4.20),
     $ ( 7.09, 2.72),(-1.89, 3.27),( 4.10,-6.70),
     $ ( 7.78,-4.06),( 4.57,-2.07),( 3.28,-3.84),
     $ (-0.79,-7.21),(-3.88,-3.30),( 3.84, 1.19)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CGESVD
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, RWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL CGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, RWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF
*
*     Print singular values.
*
      CALL PRINT_RMATRIX( 'Singular values', 1, M, S, 1 )
*
*     Print left singular vectors.
*
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)', 
     $                   M, M, U, LDU )
*
*     Print right singular vectors.
*
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)', 
     $                   M, N, VT, LDVT )
      STOP
      END
*
*     End of CGESVD Example.
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
