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
*  CGESDD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular complex matrix A using a divide and conquer method, where A is:
*
*  ( -5.40,  7.40) (  6.00,  6.38) (  9.91,  0.16) ( -5.28, -4.16)
*  (  1.09,  1.55) (  2.60,  0.07) (  3.98, -5.26) (  2.03,  1.11)
*  (  9.88,  1.91) (  4.92,  6.31) ( -2.11,  7.39) ( -9.81, -8.98)
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a complex
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. If singular vectors are desired, it uses a divide and conquer
*  algorithm. The SVD is written as
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
* CGESDD Example Program Results
* 
* Singular values
*  21.76  16.60   3.97
* 
* Left singular vectors (stored columnwise)
* (  0.55,  0.00) (  0.76,  0.00) ( -0.34,  0.00)
* ( -0.04, -0.15) (  0.27, -0.23) (  0.55, -0.74)
* (  0.81,  0.12) ( -0.52, -0.14) (  0.13, -0.11)
* 
* Right singular vectors (stored rowwise)
* (  0.23,  0.21) (  0.37,  0.39) (  0.24,  0.33) ( -0.56, -0.37)
* ( -0.58,  0.40) (  0.11,  0.17) (  0.60, -0.27) (  0.16,  0.06)
* (  0.60,  0.12) ( -0.19,  0.30) (  0.39,  0.20) (  0.45,  0.31)
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
*     IWORK dimension should be at least 8*MIN(M,N)
      INTEGER          IWORK( 8*M )
*     RWORK dimension should be at least 5*(MIN(M,N))**2 + 7*MIN(M,N))
      REAL             S( M ), RWORK( 5*M*M + 7*M )
      COMPLEX          A( LDA, N ), U( LDU, M ), VT( LDVT, N ), 
     $                 WORK( LWMAX )
      DATA             A/
     $ (-5.40, 7.40),( 1.09, 1.55),( 9.88, 1.91),
     $ ( 6.00, 6.38),( 2.60, 0.07),( 4.92, 6.31),
     $ ( 9.91, 0.16),( 3.98,-5.26),(-2.11, 7.39),
     $ (-5.28,-4.16),( 2.03, 1.11),(-9.81,-8.98)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CGESDD
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CGESDD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL CGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT, 
     $             LDVT, WORK, LWORK, RWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL CGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT, 
     $             LDVT, WORK, LWORK, RWORK, IWORK, INFO )
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
*     End of CGESDD Example.
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
