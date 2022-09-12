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
*  DGESDD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular matrix A using a divide and conquer method, where A is:
*
*    7.52  -1.10  -7.95   1.08
*   -0.76   0.62   9.34  -7.10
*    5.13   6.62  -5.66   0.87
*   -4.75   8.52   5.75   5.30
*    1.33   4.91  -5.49  -3.52
*   -2.40  -6.77   2.34   3.95
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a real
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. If singular vectors are desired, it uses a divide and conquer
*  algorithm. The SVD is written as
*
*  A = U*SIGMA*VT
*
*  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
*  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
*  is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and are
*  returned in descending order. The first min(m, n) columns of U and V are
*  the left and right singular vectors of A.
*
*  Note that the routine returns VT, not V.
*
*  Example Program Results.
*  ========================
*
* DGESDD Example Program Results
* 
* Singular values
*  18.37  13.63  10.85   4.49
* 
* Left singular vectors (stored columnwise)
*  -0.57   0.18   0.01   0.53
*   0.46  -0.11  -0.72   0.42
*  -0.45  -0.41   0.00   0.36
*   0.33  -0.69   0.49   0.19
*  -0.32  -0.31  -0.28  -0.61
*   0.21   0.46   0.39   0.09
* 
* Right singular vectors (stored rowwise)
*  -0.52  -0.12   0.85  -0.03
*   0.08  -0.99  -0.09  -0.01
*  -0.28  -0.02  -0.14   0.95
*   0.81   0.01   0.50   0.31
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N
      PARAMETER        ( M = 6, N = 4 )
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
      INTEGER          IWORK( 8*N )
      DOUBLE PRECISION A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ), 
     $                 WORK( LWMAX )
      DATA             A/
     $  7.52,-0.76, 5.13,-4.75, 1.33,-2.40,
     $ -1.10, 0.62, 6.62, 8.52, 4.91,-6.77,
     $ -7.95, 9.34,-5.66, 5.75,-5.49, 2.34,
     $  1.08,-7.10, 0.87, 5.30,-3.52, 3.95
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DGESDD
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DGESDD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT, 
     $             LDVT, WORK, LWORK, IWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESDD( 'Singular vectors', M, N, A, LDA, S, U, LDU, VT, 
     $             LDVT, WORK, LWORK, IWORK, INFO )
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
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
*
*     Print left singular vectors.
*
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)', 
     $                   M, N, U, LDU )
*
*     Print right singular vectors.
*
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)', 
     $                   N, N, VT, LDVT )
      STOP
      END
*
*     End of DGESDD Example.
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
