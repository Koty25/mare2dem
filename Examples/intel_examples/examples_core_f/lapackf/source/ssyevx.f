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
*  SSYEV Example.
*  ==============
*
*  Program computes all eigenvalues and eigenvectors of a real symmetric
*  matrix A:
*
*    1.96  -6.49  -0.47  -7.20  -0.65
*   -6.49   3.80  -6.39   1.50  -6.34
*   -0.47  -6.39   4.17  -1.51   2.67
*   -7.20   1.50  -1.51   5.70   1.80
*   -0.65  -6.34   2.67   1.80  -7.10
*
*  Description.
*  ============
*
*  The routine computes all eigenvalues and, optionally, eigenvectors of an
*  n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
*
*  A*v(j) = lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The computed eigenvectors are
*  orthonormal.
*
*  Example Program Results.
*  ========================
*
* SSYEV Example Program Results
* 
* Eigenvalues
* -11.07  -6.23   0.86   8.87  16.09
* 
* Eigenvectors (stored columnwise)
*  -0.30  -0.61   0.40  -0.37   0.49
*  -0.51  -0.29  -0.41  -0.36  -0.61
*  -0.08  -0.38  -0.66   0.50   0.40
*   0.00  -0.45   0.46   0.62  -0.46
*  -0.80   0.45   0.17   0.31   0.16
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 5 )
      INTEGER          LDA
      PARAMETER        ( LDA = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      REAL             A( LDA, N ), W( N ), WORK( LWMAX )
      DATA             A/
     $  1.96, 0.00, 0.00, 0.00, 0.00,
     $ -6.49, 3.80, 0.00, 0.00, 0.00,
     $ -0.47,-6.39, 4.17, 0.00, 0.00,
     $ -7.20, 1.50,-1.51, 5.70, 0.00,
     $ -0.65,-6.34, 2.67, 1.80,-7.10
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         SSYEV
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'SSYEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL SSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL SSYEV( 'Vectors', 'Upper', N, A, LDA, W, WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print eigenvalues.
*
      CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, N, A, 
     $                   LDA )
      STOP
      END
*
*     End of SSYEV Example.
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
