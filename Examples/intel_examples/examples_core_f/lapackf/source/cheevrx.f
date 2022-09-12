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
*  CHEEVR Example.
*  ==============
*
*  Program computes eigenvalues specified by a selected range of values
*  and corresponding eigenvectors of a complex Hermitian matrix A using the
*  Relatively Robust Representations, where A is:
*
*  ( -2.16,  0.00) ( -0.16, -4.86) ( -7.23, -9.38) ( -0.04,  6.86)
*  ( -0.16,  4.86) (  7.45,  0.00) (  4.39,  6.29) ( -8.11, -4.41)
*  ( -7.23,  9.38) (  4.39, -6.29) ( -9.03,  0.00) ( -6.89, -7.66)
*  ( -0.04, -6.86) ( -8.11,  4.41) ( -6.89,  7.66) (  7.76,  0.00)
*
*  Description.
*  ============
*
*  The routine computes selected eigenvalues and, optionally, eigenvectors of
*  an n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies
*
*  A*v(j) = lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The computed eigenvectors are
*  orthonormal.
*  Eigenvalues and eigenvectors can be selected by specifying either a range
*  of values or a range of indices for the desired eigenvalues.
*
*  Example Program Results.
*  ========================
*
* CHEEVR Example Program Results
*
* The total number of eigenvalues found: 2
* 
* Selected eigenvalues
*  -4.18   3.57
* 
* Selected eigenvectors (stored columnwise)
* (  0.68,  0.00) (  0.38,  0.00)
* (  0.03,  0.18) (  0.54, -0.57)
* ( -0.03,  0.21) ( -0.40,  0.04)
* (  0.20,  0.64) ( -0.14, -0.26)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 4 )
      INTEGER          LDA, LDZ
      PARAMETER        ( LDA = N, LDZ = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, LRWORK, LIWORK, IL, IU, M
      REAL             ABSTOL, VL, VU
*
*     .. Local Arrays ..
      INTEGER          ISUPPZ( 2*N ), IWORK( LWMAX )
      REAL             W( N ), RWORK( LWMAX )
      COMPLEX          A( LDA, N ), Z( LDZ, N ), WORK( LWMAX )
      DATA             A/
     $ (-2.16, 0.00),(-0.16, 4.86),(-7.23, 9.38),(-0.04,-6.86),
     $ ( 0.00, 0.00),( 7.45, 0.00),( 4.39,-6.29),(-8.11, 4.41),
     $ ( 0.00, 0.00),( 0.00, 0.00),(-9.03, 0.00),(-6.89, 7.66),
     $ ( 0.00, 0.00),( 0.00, 0.00),( 0.00, 0.00),( 7.76, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         CHEEVR
      EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'CHEEVR Example Program Results'
*     Negative ABSTOL means using the default value
      ABSTOL = -1.0
*     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
      VL = -5.0
      VU = 5.0
*
*     Query the optimal workspace.
*
      LWORK = -1
      LRWORK = -1
      LIWORK = -1
      CALL CHEEVR( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL CHEEVR( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK,
     $             LRWORK, IWORK, LIWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
*
*     Print the number of eigenvalues found.
*
      WRITE(*,'(/A,I2)')' The total number of eigenvalues found:', M
*
*     Print eigenvalues.
*
      CALL PRINT_RMATRIX( 'Selected eigenvalues', 1, M, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Selected eigenvectors (stored columnwise)',
     $                   N, M, Z, LDZ )
      STOP
      END
*
*     End of CHEEVR Example.
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
