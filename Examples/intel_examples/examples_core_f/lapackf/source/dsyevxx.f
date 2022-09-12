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
*  DSYEVX Example.
*  ==============
*
*  Program computes the smallest eigenvalues and the corresponding
*  eigenvectors of a real symmetric matrix A:
*
*    6.29  -0.39   0.61   1.18  -0.08
*   -0.39   7.19   0.81   1.19  -0.08
*    0.61   0.81   5.48  -3.13   0.22
*    1.18   1.19  -3.13   3.79  -0.26
*   -0.08  -0.08   0.22  -0.26   0.83
*
*  Description.
*  ============
*
*  The routine computes selected eigenvalues and, optionally, eigenvectors of
*  an n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies
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
* DSYEVX Example Program Results
*
* The total number of eigenvalues found: 3
* 
* Selected eigenvalues
*   0.71   0.82   6.58
* 
* Selected eigenvectors (stored columnwise)
*   0.22   0.09  -0.95
*   0.21   0.08  -0.04
*  -0.52  -0.22  -0.29
*  -0.73  -0.21  -0.09
*  -0.32   0.94   0.01
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N, NSELECT
      PARAMETER        ( N = 5, NSELECT = 3 )
      INTEGER          LDA, LDZ
      PARAMETER        ( LDA = N, LDZ = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK, IL, IU, M
      DOUBLE PRECISION ABSTOL, VL, VU
*
*     .. Local Arrays ..
*     IWORK dimension should be at least 5*N
      INTEGER          IFAIL( N ), IWORK( 5*N )
      DOUBLE PRECISION A( LDA, N ), W( N ), Z( LDZ, NSELECT ),
     $                 WORK( LWMAX )
      DATA             A/
     $  6.29, 0.00, 0.00, 0.00, 0.00,
     $ -0.39, 7.19, 0.00, 0.00, 0.00,
     $  0.61, 0.81, 5.48, 0.00, 0.00,
     $  1.18, 1.19,-3.13, 3.79, 0.00,
     $ -0.08,-0.08, 0.22,-0.26, 0.83
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DSYEVX
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DSYEVX Example Program Results'
*     Negative ABSTOL means using the default value
      ABSTOL = -1.0
*     Set IL, IU to compute NSELECT smallest eigenvalues
      IL = 1
      IU = NSELECT
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DSYEVX( 'Vectors', 'Indices', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL,
     $             INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL DSYEVX( 'Vectors', 'Indices', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL,
     $             INFO )
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
      CALL PRINT_MATRIX( 'Selected eigenvalues', 1, M, W, 1 )
*
*     Print eigenvectors.
*
      CALL PRINT_MATRIX( 'Selected eigenvectors (stored columnwise)',
     $                   N, M, Z, LDZ )
      STOP
      END
*
*     End of DSYEVX Example.
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
