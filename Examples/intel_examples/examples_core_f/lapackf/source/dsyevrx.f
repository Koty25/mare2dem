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
*  DSYEVR Example.
*  ==============
*
*  Program computes the smallest eigenvalues and the corresponding
*  eigenvectors of a real symmetric matrix A using the Relatively Robust
*  Representations, where A is:
*
*    0.67  -0.20   0.19  -1.06   0.46
*   -0.20   3.82  -0.13   1.06  -0.48
*    0.19  -0.13   3.27   0.11   1.10
*   -1.06   1.06   0.11   5.86  -0.98
*    0.46  -0.48   1.10  -0.98   3.54
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
* DSYEVR Example Program Results
*
* The total number of eigenvalues found: 3
* 
* Selected eigenvalues
*   0.43   2.14   3.37
* 
* Selected eigenvectors (stored columnwise)
*  -0.98  -0.01  -0.08
*   0.01   0.02  -0.93
*   0.04  -0.69  -0.07
*  -0.18   0.19   0.31
*   0.07   0.69  -0.13
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
      INTEGER          INFO, LWORK, LIWORK, IL, IU, M
      DOUBLE PRECISION ABSTOL, VL, VU
*
*     .. Local Arrays ..
      INTEGER          ISUPPZ( 2*N ), IWORK( LWMAX )
      DOUBLE PRECISION A( LDA, N ), W( N ), Z( LDZ, NSELECT ),
     $                 WORK( LWMAX )
      DATA             A/
     $  0.67, 0.00, 0.00, 0.00, 0.00,
     $ -0.20, 3.82, 0.00, 0.00, 0.00,
     $  0.19,-0.13, 3.27, 0.00, 0.00,
     $ -1.06, 1.06, 0.11, 5.86, 0.00,
     $  0.46,-0.48, 1.10,-0.98, 3.54
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DSYEVR
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DSYEVR Example Program Results'
*     Negative ABSTOL means using the default value
      ABSTOL = -1.0
*     Set IL, IU to compute NSELECT smallest eigenvalues
      IL = 1
      IU = NSELECT
*
*     Query the optimal workspace.
*
      LWORK = -1
      LIWORK = -1
      CALL DSYEVR( 'Vectors', 'Indices', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      LIWORK = MIN( LWMAX, IWORK( 1 ) )
*
*     Solve eigenproblem.
*
      CALL DSYEVR( 'Vectors', 'Indices', 'Upper', N, A, LDA, VL, VU, IL,
     $             IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
     $             LIWORK, INFO )
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
*     End of DSYEVR Example.
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
