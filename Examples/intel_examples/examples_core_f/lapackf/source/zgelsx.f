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
*  ZGELS Example.
*  ==============
*
*  Program computes the minimum norm solution to the underdetermined linear
*  system A*X = B with full rank matrix A using LQ factorization,
*  where A is the coefficient matrix:
*
*  ( -4.20, -3.44) ( -3.35,  1.52) (  1.73,  8.85) (  2.35,  0.34)
*  ( -5.43, -8.81) ( -4.53, -8.47) (  5.93,  3.75) ( -3.75, -5.66)
*  ( -5.56,  3.39) (  2.90, -9.22) (  8.03,  9.37) (  5.69, -0.47)
*
*  and B is the right-hand side matrix:
*
*  ( -7.02,  4.80) (  3.88, -2.59)
*  (  0.62, -2.40) (  1.57,  3.24)
*  (  3.10, -2.19) ( -6.93, -5.99)
*
*  Description.
*  ============
*
*  The routine solves overdetermined or underdetermined complex linear systems
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
* ZGELS Example Program Results
*
* Minimum norm solution
* ( -0.25, -0.04) ( -0.21,  0.42)
* (  0.99,  0.27) ( -0.21, -0.43)
* (  0.25,  0.43) ( -0.24, -0.13)
* ( -0.32,  0.14) ( -0.23, -0.09)
*
* Details of LQ factorization
* ( 11.40,  0.00) (  0.18, -0.14) ( -0.23, -0.52) ( -0.15,  0.01)
* (  7.73, -0.39) ( 15.32,  0.00) ( -0.22,  0.42) (  0.45,  0.17)
* (  8.60, -5.68) (  3.96,  6.46) ( 12.54,  0.00) ( -0.02, -0.47)
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          M, N, NRHS
      PARAMETER        ( M = 3, N = 4, NRHS = 2 )
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = M, LDB = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 100 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      COMPLEX*16       A( LDA, N ), B( LDB, NRHS ), WORK( LWMAX )
      DATA             A/
     $ (-4.20,-3.44),(-5.43,-8.81),(-5.56, 3.39),
     $ (-3.35, 1.52),(-4.53,-8.47),( 2.90,-9.22),
     $ ( 1.73, 8.85),( 5.93, 3.75),( 8.03, 9.37),
     $ ( 2.35, 0.34),(-3.75,-5.66),( 5.69,-0.47)
     $                  /
      DATA             B/
     $ (-7.02, 4.80),( 0.62,-2.40),( 3.10,-2.19),( 0.00, 0.00),
     $ ( 3.88,-2.59),( 1.57, 3.24),(-6.93,-5.99),( 0.00, 0.00)
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         ZGELS
      EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'ZGELS Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL ZGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK, 
     $            LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve the equations A*X = B.
*
      CALL ZGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK, 
     $            LWORK, INFO )
*
*     Check for the full rank.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
         WRITE(*,*)'factor of A is zero, so that A does not have full '
         WRITE(*,*)'rank; the minimum norm solution could not be '
         WRITE(*,*)'computed.'
         STOP
      END IF
*
*     Print minimum norm solution.
*
      CALL PRINT_MATRIX( 'Minimum norm solution', N, NRHS, B, LDB )
*
*     Print details of LQ factorization.
*
      CALL PRINT_MATRIX( 'Details of LQ factorization', M, N, A, LDA )
      STOP
      END
*
*     End of ZGELS Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
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
