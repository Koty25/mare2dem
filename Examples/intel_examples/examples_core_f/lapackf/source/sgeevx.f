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
*  SGEEV Example.
*  ==============
*
*  Program computes the eigenvalues and left and right eigenvectors of a general
*  rectangular matrix A:
*
*   -1.01   0.86  -4.60   3.31  -4.81
*    3.98   0.53  -7.04   5.29   3.55
*    3.30   8.26  -3.89   8.20  -1.51
*    4.43   4.96  -7.66  -7.33   6.18
*    7.31  -6.43  -6.16   2.47   5.58
*
*  Description.
*  ============
*
*  The routine computes for an n-by-n real nonsymmetric matrix A, the
*  eigenvalues and, optionally, the left and/or right eigenvectors. The right
*  eigenvector v(j) of A satisfies
*
*  A*v(j)= lambda(j)*v(j)
*
*  where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies
*
*  u(j)H*A = lambda(j)*u(j)H
*
*  where u(j)H denotes the conjugate transpose of u(j). The computed
*  eigenvectors are normalized to have Euclidean norm equal to 1 and
*  largest component real.
*
*  Example Program Results.
*  ========================
*
* SGEEV Example Program Results
* 
* Eigenvalues
* (  2.86, 10.76) (  2.86,-10.76) ( -0.69,  4.70) ( -0.69, -4.70) -10.46
* 
* Left eigenvectors
* (  0.04,  0.29) (  0.04, -0.29) ( -0.13, -0.33) ( -0.13,  0.33)   0.04
* (  0.62,  0.00) (  0.62,  0.00) (  0.69,  0.00) (  0.69,  0.00)   0.56
* ( -0.04, -0.58) ( -0.04,  0.58) ( -0.39, -0.07) ( -0.39,  0.07)  -0.13
* (  0.28,  0.01) (  0.28, -0.01) ( -0.02, -0.19) ( -0.02,  0.19)  -0.80
* ( -0.04,  0.34) ( -0.04, -0.34) ( -0.40,  0.22) ( -0.40, -0.22)   0.18
* 
* Right eigenvectors
* (  0.11,  0.17) (  0.11, -0.17) (  0.73,  0.00) (  0.73,  0.00)   0.46
* (  0.41, -0.26) (  0.41,  0.26) ( -0.03, -0.02) ( -0.03,  0.02)   0.34
* (  0.10, -0.51) (  0.10,  0.51) (  0.19, -0.29) (  0.19,  0.29)   0.31
* (  0.40, -0.09) (  0.40,  0.09) ( -0.08, -0.08) ( -0.08,  0.08)  -0.74
* (  0.54,  0.00) (  0.54,  0.00) ( -0.29, -0.49) ( -0.29,  0.49)   0.16
*  =============================================================================
*
*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 5 )
      INTEGER          LDA, LDVL, LDVR
      PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      REAL             A( LDA, N ), VL( LDVL, N ), VR( LDVR, N ), 
     $                 WR( N ), WI( N ), WORK( LWMAX )
      DATA             A/
     $ -1.01, 3.98, 3.30, 4.43, 7.31,
     $  0.86, 0.53, 8.26, 4.96,-6.43,
     $ -4.60,-7.04,-3.89,-7.66,-6.16,
     $  3.31, 5.29, 8.20,-7.33, 2.47,
     $ -4.81, 3.55,-1.51, 6.18, 5.58
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         SGEEV
      EXTERNAL         PRINT_EIGENVALUES, PRINT_EIGENVECTORS
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'SGEEV Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL SGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Solve eigenproblem.
*
      CALL SGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,
     $            VR, LDVR, WORK, LWORK, INFO )
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
      CALL PRINT_EIGENVALUES( 'Eigenvalues', N, WR, WI )
*
*     Print left eigenvectors.
*
      CALL PRINT_EIGENVECTORS( 'Left eigenvectors', N, WI, VL, LDVL )
*
*     Print right eigenvectors.
*
      CALL PRINT_EIGENVECTORS( 'Right eigenvectors', N, WI, VR, LDVR )
      STOP
      END
*
*     End of SGEEV Example.
*
*  =============================================================================
*
*     Auxiliary routine: printing eigenvalues.
*
      SUBROUTINE PRINT_EIGENVALUES( DESC, N, WR, WI )
      CHARACTER*(*)    DESC
      INTEGER          N
      REAL             WR( * ), WI( * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO J = 1, N
         IF( WI( J ).EQ.ZERO ) THEN
            WRITE(*,9998,ADVANCE='NO') WR( J )
         ELSE
            WRITE(*,9999,ADVANCE='NO') WR( J ), WI( J )
         END IF
      END DO
      WRITE(*,*)
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
*
*     Auxiliary routine: printing eigenvectors.
*
      SUBROUTINE PRINT_EIGENVECTORS( DESC, N, WI, V, LDV )
      CHARACTER*(*)    DESC
      INTEGER          N, LDV
      REAL             WI( * ), V( LDV, * )
*
      REAL             ZERO
      PARAMETER        ( ZERO = 0.0 )
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, N
         J = 1
         DO WHILE( J.LE.N )
            IF( WI( J ).EQ.ZERO ) THEN
               WRITE(*,9998,ADVANCE='NO') V( I, J )
               J = J + 1
            ELSE
               WRITE(*,9999,ADVANCE='NO') V( I, J ), V( I, J+1 )
               WRITE(*,9999,ADVANCE='NO') V( I, J ), -V( I, J+1 )
               J = J + 2
            END IF
         END DO
         WRITE(*,*)
      END DO
*
 9998 FORMAT( 11(:,1X,F6.2) )
 9999 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
      RETURN
      END
