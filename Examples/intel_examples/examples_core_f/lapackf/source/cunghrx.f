*     CUNGHR (F08NTE) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*
********************************************
*                                          *
* Modified by Intel Corporation, July 2017 *
*                                          *
********************************************
*
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDA, LDZ, LWORK
      PARAMETER        (NMAX=8,LDA=NMAX,LDZ=NMAX,LWORK=64*(NMAX-1))
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N
*     .. Local Arrays ..
      COMPLEX          A(LDA,NMAX), TAU(NMAX), W(NMAX), WORK(LWORK),
     +                 Z(LDZ,NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         CGEHRD, CHSEQR, CUNGHR, CLACPY
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CUNGHR Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Reduce A to upper Hessenberg form H = (Q**H)*A*Q
*
         CALL CGEHRD(N,1,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*        Copy A into Z
*
         CALL CLACPY('General',N,N,A,LDA,Z,LDZ)
*
*        Form Q explicitly, storing the result in Z
*
         CALL CUNGHR(N,1,N,Z,LDZ,TAU,WORK,LWORK,INFO)
*
*        Calculate the Schur factorization of H = Y*T*(Y**H) and form
*        Q*Y explicitly, storing the result in Z
*
*        Note that A = Z*T*(Z**H), where Z = Q*Y
*
         CALL CHSEQR('Schur form','Vectors',N,1,N,A,LDA,W,Z,LDZ,WORK,
     +               LWORK,INFO)
*
*        Print Schur form
*
         WRITE (NOUT,*)
*
         CALL PRINT_MATRIX( 'Schur form', N, N, A, LDA )
*
*        Print Schur vectors
*
         WRITE (NOUT,*)
*
         CALL PRINT_MATRIX( 'Schur vectors of A', N, N, Z, LDZ )
*
      END IF
*
      STOP
      END
*
*     End of CUNGHR Example
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
      WRITE(*,*) DESC
      WRITE(*, 9999) ( J, J = 1, N)
      DO I = 1, M
         WRITE(*,9998) I, ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( I2, ' ', 11(:,1X,'(',F7.4,',',F7.4,')') )
 9999 FORMAT( '   ', 11(:,1X,I17) )
*
      RETURN
      END
