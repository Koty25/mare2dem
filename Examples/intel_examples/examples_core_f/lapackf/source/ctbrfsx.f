*     CTBRFS (F07VVE) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
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
      INTEGER          NMAX, KDMAX, LDAB, NRHMAX, LDB, LDX
      PARAMETER        (NMAX=8,KDMAX=NMAX,LDAB=KDMAX+1,NRHMAX=NMAX,
     +                 LDB=NMAX,LDX=NMAX)
      CHARACTER        TRANS, DIAG
      PARAMETER        (TRANS='N',DIAG='N')
*     .. Local Scalars ..
      INTEGER          I, INFO, J, KD, N, NRHS
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX          AB(LDAB,NMAX), B(LDB,NRHMAX), WORK(2*NMAX),
     +                 X(LDX,NMAX)
      REAL             BERR(NRHMAX), FERR(NRHMAX), RWORK(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         CTBRFS, CTBTRS, CLACPY
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CTBRFS Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, KD, NRHS
      IF (N.LE.NMAX .AND. KD.LE.KDMAX .AND. NRHS.LE.NRHMAX) THEN
*
*        Read A and B from data file, and copy B to X
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            DO 20 I = 1, N
               READ (NIN,*) (AB(KD+1+I-J,J),J=I,MIN(N,I+KD))
   20       CONTINUE
         ELSE IF (UPLO.EQ.'L') THEN
            DO 40 I = 1, N
               READ (NIN,*) (AB(1+I-J,J),J=MAX(1,I-KD),I)
   40       CONTINUE
         END IF
         READ (NIN,*) ((B(I,J),J=1,NRHS),I=1,N)
         CALL CLACPY('General',N,NRHS,B,LDB,X,LDX)
*
*        Compute solution in the array X
*
         CALL CTBTRS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,X,LDX,INFO)
*
*        Compute backward errors and estimated bounds on the
*        forward errors
*
         CALL CTBRFS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,X,LDX,FERR,
     +               BERR,WORK,RWORK,INFO)
*
*        Print solution
*
         WRITE (NOUT,*)
         CALL PRINT_MATRIX( 'Solution(s)', N, NRHS, X, LDX )
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'Backward errors (machine-dependent)'
         WRITE (NOUT,99999) (BERR(J),J=1,NRHS)
         WRITE (NOUT,*)
     +     'Estimated forward error bounds (machine-dependent)'
         WRITE (NOUT,99999) (FERR(J),J=1,NRHS)
      END IF
*
*
      STOP
*
99999 FORMAT ((5X,1P,4(E11.1,7X)))
      END
*
*     End of CTBRFS Example
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
