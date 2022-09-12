*     CPOCON (F07FUE) Example Program Text
*     Mark 16 Release. NAG Copyright 1993.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDA
      PARAMETER        (NMAX=8,LDA=NMAX)
*     .. Local Scalars ..
      REAL             ANORM, RCOND
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX          A(LDA,NMAX), WORK(2*NMAX)
      REAL             RWORK(NMAX)
*     .. External Functions ..
      REAL             CLANHE, SLAMCH
      EXTERNAL         CLANHE, SLAMCH
*     .. External Subroutines ..
      EXTERNAL         CPOCON, CPOTRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CPOCON Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((A(I,J),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((A(I,J),J=1,I),I=1,N)
         END IF
*
*        Compute norm of A
*
         ANORM = CLANHE('1-norm',UPLO,N,A,LDA,RWORK)
*
*        Factorize A
*
         CALL CPOTRF(UPLO,N,A,LDA,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL CPOCON(UPLO,N,A,LDA,ANORM,RCOND,WORK,RWORK,INFO)
*
            IF (RCOND.GE.SLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'A is not positive-definite'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,E10.2)
      END
