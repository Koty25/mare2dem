*     ZTRCON (F07TUF) Example Program Text
*     Mark 16 Release. NAG Copyright 1993.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDA
      PARAMETER        (NMAX=8,LDA=NMAX)
      CHARACTER        NORM, DIAG
      PARAMETER        (NORM='1',DIAG='N')
*     .. Local Scalars ..
      DOUBLE PRECISION RCOND
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX*16       A(LDA,NMAX), WORK(2*NMAX)
      DOUBLE PRECISION RWORK(NMAX)
*     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL         DLAMCH
*     .. External Subroutines ..
      EXTERNAL         ZTRCON
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZTRCON Example Program Results'
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
*        Estimate condition number
*
         CALL ZTRCON(NORM,UPLO,DIAG,N,A,LDA,RCOND,WORK,RWORK,INFO)
*
         WRITE (NOUT,*)
         IF (RCOND.GE.DLAMCH('E')) THEN
            WRITE (NOUT,99999) 'Estimate of condition number =',
     +        1.0D0/RCOND
         ELSE
            WRITE (NOUT,*) 'A is singular to working precision'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,D10.2)
      END
