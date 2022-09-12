*     CHPCON (F07PUE) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
*     .. Local Scalars ..
      REAL             ANORM, RCOND
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX          AP(NMAX*(NMAX+1)/2), WORK(2*NMAX)
      REAL             RWORK(NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Functions ..
      REAL             CLANHP, SLAMCH
      EXTERNAL         CLANHP, SLAMCH
*     .. External Subroutines ..
      EXTERNAL         CHPCON, CHPTRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CHPCON Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((AP(I+J*(J-1)/2),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((AP(I+(2*N-J)*(J-1)/2),J=1,I),I=1,N)
         END IF
*
*        Compute norm of A
*
         ANORM = CLANHP('1-norm',UPLO,N,AP,RWORK)
*
*        Factorize A
*
         CALL CHPTRF(UPLO,N,AP,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL CHPCON(UPLO,N,AP,IPIV,ANORM,RCOND,WORK,INFO)
*
            IF (RCOND.GE.SLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'The factor D is singular'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,E10.2)
      END
