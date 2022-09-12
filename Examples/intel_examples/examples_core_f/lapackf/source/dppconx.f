*     DPPCON (F07GGF) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
*     .. Local Scalars ..
      DOUBLE PRECISION ANORM, RCOND
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      DOUBLE PRECISION AP(NMAX*(NMAX+1)/2), WORK(3*NMAX)
      INTEGER          IWORK(NMAX)
*     .. External Functions ..
      DOUBLE PRECISION DLANSP, DLAMCH
      EXTERNAL         DLANSP, DLAMCH
*     .. External Subroutines ..
      EXTERNAL         DPPCON, DPPTRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'DPPCON Example Program Results'
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
         ANORM = DLANSP('1-norm',UPLO,N,AP,WORK)
*
*        Factorize A
*
         CALL DPPTRF(UPLO,N,AP,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL DPPCON(UPLO,N,AP,ANORM,RCOND,WORK,IWORK,INFO)
*
            IF (RCOND.GE.DLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0D0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'A is not positive-definite'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,D10.2)
      END
