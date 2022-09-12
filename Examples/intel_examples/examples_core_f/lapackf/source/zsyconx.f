*     ZSYCON (F07NUF) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDA, LWORK
      PARAMETER        (NMAX=8,LDA=NMAX,LWORK=64*NMAX)
*     .. Local Scalars ..
      DOUBLE PRECISION ANORM, RCOND
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX*16       A(LDA,NMAX), WORK(LWORK)
      DOUBLE PRECISION RWORK(NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Functions ..
      DOUBLE PRECISION ZLANSY, DLAMCH
      EXTERNAL         ZLANSY, DLAMCH
*     .. External Subroutines ..
      EXTERNAL         ZSYCON, ZSYTRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZSYCON Example Program Results'
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
         ANORM = ZLANSY('1-norm',UPLO,N,A,LDA,RWORK)
*
*        Factorize A
*
         CALL ZSYTRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL ZSYCON(UPLO,N,A,LDA,IPIV,ANORM,RCOND,WORK,INFO)
*
            IF (RCOND.GE.DLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0D0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'The factor D is singular'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,D10.2)
      END
