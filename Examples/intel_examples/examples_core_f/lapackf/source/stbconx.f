*     STBCON (F07VGE) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, KDMAX, LDAB
      PARAMETER        (NMAX=8,KDMAX=NMAX,LDAB=KDMAX+1)
      CHARACTER        NORM, DIAG
      PARAMETER        (NORM='1',DIAG='N')
*     .. Local Scalars ..
      REAL             RCOND
      INTEGER          I, INFO, J, KD, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      REAL             AB(LDAB,NMAX), WORK(3*NMAX)
      INTEGER          IWORK(NMAX)
*     .. External Functions ..
      REAL             SLAMCH
      EXTERNAL         SLAMCH
*     .. External Subroutines ..
      EXTERNAL         STBCON
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'STBCON Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, KD
      IF (N.LE.NMAX .AND. KD.LE.KDMAX) THEN
*
*        Read A from data file
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
*
*        Estimate condition number
*
         CALL STBCON(NORM,UPLO,DIAG,N,KD,AB,LDAB,RCOND,WORK,IWORK,INFO)
*
         WRITE (NOUT,*)
         IF (RCOND.GE.SLAMCH('E')) THEN
            WRITE (NOUT,99999) 'Estimate of condition number =',
     +        1.0/RCOND
         ELSE
            WRITE (NOUT,*) 'A is singular to working precision'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,E10.2)
      END
