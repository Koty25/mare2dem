*     SGBCON (F07BGE) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, KLMAX, KUMAX, LDAB
      PARAMETER        (NMAX=8,KLMAX=8,KUMAX=8,LDAB=2*KLMAX+KUMAX+1)
      CHARACTER        NORM
      PARAMETER        (NORM='1')
*     .. Local Scalars ..
      REAL             ANORM, RCOND
      INTEGER          I, INFO, J, K, KL, KU, N
*     .. Local Arrays ..
      REAL             AB(LDAB,NMAX), WORK(3*NMAX)
      INTEGER          IPIV(NMAX), IWORK(NMAX)
*     .. External Functions ..
      REAL             SLANGB, SLAMCH
      EXTERNAL         SLANGB, SLAMCH
*     .. External Subroutines ..
      EXTERNAL         SGBCON, SGBTRF
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SGBCON Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, KL, KU
      IF (N.LE.NMAX .AND. KL.LE.KLMAX .AND. KU.LE.KUMAX) THEN
*
*        Read A from data file
*
         K = KL + KU + 1
         READ (NIN,*) ((AB(K+I-J,J),J=MAX(I-KL,1),MIN(I+KU,N)),I=1,
     +     N)
*
*        Compute norm of A
*
         ANORM = SLANGB(NORM,N,KL,KU,AB(KL+1,1),LDAB,WORK)
*
*        Factorize A
*
         CALL SGBTRF(N,N,KL,KU,AB,LDAB,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL SGBCON(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,WORK,
     +                  IWORK,INFO)
*
            IF (RCOND.GE.SLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,E10.2)
      END
