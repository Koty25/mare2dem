*     CGBCON (F07BUE) Example Program Text
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
      COMPLEX          AB(LDAB,NMAX), WORK(2*NMAX)
      REAL             RWORK(NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Functions ..
      REAL             CLANGB, SLAMCH
      EXTERNAL         CLANGB, SLAMCH
*     .. External Subroutines ..
      EXTERNAL         CGBCON, CGBTRF
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CGBCON Example Program Results'
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
         ANORM = CLANGB(NORM,N,KL,KU,AB(KL+1,1),LDAB,RWORK)
*
*        Factorize A
*
         CALL CGBTRF(N,N,KL,KU,AB,LDAB,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL CGBCON(NORM,N,KL,KU,AB,LDAB,IPIV,ANORM,RCOND,WORK,
     +                  RWORK,INFO)
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
