*     DGECON (F07AGF) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDA
      PARAMETER        (NMAX=8,LDA=NMAX)
      CHARACTER        NORM
      PARAMETER        (NORM='1')
*     .. Local Scalars ..
      DOUBLE PRECISION ANORM, RCOND
      INTEGER          I, INFO, J, N
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), WORK(4*NMAX)
      INTEGER          IPIV(NMAX), IWORK(NMAX)
*     .. External Functions ..
      DOUBLE PRECISION DLANGE, DLAMCH
      EXTERNAL         DLANGE, DLAMCH
*     .. External Subroutines ..
      EXTERNAL         DGECON, DGETRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'DGECON Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Compute norm of A
*
         ANORM = DLANGE(NORM,N,N,A,LDA,WORK)
*
*        Factorize A
*
         CALL DGETRF(N,N,A,LDA,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Estimate condition number
*
            CALL DGECON(NORM,N,A,LDA,ANORM,RCOND,WORK,IWORK,INFO)
*
            IF (RCOND.GE.DLAMCH('E')) THEN
               WRITE (NOUT,99999) 'Estimate of condition number =',
     +           1.0D0/RCOND
            ELSE
               WRITE (NOUT,*) 'A is singular to working precision'
            END IF
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
         END IF
      END IF
      STOP
*
99999 FORMAT (1X,A,1P,D10.2)
      END
