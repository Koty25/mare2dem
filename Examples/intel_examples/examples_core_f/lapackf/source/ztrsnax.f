*     ZTRSNA (F08QYF) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDT, LDWORK, LDVL, LDVR
      PARAMETER        (NMAX=8,LDT=NMAX,LDWORK=NMAX,LDVL=NMAX,LDVR=NMAX)
*     .. Local Scalars ..
      DOUBLE PRECISION EPS, TNORM
      INTEGER          I, INFO, J, M, N
*     .. Local Arrays ..
      COMPLEX*16       T(LDT,NMAX), VL(LDVL,NMAX), VR(LDVR,NMAX),
     +                 WORK(LDWORK,NMAX+6)
      DOUBLE PRECISION RWORK(NMAX), S(NMAX), SEP(NMAX)
      LOGICAL          SELECT(1)
*     .. External Functions ..
      DOUBLE PRECISION ZLANGE, DLAMCH
      EXTERNAL         ZLANGE, DLAMCH
*     .. External Subroutines ..
      EXTERNAL         ZTREVC, ZTRSNA
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZTRSNA Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read T from data file
*
         READ (NIN,*) ((T(I,J),J=1,N),I=1,N)
*
*        Calculate the left and right eigenvectors of T
*
         CALL ZTREVC('Both','All',SELECT,N,T,LDT,VL,LDVL,VR,LDVR,N,M,
     +               WORK,RWORK,INFO)
*
*        Estimate condition numbers for all the eigenvalues and right
*        eigenvectors of T
*
         CALL ZTRSNA('Both','All',SELECT,N,T,LDT,VL,LDVL,VR,LDVR,S,SEP,
     +               N,M,WORK,LDWORK,RWORK,INFO)
*
*        Print condition numbers of eigenvalues and right eigenvectors
*
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'S'
         WRITE (NOUT,99999) (S(I),I=1,M)
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'SEP'
         WRITE (NOUT,99999) (SEP(I),I=1,M)
*
*        Calculate approximate error estimates (using the 1-norm)
*
         EPS = DLAMCH('E')
         TNORM = ZLANGE('1-norm',N,N,T,LDT,RWORK)
         WRITE (NOUT,*)
         WRITE (NOUT,*)
     +     'Approximate error estimates for eigenvalues ',
     +     'of T (machine-dependent)'
         WRITE (NOUT,99999) (EPS*TNORM/S(I),I=1,M)
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'Approximate error estimates for right ',
     +     'eigenvectors of T (machine-dependent)'
         WRITE (NOUT,99999) (EPS*TNORM/SEP(I),I=1,M)
      END IF
      STOP
*
99999 FORMAT ((3X,1P,7D11.1))
      END
