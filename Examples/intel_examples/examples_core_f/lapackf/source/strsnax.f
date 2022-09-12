*     STRSNA (F08QLE) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX, LDT, LDWORK, LDVL, LDVR
      PARAMETER        (NMAX=8,LDT=NMAX,LDWORK=NMAX,LDVL=NMAX,LDVR=NMAX)
*     .. Local Scalars ..
      REAL             EPS, TNORM
      INTEGER          I, INFO, J, M, N
*     .. Local Arrays ..
      REAL             S(NMAX), SEP(NMAX), T(LDT,NMAX), VL(LDVL,NMAX),
     +                 VR(LDVR,NMAX), WORK(LDWORK,NMAX+6)
      INTEGER          IWORK(2*NMAX-2)
      LOGICAL          SELECT(1)
*     .. External Functions ..
      REAL             SLANGE, SLAMCH
      EXTERNAL         SLANGE, SLAMCH
*     .. External Subroutines ..
      EXTERNAL         STREVC, STRSNA
*     .. Executable Statements ..
      WRITE (NOUT,*) 'STRSNA Example Program Results'
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
         CALL STREVC('Both','All',SELECT,N,T,LDT,VL,LDVL,VR,LDVR,N,M,
     +               WORK,INFO)
*
*        Estimate condition numbers for all the eigenvalues and right
*        eigenvectors of T
*
         CALL STRSNA('Both','All',SELECT,N,T,LDT,VL,LDVL,VR,LDVR,S,SEP,
     +               N,M,WORK,LDWORK,IWORK,INFO)
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
         EPS = SLAMCH('E')
         TNORM = SLANGE('1-norm',N,N,T,LDT,WORK)
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
99999 FORMAT ((3X,1P,7E11.1))
      END
