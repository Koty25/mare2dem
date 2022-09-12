*     SGBBRD (F08LEE) Example Program Text
*     Mark 19 Release. NAG Copyright 1998.
C     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          MMAX, NMAX, NCCMAX, KLMAX, KUMAX, LDAB, LDQ,
     *                 LDPT, LDC
      PARAMETER        (MMAX=8,NMAX=8,NCCMAX=8,KLMAX=8,KUMAX=8,
     *                 LDAB=KLMAX+KUMAX+1,LDQ=MMAX,LDPT=NMAX,LDC=MMAX)
      CHARACTER        VECT
      PARAMETER        (VECT='N')
C     .. Local Scalars ..
      INTEGER          I, INFO, J, KL, KU, M, N, NCC
C     .. Local Arrays ..
      REAL             AB(LDAB,NMAX), C(MMAX,NCCMAX), D(NMAX),
     *                 E(NMAX-1), PT(LDPT,NMAX), Q(LDQ,MMAX),
     *                 WORK(2*MMAX+2*NMAX)
C     .. External Subroutines ..
      EXTERNAL         SGBBRD
C     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
C     .. Executable Statements ..
      WRITE (NOUT,FMT=*) 'SGBBRD Example Program Results'
*     Skip heading in data file
      READ (NIN,FMT=*)
      READ (NIN,FMT=*) M, N, KL, KU, NCC
      IF (M.LE.MMAX .AND. N.LE.NMAX .AND. KL.LE.KLMAX .AND. KU.LE.
     *    KUMAX .AND. NCC.LE.NCCMAX) THEN
*
*        Read A from data file
*
         READ (NIN,FMT=*) ((AB(KU+1+I-J,J),J=MAX(I-KL,1),MIN(I+KU,N)),
     *     I=1,M)
*
*        Reduce A to upper bidiagonal form
*
         CALL SGBBRD(VECT,M,N,NCC,KL,KU,AB,LDAB,D,E,Q,LDQ,PT,LDPT,C,LDC,
     *               WORK,INFO)
*
*        Print bidiagonal form
*
         WRITE (NOUT,FMT=*)
         WRITE (NOUT,FMT=*) 'Diagonal'
         WRITE (NOUT,FMT=99999) (D(I),I=1,MIN(M,N))
         WRITE (NOUT,FMT=*) 'Super-diagonal'
         WRITE (NOUT,FMT=99999) (E(I),I=1,MIN(M,N)-1)
      END IF
      STOP
*
99999 FORMAT (1X,8F9.4)
      END
