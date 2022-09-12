*     DSPTRD (F08GEF) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      DOUBLE PRECISION AP(NMAX*(NMAX+1)/2), D(NMAX), E(NMAX-1),
     +                 TAU(NMAX-1)
*     .. External Subroutines ..
      EXTERNAL         DSPTRD
*     .. Executable Statements ..
      WRITE (NOUT,*) 'DSPTRD Example Program Results'
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
*        Reduce A to tridiagonal form
*
         CALL DSPTRD(UPLO,N,AP,D,E,TAU,INFO)
*
*        Print tridiagonal form
*
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'Diagonal'
         WRITE (NOUT,99999) (D(I),I=1,N)
         WRITE (NOUT,*) 'Off-diagonal'
         WRITE (NOUT,99999) (E(I),I=1,N-1)
      END IF
      STOP
*
99999 FORMAT (1X,8F9.4)
      END
