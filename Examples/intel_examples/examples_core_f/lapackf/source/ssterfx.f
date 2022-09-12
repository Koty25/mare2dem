*     SSTERF (F08JFE) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
*     .. Local Scalars ..
      INTEGER          I, INFO, N
*     .. Local Arrays ..
      REAL             D(NMAX), E(NMAX-1)
*     .. External Subroutines ..
      EXTERNAL         SSTERF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SSTERF Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read T from data file
*
         READ (NIN,*) (D(I),I=1,N)
         READ (NIN,*) (E(I),I=1,N-1)
*
*        Calculate the eigenvalues of T
*
         CALL SSTERF(N,D,E,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.GT.0) THEN
            WRITE (NOUT,*) 'Failure to converge.'
         ELSE
            WRITE (NOUT,*) 'Eigenvalues'
            WRITE (NOUT,99999) (D(I),I=1,N)
         END IF
      END IF
      STOP
*
99999 FORMAT (3X,(9F8.4))
      END
