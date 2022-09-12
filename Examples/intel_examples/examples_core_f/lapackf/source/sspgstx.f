*     SSPGST (F08TEE) Example Program Text
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
      REAL             AP(NMAX*(NMAX+1)/2), BP(NMAX*(NMAX+1)/2),
     +                 D(NMAX), E(NMAX-1), TAU(NMAX)
*     .. External Subroutines ..
      EXTERNAL         SPPTRF, SSPGST, SSPTRD, SSTERF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SSPGST Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A and B from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((AP(I+J*(J-1)/2),J=I,N),I=1,N)
            READ (NIN,*) ((BP(I+J*(J-1)/2),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((AP(I+(2*N-J)*(J-1)/2),J=1,I),I=1,N)
            READ (NIN,*) ((BP(I+(2*N-J)*(J-1)/2),J=1,I),I=1,N)
         END IF
*
*        Compute the Cholesky factorization of B
*
         CALL SPPTRF(UPLO,N,BP,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.GT.0) THEN
            WRITE (NOUT,*) 'B is not positive-definite.'
         ELSE
*
*           Reduce the problem to standard form C*y = lambda*y, storing
*           the result in A
*
            CALL SSPGST(1,UPLO,N,AP,BP,INFO)
*
*           Reduce C to tridiagonal form T = (Q**T)*C*Q
*
            CALL SSPTRD(UPLO,N,AP,D,E,TAU,INFO)
*
*           Calculate the eigenvalues of T (same as C)
*
            CALL SSTERF(N,D,E,INFO)
*
            IF (INFO.GT.0) THEN
               WRITE (NOUT,*) 'Failure to converge in SSTERF.'
            ELSE
*
*              Print eigenvalues
*
               WRITE (NOUT,*) 'Eigenvalues'
               WRITE (NOUT,99999) (D(I),I=1,N)
            END IF
         END IF
      END IF
      STOP
*
99999 FORMAT (3X,(9F8.4))
      END
