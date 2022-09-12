*     ZUNGBR (F08KTF) Example Program Text
*     Mark 16 Release. NAG Copyright 1992.
*
********************************************
*                                          *
* Modified by Intel Corporation, July 2017 *
*                                          *
********************************************
*
*     .. Parameters ..
      INTEGER          NIN, NOUT
      PARAMETER        (NIN=5,NOUT=6)
      INTEGER          MMAX, NMAX, LDA, LDVT, LDU, LDC, LWORK
      PARAMETER        (MMAX=8,NMAX=8,LDA=MMAX,LDVT=NMAX,LDU=MMAX,
     +                 LDC=NMAX,LWORK=64*(MMAX+NMAX))
*     .. Local Scalars ..
      INTEGER          I, IC, INFO, J, M, N
*     .. Local Arrays ..
      COMPLEX*16       A(LDA,NMAX), C(LDC,NMAX), TAUP(NMAX), TAUQ(NMAX),
     +                 U(LDU,NMAX), VT(LDVT,NMAX), WORK(LWORK)
      DOUBLE PRECISION D(NMAX), E(NMAX-1), RWORK(4*NMAX-4)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         ZLACPY, ZBDSQR, ZGEBRD, ZUNGBR
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZUNGBR Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      DO 20 IC = 1, 2
         READ (NIN,*) M, N
         IF (M.LE.MMAX .AND. N.LE.NMAX) THEN
*
*           Read A from data file
*
            READ (NIN,*) ((A(I,J),J=1,N),I=1,M)
*
*           Reduce A to bidiagonal form
*
            CALL ZGEBRD(M,N,A,LDA,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
*
            IF (M.GE.N) THEN
*
*              Copy A to VT and U
*
               CALL ZLACPY('Upper',N,N,A,LDA,VT,LDVT)
*
               CALL ZLACPY('Lower',M,N,A,LDA,U,LDU)
*
*              Form P**H explicitly, storing the result in VT
*
               CALL ZUNGBR('P',N,N,M,VT,LDVT,TAUP,WORK,LWORK,INFO)
*
*              Form Q explicitly, storing the result in U
*
               CALL ZUNGBR('Q',M,N,N,U,LDU,TAUQ,WORK,LWORK,INFO)
*
*              Compute the SVD of A
*
               CALL ZBDSQR('Upper',N,N,M,0,D,E,VT,LDVT,U,LDU,C,LDC,
     +                     RWORK,INFO)
*
*              Print singular values, left & right singular vectors
*
               WRITE (NOUT,*)
               WRITE (NOUT,*) 'Example 1: singular values'
               WRITE (NOUT,99999) (D(I),I=1,N)
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Right singular vectors, by row', N,
     +                            N, VT, LDVT )

*
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Left singular vectors, by column',
     +                            M, N, U, LDU )

*
            ELSE
*
*              Copy A to VT and U
*
               CALL ZLACPY('Upper',M,N,A,LDA,VT,LDVT)
*
               CALL ZLACPY('Lower',M,M,A,LDA,U,LDU)
*
*              Form P**H explicitly, storing the result in VT
*
               CALL ZUNGBR('P',M,N,M,VT,LDVT,TAUP,WORK,LWORK,INFO)
*
*              Form Q explicitly, storing the result in U
*
               CALL ZUNGBR('Q',M,M,N,U,LDU,TAUQ,WORK,LWORK,INFO)
*
*              Compute the SVD of A
*
               CALL ZBDSQR('Lower',M,N,M,0,D,E,VT,LDVT,U,LDU,C,LDC,
     +                     RWORK,INFO)
*
*              Print singular values, left & right singular vectors
*
               WRITE (NOUT,*)
               WRITE (NOUT,*) 'Example 2: singular values'
               WRITE (NOUT,99999) (D(I),I=1,M)
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Right singular vectors, by row', M,
     +                            N, VT, LDVT )

*
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Left singular vectors, by column',
     +                            M, M, U, LDU )

*
            END IF
         END IF
   20 CONTINUE
*
*
      STOP
*
99999 FORMAT (8X,4(F7.4,11X,:))
      END
*
*     End of ZUNGBR Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX*16       A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*) DESC
      WRITE(*, 9999) ( J, J = 1, N)
      DO I = 1, M
         WRITE(*,9998) I, ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( I2, ' ', 11(:,1X,'(',F7.4,',',F7.4,')') )
 9999 FORMAT( '   ', 11(:,1X,I17) )
*
      RETURN
      END
