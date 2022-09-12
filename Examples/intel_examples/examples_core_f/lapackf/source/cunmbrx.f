*     CUNMBR (F08KUE) Example Program Text
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
      INTEGER          MMAX, NMAX, LDA, LDPH, LDU, LWORK
      PARAMETER        (MMAX=8,NMAX=8,LDA=MMAX,LDPH=NMAX,LDU=MMAX,
     +                 LWORK=64*(MMAX+NMAX))
      COMPLEX          ZERO
      PARAMETER        (ZERO=(0.0,0.0))
*     .. Local Scalars ..
      INTEGER          I, IC, INFO, J, M, N
*     .. Local Arrays ..
      COMPLEX          A(LDA,NMAX), PH(LDPH,NMAX), TAU(NMAX),
     +                 TAUP(NMAX), TAUQ(NMAX), U(LDU,NMAX), WORK(LWORK)
      REAL             D(NMAX), E(NMAX-1)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         CGEBRD, CGELQF, CGEQRF, CUNGLQ, CUNGQR, CUNMBR,
     +                 CLACPY, CLASET
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CUNMBR Example Program Results'
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
            IF (M.GE.N) THEN
*
*              Compute the QR factorization of A
*
               CALL CGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*              Copy A to U
*
               CALL CLACPY('Lower',M,N,A,LDA,U,LDU)
*
*              Form Q explicitly, storing the result in U
*
               CALL CUNGQR(M,M,N,U,LDU,TAU,WORK,LWORK,INFO)
*
*              Copy R to PH (used as workspace)
*
               CALL CLACPY('Upper',N,N,A,LDA,PH,LDPH)
*
*              Set the strictly lower triangular part of R to zero
*
               CALL CLASET('Lower',N-1,N-1,ZERO,ZERO,PH(2,1),LDPH)
*
*              Bidiagonalize R
*
               CALL CGEBRD(N,N,PH,LDPH,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
*
*              Update Q, storing the result in U
*
               CALL CUNMBR('Q','Right','No transpose',M,N,N,PH,LDPH,
     +                     TAUQ,U,LDU,WORK,LWORK,INFO)
*
*              Print bidiagonal form and matrix Q
*
               WRITE (NOUT,*)
               WRITE (NOUT,*) 'Example 1: bidiagonal matrix B'
               WRITE (NOUT,*) 'Diagonal'
               WRITE (NOUT,99999) (D(I),I=1,N)
               WRITE (NOUT,*) 'Super-diagonal'
               WRITE (NOUT,99999) (E(I),I=1,N-1)
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Example 1: matrix Q', M, N, U, LDU )
*
            ELSE
*
*              Compute the LQ factorization of A
*
               CALL CGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*              Copy A to PH
*
               CALL CLACPY('Upper',M,N,A,LDA,PH,LDPH)
*
*              Form Q explicitly, storing the result in PH
*
               CALL CUNGLQ(N,N,M,PH,LDPH,TAU,WORK,LWORK,INFO)
*
*              Copy L to U (used as workspace)
*
               CALL CLACPY('Lower',M,M,A,LDA,U,LDU)
*
*              Set the strictly upper triangular part of L to zero
*
               CALL CLASET('Upper',M-1,M-1,ZERO,ZERO,U(1,2),LDU)
*
*              Bidiagonalize L
*
               CALL CGEBRD(M,M,U,LDU,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
*
*              Update P**H, storing the result in PH
*
               CALL CUNMBR('P','Left','Conjugate transpose',M,N,M,U,LDU,
     +                     TAUP,PH,LDPH,WORK,LWORK,INFO)
*
*              Print bidiagonal form and matrix P**H
*
               WRITE (NOUT,*)
               WRITE (NOUT,*) 'Example 2: bidiagonal matrix B'
               WRITE (NOUT,*) 'Diagonal'
               WRITE (NOUT,99999) (D(I),I=1,M)
               WRITE (NOUT,*) 'Super-diagonal'
               WRITE (NOUT,99999) (E(I),I=1,M-1)
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Example 2: matrix P**H', M, N, PH,
     +                            LDPH )

*
            END IF
         END IF
   20 CONTINUE
*
*
      STOP
*
99999 FORMAT (3X,(8F8.4))
      END
*
*     End of CUNMBR Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      COMPLEX          A( LDA, * )
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
