*     DORMBR (F08KGF) Example Program Text
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
      INTEGER          MMAX, NMAX, LDA, LDPT, LDU, LWORK
      PARAMETER        (MMAX=8,NMAX=8,LDA=MMAX,LDPT=NMAX,LDU=MMAX,
     +                 LWORK=64*(MMAX+NMAX))
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
*     .. Local Scalars ..
      INTEGER          I, IC, INFO, J, M, N
*     .. Local Arrays ..
      DOUBLE PRECISION A(LDA,NMAX), D(NMAX), E(NMAX-1), PT(LDPT,NMAX),
     +                 TAU(NMAX), TAUP(NMAX), TAUQ(NMAX), U(LDU,NMAX),
     +                 WORK(LWORK)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         DGEBRD, DGELQF, DGEQRF, DORGLQ, DORGQR, DORMBR,
     +                 DLACPY, DLASET
*     .. Executable Statements ..
      WRITE (NOUT,*) 'DORMBR Example Program Results'
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
               CALL DGEQRF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*              Copy A to U
*
               CALL DLACPY('Lower',M,N,A,LDA,U,LDU)
*
*              Form Q explicitly, storing the result in U
*
               CALL DORGQR(M,M,N,U,LDU,TAU,WORK,LWORK,INFO)
*
*              Copy R to PT (used as workspace)
*
               CALL DLACPY('Upper',N,N,A,LDA,PT,LDPT)
*
*              Set the strictly lower triangular part of R to zero
*
               CALL DLASET('Lower',N-1,N-1,ZERO,ZERO,PT(2,1),LDPT)
*
*              Bidiagonalize R
*
               CALL DGEBRD(N,N,PT,LDPT,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
*
*              Update Q, storing the result in U
*
               CALL DORMBR('Q','Right','No transpose',M,N,N,PT,LDPT,
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
            ELSE
*
*              Compute the LQ factorization of A
*
               CALL DGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*              Copy A to PT
*
               CALL DLACPY('Upper',M,N,A,LDA,PT,LDPT)
*
*              Form Q explicitly, storing the result in PT
*
               CALL DORGLQ(N,N,M,PT,LDPT,TAU,WORK,LWORK,INFO)
*
*              Copy L to U (used as workspace)
*
               CALL DLACPY('Lower',M,M,A,LDA,U,LDU)
*
*              Set the strictly upper triangular part of L to zero
*
               CALL DLASET('Upper',M-1,M-1,ZERO,ZERO,U(1,2),LDU)
*
*              Bidiagonalize L
*
               CALL DGEBRD(M,M,U,LDU,D,E,TAUQ,TAUP,WORK,LWORK,INFO)
*
*              Update P**T, storing the result in PT
*
               CALL DORMBR('P','Left','Transpose',M,N,M,U,LDU,TAUP,PT,
     +                     LDPT,WORK,LWORK,INFO)
*
*              Print bidiagonal form and matrix P**T
*
               WRITE (NOUT,*)
               WRITE (NOUT,*) 'Example 2: bidiagonal matrix B'
               WRITE (NOUT,*) 'Diagonal'
               WRITE (NOUT,99999) (D(I),I=1,M)
               WRITE (NOUT,*) 'Super-diagonal'
               WRITE (NOUT,99999) (E(I),I=1,M-1)
               WRITE (NOUT,*)
*
               CALL PRINT_MATRIX( 'Example 2: matrix P**T', M, N, PT,
     +                            LDPT )

            END IF
         END IF
   20 CONTINUE
*
99999 FORMAT (3X,(8F8.4))
*
      STOP
      END
*
*     End of DORMBR Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*) DESC
      WRITE(*, 9999) ( J, J = 1, N)
      DO I = 1, M
         WRITE(*, 9998) I, ( A( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( I2, ' ', 11(:,1X,F10.4) )
 9999 FORMAT( '   ', 11(:,1X,I10) )
*
      RETURN
      END