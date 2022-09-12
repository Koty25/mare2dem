*     CGEBAL (F08NVE) Example Program Text
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
      INTEGER          NMAX, LDA, LDH, LWORK, LDVL, LDVR
      PARAMETER        (NMAX=8,LDA=NMAX,LDH=NMAX,LWORK=64*NMAX,LDVL=1,
     +                 LDVR=NMAX)
*     .. Local Scalars ..
      INTEGER          I, IHI, ILO, INFO, J, M, N
*     .. Local Arrays ..
      COMPLEX          A(LDA,NMAX), H(LDH,NMAX), TAU(NMAX), VL(LDVL,1),
     +                 VR(LDVR,NMAX), W(NMAX), WORK(LWORK)
      REAL             RWORK(NMAX), SCALE(NMAX)
      LOGICAL          SELECT(1)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         CGEBAK, CGEBAL, CGEHRD, CHSEQR, CTREVC, CUNGHR,
     +                 CLACPY
*     .. Intrinsic Functions ..
      INTRINSIC        AIMAG, REAL
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CGEBAL Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Balance A
*
         CALL CGEBAL('Both',N,A,LDA,ILO,IHI,SCALE,INFO)
*
*        Reduce A to upper Hessenberg form H = (Q**H)*A*Q
*
         CALL CGEHRD(N,ILO,IHI,A,LDA,TAU,WORK,LWORK,INFO)
*
*        Copy A to H
*
         CALL CLACPY('General',N,N,A,LDA,H,LDH)
*
*        Copy A into VR
*
         CALL CLACPY('General',N,N,A,LDA,VR,LDVR)
*
*        Form Q explicitly, storing the result in VR
*
         CALL CUNGHR(N,1,N,VR,LDVR,TAU,WORK,LWORK,INFO)
*
*        Calculate the eigenvalues and Schur factorization of A
*
         CALL CHSEQR('Schur form','Vectors',N,ILO,IHI,H,LDH,W,VR,LDVR,
     +               WORK,LWORK,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.GT.0) THEN
            WRITE (NOUT,*) 'Failure to converge.'
         ELSE
            WRITE (NOUT,*) 'Eigenvalues'
*            WRITE (NOUT,99999) (' (',REAL(W(I)),',',AIMAG(W(I)),')',
*     *             I=1,N)
             WRITE (NOUT,99999) 
     *             (' (',REAL(W(I)),',',AIMAG(W(I)),')', I=1,N)
*
*           Calculate the eigenvectors of A, storing the result in VR
*
            CALL CTREVC('Right','Back',SELECT,N,H,LDH,VL,LDVL,VR,
     +                  LDVR,N,M,WORK,RWORK,INFO)
*
*           Backtransform eigenvectors
*
            CALL CGEBAK('Both','Right',N,ILO,IHI,SCALE,M,VR,LDVR,INFO)
*
*           Print eigenvectors
*
            WRITE (NOUT,*)
*
            CALL PRINT_MATRIX( 'Contents of array VR', N, M, VR, LDVR )
*
         END IF
      END IF
      STOP
*
99999 FORMAT ((3X,4(A,F7.4,A,F7.4,A,:)))
*
      END
*
*     End of CGEBAL Example
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
