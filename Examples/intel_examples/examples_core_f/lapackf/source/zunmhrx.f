*     ZUNMHR (F08NUF) Example Program Text
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
      INTEGER          NMAX, LDA, LDH, LDZ, LWORK, LDVL, LDVR
      PARAMETER        (NMAX=8,LDA=NMAX,LDH=NMAX,LDZ=1,LWORK=64*NMAX,
     +                 LDVL=NMAX,LDVR=NMAX)
*     .. Local Scalars ..
      DOUBLE PRECISION THRESH
      INTEGER          I, INFO, J, M, N
*     .. Local Arrays ..
      COMPLEX*16       A(LDA,NMAX), H(LDH,NMAX), TAU(NMAX),
     +                 VL(LDVL,NMAX), VR(LDVR,NMAX), W(NMAX),
     +                 WORK(LWORK), Z(LDZ,1)
      DOUBLE PRECISION RWORK(NMAX)
      INTEGER          IFAILL(NMAX), IFAILR(NMAX)
      LOGICAL          SELECT(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         ZLACPY, ZGEHRD, ZHSEIN, ZHSEQR, ZUNMHR
*     .. Intrinsic Functions ..
      INTRINSIC        DBLE, DIMAG
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZUNMHR Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
         READ (NIN,*) THRESH
*
*        Reduce A to upper Hessenberg form H = (Q**H)*A*Q
*
         CALL ZGEHRD(N,1,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*        Copy A to H
*
         CALL ZLACPY('General',N,N,A,LDA,H,LDH)
*
*        Calculate the eigenvalues of H (same as A)
*
         CALL ZHSEQR('Eigenvalues','No vectors',N,1,N,H,LDH,W,Z,LDZ,
     +               WORK,LWORK,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.GT.0) THEN
            WRITE (NOUT,*) 'Failure to converge.'
         ELSE
            WRITE (NOUT,*) 'Eigenvalues'
*            WRITE (NOUT,99999) (' (',DBLE(W(I)),',',DIMAG(W(I)),')',
*     +        I=1,N)
            WRITE (NOUT,99999) 
     *            (' (',DBLE(W(I)),',',DIMAG(W(I)),')', I=1,N)
*
            DO 20 I = 1, N
               SELECT(I) = DBLE(W(I)) .LT. THRESH
   20       CONTINUE
*
*           Calculate the eigenvectors of H (as specified by SELECT),
*           storing the result in VR
*
            CALL ZHSEIN('Right','QR','No initial vectors',SELECT,N,A,
     +                  LDA,W,VL,LDVL,VR,LDVR,N,M,WORK,RWORK,IFAILL,
     +                  IFAILR,INFO)
*
*           Calculate the eigenvectors of A = Q * (eigenvectors of H)
*
            CALL ZUNMHR('Left','No transpose',N,M,1,N,A,LDA,TAU,VR,LDVR,
     +                  WORK,LWORK,INFO)
*
*           Print eigenvectors
*
            WRITE (NOUT,*)
*
            CALL PRINT_MATRIX( 'Contents of array VR', N, M, VR, LDVR )
*
         END IF
      END IF
*
*
      STOP
*
99999 FORMAT ((3X,4(A,F7.4,A,F7.4,A,:)))
      END
*
*     End of ZUNMHR Example
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