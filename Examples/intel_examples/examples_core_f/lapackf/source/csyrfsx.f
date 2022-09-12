*     CSYRFS (F07NVE) Example Program Text
*     Mark 15 Release. NAG Copyright 1991.
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
      INTEGER          NMAX, NRHMAX, LDA, LWORK, LDAF, LDB, LDX
      PARAMETER        (NMAX=8,NRHMAX=NMAX,LDA=NMAX,LWORK=64*NMAX,
     +                 LDAF=NMAX,LDB=NMAX,LDX=NMAX)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N, NRHS
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX          A(LDA,NMAX), AF(LDAF,NMAX), B(LDB,NRHMAX),
     +                 WORK(LWORK), X(LDX,NMAX)
      REAL             BERR(NRHMAX), FERR(NRHMAX), RWORK(NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         CSYRFS, CSYTRF, CSYTRS, CLACPY
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CSYRFS Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, NRHS
      IF (N.LE.NMAX .AND. NRHS.LE.NRHMAX) THEN
*
*        Read A and B from data file, and copy A to AF and B to X
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((A(I,J),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((A(I,J),J=1,I),I=1,N)
         END IF
         READ (NIN,*) ((B(I,J),J=1,NRHS),I=1,N)
         CALL CLACPY(UPLO,N,N,A,LDA,AF,LDAF)
         CALL CLACPY('General',N,NRHS,B,LDB,X,LDX)
*
*        Factorize A in the array AF
*
         CALL CSYTRF(UPLO,N,AF,LDAF,IPIV,WORK,LWORK,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Compute solution in the array X
*
            CALL CSYTRS(UPLO,N,NRHS,AF,LDAF,IPIV,X,LDX,INFO)
*
*           Improve solution, and compute backward errors and
*           estimated bounds on the forward errors
*
            CALL CSYRFS(UPLO,N,NRHS,A,LDA,AF,LDAF,IPIV,B,LDB,X,LDX,FERR,
     +                  BERR,WORK,RWORK,INFO)
*
*           Print solution
*
            CALL PRINT_MATRIX( 'Solution(s)', N, NRHS, X, LDX )
            WRITE (NOUT,*)
            WRITE (NOUT,*) 'Backward errors (machine-dependent)'
            WRITE (NOUT,99999) (BERR(J),J=1,NRHS)
            WRITE (NOUT,*)
     +        'Estimated forward error bounds (machine-dependent)'
            WRITE (NOUT,99999) (FERR(J),J=1,NRHS)
         ELSE
            WRITE (NOUT,*) 'The factor D is singular'
         END IF
      END IF
*
*
      STOP
*
99999 FORMAT ((5X,1P,4(E11.1,7X)))
      END
*
*     End of CSYRFS Example
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
