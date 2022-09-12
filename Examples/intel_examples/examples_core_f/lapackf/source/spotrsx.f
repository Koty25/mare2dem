*     SPOTRS (F07FEE) Example Program Text
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
      INTEGER          NMAX, LDA, NRHMAX, LDB
      PARAMETER        (NMAX=8,LDA=NMAX,NRHMAX=NMAX,LDB=NMAX)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N, NRHS
      CHARACTER        UPLO
*     .. Local Arrays ..
      REAL             A(LDA,NMAX), B(LDB,NRHMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         SPOTRF, SPOTRS
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SPOTRS Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, NRHS
      IF (N.LE.NMAX .AND. NRHS.LE.NRHMAX) THEN
*
*        Read A and B from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((A(I,J),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((A(I,J),J=1,I),I=1,N)
         END IF
         READ (NIN,*) ((B(I,J),J=1,NRHS),I=1,N)
*
*        Factorize A
*
         CALL SPOTRF(UPLO,N,A,LDA,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Compute solution
*
            CALL SPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO)
*
*           Print solution
*
            CALL PRINT_MATRIX( 'Solution(s)', N, NRHS, B, LDB )
         ELSE
            WRITE (NOUT,*) 'A is not positive-definite'
         END IF
      END IF
*
      STOP
      END
*
*     End of SPOTRS Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      REAL             A( LDA, * )
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