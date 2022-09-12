*     SGETRI (F07AJE) Example Program Text
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
      INTEGER          NMAX, LDA, LWORK
      PARAMETER        (NMAX=8,LDA=NMAX,LWORK=64*NMAX)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N
*     .. Local Arrays ..
      REAL             A(LDA,NMAX), WORK(LWORK)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         SGETRF, SGETRI
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SGETRI Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,N)
*
*        Factorize A
*
         CALL SGETRF(N,N,A,LDA,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Compute inverse of A
*
            CALL SGETRI(N,A,LDA,IPIV,WORK,LWORK,INFO)
*
*           Print inverse
*
            CALL PRINT_MATRIX( 'Inverse', N, N, A, LDA )
         ELSE
            WRITE (NOUT,*) 'The factor U is singular'
         END IF
      END IF
*
      STOP
      END
*
*     End of SGETRI Example
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
