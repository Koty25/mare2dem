*     ZUNGLQ (F08AWF) Example Program Text
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
      INTEGER          MMAX, NMAX, LDA, LWORK
      PARAMETER        (MMAX=8,NMAX=8,LDA=MMAX,LWORK=64*MMAX)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, M, N
      CHARACTER*30     TITLE
*     .. Local Arrays ..
      COMPLEX*16       A(LDA,NMAX), TAU(NMAX), WORK(LWORK)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         ZGELQF, ZUNGLQ
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZUNGLQ Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) M, N
      IF (M.LE.MMAX .AND. N.LE.NMAX .AND. M.LE.N) THEN
*
*        Read A from data file
*
         READ (NIN,*) ((A(I,J),J=1,N),I=1,M)
*
*        Compute the LQ factorization of A
*
         CALL ZGELQF(M,N,A,LDA,TAU,WORK,LWORK,INFO)
*
*        Form the leading M rows of Q explicitly
*
         CALL ZUNGLQ(M,N,M,A,LDA,TAU,WORK,LWORK,INFO)
*
*        Print the leading M rows of Q only
*
         WRITE (NOUT,*)
         WRITE (TITLE,99999) M
*
         CALL PRINT_MATRIX( TITLE, M, N, A, LDA )
      END IF
*
*
      STOP
*
99999 FORMAT ('The leading ',I2,' rows of Q')
      END
*
*     End of ZUNGLQ Example
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