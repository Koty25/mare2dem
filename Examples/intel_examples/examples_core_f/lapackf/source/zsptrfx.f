*     ZSPTRF (F07QRF) Example Program Text
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
      INTEGER          NMAX
      PARAMETER        (NMAX=8)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      COMPLEX*16       AP(NMAX*(NMAX+1)/2)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_TRI_MATRIX
      EXTERNAL         ZSPTRF
*     .. Executable Statements ..
      WRITE (NOUT,*) 'ZSPTRF Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            READ (NIN,*) ((AP(I+J*(J-1)/2),J=I,N),I=1,N)
         ELSE IF (UPLO.EQ.'L') THEN
            READ (NIN,*) ((AP(I+(2*N-J)*(J-1)/2),J=1,I),I=1,N)
         END IF
*
*        Factorize A
*
         CALL ZSPTRF(UPLO,N,AP,IPIV,INFO)
*
         WRITE (NOUT,*)
*
*        Print details of factorization
*
         CALL PRINT_TRI_MATRIX( 'Details of factorization', UPLO, N,
     =                          AP )
*
*        Print pivot indices
*
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'IPIV'
         WRITE (NOUT,99999) (IPIV(I),I=1,N)
*
         IF (INFO.NE.0) WRITE (NOUT,*) 'The factor D is singular'
*
      END IF
*
*
      STOP
*
99999 FORMAT ((1X,I12,3I18))
      END
*
*     End of ZSPTRF Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_TRI_MATRIX( DESC, UPLO, N, A )
      CHARACTER*(*)    DESC
      CHARACTER        UPLO
      INTEGER          N
      COMPLEX*16       A( * )
      COMPLEX*16       A_TRI( N, N )
*
      COMPLEX*16       ZERO
      PARAMETER        (ZERO=(0.0D0,0.0D0))
*
      EXTERNAL         ZLASET
      INTEGER          I, J, K
*
      IF (UPLO.EQ.'U') THEN
         DO J = 1, N
            DO I = 1, J
               K = J * (J-1) / 2
               A_TRI( I, J ) = A( I + K )
            END DO
         END DO
         CALL ZLASET( 'L', N-1, N-1, ZERO, ZERO,
     +                A_TRI(2,1), N)
      ELSE IF (UPLO.EQ.'L') THEN
         DO I = 1, N
            DO J = 1, I
               K = (2*N-J)*(J-1)/2
               A_TRI( I, J ) = A( I + K )
            END DO
         END DO
         CALL ZLASET( 'U', N-1, N-1, ZERO, ZERO,
     +                A_TRI(1,2), N)
      END IF
*
      WRITE(*,*) DESC
      WRITE(*, 9999) ( J, J = 1, N)
      DO I = 1, N
         WRITE(*,9998) I, ( A_TRI( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( I2, ' ', 11(:,1X,'(',F7.4,',',F7.4,')') )
 9999 FORMAT( '   ', 11(:,1X,I17) )
*
      RETURN
      END
