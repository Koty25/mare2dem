*     CSPTRI (F07QWE) Example Program Text
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
      COMPLEX          AP(NMAX*(NMAX+1)/2), WORK(NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_TRI_MATRIX
      EXTERNAL         CSPTRF, CSPTRI
*     .. Executable Statements ..
      WRITE (NOUT,*) 'CSPTRI Example Program Results'
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
         CALL CSPTRF(UPLO,N,AP,IPIV,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Compute inverse of A
*
            CALL CSPTRI(UPLO,N,AP,IPIV,WORK,INFO)
*
*           Print inverse
*
            CALL PRINT_TRI_MATRIX( 'Inverse', UPLO, N, AP )
         ELSE
            WRITE (NOUT,*) 'The factor D is singular'
         END IF
      END IF
*
      STOP
      END
*
*     End of CSPTRI Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_TRI_MATRIX( DESC, UPLO, N, A )
      CHARACTER*(*)    DESC
      CHARACTER        UPLO
      INTEGER          N
      COMPLEX          A( * )
      COMPLEX          A_TRI( N, N )
*
      COMPLEX          ZERO
      PARAMETER        (ZERO=(0.0,0.0))
*
      EXTERNAL         CLASET
      INTEGER          I, J, K
*
      IF (UPLO.EQ.'U') THEN
         DO J = 1, N
            DO I = 1, J
               K = J * (J-1) / 2
               A_TRI( I, J ) = A( I + K )
            END DO
         END DO
         CALL CLASET( 'L', N-1, N-1, ZERO, ZERO,
     +                A_TRI(2,1), N)
      ELSE IF (UPLO.EQ.'L') THEN
         DO I = 1, N
            DO J = 1, I
               K = (2*N-J)*(J-1)/2
               A_TRI( I, J ) = A( I + K )
            END DO
         END DO
         CALL CLASET( 'U', N-1, N-1, ZERO, ZERO,
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
