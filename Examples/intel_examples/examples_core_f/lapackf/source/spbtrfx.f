*     SPBTRF (F07HDE) Example Program Text
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
      INTEGER          NMAX, KMAX, LDAB
      PARAMETER        (NMAX=8,KMAX=8,LDAB=KMAX+1)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, KD, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      REAL             AB(LDAB,NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_BAND_MATRIX
      EXTERNAL         SPBTRF
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SPBTRF Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N, KD
      IF (N.LE.NMAX .AND. KD.LE.KMAX) THEN
*
*        Read A from data file
*
         READ (NIN,*) UPLO
         IF (UPLO.EQ.'U') THEN
            DO 20 I = 1, N
               READ (NIN,*) (AB(KD+1+I-J,J),J=I,MIN(N,I+KD))
   20       CONTINUE
         ELSE IF (UPLO.EQ.'L') THEN
            DO 40 I = 1, N
               READ (NIN,*) (AB(1+I-J,J),J=MAX(1,I-KD),I)
   40       CONTINUE
         END IF
*
*        Factorize A
*
         CALL SPBTRF(UPLO,N,KD,AB,LDAB,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.EQ.0) THEN
*
*           Print factor
*
*
            IF (UPLO.EQ.'U') THEN
*
               CALL PRINT_BAND_MATRIX( 'Factor', N, N, 0, KD, AB,
     +                                  LDAB )
*
            ELSE IF (UPLO.EQ.'L') THEN
*
               CALL PRINT_BAND_MATRIX( 'Factor', N, N, KD, 0, AB,
     +                                  LDAB )
*
            END IF
*
         ELSE
            WRITE (NOUT,*) 'A is not positive-definite'
         END IF
      END IF
      STOP
*
      END
*
*     End of SPBTRF Example
*
*  =============================================================================
*
*     Auxiliary routine: printing a banded matrix stored in packed form.
*
      SUBROUTINE PRINT_BAND_MATRIX( DESC, M, N, KL, KU, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, KL, KU, LDA
      REAL             A( LDA, * )
      REAL             A_B( M, N )
      REAL             ZERO
      PARAMETER        (ZERO=0.0)
*
      INTEGER          I, J, LB, UB
*
      DO J = 1, N
         LB = MAX(1, J - KU)
         UB = MIN(M, J + KL)
         DO I = 1, M
            IF ((LB.LE.I).AND.(I.LE.UB)) THEN
               A_B( I, J ) = A( KU + 1 + I - J, J )
            ELSE
               A_B( I, J) = ZERO
            END IF
         END DO
      END DO

      WRITE(*,*) DESC
      WRITE(*, 9999) ( J, J = 1, N)
      DO I = 1, M
         WRITE(*, 9998) I, ( A_B( I, J ), J = 1, N )
      END DO
*
 9998 FORMAT( I2, ' ', 11(:,1X,F10.4) )
 9999 FORMAT( '   ', 11(:,1X,I10) )
*
      RETURN
      END
