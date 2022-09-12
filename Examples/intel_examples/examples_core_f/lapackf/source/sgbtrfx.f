*     SGBTRF (F07BDE) Example Program Text
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
      INTEGER          MMAX, NMAX, KLMAX, KUMAX, LDAB
      PARAMETER        (MMAX=8,NMAX=8,KLMAX=8,KUMAX=8,
     +                 LDAB=2*KLMAX+KUMAX+1)
*     .. Local Scalars ..
      INTEGER          I, INFO, J, K, KL, KU, M, N
*     .. Local Arrays ..
      REAL             AB(LDAB,NMAX)
      INTEGER          IPIV(NMAX)
*     .. External Subroutines ..
      EXTERNAL         PRINT_BAND_MATRIX
      EXTERNAL         SGBTRF
*     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SGBTRF Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) M, N, KL, KU
      IF (M.LE.MMAX .AND. N.LE.NMAX .AND. KL.LE.KLMAX .AND. KU.LE.KUMAX)
     +    THEN
*
*        Read A from data file
*
         K = KL + KU + 1
         READ (NIN,*) ((AB(K+I-J,J),J=MAX(I-KL,1),MIN(I+KU,N)),I=1,
     +     M)
*
*        Factorize A
*
         CALL SGBTRF(M,N,KL,KU,AB,LDAB,IPIV,INFO)
*
*        Print details of factorization
*
         WRITE (NOUT,*)
*
         CALL PRINT_BAND_MATRIX( 'Details of factorization', M, N, KL,
     +                           KL+KU, AB, LDAB )
*
*        Print pivot indices
*
         WRITE (NOUT,*)
         WRITE (NOUT,*) 'IPIV'
         WRITE (NOUT,99999) (IPIV(I),I=1,MIN(M,N))
*
         IF (INFO.NE.0) WRITE (NOUT,*) 'The factor U is singular'
*
      END IF
      STOP
*
99999 FORMAT ((3X,7I11))
      END
*
*     End of SGBTRF Example
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