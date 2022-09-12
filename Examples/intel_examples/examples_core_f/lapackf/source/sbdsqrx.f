*     SBDSQR (F08MEE) Example Program Text
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
      INTEGER          NMAX, LDVT, LDU, LDC
      PARAMETER        (NMAX=8,LDVT=NMAX,LDU=NMAX,LDC=1)
      REAL             ZERO, ONE
      PARAMETER        (ZERO=0.0,ONE=1.0)
*     .. Local Scalars ..
      INTEGER          I, INFO, N
      CHARACTER        UPLO
*     .. Local Arrays ..
      REAL             C(LDC,1), D(NMAX), E(NMAX-1), U(LDU,NMAX),
     +                 VT(LDVT,NMAX), WORK(4*NMAX-4)
*     .. External Subroutines ..
      EXTERNAL         PRINT_MATRIX
      EXTERNAL         SLASET, SBDSQR
*     .. Executable Statements ..
      WRITE (NOUT,*) 'SBDSQR Example Program Results'
*     Skip heading in data file
      READ (NIN,*)
      READ (NIN,*) N
      IF (N.LE.NMAX) THEN
*
*        Read B from data file
*
         READ (NIN,*) (D(I),I=1,N)
         READ (NIN,*) (E(I),I=1,N-1)
*
         READ (NIN,*) UPLO
*
*        Initialise U and VT to be the unit matrix
*
         CALL SLASET('General',N,N,ZERO,ONE,U,LDU)
*
         CALL SLASET('General',N,N,ZERO,ONE,VT,LDVT)
*
*        Calculate the SVD of B
*
         CALL SBDSQR(UPLO,N,N,N,0,D,E,VT,LDVT,U,LDU,C,LDC,WORK,INFO)
*
         WRITE (NOUT,*)
         IF (INFO.GT.0) THEN
            WRITE (NOUT,*) 'Failure to converge.'
         ELSE
*
*           Print singular values, left & right singular vectors
*
            WRITE (NOUT,*) 'Singular values'
            WRITE (NOUT,99999) (D(I),I=1,N)
            WRITE (NOUT,*)
*
            CALL PRINT_MATRIX( 'Right singular vectors, by row', N,
     +                         N, VT, LDVT )

            WRITE (NOUT,*)
*
            CALL PRINT_MATRIX( 'Left singular vectors, by column', N,
     +                         N, U, LDU )
         END IF
      END IF
      STOP
*
99999 FORMAT (3X,(8F8.4))
      END
*
*     End of SBDSQR Example
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
