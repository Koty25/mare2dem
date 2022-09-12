!===============================================================================
! Copyright 2003-2020 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!    sconf  Example Program Text
!*******************************************************************************

      EXTERNAL   SCONF
      INTEGER  INCH, INCX, INCY
      INTEGER  NH,   NX,   NY
      INTEGER(KIND=4)  IY0
      REAL(KIND=4)     H(3), X(15), Y(24)
      REAL(KIND=4)     AUX1(32), AUX2(32)
      INTEGER(KIND=4)  INIT, M,   NAUX1, NAUX2
      INTEGER(KIND=4)  INC2X, INC2Y
      INTEGER(KIND=4)  I
      REAL(KIND=8)     R10
      INTEGER(KIND=4)  R1(6)
      INTEGER(KIND=4)  R2(3)
      DATA       R1 /1, 4, 7, 10, 13, 10 /
!************* Initialize data *****/
      R10 = 2.0E-5
      INCH = 2
      INCX = 3
      INCY = 4
      NH   = 2
      NX   = 5
      NY   = 6
      NAUX1 = 32
      NAUX2 = 32
      DO I = 1, NH
         H(1+(I-1)*INCH) = I
      END DO
      PRINT *, '  '
       DO I = 1, NX
         X(1+(I-1)*INCX) = I
      END DO
      M     = 1
      INC2X = 1
      INC2Y = 1
      IY0   = 0
!************* 1-st step of Sample( preliminary ) **********/

      INIT = 1

!*********** Call sconf *****/

      CALL SCONF(INIT, H, INCH, X, INCX, INC2X, Y, INCY, INC2Y, NH, NX, &
     &  M, IY0, NY, AUX1, NAUX1, AUX2, NAUX2 )

!************ 2-nd step of Sample( main ) **********/

      INIT = 0

!************ Call sconf *****/

      CALL SCONF( INIT, H, INCH, X, INCX, INC2X, Y, INCY, INC2Y, NH, NX,&
     &  M, IY0, NY, AUX1, NAUX1, AUX2, NAUX2 )

!************* Printing results *****/

      PRINT *,' SAMPLE OF SCONF.'
      PRINT *,'----------------------'
      PRINT *,' PARAMETERS:  '
      PRINT 10, '    INCH = ',INCH,'   INCX = ',INCX,'   INCY = ',INCY
      PRINT 10, '    NH   = ',NH,'    NX   = ',NX,'   NY   = ',NY
      PRINT 10, '    IY0  = ',IY0
      DO I = 1, NH
        PRINT 11, '  H( ',1+(I-1)*INCH,') = ',H(1+(I-1)*INCH)
      END DO
        PRINT *,'  '
      DO I = 1, NX
        PRINT 11, ' X( ',1+(I-1)*INCX,') = ',X(1+(I-1)*INCX)
      END DO
      PRINT *, ' '
      PRINT *, 'RESULTS:'
      PRINT *, '---------------------------'
      DO I = 1, NY
         PRINT 11, ' Y( ',1+(I-1)*INCY,') = ',Y(1+(I-1)*INCY)
      END DO
      DO I = 1, 6
         IF(ABS(Y(1+(I-1)*INCY)-R1(I)) .GT. R10) THEN
         PRINT *, 'ERROR: wrong result: Y(',1+(I-1)*INCY,')'
         PRINT *, '---------------------------'
         PRINT 10, ' TEST FAILED '
         PRINT *, '---------------------------'
         STOP 1
         ENDIF
      END DO
      PRINT *, ' '
      PRINT *, '---------------------------'
10    FORMAT(A,I4,A,I4,A,I4)
11    FORMAT(A,I2,A,F5.2)
      PRINT 10, ' TEST PASSED '
      PRINT *, '---------------------------'
      STOP
      END
