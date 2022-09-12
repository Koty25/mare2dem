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
!    sddcor  Example Program Text
!*******************************************************************************

      EXTERNAL   SDDCOR
      INTEGER  INCH, INCX, INCY
      INTEGER  NH,   NX,   NY
      INTEGER(KIND=4)  IY0,  ID
      REAL(KIND=4)     H(3), X(15), Y(24)
      INTEGER(KIND=4)  I
      REAL(KIND=8)     R10
      INTEGER(KIND=4)  R1(6)
      INTEGER(KIND=4)  R2(3)
      DATA       R1 /2, 5, 8, 11, 14, 5 /
      DATA       R2 / 5, 11, 5 /
!************* Initialize data *****/
      R10 = 2.0E-5
      INCH = 2
      INCX = 3
      INCY = 4
      NH   = 2
      NX   = 5
      NY   = 6
      DO I = 1, NH
         H(1+(I-1)*INCH) = I
      END DO
      PRINT *, '   '
      DO I = 1, NX
         X(1+(I-1)*INCX) = I
      END DO

!************* 1-st Sample **********/

      IY0 = -1
      ID  = 1

!************* Call sddcon ***********/

      CALL SDDCOR( H, INCH, X, INCX, Y, INCY, NH, NX, IY0, NY, ID )

!************* Printing results *****/

      PRINT *, ' 1-ST SAMPLE OF SDDCOR.'
      PRINT *, '----------------------'
      PRINT *, 'PARAMETERS:'
      PRINT 10, '    INCH = ',INCH,'    INCX = ',INCX,'    INCY = ',INCY
      PRINT 10, '    NH   = ',NH,  '    NX   = ',NX,  '    NY   = ',NY
      PRINT 10, '    IY0  = ',IY0, '    ID   = ',ID
      DO I = 1, NH
         PRINT 11, ' H( ',1+(I-1)*INCH,') = ',H(1+(I-1)*INCH)
      END DO
      PRINT *, ' '
      DO I = 1, NX
         PRINT 11, ' X( ',1+(I-1)*INCX,') = ',X(1+(I-1)*INCX)
      END DO
      PRINT *, ' '
      PRINT *, 'RESULTS:'
      PRINT *, '---------------------------'
      DO I = 1, NY
         PRINT 11, ' Y( ',1+(I-1)*INCY,') = ',Y(1+(I-1)*INCY)
      END DO
      PRINT *, ' '
       DO I = 1, 6
         IF(ABS(Y(1+(I-1)*INCY)-R1(I)) .GT. R10) THEN
         PRINT *, 'ERROR: wrong result: Y(',1+(I-1)*INCY,')'
         PRINT *, '---------------------------'
         PRINT 10, ' TEST FAILED '
         PRINT *, '---------------------------'
         STOP 1
         ENDIF
      END DO
!**************** 2-nd Sample **********/

      IY0 = 0
      ID  = 2
      NY  = 3

!************* Call sddcon *****/

      CALL SDDCOR( H, INCH, X, INCX, Y, INCY, NH, NX, IY0, NY, ID )

!************* Printing results *****/

      PRINT *, ' 2-ND SAMPLE OF SDDCOR.'
      PRINT *, '----------------------'
      PRINT *, 'PARAMETERS:'
      PRINT 10, '    INCH = ',INCH,'    INCX = ',INCX,'    INCY = ',INCY
      PRINT 10, '    NH   = ',NH,  '    NX   = ',NX,  '    NY   = ',NY
      PRINT 10, '    IY0  = ',IY0, '    ID   = ',ID
      DO I = 1, NH
         PRINT 11, ' H( ',1+(I-1)*INCH,') = ',H(1+(I-1)*INCH)
      END DO
      PRINT *, ' '
      DO I = 1, NX
         PRINT 11, ' X( ',1+(I-1)*INCX,') = ',X(1+(I-1)*INCX)
      END DO
      PRINT *, ' '
      PRINT *, 'RESULTS:'
      PRINT *, '---------------------------'
      DO I = 1, NY
         PRINT 11, ' Y( ',1+(I-1)*INCY,') = ',Y(1+(I-1)*INCY)
      END DO
       DO I = 1, 3
         IF(ABS(Y(1+(I-1)*INCY)-R2(I)) .GT. R10) THEN
         PRINT *, 'ERROR: wrong result: Y(',1+(I-1)*INCY,')'
         PRINT *, '---------------------------'
         PRINT 10, ' TEST FAILED '
         PRINT *, '---------------------------'
         STOP 1
         ENDIF
      END DO
      PRINT *, ' '
      PRINT *, '---------------------------'
      PRINT 10, ' TEST PASSED '
      PRINT *, '---------------------------'
10    FORMAT(A,I4,A,I4,A,I4)
11    FORMAT(A,I2,A,F5.2)

      END
