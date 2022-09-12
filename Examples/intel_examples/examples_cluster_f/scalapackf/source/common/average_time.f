!===============================================================================
! Copyright 2015-2020 Intel Corporation.
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

      SUBROUTINE CALC_AVERAGE_TIME( TIME_ITER, ITERATIONS, TIME )
*     ==== Subroutine arguments ========================================
      INTEGER          ITERATIONS
      DOUBLE PRECISION TIME_ITER( ITERATIONS )
      DOUBLE PRECISION TIME
*
*     ==== External functions ========================================== 
      EXTERNAL SORT_ARRAY
*      
*     This subroutine computes average value of time array, excluding 
*     quarter of best times and quarter of worst times
*
*     ==== Local scalars ===============================================
      DOUBLE PRECISION FULLTIME
      INTEGER          I, LOWBND
*
*     ==== Executable statements =======================================
*
*     Sort array by increase index
      CALL SORT_ARRAY( TIME_ITER, ITERATIONS )
*
*     Calculate average time, excluding upper and lower quarters of array
      LOWBND = ITERATIONS / 4
      FULLTIME = 0.0D+0
      DO 10 I = 1 + LOWBND, ITERATIONS - LOWBND
         FULLTIME = FULLTIME + TIME_ITER( I )
   10 CONTINUE
      TIME = FULLTIME / DBLE( ITERATIONS - 2 * LOWBND )
*
*     ==== End of CALC_AVERAGE_TIME ====================================
      END
*
*=======================================================================
*
      SUBROUTINE SORT_ARRAY( ARRAY, DIM )
*     ==== Subroutine parameters =======================================
      INTEGER DIM
      DOUBLE PRECISION ARRAY( DIM )
*
*     This is subroutine for array sorting
*
*     ==== Local scalars ===============================================
      INTEGER I, J
      DOUBLE PRECISION TMP
*
*     ==== Executable statements =======================================
      DO 20 I = 2, DIM
         TMP = ARRAY( I )
         J = I - 1
         DO WHILE ( J.GE.1 .AND. ARRAY( J ).GT.TMP )
            ARRAY( J + 1 ) = ARRAY( J )
            J = J - 1
         ENDDO
         ARRAY( J + 1 ) = TMP
*
   20 CONTINUE
*
*     ==== End of SORT_ARRAY ===========================================
      END
