!===============================================================================
! Copyright 2004-2020 Intel Corporation.
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

*   Content : TR Solver Fortran-77 example
*
********************************************************************************

!** NONLINEAR LEAST SQUARE PROBLEM WITH BOUNDARY CONSTRAINTS
      PROGRAM EXAMPLE_DTRNSPBC_POWELL
        IMPLICIT NONE
!** HEADER-FILE WITH DEFINITIONS (CONSTANTS, EXTERNALS)
        INCLUDE "mkl_rci.fi"
!** USER�S OBJECTIVE FUNCTION
        EXTERNAL EXTENDED_POWELL
!** N - NUMBER OF FUNCTION VARIABLES
        INTEGER             N
        PARAMETER           (N = 4)
!** M - DIMENSION OF FUNCTION VALUE
        INTEGER             M
        PARAMETER           (M = 4)
!** SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
        REAL*4              X (N)
!** PRECISIONS FOR STOP-CRITERIA (SEE MANUAL FOR MORE DETAILS)
        REAL*4              EPS (6)
!** JACOBI CALCULATION PRECISION
        REAL*4              JAC_EPS
!** LOWER AND UPPER BOUNDS
        REAL*4              LW (N), UP (N)
!** REVERSE COMMUNICATION INTERFACE PARAMETER
        INTEGER             RCI_REQUEST
!** FUNCTION (F(X)) VALUE VECTOR
        REAL*4              FVEC (M)
!** JACOBI MATRIX
        REAL*4              FJAC (M, N)
!** NUMBER OF ITERATIONS
        INTEGER             ITER
!** NUMBER OF STOP-CRITERION
        INTEGER             ST_CR
!** CONTROLS OF RCI CYCLE
        INTEGER             SUCCESSFUL
!** MAXIMUM NUMBER OF ITERATIONS
        INTEGER             ITER1
!** MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
        INTEGER             ITER2
!** INITIAL STEP BOUND
        REAL*4              RS
!** INITIAL AND FINAL RESIDUALS
        REAL*4              R1, R2
!** TR SOLVER HANDLE
        INTEGER*8           HANDLE
!** CYCLE�S COUNTERS
        INTEGER             I, J
!** RESULTS OF INPUT PARAMETER CHECKING
        INTEGER             INFO(6)
!** SET PRECISIONS FOR STOP-CRITERIA
        DO I = 1, 6
            EPS (I) = 1.D-3
        END DO
!** SET MAXIMUM NUMBER OF ITERATIONS
        ITER1 = 1000
!** SET MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
        ITER2 = 100
!** SET INITIAL STEP BOUND
        RS = 100.D0
!** PRECISIONS FOR JACOBI CALCULATION
        JAC_EPS = 1.D-6
!** SET THE INITIAL GUESS
        DO I = 1, N/4
            X (4*I - 3) =  3.D0
            X (4*I - 2) = -1.D0
            X (4*I - 1) =  0.D0
            X (4*I)     =  1.D0
        END DO
!** SET LOWER AND UPPER BOUNDS
        DO I = 1, N/4
            LW(4*I-3) =  0.1D0
            LW(4*I-2) = -20.D0
            LW(4*I-1) =  -1.D0
            LW(4*I)   =  -1.D0

            UP(4*I-3) = 100.D0
            UP(4*I-2) =  20.D0
            UP(4*I-1) =   1.D0
            UP(4*I)   =  50.D0
        END DO
!** SET INITIAL VALUES
        DO I = 1, M
            FVEC (I) = 0.D0
            DO J = 1, N
                FJAC (I, J) = 0.D0
            END DO
        END DO
!** INITIALIZE SOLVER (ALLOCATE MEMORY, SET INITIAL VALUES)
!**   HANDLE    IN/OUT: TR SOLVER HANDLE
!**   N         IN:     NUMBER OF FUNCTION VARIABLES
!**   M         IN:     DIMENSION OF FUNCTION VALUE
!**   X         IN:     SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
!**   LW        IN:     LOWER BOUND
!**   UP        IN:     UPPER BOUND
!**   EPS       IN:     PRECISIONS FOR STOP-CRITERIA
!**   ITER1     IN:     MAXIMUM NUMBER OF ITERATIONS
!**   ITER2     IN:     MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
!**   RS        IN:     INITIAL STEP BOUND
        IF (STRNLSPBC_INIT (HANDLE, N, M, X, LW, UP, EPS, ITER1, ITER2
     &          , RS) .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            PRINT *, '| ERROR IN DTRNLSPBC_INIT'
!** RELEASE INTERNAL Intel(R) Math Kernel Library (Intel(R) MKL) MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
!** AND STOP
            STOP 1
        END IF
!** CHECKS THE CORRECTNESS OF HANDLE AND ARRAYS CONTAINING JACOBIAN MATRIX, 
!** OBJECTIVE FUNCTION, LOWER AND UPPER BOUNDS, AND STOPPING CRITERIA.
        IF (STRNLSPBC_CHECK (HANDLE, N, M, FJAC, FVEC, LW, UP, EPS, 
     &         INFO) .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            PRINT *, '| ERROR IN DTRNLSPBC_INIT'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
!** AND STOP
            STOP 1
        ELSE
!** THE HANDLE IS NOT VALID.
            IF( INFO(1) .NE. 0 .OR. 
!** THE FJAC ARRAY IS NOT VALID.
     &              INFO(2) .NE. 0 .OR. 
!** THE FVEC ARRAY IS NOT VALID.
     &              INFO(3) .NE. 0 .OR. 
!** THE LW ARRAY IS NOT VALID.
     &              INFO(4) .NE. 0 .OR. 
!** THE UP ARRAY IS NOT VALID.
     &              INFO(5) .NE. 0 .OR. 
!** THE EPS ARRAY IS NOT VALID.
     &              INFO(6) .NE. 0 ) THEN
                PRINT *, '| INPUT PARAMETERS ARE NOT VALID'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
                CALL MKL_FREE_BUFFERS
!** AND STOP
                STOP 1
            END IF
        END IF
!** SET INITIAL RCI CYCLE VARIABLES
        RCI_REQUEST = 0
        SUCCESSFUL = 0
!** RCI CYCLE
        DO WHILE (SUCCESSFUL == 0)
!** CALL TR SOLVER
!**   HANDLE        IN/OUT: TR SOLVER HANDLE
!**   FVEC          IN:     VECTOR
!**   FJAC          IN:     JACOBI MATRIX
!**   RCI_REQUEST   IN/OUT: RETURN NUMBER WHICH DENOTE NEXT STEP FOR PERFORMING
            IF (STRNLSPBC_SOLVE (HANDLE, FVEC, FJAC, RCI_REQUEST)
     &              .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
                PRINT *, '| ERROR IN DTRNLSPBC_SOLVE'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
                CALL MKL_FREE_BUFFERS
!** AND STOP
                STOP 1
            END IF
!** ACCORDING WITH RCI_REQUEST VALUE WE DO NEXT STEP
            SELECT CASE (RCI_REQUEST)
            CASE (-1, -2, -3, -4, -5, -6)
!**   STOP RCI CYCLE
                SUCCESSFUL = 1
            CASE (1)
!**   RECALCULATE FUNCTION VALUE
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     X               IN:     SOLUTION VECTOR
!**     FVEC            OUT:    FUNCTION VALUE F(X)
                CALL EXTENDED_POWELL (M, N, X, FVEC)
            CASE (2)
!**   COMPUTE JACOBI MATRIX
!**     EXTENDED_POWELL IN:     EXTERNAL OBJECTIVE FUNCTION
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     FJAC            OUT:    JACOBI MATRIX
!**     X               IN:     SOLUTION VECTOR
!**     JAC_EPS         IN:     JACOBI CALCULATION PRECISION
                IF (SJACOBI (EXTENDED_POWELL, N, M, FJAC, X, JAC_EPS)
     &                  .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
                    PRINT *, '| ERROR IN DJACOBI'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
                    CALL MKL_FREE_BUFFERS
!** AND STOP
                    STOP 1
                END IF
            ENDSELECT
        END DO
!** GET SOLUTION STATUSES
!**   HANDLE            IN:
!**   ITER              OUT: NUMBER OF ITERATIONS
!**   ST_CR             OUT: NUMBER OF STOP CRITERION
!**   R1                OUT: INITIAL RESIDUALS
!**   R2                OUT: FINAL RESIDUALS
        IF (STRNLSPBC_GET (HANDLE, ITER, ST_CR, R1, R2)
     &          .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            PRINT *, '| ERROR IN DTRNLSPBC_GET'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
!** AND STOP
            STOP 1
        END IF
!** FREE HANDLE MEMORY
        IF (STRNLSPBC_DELETE (HANDLE) .NE. TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            PRINT *, '| ERROR IN DTRNLSPBC_DELETE'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            CALL MKL_FREE_BUFFERS
!** AND STOP
            STOP 1
        END IF
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
        CALL MKL_FREE_BUFFERS
!** IF FINAL RESIDUAL LESS THEN REQUIRED PRECISION THEN PRINT PASS
        IF (R2 .LT. 0.1D0) THEN
            PRINT *, '|         DTRNLSPBC POWELL..........PASS'
            STOP 0
!** ELSE PRINT FAILED
        ELSE
            PRINT *, '|         DTRNLSPBC POWELL..........FAILED'
            STOP 1
        END IF
      END PROGRAM EXAMPLE_DTRNSPBC_POWELL

!** ROUTINE FOR EXTENDED POWELL FUNCTION CALCULATION
!**   M     IN:     DIMENSION OF FUNCTION VALUE
!**   N     IN:     NUMBER OF FUNCTION VARIABLES
!**   X     IN:     VECTOR FOR FUNCTION CALCULATING
!**   F     OUT:    FUNCTION VALUE F(X)
      SUBROUTINE EXTENDED_POWELL (M, N, X, F)
        IMPLICIT NONE
        INTEGER M, N
        REAL*4 X (*), F (*)
        INTEGER I

        DO I = 1, N/4
            F (4*I-3) = X(4*I - 3) + 10.D0 * X(4*I - 2)
            F (4*I-2) = 2.2D0 * (X(4*I-1) - X(4*I))
            F (4*I-1) = ( X(4*I-2) - 2.D0*X(4*I-1) )**2
            F (4*I)   = 3.1D0 * (X(4*I-3) - X(4*I))**2
        END DO
      END SUBROUTINE EXTENDED_POWELL
