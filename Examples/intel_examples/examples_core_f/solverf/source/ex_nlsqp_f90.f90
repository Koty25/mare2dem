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

!   Content : TR Solver Fortran-77 example
!
!*******************************************************************************

!** NONLINEAR LEAST SQUARE PROBLEM WITHOUT BOUNDARY CONSTRAINTS
    include 'mkl_rci.f90'
program EXAMPLE_DTRNLSP_POWELL
    use MKL_RCI
    use MKL_RCI_TYPE
    implicit none
!** HEADER-FILE WITH DEFINITIONS (CONSTANTS, EXTERNALS)

!** USER’S OBJECTIVE FUNCTION
    external            EXTENDED_POWELL
!** N - NUMBER OF FUNCTION VARIABLES
    integer             N
    parameter           (N = 40)
!** M - DIMENSION OF FUNCTION VALUE
    integer             M
    parameter           (M = 40)
!** SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
    double precision    X (N)
!** PRECISIONS FOR STOP-CRITERIA (SEE MANUAL FOR MORE DETAILS)
    double precision    EPS (6)
!** JACOBI CALCULATION PRECISION
    double precision    JAC_EPS
!** REVERSE COMMUNICATION INTERFACE PARAMETER
    integer             RCI_REQUEST
!** FUNCTION (F(X)) VALUE VECTOR
    double precision    FVEC (M)
!** JACOBI MATRIX
    double precision    FJAC (M, N)
!** NUMBER OF ITERATIONS
    integer             ITER
!** NUMBER OF STOP-CRITERION
    integer             ST_CR
!** CONTROLS OF RCI CYCLE
    integer             SUCCESSFUL
!** MAXIMUM NUMBER OF ITERATIONS
    integer             ITER1
!** MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
    integer             ITER2
!** INITIAL STEP BOUND
    double precision    RS
!** INITIAL AND FINAL RESIDUALS
    double precision    R1, R2
!** TR SOLVER HANDLE
    TYPE(HANDLE_TR) :: HANDLE
!** CYCLES COUNTERS
    integer             I, J
!** RESULTS OF INPUT PARAMETER CHECKING
    integer INFO(6)
!** SET PRECISIONS FOR STOP-CRITERIA
    do I = 1, 6
        EPS (I) = 1.D-5
    end do
!** SET MAXIMUM NUMBER OF ITERATIONS
    ITER1 = 1000
!** SET MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
    ITER2 = 100
!** SET INITIAL STEP BOUND
    RS = 100.D0
!** PRECISIONS FOR JACOBI CALCULATION
    JAC_EPS = 1.D-8
!** SET THE INITIAL GUESS
    do I = 1, N/4
        X (4*I - 3) =  3.D0
        X (4*I - 2) = -1.D0
        X (4*I - 1) =  0.D0
        X (4*I)     =  1.D0
    end do
!** SET INITIAL VALUES
    do I = 1, M
        FVEC (I) = 0.D0
        do J = 1, N
            FJAC (I, J) = 0.D0
        end do
    end do
!** INITIALIZE SOLVER (ALLOCATE MEMORY, SET INITIAL VALUES)
!**     HANDLE    IN/OUT: TR SOLVER HANDLE
!**     N         IN:     NUMBER OF FUNCTION VARIABLES
!**     M         IN:     DIMENSION OF FUNCTION VALUE
!**     X         IN:     SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
!**     EPS       IN:     PRECISIONS FOR STOP-CRITERIA
!**     ITER1     IN:     MAXIMUM NUMBER OF ITERATIONS
!**     ITER2     IN:     MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
!**     RS        IN:     INITIAL STEP BOUND
    if (DTRNLSP_INIT (HANDLE, N, M, X, EPS, ITER1, ITER2, RS) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_INIT'
!** RELEASE INTERNAL Intel(R) Math Kernel Library (Intel(R) MKL) MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if
!** CHECKS THE CORRECTNESS OF HANDLE AND ARRAYS CONTAINING JACOBIAN MATRIX, 
!** OBJECTIVE FUNCTION, LOWER AND UPPER BOUNDS, AND STOPPING CRITERIA.
    if (DTRNLSP_CHECK (HANDLE, N, M, FJAC, FVEC, EPS, INFO) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSPBC_INIT'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    else
!** THE HANDLE IS NOT VALID.
        if ( INFO(1) /= 0 .or. INFO(2) /= 0 .or. INFO(3) /= 0 .or. INFO(4) /= 0 ) THEN
!** THE FJAC ARRAY IS NOT VALID.
!     +     INFO(2) /= 0 .or. 
!** THE FVEC ARRAY IS NOT VALID.
!     +     INFO(3) /= 0 .or. 
!** THE EPS ARRAY IS NOT VALID.
!     +     INFO(4) /= 0 ) THEN
            print *, '| INPUT PARAMETERS ARE NOT VALID'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            call MKL_FREE_BUFFERS
!** AND STOP
            stop 1
        end if
    end if
!** SET INITIAL RCI CYCLE VARIABLES
    RCI_REQUEST = 0
    SUCCESSFUL = 0
!** RCI CYCLE
    do while (SUCCESSFUL == 0)
!** CALL TR SOLVER
!**   HANDLE        IN/OUT: TR SOLVER HANDLE
!**   FVEC          IN:     VECTOR
!**   FJAC          IN:     JACOBI MATRIX
!**   RCI_REQUEST   IN/OUT: RETURN NUMBER WHICH DENOTE NEXT STEP FOR PERFORMING
        if (DTRNLSP_SOLVE (HANDLE, FVEC, FJAC, RCI_REQUEST) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
            print *, '| ERROR IN DTRNLSP_SOLVE'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
            call MKL_FREE_BUFFERS
!** AND STOP
            stop 1
        end if
!** ACCORDING WITH RCI_REQUEST VALUE WE DO NEXT STEP
        select case (RCI_REQUEST)
        case (-1, -2, -3, -4, -5, -6)
!**   STOP RCI CYCLE
            SUCCESSFUL = 1
        case (1)
!**   RECALCULATE FUNCTION VALUE
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     X               IN:     SOLUTION VECTOR
!**     FVEC            OUT:    FUNCTION VALUE F(X)
            call EXTENDED_POWELL (M, N, X, FVEC)
        case (2)
!**   COMPUTE JACOBI MATRIX
!**     EXTENDED_POWELL IN:     EXTERNAL OBJECTIVE FUNCTION
!**     N               IN:     NUMBER OF FUNCTION VARIABLES
!**     M               IN:     DIMENSION OF FUNCTION VALUE
!**     FJAC            OUT:    JACOBI MATRIX
!**     X               IN:     SOLUTION VECTOR
!**     JAC_EPS         IN:     JACOBI CALCULATION PRECISION
            if (DJACOBI (EXTENDED_POWELL, N, M, FJAC, X, JAC_EPS) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
                print *, '| ERROR IN DJACOBI'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
                call MKL_FREE_BUFFERS
!** AND STOP
                stop 1
            end if
        end select
    end do
!** GET SOLUTION STATUSES
!**   HANDLE            IN: TR SOLVER HANDLE
!**   ITER              OUT: NUMBER OF ITERATIONS
!**   ST_CR             OUT: NUMBER OF STOP CRITERION
!**   R1                OUT: INITIAL RESIDUALS
!**   R2                OUT: FINAL RESIDUALS
    if (DTRNLSP_GET (HANDLE, ITER, ST_CR, R1, R2) /= TR_SUCCESS) THEN
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_GET'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if
!** FREE HANDLE MEMORY
    if (DTRNLSP_DELETE (HANDLE) /= TR_SUCCESS) then
!** IF FUNCTION DOES NOT COMPLETE SUCCESSFULLY THEN PRINT ERROR MESSAGE
        print *, '| ERROR IN DTRNLSP_DELETE'
!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
        call MKL_FREE_BUFFERS
!** AND STOP
        stop 1
    end if

!** RELEASE INTERNAL Intel MKL MEMORY THAT MIGHT BE USED FOR COMPUTATIONS.
!** NOTE: IT IS IMPORTANT TO CALL THE ROUTINE BELOW TO AVOID MEMORY LEAKS
!** UNLESS YOU DISABLE Intel MKL MEMORY MANAGER
    call MKL_FREE_BUFFERS
!** IF FINAL RESIDUAL LESS THEN REQUIRED PRECISION THEN PRINT PASS
    if (R2 < 1.D-5) then
        print *, '|         DTRNLSP POWELL............PASS'
        stop 0
!** ELSE PRINT FAILED
    else
        print *, '|         DTRNLSP POWELL............FAILED'
        stop 1
    end if
end program EXAMPLE_DTRNLSP_POWELL

!** ROUTINE FOR EXTENDED POWELL FUNCTION CALCULATION
!**   M     IN:     DIMENSION OF FUNCTION VALUE
!**   N     IN:     NUMBER OF FUNCTION VARIABLES
!**   X     IN:     VECTOR FOR FUNCTION CALCULATING
!**   F     OUT:    FUNCTION VALUE F(X)
subroutine EXTENDED_POWELL (M, N, X, F)
    implicit none
    integer M, N
    double precision X (*), F (*)
    integer I

    do I = 1, N/4
        F (4*I-3) = X(4*I - 3) + 10.D0 * X(4*I - 2)
        F (4*I-2) = 2.2360679774998D0 * (X(4*I-1) - X(4*I))
        F (4*I-1) = ( X(4*I-2) - 2.D0*X(4*I-1) )**2
        F (4*I)   = 3.1622776601684D0 * (X(4*I-3) - X(4*I))**2
    end do
end subroutine EXTENDED_POWELL
