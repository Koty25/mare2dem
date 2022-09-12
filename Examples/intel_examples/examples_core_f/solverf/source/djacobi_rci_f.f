!===============================================================================
! Copyright 2009-2020 Intel Corporation.
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

!  Content: DJACOBI RCI Example
!
!  The program computes the Jacobi matrix of the function on the basis of RCI
!  using the central difference.
!*******************************************************************************

      PROGRAM JACOBI_MATRIX 
        IMPLICIT NONE 
C**
        INCLUDE 'mkl_rci.fi'
C**
        EXTERNAL EXTENDED_POWELL 
C**
C** N - Number of function variables 
C** M - Dimension of function value 
        INTEGER N, M, I 
        PARAMETER (N = 4)
        PARAMETER (M = 4)
C**
C** Jacobi matrix 
        DOUBLE PRECISION A (M,N)
C** Solution vector. contains values x for f(x) 
C** Temporary arrays f1 & f2 which contains f1 = f(x+eps) | f2 = f(x-eps) 
        DOUBLE PRECISION F1(M), F2(M), X(N)
C** Precisions for jacobi_matrix calculation 
        DOUBLE PRECISION EPS 
C**
C** Jacobi-matrix solver handle 
        INTEGER*8 HANDLE 
C** Controls of rci cycle 
        INTEGER SUCCESSFUL, RCI_REQUEST 
        INTEGER RESULT
C**
C** Set the x values 
C** X   = 10.D0
        DO I = 1, N
            X(I) = 10.0D0
        END DO
C**
        EPS = 1.D-6
        PRINT *, 'START TESTING ...'
C** Initialize solver (allocate memory, set initial values) 
        RESULT = DJACOBI_INIT (HANDLE, N, M, X, A, EPS)
        IF (RESULT .NE. TR_SUCCESS) THEN
C** If function does not complete successfully then print error message 
            PRINT *, '#FAIL: ERROR IN DJACOBI_INIT' 
            CALL MKL_FREE_BUFFERS
            STOP 1
        END IF
C** Set initial rci cycle variables
        RCI_REQUEST = 0
        SUCCESSFUL  = 0
C** Rci cycle
        DO WHILE (SUCCESSFUL.EQ.0)
C** Call solver
            IF (DJACOBI_SOLVE (HANDLE, F1, F2, RCI_REQUEST) .NE. 
     &          TR_SUCCESS) THEN
C** If function does not complete successfully then print error message
                PRINT *, '#FAIL: ERROR IN DJACOBI_SOLVE'
                CALL MKL_FREE_BUFFERS
                STOP 1
            END IF
            IF (RCI_REQUEST .EQ. 1) THEN
C** Calculate the function value f1 = f(x+eps) 
                CALL EXTENDED_POWELL (M, N, X, F1)
            ELSE IF (RCI_REQUEST .EQ. 2) THEN
C** Calculate the function value f2 = f(x-eps) 
                CALL EXTENDED_POWELL (M, N, X, F2) 
            ELSE IF (RCI_REQUEST .EQ. 0) THEN
C** Exit rci cycle 
                SUCCESSFUL = 1 
            END IF
        END DO
C** Free handle memory
        IF (DJACOBI_DELETE (HANDLE) .NE. TR_SUCCESS) THEN
C** If function does not complete successfully then print error message 
            PRINT *, '#FAIL: ERROR IN DJACOBI_DELETE' 
            CALL MKL_FREE_BUFFERS
            STOP 1
        END IF
        PRINT *, '#PASS'
      END PROGRAM JACOBI_MATRIX
C**
C** Routine for extended Powell function calculation 
C** M in: dimension of function value 
C** N in: number of function variables 
C** X in: vector for function calculation 
C** F out: function value f(x) 
C**
      SUBROUTINE EXTENDED_POWELL (M, N, X, F) 
        IMPLICIT NONE 
        INTEGER M, N 
        DOUBLE PRECISION X (*), F (*) 
        INTEGER I 
        DO I = 1, N/4 
            F (4*I-3) = X(4*I - 3) + 10.D0 * X(4*I - 2) 
            F (4*I-2) = DSQRT(5.D0) * (X(4*I-1) - X(4*I)) 
            F (4*I-1) = ( X(4*I-2) - 2.D0*X(4*I-1) )**2 
            F (4*I)   = DSQRT(10.D0) * (X(4*I-3) - X(4*I))**2 
        END DO 
      END SUBROUTINE EXTENDED_POWELL
