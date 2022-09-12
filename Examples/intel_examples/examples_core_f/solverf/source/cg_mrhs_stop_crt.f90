!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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

!  Content: Intel(R) Math Kernel Library (Intel(R) MKL) RCI (P)CG Fortran example
!
!*******************************************************************************
!---------------------------------------------------------------------------
!  Example program for solving symmetric positive definite system of equations.
!  Simplest case: no preconditioning and the user-defined stopping tests.
!---------------------------------------------------------------------------
      PROGRAM rci_pcgmrhs_f_stop_crt
        USE MKL_SPBLAS
        IMPLICIT NONE

        INCLUDE 'mkl_rci.fi'
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix and rhs vector
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER N, nRhs
        PARAMETER (N=8,nRhs=2)
        INTEGER RCI_request, itercount(nRhs), i
        DOUBLE PRECISION  rhs(N,nRhs)

! Matrix
        INTEGER IA(9)
        INTEGER JA(18)
        DOUBLE PRECISION A(18)
! Fill all arrays containing matrix data.
        DATA IA /1,5,8,10,12,15,17,18,19/
        DATA JA &
&       /1,  3,    6,7,   &
&          2,3,  5,       &
&            3,        8, &
&              4,    7,   &
&                5,6,7,   &
&                  6,  8, &
&                    7,   &
&                      8/
        DATA A &
&       /7.D0,       1.D0,             2.D0, 7.D0,       &
&             -4.D0, 8.D0,       2.D0,                  &
&                    1.D0,                         5.D0,&
&                          7.D0,             9.D0,      &
&                                5.D0, 1.D0, 5.D0,      &
&                                     -1.D0,       5.D0,&
&                                           11.D0,      &
&                                                  5.D0/

!---------------------------------------------------------------------------
! Allocate storage for the solver ?par and the initial solution vector
!---------------------------------------------------------------------------
        INTEGER length, method, info
        PARAMETER (length=128, method=1)
        DOUBLE PRECISION expected_sol(N,nRhs)
        DOUBLE PRECISION solution(N,nRhs)

        INTEGER ipar(length + 2*nRhs)
        DOUBLE PRECISION dpar(length + 2*nRhs),TMP(N,3 + nRhs)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)CG solver
!---------------------------------------------------------------------------
        DATA expected_sol/1.D0, 0.D0, 1.D0, 0.D0, 1.D0, 0.D0, 1.D0, &
     &  0.D0, 0.D0, 2.D0, 0.D0, 2.D0, 0.D0, 2.D0, 0.D0, 2.D0/

        DOUBLE PRECISION DNRM2, Euclidean_norm, temp(N)
        EXTERNAL DNRM2
        DOUBLE PRECISION alpha, beta
!   Matrix descriptor
        TYPE(MATRIX_DESCR) descrA
!   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) csrA
        !   Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
        descrA % MODE = SPARSE_FILL_MODE_UPPER
        descrA % DIAG = SPARSE_DIAG_NON_UNIT
        alpha = 1.0
        beta  = 0.0
        info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,N,N,IA,IA(2),JA,A)
        DO I = 1, N
            rhs(I,1)     = 0.D0
            rhs(I,2)     = 0.D0
            temp(I)      = 0.D0
        END DO
!---------------------------------------------------------------------------
! Initialize the right hand side through matrix-vector product
!---------------------------------------------------------------------------
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,expected_sol,beta,rhs)
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,expected_sol(1,2),beta,rhs(1,2))
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
        DO I = 1, length + 2 * nRhs
            ipar(i) = 0;
            dpar(i) = 0;
        END DO
!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
        DO I = 1, N
            solution(I,1) = 1.D0
            solution(I,2) = 2.D0
        END DO

        CALL dcgmrhs_init(N,solution,nRhs,rhs,method,RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0 ) GO TO 999
!---------------------------------------------------------------------------
! Set the desired parameters:
! LOGICAL parameters:
! do residual stopping test
! request for the user defined stopping test
!---------------------------------------------------------------------------
        ipar(5)  = 100
        ipar(10) = 1
        dpar(1) = 1.D-5
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)CG solver without preconditioning
! Reverse Communications starts here
!---------------------------------------------------------------------------
1       CALL dcgmrhs(N,solution,nRhs,rhs,RCI_request,ipar,dpar,TMP)
!---------------------------------------------------------------------------
! If RCI_request=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
        IF (RCI_request .EQ. 0) THEN
            GO TO 700
!---------------------------------------------------------------------------
! If RCI_request=1, then compute the vector A*TMP(:,1)
! and put the result in vector TMP(:,2)
!---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 1) THEN
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,TMP,beta,TMP(1,2))
            GO TO 1
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test: compute the
! Euclidean norm of the actual residual using Intel MKL routines and check if
! it is less than 1.D-8
!---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 2) THEN
            info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,solution(1,ipar(3)),beta,temp)
            CALL DAXPY(N,-1.D0,rhs(1,ipar(3)),1,temp,1)
            Euclidean_norm = DNRM2(N,temp,1)
            IF (Euclidean_norm .GT. 1.D-4) THEN
!---------------------------------------------------------------------------
! The solution has not been found yet according to the user-defined stopping
! test. Continue RCI (P)CG iterations.
!---------------------------------------------------------------------------
                GO TO 1
            ELSE
                GO TO 700
            END IF
        ELSE
            GO TO 999
        END IF

!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number
!---------------------------------------------------------------------------
700     CALL dcgmrhs_get(N,solution,nRhs,rhs,RCI_request,ipar,dpar,TMP,itercount)
!---------------------------------------------------------------------------
! Print solution vector: solution(N) and number of iterations: itercount
!---------------------------------------------------------------------------
        WRITE(*, *) ' The system has been solved '
        WRITE(*, *) ' The following solution obtained '
        WRITE(*,800) (solution(i,1), i=1,N)
        WRITE(*,800) (solution(i,2), i=1,N)
        WRITE(*, *) ' expected solution '
        WRITE(*,800)(expected_sol(i,1), i=1,N)
        WRITE(*,800)(expected_sol(i,2), i=1,N)
800     FORMAT(4(F10.3))
        DO i = 1, N
            expected_sol(i,1) = expected_sol(i,1) - solution(i,1)
            expected_sol(i,2) = expected_sol(i,2) - solution(i,2)
        END DO

        Euclidean_norm=DNRM2(2*N,expected_sol,1)

!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
        CALL MKL_FREE_BUFFERS

        IF (Euclidean_norm .LE. 1.0D-12) THEN
            WRITE( *,'(A,A)') 'This example has successfully PASSED', &
     &  ' through all steps of computation!'
            STOP 0
        ELSE
            WRITE( *,'(A,A,A,E12.5,A)') 'This example may have', &
     &  ' FAILED as the computed solution differs much from the', &
     &  ' expected solution (Euclidean norm is',Euclidean_norm,').'
            STOP 1
        END IF
!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
999     WRITE( *,'(A,A)') 'This example FAILED as the solver has',&
     &  ' returned the ERROR code', RCI_request
        info = MKL_SPARSE_DESTROY(csrA)
        CALL MKL_FREE_BUFFERS
        STOP 1

      END PROGRAM rci_pcgmrhs_f_stop_crt
