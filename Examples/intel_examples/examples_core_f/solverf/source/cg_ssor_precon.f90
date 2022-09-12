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
! Example program for solving symmetric positive definite system of
! equations A*x =f where A is a symmetric tridiagonal matrix of order N whose
! nonzero elements are defined as follows
!
!         a(i,  i  ) = 1,     i=1, N;
!         a(i,  i+1) = 1/2,   i=1, N-1;
!         a(i+1, i ) = 1/2,   i=1, N-1;
!
! The system is solved with the help of the conjugate gradient method where
! the Symmetric Successive Over Relaxation (SSOR) stationary iterative solver
! with the given number of iterations is used as the preconditioner. The relaxation
! parameter is set to 1/2, the number of iterations for the SSOR preconditioner
! is equal to 20. For simplicity, we don't check the convergence while solving
! equations with the help of SSOR. Let us recall that the scheme of SSOR is
! the following:
!
!         (D+w*U^{T}) x_{k+1} = (-w*U+(1-w)*D)*x_{k}+w*f
!
! where w is the relaxation parameter, f is the right hand side, D is the diagonal
! of A, U its strict upper part so that
!
!                  A   =  U^{T} + D + U
!
! The compressed sparse row format is used for storing nonzeros of A and
! since D is the identity matrix, we don't need to store diagonal elements. So
! the code given below only uses the sparse representation of U. The Intel
! Intel MKL Sparse BLAS is designed so that it allows the user to perform large
! variety of operations with one sparse representation setting appropriate
! values of the descriptor array (see Intel MKL User Guide for the further details).
!
! Full case: full functionality of RCI (P)CG is used.
!---------------------------------------------------------------------------
      PROGRAM rci_pcg_f_test_ssor
        USE MKL_SPBLAS
        IMPLICIT NONE
        INCLUDE 'mkl_rci.fi'
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix and
! preconditioner as well as an array for rhs vector
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER N, RCI_request, itercount, expected_itercount, i, info
        PARAMETER (N=100)
        PARAMETER (expected_itercount=19)
        DOUBLE PRECISION  rhs(N)
        INTEGER IA(N+1)
        INTEGER JA(N-1)
        DOUBLE PRECISION A(N-1), A1(N-1)
!---------------------------------------------------------------------------
! Allocate storage for the solver ?par and the initial solution vector
!---------------------------------------------------------------------------
        INTEGER length
        PARAMETER (length=128)
        INTEGER ipar(length)
        DOUBLE PRECISION dpar(length),TMP(N,4)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)CG solver
! OMEGA is the relaxation parameter, NITER_SSOR is the maximum number of
! iterations for the SSOR preconditioner
!---------------------------------------------------------------------------
        DOUBLE PRECISION solution(N)
        DOUBLE PRECISION expected_sol(N)
        DOUBLE PRECISION OMEGA, ONE, ZERO
        DATA  OMEGA/0.5D0/, ONE/1.D0/, ZERO/0.D0/
        DOUBLE PRECISION DNRM2, Euclidean_norm, temp(N)
        INTEGER   NITER_SSOR
        DATA      NITER_SSOR/20/

        EXTERNAL DNRM2
        DOUBLE PRECISION alpha, beta
!   Matrix descriptor
        TYPE(MATRIX_DESCR) descrA
!   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) csrA, csrA1
        !   Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
        descrA % MODE = SPARSE_FILL_MODE_UPPER
        descrA % DIAG = SPARSE_DIAG_UNIT
        alpha = 1.0
        beta  = 0.0

!---------------------------------------------------------------------------
! Initialize the coefficient matrix and expected solution
!---------------------------------------------------------------------------
        DO I = 1, N
            expected_sol(I) = 1.D0
        END DO
        DO I = 1, N-1
            JA(I) = I + 1
            IA(I) = I
            A(I) = 0.5D0
            A1(I) = OMEGA * A(I)
        END DO
        IA(N) = N
        IA(N+1) = IA(N)
        info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,N,N,IA,IA(2),JA,A)
        info = MKL_SPARSE_D_CREATE_CSR(csrA1,SPARSE_INDEX_BASE_ONE,N,N,IA,IA(2),JA,A1)
!---------------------------------------------------------------------------
! Initialize vectors rhs, TEMP, and TMP(:,2) with zeros as MKL_DCSRMV
! routine does not set NAN to zero. Thus, if any of the values in the
! vectors above accidentally happens to be NAN, the example will fail
! to complete.
! Initialize the right hand side through matrix-vector product
!---------------------------------------------------------------------------
        DO I = 1, N
            rhs(I) = ZERO
            TEMP(I) = ZERO
            TMP(I,2) = ZERO
        END DO
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,expected_sol,beta,rhs)

!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
        DO I = 1, N
            solution(I) = ZERO
        END DO
!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
        CALL dcg_init(N,solution,rhs,RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0) GO TO 999
!---------------------------------------------------------------------------
! Set the desired parameters:
! INTEGER parameters:
! set the maximal number of iterations to 100
! LOGICAL parameters:
! run the Preconditioned version of RCI (P)CG with preconditioner C_inverse
! DOUBLE PRECISION parameters
! -
!---------------------------------------------------------------------------
        ipar(5) = 100
        ipar(11) = 1
!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
        CALL dcg_check(N,solution,rhs,RCI_request,ipar,dpar,TMP)
        IF (RCI_request .NE. 0 .AND. RCI_request .NE. -1001) GO TO 999
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)CG solver
! Reverse Communications starts here
!---------------------------------------------------------------------------
1       CALL dcg(N,solution,rhs,RCI_request,ipar,dpar,TMP)
!---------------------------------------------------------------------------
! If RCI_request=0, then the solution was found according to the requested
! stopping tests. In this case, this means that it was found after 100
! iterations.
!---------------------------------------------------------------------------
        IF (RCI_request .EQ. 0) THEN
            GO TO 700
!---------------------------------------------------------------------------
! If RCI_request=1, then compute the vector A*TMP(:,1)
! and put the result in vector TMP(:,2)
!---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 1) THEN
            descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,TMP,beta,TMP(1,2))
            GO TO 1
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test: compute the
! Euclidean norm of the actual residual using Intel MKL routines and check if
! it is less than 1.D-8
!---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 2) THEN
            descrA % TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
            info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,SOLUTION,beta,TEMP)
           CALL DAXPY(N,-1.D0,rhs,1,TEMP,1)
           Euclidean_norm = DNRM2(N,TEMP,1)
           IF (Euclidean_norm .GT. 1.D-6) THEN
!---------------------------------------------------------------------------
! The solution has not been found yet according to the user-defined stopping
! test. Continue RCI (P)CG iterations.
!---------------------------------------------------------------------------
               GO TO 1
            ELSE
!---------------------------------------------------------------------------
! The solution has been found according to the user-defined stopping test
!---------------------------------------------------------------------------
              GO TO 700
            END IF
!---------------------------------------------------------------------------
! If RCI_request=3, then  apply the simplest SSOR preconditioning
! on vector TMP(:,3) and put the result in vector TMP(:,4)
!---------------------------------------------------------------------------
        ELSE IF (RCI_request .EQ. 3) THEN
            CALL DCOPY(N, TMP(1,3), 1, TMP(1, 4), 1)
            descrA % TYPE = SPARSE_MATRIX_TYPE_TRIANGULAR
            DO I = 1, NITER_SSOR
                CALL DCOPY(N, TMP(1, 3), 1, TEMP, 1)
                descrA % DIAG = SPARSE_DIAG_NON_UNIT
                info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,-one,csrA1,descrA,TMP(1, 4),OMEGA,temp)
                CALL DAXPY(N, ONE-OMEGA, TMP(1,4), 1, TEMP, 1)
                descrA % DIAG = SPARSE_DIAG_UNIT
                info = MKL_SPARSE_D_TRSV(SPARSE_OPERATION_TRANSPOSE,alpha,csrA1,descrA,TEMP,TMP(1,4))
            END DO
            GO TO 1
        ELSE
!---------------------------------------------------------------------------
! If RCI_request=anything else, then dcg subroutine failed
! to compute the solution vector: solution(N)
!---------------------------------------------------------------------------
            GO TO 999
        END IF
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number
!---------------------------------------------------------------------------
700     CALL dcg_get(N,solution,rhs,RCI_request,ipar,dpar,TMP, itercount)
!---------------------------------------------------------------------------
! Print solution vector: solution(N) and number of iterations: itercount
!---------------------------------------------------------------------------
        WRITE(*, *) ' The system has been solved '
        WRITE(*, *) ' The following solution obtained '
        WRITE(*,800) (solution(i), i=1,N)
        WRITE(*, *) ' expected solution '
        WRITE(*,800)(expected_sol(i), i=1,N)
800     FORMAT(4(F10.3))
        WRITE(*,900)(itercount)
900     FORMAT(' Number of iterations: ', I5)
        DO I = 1, N
            expected_sol(I) = expected_sol(I) - solution(I)
        END DO

       Euclidean_norm = DNRM2(N,expected_sol,1)

!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
        CALL MKL_FREE_BUFFERS

        IF (itercount .EQ. expected_itercount .AND. Euclidean_norm .LE. 1.0D-4) THEN
            WRITE( *,'(A,A)') 'This example has successfully PASSED', &
     &  ' through all steps of computation!'
            STOP 0
        ELSE
            WRITE( *,'(A,A,A,I5,A,A,A,E12.5,A)') 'This example may',&
     &  ' have FAILED as either the number of iterations differs from',&
     &  ' the expected number of iterations ',expected_itercount,' or',&
     &  ' the computed solution differs much from the expected',&
     &  ' solution (Euclidean norm is ',Euclidean_norm,'), or both.'
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
        info = MKL_SPARSE_DESTROY(csrA1)
        CALL MKL_FREE_BUFFERS
        STOP 1

      END PROGRAM rci_pcg_f_test_ssor
