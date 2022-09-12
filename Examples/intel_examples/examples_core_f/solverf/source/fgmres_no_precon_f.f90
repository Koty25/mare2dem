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

!  Content:
!  Intel(R) Math Kernel Library (Intel(R) MKL) RCI (P)FGMRES ((Preconditioned)
!                          Flexible Generalized Minimal RESidual method) example
!*******************************************************************************

!---------------------------------------------------------------------------
!  Example program for solving non-symmetric indefinite system of equations
!  Simplest case: no preconditioning and no user-defined stopping tests
!---------------------------------------------------------------------------

      PROGRAM FGMRES_NO_PRECON_F
        USE MKL_SPBLAS
        INCLUDE "mkl_rci.fi"

        INTEGER N
        PARAMETER (N=5)
        INTEGER SIZE
        PARAMETER (SIZE=128)
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER IA(6)
        DATA IA /1,3,6,9,12,14/
        INTEGER JA(13)
        DATA JA /1,  3,    &
&                1,2,  4,  &
&                  2,3,  5,&
&                    3,4,5,&
&                      4,5/
        DOUBLE PRECISION A(13)
        DATA A  / 1.0,     -1.0,          &
&                -1.0, 1.0,     -1.0,     &
&                      1.0,-2.0,      1.0,&
&                          -1.0, 2.0,-1.0,&
&                               -1.0,-3.0/
!---------------------------------------------------------------------------
! Allocate storage for the ?par parameters and the solution/rhs vectors
!---------------------------------------------------------------------------
        INTEGER IPAR(SIZE)
        DOUBLE PRECISION DPAR(SIZE), TMP(N*(2*N+1)+(N*(N+9))/2+1)
        DOUBLE PRECISION EXPECTED_SOLUTION(N)
        DATA EXPECTED_SOLUTION /-1.0,1.0,0.0,1.0,-1.0/
        DOUBLE PRECISION RHS(N)
        DOUBLE PRECISION COMPUTED_SOLUTION(N)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
        INTEGER ITERCOUNT, EXPECTED_ITERCOUNT
        PARAMETER (EXPECTED_ITERCOUNT=5)
        INTEGER RCI_REQUEST, I, info
        DOUBLE PRECISION DVAR
      
        DOUBLE PRECISION DNRM2
        EXTERNAL DNRM2
        DOUBLE PRECISION alpha, beta
!   Matrix descriptor
        TYPE(MATRIX_DESCR) descrA
!   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) csrA
        !   Create matrix descriptor
        descrA % TYPE = SPARSE_MATRIX_TYPE_GENERAL
        descrA % MODE = SPARSE_FILL_MODE_UPPER
        descrA % DIAG = SPARSE_DIAG_NON_UNIT
        alpha = 1.0
        beta  = 0.0
        info = MKL_SPARSE_D_CREATE_CSR(csrA,SPARSE_INDEX_BASE_ONE,N,N,IA,IA(2),JA,A)
        DO I = 1, N
            RHS(I)     = 0.D0
        END DO
        DO I = 1, N*(2*N+1)+(N*(N+9))/2+1
            TMP(I)     = 0.D0
        END DO
        PRINT *,'--------------------------------------------------'
        PRINT *,'The SIMPLEST example of usage of RCI FGMRES solver'
        PRINT *,'to solve a non-symmetric indefinite non-degenerate'
        PRINT *,'       algebraic system of linear equations'
        PRINT *,'--------------------------------------------------'
!---------------------------------------------------------------------------
! Initialize variables and the right hand side through matrix-vector product
!---------------------------------------------------------------------------
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,EXPECTED_SOLUTION,beta,RHS)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
        DO I = 1, N
            COMPUTED_SOLUTION(I) = 1.0D0
        END DO
!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
        CALL DFGMRES_INIT(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, &
&       DPAR, TMP)
        IF (RCI_REQUEST .NE. 0) GO TO 999
!---------------------------------------------------------------------------
! Set the desired parameters:
! LOGICAL parameters:
! do residual stopping test
! do not request for the user defined stopping test
! do the check of the norm of the next generated vector automatically
! DOUBLE PRECISION parameters
! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
!---------------------------------------------------------------------------
        IPAR(9) = 1
        IPAR(10) = 0
        IPAR(12) = 1
        DPAR(1) = 1.0D-3
!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
        CALL DFGMRES_CHECK(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, &
&       IPAR, DPAR, TMP)
        IF (RCI_REQUEST .NE. 0 .AND. RCI_REQUEST .NE. -1001) GO TO 999
!---------------------------------------------------------------------------
! Print the info about the RCI FGMRES method
!---------------------------------------------------------------------------
        WRITE( *,'(A)') ' '
        WRITE( *,'(A,A)') 'Some info about the current run of', &
&       ' RCI FGMRES method:'
        WRITE( *,'(A)') ' '
        IF (IPAR(8) .NE. 0) THEN
            WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the', &
&       ' automatic test for the maximal number of iterations will', &
&       ' be performed'
        ELSE
            WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the', &
&       ' automatic test for the maximal number of iterations will be', &
&       ' skipped'
        END IF
        WRITE( *,'(A)') '+++'
        IF (IPAR(9) .NE. 0) THEN
            WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the', &
&       ' automatic residual test will be performed'
        ELSE
            WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the', &
&       ' automatic residual test will be skipped'
        END IF
        WRITE( *,'(A)') '+++'
        IF (IPAR(10) .NE. 0) THEN
            WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the', &
&       ' user-defined stopping test will be requested via', &
&       ' RCI_REQUEST=2'
        ELSE
            WRITE(*,'(A,I1,A,A,A)') 'As IPAR(10)=',IPAR(10),', the', &
&       ' user-defined stopping test will not be requested, thus,', &
&       ' RCI_REQUEST will not take the value 2'
        END IF
        WRITE( *,'(A)') '+++'
        IF (IPAR(11) .NE. 0) THEN
            WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the', &
&       ' Preconditioned FGMRES iterations will be performed, thus,'
        WRITE(*,'(A,A)') 'the preconditioner action will be requested', &
&       ' via RCI_REQUEST=3'
        ELSE
            WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the', &
&       ' Preconditioned FGMRES iterations will not be performed,'
            WRITE( *,'(A)') 'thus, RCI_REQUEST will not take the', &
&       ' value 3'
        END IF
        WRITE( *,'(A)') '+++'
        IF (IPAR(12) .NE. 0) THEN
            WRITE(*,'(A,I1,A,A)')'As IPAR(12)=',IPAR(12),', the', &
&       ' automatic test for the norm of the next generated vector', &
&       ' is not'
            WRITE( *,'(A,A)') ' equal to zero up to rounding and', &
&       ' computational errors will be performed,'
          WRITE( *,'(A)') 'thus, RCI_REQUEST will not take the value 4'
        ELSE
            WRITE(*,'(A,I1,A,A)')'As IPAR(12)=',IPAR(12),', the', &
&       ' automatic test for the norm of the next generated vector is'
          WRITE(*,'(A,A)') 'not equal to zero up to rounding and', &
&       ' computational errors will be skipped,'
          WRITE(*,'(A,A)') 'thus, the user-defined test will be', &
&       ' requested via RCI_REQUEST=4'
        END IF
        WRITE( *,'(A)') '+++'
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)FGMRES solver without preconditioning
! Reverse Communication starts here
!---------------------------------------------------------------------------
1       CALL DFGMRES(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, &
&       DPAR, TMP)
!---------------------------------------------------------------------------
! If RCI_REQUEST=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
        IF (RCI_REQUEST .EQ. 0) GO TO 3
!---------------------------------------------------------------------------
! If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
! and put the result in vector TMP(IPAR(23))
!---------------------------------------------------------------------------
        IF (RCI_REQUEST .EQ. 1) THEN
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,&
&                              descrA,TMP(IPAR(22)),beta,TMP(IPAR(23)))
            GO TO 1
!---------------------------------------------------------------------------
! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
! to compute the solution vector: COMPUTED_SOLUTION(N)
!---------------------------------------------------------------------------
        ELSE
            GO TO 999
        END IF
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution (DO NOT FORGET to
! call DFGMRES_GET routine as COMPUTED_SOLUTION is still containing
! the initial guess!)
!---------------------------------------------------------------------------
3        CALL DFGMRES_GET(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, &
&       DPAR, TMP, ITERCOUNT)
!---------------------------------------------------------------------------
! Print solution vector: COMPUTED_SOLUTION(N) and
! the number of iterations: ITERCOUNT
!---------------------------------------------------------------------------
        WRITE( *,'(A)') ' '
        WRITE( *,'(A)') 'The system has been solved'
        WRITE( *,'(A)') ' '
        WRITE( *,'(A)') 'The following solution has been obtained:'
        DO I = 1, N
            WRITE(*,'(A18,I1,A2,E10.3)') 'COMPUTED_SOLUTION(',I,')=', &
&       COMPUTED_SOLUTION(I)
        END DO
        WRITE( *,'(A)') ' '
        WRITE( *,'(A)') 'The expected solution is:'
        DO I = 1, N
            WRITE(*,'(A18,I1,A2,E10.3)') 'EXPECTED_SOLUTION(',I,')=', &
&       EXPECTED_SOLUTION(I)
            EXPECTED_SOLUTION(I) = EXPECTED_SOLUTION(I) - COMPUTED_SOLUTION(I)
        END DO
        WRITE( *,'(A)') ' '
        WRITE( *,'(A,I2)') 'Number of iterations: ',ITERCOUNT
        WRITE( *,'(A)') ' '

!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
        CALL MKL_FREE_BUFFERS

        DVAR = DNRM2(N,EXPECTED_SOLUTION,1)
        IF (ITERCOUNT .EQ. EXPECTED_ITERCOUNT .AND. DVAR .LE. 1.0D-14) THEN
            WRITE( *,'(A,A)') 'This example has successfully PASSED',' through all steps of computation!'
            STOP 0
        ELSE
            WRITE( *,'(A,A,A,I5,A,A,A,E12.5,A)') 'This example may', &
     &  ' have FAILED as either the number of iterations differs from', &
     &  ' the expected number of iterations ',EXPECTED_ITERCOUNT,' or', &
     &  ' the computed solution differs much from the expected', &
     &  ' solution (Euclidean norm is ',DVAR,'), or both.'
            STOP 1
        END IF

!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
999     WRITE( *,'(A,A,I5)') 'This example FAILED as the solver has', &
     &  ' returned the ERROR code', RCI_REQUEST
        info = MKL_SPARSE_DESTROY(csrA)
        CALL MKL_FREE_BUFFERS
        STOP 1

      END PROGRAM FGMRES_NO_PRECON_F
