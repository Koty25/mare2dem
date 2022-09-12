!===============================================================================
! Copyright 2006-2020 Intel Corporation.
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
!  Intel(R) Math Kernel Library (Intel(R) MKL) example of RCI Flexible
!  Generalized Minimal RESidual method with ILU0 Preconditioner
!*******************************************************************************

!---------------------------------------------------------------------------
!  Example program for solving non-degenerate system of equations.
!  Full functionality of RCI FGMRES solver is exploited. Example shows how
!  ILU0 preconditioner accelerates the solver by reducing the number of
!  iterations.
!---------------------------------------------------------------------------

      PROGRAM FGMRES_FULL_FUNCT_F
        USE MKL_SPBLAS
        IMPLICIT NONE
        INCLUDE "mkl_rci.fi"

        INTEGER N
        PARAMETER(N=4)
        INTEGER SIZE
        PARAMETER (SIZE=128)
!---------------------------------------------------------------------------
! Define arrays for the upper triangle of the coefficient matrix
! Compressed sparse row storage is used for sparse representation
!---------------------------------------------------------------------------
        INTEGER IA(5)
        DATA IA /1,4,7,10,13/
        INTEGER JA(12)
        DATA JA /1,2,3,1,2,4,1,3,4,2,3,4/
        DOUBLE PRECISION A(12),BILU0(12),TRVEC(N)
        DATA A / 4.,-1.,-1.,-1.,4.,-1.,-1.,4.,-1.,-1.,-1.,4./
!---------------------------------------------------------------------------
! Allocate storage for the ?par parameters and the solution/rhs/residual vectors
!---------------------------------------------------------------------------
        INTEGER IPAR(SIZE),IERR
        DOUBLE PRECISION DPAR(SIZE), TMP(N*(2*N+1)+(N*(N+9))/2+1)
        DOUBLE PRECISION EXPECTED_SOLUTION(N)
        DATA EXPECTED_SOLUTION /1.0,1.0,1.0,1.0/
        DOUBLE PRECISION RHS(N), B(N)
        DOUBLE PRECISION COMPUTED_SOLUTION(N)
        DOUBLE PRECISION RESIDUAL(N)

        INTEGER MATSIZE, INCX, REF_NIT
        DOUBLE PRECISION REF_NORM2, NRM2
        PARAMETER (MATSIZE=12, INCX=1, REF_NIT=2, REF_NORM2=7.772387D0)
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
        INTEGER ITERCOUNT
        INTEGER RCI_REQUEST, I, info
        DOUBLE PRECISION DVAR
!---------------------------------------------------------------------------
! An external BLAS function is taken from Intel MKL BLAS to use
! with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
        DOUBLE PRECISION DNRM2
        EXTERNAL DNRM2
!---------------------------------------------------------------------------
! Some additional variables to use for ILUT preconditioner call
!---------------------------------------------------------------------------
        INTEGER MAXFIL
        DOUBLE PRECISION TOL
        DOUBLE PRECISION alpha, beta
!   Matrix descriptor
        TYPE(MATRIX_DESCR) descrA, descrL
!   CSR matrix representation
        TYPE(SPARSE_MATRIX_T) csrA, csrL
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

        WRITE( *,'(A,A)') '-----------------------------------------', &
&       '--------------------------'
        WRITE(*,'(A,A)') 'The FULLY ADVANCED example RCI FGMRES with', &
&       ' ILU0 preconditioner'
        WRITE(*,'(A,A)') 'to solve the non-degenerate algebraic system', &
&       ' system of linear equations'
        WRITE( *,'(A,A)') '-----------------------------------------', &
&       '--------------------------'

!---------------------------------------------------------------------------
! Initialize variables and the right hand side through matrix-vector product
!---------------------------------------------------------------------------
        info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrA,descrA,EXPECTED_SOLUTION,beta,RHS)
!---------------------------------------------------------------------------
! Save the right-hand side in vector B for future use
!---------------------------------------------------------------------------
        CALL DCOPY(N, RHS, 1, B, 1)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------
        DO I = 1, N
            COMPUTED_SOLUTION(I) = 0.D0
        END DO
        COMPUTED_SOLUTION(1) = 100.D0

!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
        CALL DFGMRES_INIT(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, &
&       DPAR, TMP)
        IF (RCI_REQUEST .NE. 0) GO TO 999

!---------------------------------------------------------------------------
! Calculate ILU0 preconditioner.
!                      !ATTENTION!
! DCSRILU0 routine uses some IPAR, DPAR set by DFGMRES_INIT routine.
! Important for DCSRILU0 default entries set by DFGMRES_INIT are
! ipar(2) = 6 - output of error messages to the screen,
! ipar(6) = 1 - allow output of error messages,
! ipar(31)= 0 - abort DCSRILU0 calculations if routine meets zero diagonal element.
!
! If ILU0 is going to be used out of Intel MKL FGMRES context, than the values
! of ipar(2), ipar(6), ipar(31), dpar(31), and dpar(32) should be user
! provided before the DCSRILU0 routine call.
!
! In this example, specific for DCSRILU0 entries are set in turn:
! ipar(31)= 1 - change small diagonal value to that given by dpar(32),
! dpar(31)= 1.D-20 instead of the default value set by DFGMRES_INIT.
!                  It is a small value to compare a diagonal entry with it.
! dpar(32)= 1.D-16 instead of the default value set by DFGMRES_INIT.
!                  It is the target value of the diagonal value if it is
!                  small as compared to dpar(31) and the routine should change
!                  it rather than abort DCSRILU0 calculations.
!---------------------------------------------------------------------------

          IPAR(31) = 1
          DPAR(31) = 1.D-20
          DPAR(32) = 1.D-16
          CALL DCSRILU0(N, A, IA, JA, BILU0, IPAR, DPAR, IERR)
        info = MKL_SPARSE_D_CREATE_CSR(csrL,SPARSE_INDEX_BASE_ONE,N,N,IA,&
&                                      IA(2),JA,BILU0)
          NRM2 = DNRM2(MATSIZE, BILU0, INCX)

        IF(IERR .NE. 0) THEN
            WRITE(*,'(A,A,I1)') ' Error after calculation of the',&
&       ' preconditioner DCSRILU0', IERR
            GO TO 998
        END IF

!---------------------------------------------------------------------------
! Set the desired parameters:
! do the restart after 2 iterations
! LOGICAL parameters:
! do not do the stopping test for the maximal number of iterations
! do the Preconditioned iterations of FGMRES method
! Set parameter IPAR(11) for preconditioner call. For this example,
! it reduces the number of iterations.
! DOUBLE PRECISION parameters
! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
! NOTE. Preconditioner may increase the number of iterations for an
! arbitrary case of the system and initial guess and even ruin the
! convergence. It is user's responsibility to use a suitable preconditioner
! and to apply it skillfully.
!---------------------------------------------------------------------------
        IPAR(15) = 2
        IPAR(8) = 0
        IPAR(11) = 1
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
! Compute the solution by RCI (P)FGMRES solver with preconditioning
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
        END IF
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test
! The residual stopping test for the computed solution is performed here
!---------------------------------------------------------------------------
! NOTE: from this point vector B(N) is no longer containing the right-hand
! side of the problem! It contains the current FGMRES approximation to the
! solution. If you need to keep the right-hand side, save it in some other
! vector before the call to DFGMRES routine. Here we saved it in vector
! RHS(N). The vector B is used instead of RHS to preserve the original
! right-hand side of the problem and guarantee the proper restart of FGMRES
! method. Vector B will be altered when computing the residual stopping
! criterion!
!---------------------------------------------------------------------------
        IF (RCI_REQUEST .EQ. 2) THEN
! Request to the DFGMRES_GET routine to put the solution into B(N) via IPAR(13)
            IPAR(13)=1
! Get the current FGMRES solution in the vector B(N)
            CALL DFGMRES_GET(N, COMPUTED_SOLUTION, B, RCI_REQUEST,  &
&       IPAR, DPAR, TMP, ITERCOUNT)
! Compute the current true residual via Intel MKL (Sparse) BLAS routines
            info = MKL_SPARSE_D_MV(SPARSE_OPERATION_NON_TRANSPOSE,alpha, &
&                                  csrA,descrA,B,beta,RESIDUAL)
            CALL DAXPY(N, -1.0D0, RHS, 1, RESIDUAL, 1)
            DVAR = DNRM2(N, RESIDUAL, 1)
            IF (DVAR .LT. 1.0E-3) THEN
                GO TO 3
            ELSE
                GO TO 1
            END IF
        END IF
!---------------------------------------------------------------------------
! If RCI_REQUEST=3, then apply the preconditioner on the vector
! TMP(IPAR(22)) and put the result in vector TMP(IPAR(23))
! Here is the recommended usage of the result produced by ILU0 routine
! via standard Intel MKL Sparse Blas solver routine mkl_dcsrtrsv.
!---------------------------------------------------------------------------
        IF (RCI_REQUEST .EQ. 3) THEN
             descrL % TYPE = SPARSE_MATRIX_TYPE_TRIANGULAR
             descrL % MODE = SPARSE_FILL_MODE_LOWER
             descrL % DIAG = SPARSE_DIAG_UNIT
             info = MKL_SPARSE_D_TRSV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrL,&
&                                     descrL,TMP(IPAR(22)),trvec)

             descrL % MODE = SPARSE_FILL_MODE_UPPER
             descrL % DIAG = SPARSE_DIAG_NON_UNIT
             info = MKL_SPARSE_D_TRSV(SPARSE_OPERATION_NON_TRANSPOSE,alpha,csrL,&
&                                     descrL,trvec,TMP(IPAR(23)))
            GO TO 1
        END IF
!---------------------------------------------------------------------------
! If RCI_REQUEST=4, then check if the norm of the next generated vector is
! not zero up to rounding and computational errors. The norm is contained
! in DPAR(7) parameter
!---------------------------------------------------------------------------
        IF (RCI_REQUEST .EQ. 4) THEN
        IF (DPAR(7) .LT. 1.0D-12) THEN
            GO TO 3
        ELSE
            GO TO 1
        END IF
!---------------------------------------------------------------------------
! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
! to compute the solution vector: COMPUTED_SOLUTION(N)
!---------------------------------------------------------------------------
        ELSE
            GO TO 999
        END IF
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector COMPUTED_SOLUTION(N) via IPAR(13)
!---------------------------------------------------------------------------
3       IPAR(13) = 0
        CALL DFGMRES_GET(N, COMPUTED_SOLUTION, RHS, RCI_REQUEST, IPAR, &
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

        IF (ITERCOUNT .EQ. REF_NIT .AND.  &
&               DABS(REF_NORM2-NRM2) .LT. 1.D-6) THEN
            WRITE( *,'(A)') ' '
            WRITE( *,'(A)') '---------------------------------------', &
&       '------'
            WRITE( *,'(A,A)') 'Fortran example of FGMRES with ILU0', &
&       ' preconditioner '
            WRITE( *,'(A,A)') 'has successfully PASSED all stages of', &
&       ' computations'
            WRITE( *,'(A)') '---------------------------------------', &
&       '------'
            WRITE( *,'(A)') ' '
            STOP 0
        ELSE
            WRITE( *,'(A,A)') 'Probably, the preconditioner was', &
&       ' computed incorrectly:'
            WRITE( *,'(A,F9.6,A,F9.6)')  &
&       'Either the preconditioner norm',NRM2, &
&       ' differs from the expected norm',REF_NORM2
            WRITE( *,'(A,I2,A,I2)'), &
&       'and/or the number of iterations ', ITERCOUNT, &
&       ' differs from the expected number ', REF_NIT
            WRITE( *,'(A)') ' '
            WRITE( *,'(A,A)') '-------------------------------------', &
&       '---------------------------'
            WRITE( *,'(A,A)') 'Unfortunately, FGMRES+ILU0 Fortran', &
&       ' example has FAILED'
            WRITE( *,'(A,A)') '-------------------------------------', &
&        '---------------------------'
            WRITE( *,'(A)') ' '
            STOP 1
        END IF
999     WRITE( *,'(A,I2)') 'The solver has returned the ERROR code ', &
&       RCI_REQUEST
!---------------------------------------------------------------------------
! Release internal Intel MKL memory that might be used for computations
! NOTE: It is important to call the routine below to avoid memory leaks
! unless you disable Intel MKL Memory Manager
!---------------------------------------------------------------------------
998     WRITE( *,'(A)') ' '
        WRITE( *,'(A,A)') '-----------------------------------------', &
&       '--------------------------'
        WRITE( *,'(A,A)') 'Unfortunately, FGMRES+ILU0 Fortran example', &
&       ' has FAILED'
        WRITE( *,'(A,A)') '-----------------------------------------', &
&       '--------------------------'
        WRITE( *,'(A)') ' '
        info = MKL_SPARSE_DESTROY(csrA)
        info = MKL_SPARSE_DESTROY(csrL)
        CALL MKL_FREE_BUFFERS
        STOP 1

      END PROGRAM FGMRES_FULL_FUNCT_F
