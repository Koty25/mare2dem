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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) DSS Fortran-90 example
!
!*******************************************************************************
!--------------------------------------------------------------------------
!
! Example program for solving a symmetric positive definite system of
! equations.
!
!--------------------------------------------------------------------------
INCLUDE 'mkl_dss.f90' ! Include the standard DSS "header file."
PROGRAM solver_f90_test
    use mkl_dss
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(1.0D0)
    INTEGER :: error
    INTEGER :: i
    INTEGER, PARAMETER :: bufLen = 20
! Define the data arrays and the solution and rhs vectors.
    INTEGER, ALLOCATABLE :: columns( : )
    INTEGER :: nCols
    INTEGER :: nNonZeros
    INTEGER :: nRhs
    INTEGER :: nRows
    REAL(KIND=DP), ALLOCATABLE :: rhs( : )
    INTEGER, ALLOCATABLE :: rowIndex( : )
    REAL(KIND=DP), ALLOCATABLE :: solution( : )
    REAL(KIND=DP), ALLOCATABLE :: values( : )
    TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
    REAL(KIND=DP),ALLOCATABLE::statOUt( : )
    CHARACTER*15 statIn
    INTEGER perm(1)
    INTEGER buff(bufLen)
! Set the problem to be solved.
    nRows = 5
    nCols = 5
    nNonZeros = 9
    nRhs = 1
    perm(1) = 0
    ALLOCATE(rowIndex(nRows + 1))
    rowIndex = (/ 1, 6, 7, 8, 9, 10 /)
    ALLOCATE( columns( nNonZeros ) )
    columns = (/ 1, 2, 3, 4, 5, 2, 3, 4, 5 /)
    ALLOCATE(values(nNonZeros))
    values = (/ 9.0_DP, 1.5_DP, 6.0_DP, 0.75_DP, 3.0_DP, 0.5_DP, 12.0_DP, &
              0.625_DP, 16.0_DP /)
    ALLOCATE(rhs(nRows))
    rhs = (/ 1.0_DP, 2.0_DP, 3.0_DP, 4.0_DP, 5.0_DP /)
! Initialize the solver.
    error = DSS_CREATE(handle, MKL_DSS_DEFAULTS)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
! Define the non-zero structure of the matrix.
    error = DSS_DEFINE_STRUCTURE(handle, MKL_DSS_SYMMETRIC, rowIndex, nRows, &
                                 nCols, columns, nNonZeros)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
! Reorder the matrix.
    error = DSS_REORDER(handle, MKL_DSS_DEFAULTS, perm)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
! Factor the matrix.
    error = DSS_FACTOR_REAL(handle, MKL_DSS_DEFAULTS, values)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
! Allocate the solution vector and solve the problem.
    ALLOCATE(solution(nRows))
    error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, rhs, nRhs, solution)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
! Print Out the determinant of the matrix (no statistics for a diagonal matrix)
    IF(nRows < nNonZeros) THEN
       ALLOCATE(statOut(5))
       statIn = 'determinant'
       call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
       error = DSS_STATISTICS(handle, MKL_DSS_DEFAULTS, buff, statOut)
       IF (error /= MKL_DSS_SUCCESS) GOTO 999
       WRITE(*,"('pow of determinant is '(5F10.3))") (statOut(1))
       WRITE(*,"('base of determinant is '(5F10.3))") (statOut(2))
       WRITE(*,"('Determinant is '(5F10.3))") ((10**statOut(1))*statOut(2))
    END IF
! Deallocate solver storage and various local arrays.
    error = DSS_DELETE(handle, MKL_DSS_DEFAULTS)
    IF (error /= MKL_DSS_SUCCESS) GOTO 999
    IF (ALLOCATED(rowIndex)) DEALLOCATE(rowIndex)
    IF (ALLOCATED(columns )) DEALLOCATE(columns)
    IF (ALLOCATED(values )) DEALLOCATE(values)
    IF (ALLOCATED(rhs)) DEALLOCATE(rhs)
    IF (ALLOCATED(statOut)) DEALLOCATE(statOut)
! Print the solution vector, deallocate it and exit
    WRITE(*,"('Solution Array: '(5F10.3))") (solution(i), i = 1, nCols)
    IF (ALLOCATED(solution)) DEALLOCATE(solution)
    GOTO 1000
! Print an error message and exit
999 WRITE(*,*) "Solver returned error code ", error
    STOP 1
1000 CONTINUE
END PROGRAM solver_f90_test
