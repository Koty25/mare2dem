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

*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) DSS Fortran-77 example
*
********************************************************************************
C---------------------------------------------------------------------------
C Example program for solving symmetric positive definite system of
C equations.
C---------------------------------------------------------------------------
      PROGRAM solver_f_test
        IMPLICIT NONE
        INCLUDE 'mkl_dss.fi'
C---------------------------------------------------------------------------
C Define the array and rhs vectors
C---------------------------------------------------------------------------
        INTEGER nRows, nCols, nNonZeros, i, nRhs
        PARAMETER (nRows = 5,
     1 nCols = 5,
     2 nNonZeros = 9,
     3 nRhs = 1)
        INTEGER rowIndex(nRows + 1), columns(nNonZeros), idum(1)
        DOUBLE PRECISION values(nNonZeros), rhs(nRows)
        DATA rowIndex / 1, 6, 7, 8, 9, 10 /
        DATA columns / 1, 2, 3, 4, 5, 2, 3, 4, 5 /
        DATA values / 9, 1.5, 6, .75, 3, 0.5, 12, .625, 16 /
        DATA rhs / 1, 2, 3, 4, 5 /
C---------------------------------------------------------------------------
C Allocate storage for the solver handle and the solution vector
C---------------------------------------------------------------------------
        DOUBLE PRECISION solution(nRows)
        INTEGER*8 handle
        INTEGER error
        CHARACTER*15 statIn
        DOUBLE PRECISION statOut(5)
        INTEGER bufLen
        PARAMETER(bufLen = 20)
        INTEGER buff(bufLen)
C---------------------------------------------------------------------------
C Initialize the solver
C---------------------------------------------------------------------------
        error = dss_create(handle, MKL_DSS_DEFAULTS)
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Define the non-zero structure of the matrix
C---------------------------------------------------------------------------
        error = dss_define_structure( handle, MKL_DSS_SYMMETRIC,
     &  rowIndex, nRows, nCols, columns, nNonZeros )
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Reorder the matrix
C---------------------------------------------------------------------------
        error = dss_reorder( handle, MKL_DSS_DEFAULTS, idum)
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Factor the matrix
C---------------------------------------------------------------------------
        error = dss_factor_real( handle, MKL_DSS_DEFAULTS, VALUES)
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Get the solution vector
C---------------------------------------------------------------------------
        error = dss_solve_real( handle, MKL_DSS_DEFAULTS, rhs, nRhs,
     &  solution)
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Print Determinant of the matrix (no statistics for a diagonal matrix)
C---------------------------------------------------------------------------
        IF( nRows .LT. nNonZeros ) THEN
            statIn = 'determinant'
            call mkl_cvt_to_null_terminated_str(buff,bufLen,statIn)
            error = dss_statistics(handle, MKL_DSS_DEFAULTS, buff,
     &  statOut)
            WRITE(*,"(' pow of determinant is ', 5(F10.3))") statOut(1)
            WRITE(*,"(' base of determinant is ', 5(F10.3))") statOut(2)
            WRITE(*,"(' Determinant is ', 5(F10.3))")(10**statOut(1))*
     &  statOut(2)
        END IF
C---------------------------------------------------------------------------
C Deallocate solver storage
C---------------------------------------------------------------------------
        error = dss_delete( handle, MKL_DSS_DEFAULTS )
        IF (error .NE. MKL_DSS_SUCCESS ) GO TO 999
C---------------------------------------------------------------------------
C Print solution vector
C---------------------------------------------------------------------------
        WRITE(*,900) (solution(i), i = 1, nCols)
  900   FORMAT(' Solution Array: ',5(F10.3))
        GO TO 1000
  999   WRITE(*,*) "Solver returned error code ", error
        STOP 1
 1000 END PROGRAM solver_f_test
