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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) PARDISO Fortran-77
!             example
!
!*******************************************************************************
!----------------------------------------------------------------------
! Example program to show the use of the "MATRIX_CHECK" routine
! for symmetric linear systems
!---------------------------------------------------------------------

      PROGRAM MATRIX_CHECK

        USE, INTRINSIC :: ISO_C_BINDING ! required for 'mkl_sparse_handle.fi'
        IMPLICIT NONE

        INCLUDE 'mkl_sparse_handle.fi'

        INTEGER N, NNZ, CHECK_RESULT_CODE, ERROR
        INTEGER, ALLOCATABLE :: IA(:), JA(:)
        TYPE(SPARSE_STRUCT) PT

        N   = 8
        NNZ = 18
        ERROR = 0
        ALLOCATE(IA(N + 1))
        IA = (/ 1, 5, 8, 10, 12, 15, 17, 18, 19 /)
        ALLOCATE(JA(NNZ))
        JA = (/ 1,    3,       6, 7,
     .             2, 3,    5,
     .                3,             8,
     .                3,          7,
     .                      5, 6, 7,
     .                         6,    8,
     .                            7,
     .                               8 /)

        CALL SPARSE_MATRIX_CHECKER_INIT(PT)

        PT % N = N
        PT % CSR_IA = LOC(IA)
        PT % CSR_JA = LOC(JA)
        PT % INDEXING         = MKL_ONE_BASED
        PT % MATRIX_STRUCTURE = MKL_UPPER_TRIANGULAR
        PT % PRINT_STYLE      = MKL_FORTRAN_STYLE
        PT % MESSAGE_LEVEL    = MKL_PRINT

        CHECK_RESULT_CODE = SPARSE_MATRIX_CHECKER(PT)


        IF (CHECK_RESULT_CODE .EQ. MKL_SPARSE_CHECKER_NONTRIANGULAR)
     *      THEN
            WRITE(*,*)'Matrix check result:
     *                 MKL_SPARSE_CHECKER_NONTRIANGULAR'
            ERROR = 0
        ELSE
            IF (CHECK_RESULT_CODE .EQ. MKL_SPARSE_CHECKER_NON_MONOTONIC)
     *      THEN
            WRITE(*,*)'Matrix check result:
     *                     MKL_SPARSE_CHECKER_NON_MONOTONIC'
            ENDIF
            IF (CHECK_RESULT_CODE .EQ. MKL_SPARSE_CHECKER_SUCCESS)
     *      THEN
            WRITE(*,*)'Matrix check result:
     *                     MKL_SPARSE_CHECKER_SUCCESS'
            ENDIF
            IF (CHECK_RESULT_CODE .EQ. MKL_SPARSE_CHECKER_OUT_OF_RANGE)
     *      THEN
            WRITE(*,*)'Matrix check result:
     *                     MKL_SPARSE_CHECKER_OUT_OF_RANGE'
            ENDIF
            IF (CHECK_RESULT_CODE .EQ. MKL_SPARSE_CHECKER_NONORDERED)
     *      THEN
            WRITE(*,*)'Matrix check result:
     *                     MKL_SPARSE_CHECKER_NONORDERED'
            ENDIF
            ERROR = 1
        ENDIF

        WRITE(*,
     .      '("Matrix check details: (",i0,", ",i0,", ",i0,")")')
     .      PT % CHECK_RESULT


        IF (ERROR .NE. 0) STOP 1

      END PROGRAM MATRIX_CHECK
