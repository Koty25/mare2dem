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

*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) PARDISO Fortran example
*
********************************************************************************
C----------------------------------------------------------------------
C       Example program to show the use of the "PARDISO" routine
C       for complex unsymmetric linear systems with
C       storing/restoring PARDISO handle between steps
C---------------------------------------------------------------------
        PROGRAM pardiso_unsym_store_restore
        IMPLICIT NONE
        include 'mkl_pardiso.fi'
C..     Internal solver memory pointer for 64-bit architectures
C..     INTEGER*8 pt(64)
C..     Internal solver memory pointer for 32-bit architectures
C..     INTEGER*4 pt(64)
C..     This is OK in both cases.
        TYPE(MKL_PARDISO_HANDLE) pt_11(64), pt_22(64), pt_33(64)

C..     All other variables

        INTEGER     maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        INTEGER     iparm(64)
        INTEGER     ia(9)
        INTEGER     ja(20)
        COMPLEX*16  a(20)
        COMPLEX*16  b(8)
        COMPLEX*16  x(8)

        INTEGER i, idum(1)
        COMPLEX*16 ddum(1)
        REAL*8  waltime1, waltime2
        CHARACTER*15 pathToFile
        INTEGER bufLen
        PARAMETER(bufLen = 20)
        CHARACTER buff(bufLen)


C.. Fill all arrays containing matrix data.

        DATA n /8/, nrhs /1/, maxfct /1/, mnum /1/

        DATA ia /1,5,8,10,12,13,16,18,21/

        DATA ja
     1        /1,          3,                 6,    7,
     2               2,    3,          5,
     3                     3,                             8,
     4                          4,                  7,
     5               2,
     6                     3,                 6,          8,
     7               2,                             7,
     8                     3,                       7,    8/

        DATA a
     1     /(7.d0, 1.d0), (1.d0,1.d0), (2.d0,1.d0), (7.d0,1.d0),
     2      (-4.d0,0.d0), (8.d0,1.d0), (2.d0,1.d0),
     3      (1.d0,1.d0),  (5.d0,1.d0),
     4      (7.d0,0.d0),  (9.d0,1.d0),
     5      (-4d0,1.d0),
     6      (7.d0,1.d0),  (3.d0,1.d0), (8.d0,0.d0),
     7      (1.d0,1.d0),  (11.d0,1.d0),
     8      (-3.d0,1.d0), (2.d0,1.d0), (5.d0,0.d0)/


C
C  .. Setup Pardiso control parameters and initialize the solvers
C     internal address pointers. This is only necessary for the FIRST
C     call of the PARDISO solver.
C
        DO i = 1, 64
            iparm(i) = 0
        END DO

        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(3) = 1 ! numbers of processors, value of OMP_NUM_THREADS
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(7) = 0 ! not in use
        iparm(8) = 2 ! numbers of iterative refinement steps
        iparm(9) = 0 ! not in use
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(12) = 0 ! not in use
        iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric).
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(15) = 0 ! not in use
        iparm(16) = 0 ! not in use
        iparm(17) = 0 ! not in use
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations
        error = 0  ! initialize error flag
        msglvl = 1 ! print statistical information
        mtype = 13 ! complex unsymmetric matrix

C.. Initialize the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
        DO i = 1, 64
            pt_11(i)%DUMMY = 0
            pt_22(i)%DUMMY = 0
            pt_33(i)%DUMMY = 0
        END DO

C..   Reordering and Symbolic Factorization, This step also allocates
C     all memory that is necessary for the factorization
C
        phase     = 11      ! only reordering and symbolic factorization
        CALL pardiso (pt_11, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1              idum, nrhs, iparm, msglvl, ddum, ddum, error)

        WRITE(*,*) 'Reordering completed ... '

        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF

        WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
        WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

C..   Storing current handle in work folder, deallocate it and        
C     restoring new handle                                            

        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_store(pt_11,buff,error)
        IF (error .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           STOP 1
        ENDIF
      
        phase = -1 ! release internal memory
        CALL pardiso (pt_11, maxfct, mnum, mtype, phase, n, ddum, idum, 
     1         idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
      
        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_restore(pt_22,buff,error)
        IF (error .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           STOP 1
        ENDIF
      
        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_delete(buff,error)
        IF (error .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           STOP 1
        ENDIF

C.. Factorization.
        phase     = 22  ! only factorization
        CALL pardiso (pt_22, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1              idum, nrhs, iparm, msglvl, ddum, ddum, error)

        WRITE(*,*) 'Factorization completed ...  '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF

C..   Storing current handle in work folder, deallocate it and        
C     restoring new handle                                            

        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_store(pt_22,buff,error)
        IF (error .NE. 0) THEN
          WRITE(*,*) 'The following ERROR was detected: ', error
          STOP 1
        ENDIF

        phase = -1 ! release internal memory
        CALL pardiso (pt_22, maxfct, mnum, mtype, phase, n, ddum, idum, 
     1        idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)

        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_restore(pt_33,buff,error)
        IF (error .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           STOP 1
        ENDIF

        pathToFile = ''
        CALL mkl_cvt_to_null_terminated_str(buff,bufLen,pathToFile)
        CALL pardiso_handle_delete(buff,error)
        IF (error .NE. 0) THEN
           WRITE(*,*) 'The following ERROR was detected: ', error
           STOP 1
        ENDIF

C.. Back substitution and iterative refinement
        phase     = 33  ! only factorization
        iparm(8)  = 3   ! max numbers of iterative refinement steps
        DO i = 1, n
            b(i) = (1.d0,1.d0)
        END DO
  
        CALL pardiso (pt_33, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1              idum, nrhs, iparm, msglvl, b, x, error)
        WRITE(*,*) 'Solve completed ... '
  
        WRITE(*,*) 'The solution of the system is '
        DO i = 1, n
            WRITE(*,*) ' x(',i,') = ', x(i)
        END DO

C.. Termination and release of memory
        phase     = -1           ! release internal memory
        CALL pardiso (pt_33, maxfct, mnum, mtype, phase, n, ddum, idum,
     1            idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)

        STOP 0
      END PROGRAM pardiso_unsym_store_restore



