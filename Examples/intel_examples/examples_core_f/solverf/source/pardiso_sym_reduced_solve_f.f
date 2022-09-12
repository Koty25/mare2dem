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

*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) PARDISO Fortran-77
*             example
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the partial solve for sparse right-hand 
c sides and sparse solution for symmetric linear systems.
c The feature can used either a few components of the solution vector are needed 
c or the user wants to reduce computation cost at solver step.  
c---------------------------------------------------------------------
      PROGRAM pardiso_sym
        IMPLICIT NONE
        include 'mkl_pardiso.fi'
C.. Internal solver memory pointer for 64-bit architectures
C.. INTEGER*8 pt(64)
C.. Internal solver memory pointer for 32-bit architectures
C.. INTEGER*4 pt(64)
C.. This is OK in both cases
        TYPE(MKL_PARDISO_HANDLE) pt(64)
C.. All other variables
        INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        INTEGER iparm(64), perm(8)
        INTEGER ia(9)
        INTEGER ja(18)
        REAL*8 a(18)
        REAL*8 b(8)
        REAL*8 x(8), x_reduced(8)
        INTEGER i, j, idum(1)
        REAL*8 waltime1, waltime2, ddum(1), eps
C.. Fill all arrays containing matrix data.
        DATA n /8/, nrhs /1/, maxfct /1/, mnum /1/
        DATA ia /1,5,8,10,12,15,17,18,19/
        DATA ja
     1 /1,  3,    6,7,
     2    2,3,  5,
     3      3,        8,
     4        4,    7,
     5          5,6,7,
     6            6,  8,
     7              7,
     8                8/
        DATA a
     1 /7.d0,       1.d0,             2.d0, 7.d0,
     2       -4.d0, 8.d0,       2.d0,
     3              1.d0,                         5.d0,
     4                    7.d0,             9.d0,
     5                          5.d0, 1.d0, 5.d0,
     6                               -1.d0,       5.d0,
     7                                     11.d0,
     8                                            5.d0/
C..
C.. Set up PARDISO control parameters
C..
        DO i = 1, 64
            iparm(i) = 0
        END DO
        iparm(1) = 1 ! no solver default
        iparm(2) = 2 ! fill-in reordering from METIS
        iparm(3) = 0 ! not used in Intel MKL PARDISO
        iparm(4) = 0 ! no iterative-direct algorithm
        iparm(5) = 0 ! no user fill-in reducing permutation
        iparm(6) = 0 ! =0 solution on the first n components of x
        iparm(7) = 0 ! not in use
        iparm(8) = 0 ! numbers of iterative refinement steps
        iparm(9) = 0 ! not in use
        iparm(10) = 13 ! perturb the pivot elements with 1E-13
        iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
        iparm(12) = 0 ! not in use
        iparm(13) = 0 ! maximum weighted matching algorithm is switched-off 
c                     (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
        iparm(14) = 0 ! Output: number of perturbed pivots
        iparm(15) = 0 ! not in use
        iparm(16) = 0 ! not in use
        iparm(17) = 0 ! not in use
        iparm(18) = -1 ! Output: number of nonzeros in the factor LU
        iparm(19) = -1 ! Output: Mflops for LU factorization
        iparm(20) = 0 ! Output: Numbers of CG Iterations
        error = 0 ! initialize error flag
        msglvl = 0 ! no statistical information
        mtype = -2 ! symmetric, indefinite
C.. Initialize the internal solver memory pointer. This is only
C necessary for the FIRST call of the PARDISO solver.
        DO i = 1, 64
            pt(i)%DUMMY = 0
        END DO
C.. Initialize the right hand side. 
        DO i = 1, n
            b(i) = 0.0D0
        END DO
        b(1) = 1.0D0
        b(8) = 1.0D0
        phase = 13 ! get the solution of the full system
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     &  idum, nrhs, iparm, msglvl, b, x, error)

        WRITE(*,*) ' Regular solve completed ... '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
        WRITE(*,*) 'Solve completed ... '
        WRITE(*,*) 'The solution of the system is '
        DO i = 1, n
            WRITE(*,*) ' x(',i,') = ', x(i)
        END DO
C.. Termination and release of memory
        phase = -1 ! release internal memory
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, 
     &  idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
c
c     Compute the first and last components of the solution only
c     To use this option the user has to define  the input vector PERM 
c     so that PERM(i)=1 means that the i-th component in the right-hand 
c     side is nonzero. 
c     In this case, PERM(i)=1 also means that  the i-th component in  
c     the solution vector should be computed. 
        iparm(31) = 1
        DO i = 2, n-1
            perm(i) = 0
        END DO
        perm(1) = 1
        perm(n) = 1
        phase = 13 ! compute the first and last components if the solution
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     &  perm, nrhs, iparm, msglvl, b, x_reduced, error)
        WRITE(*,*) ' Finding a few components of the colution completed'
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
c    Check accuracy of the first and last components of the solution vector
        eps = 0.0D0
        DO i = n-1, n
            j = perm(i)
            eps = eps+(x(j) - x_reduced(j)) * (x(j) - x_reduced(j))
            WRITE(*,*) ' x(',j,') = ', x(j), x_reduced(j)
        END DO
        eps = dsqrt(eps)
        WRITE(*,*) ' Accuracy of the first and last components', eps
c    Setting IPARM(31)=2 allows reducing computational cost   
c    at the solver step. All components of the solution vector are computed.
c    To use IPARM(31)=2, define the components of the permutation vector   
c    PERM so that PERM(i)=1 means that the i-th component in the  
c    right-hand side is nonzero. Please note that if IPERM(i) is not   
c    equal to 1,  the i-th component of the right hand side must 
c    be set to zero explicitly.
        iparm(31) = 2
        DO i = 2, n-1
            b(i) = 0.0D0
        END DO

        phase = 33 ! find all components of the solution vector 
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     &  perm, nrhs, iparm, msglvl, b, x_reduced, error) 
        WRITE(*,*) 'Reduced forward solver step is used for computing '
        WRITE(*,*) 'all components of the solution vector '
        IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
        END IF
c    Check accuracy for the reduced forward substitution 
        eps = 0.0D0
        DO i = 1, n
            eps = eps+(x(i) - x_reduced(i)) * (x(i) - x_reduced(i))
            WRITE(*,*) ' x(',i,') = ', x(i), x_reduced(i)
        END DO
        eps = dsqrt(eps)
        WRITE(*,*) 'Accuracy of the solution vector found with the help'
        WRITE(*,*) 'of reduced forward step is ', eps

C.. Termination and release of memory
        phase = -1 ! release internal memory
        CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, 
     &  idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
      END PROGRAM pardiso_sym
