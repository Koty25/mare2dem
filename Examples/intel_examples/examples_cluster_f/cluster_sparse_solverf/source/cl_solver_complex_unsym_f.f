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

*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Cluster Sparse Solver
*             Fortran example complex, double precision, unsymmetric matrix
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the "CLUSTER_SPARSE_SOLVER" routine
C for solving complex unsymmetric linear systems.
C---------------------------------------------------------------------
      program cluster_sparse_solver_complex_unsym
      implicit none
      include 'mkl_cluster_sparse_solver.fi'
      include 'mpif.h'
C.. Internal solver memory pointer for 64-bit architectures
      TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE) pt(64)
C.. All other variables
      INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
      INTEGER*4 rank, mpi_stat
      INTEGER iparm(64)
      INTEGER ia(9)
      INTEGER ja(20)
      COMPLEX*16 a(20)
      COMPLEX*16 b(8)
      COMPLEX*16 bs(8)
      COMPLEX*16 x(8)
      INTEGER i, idum(1)
      COMPLEX*16 ddum(1)
      COMPLEX*16 res, res0
      external mkl_zcsrgemv
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
C..
C.. Initialize MPI.   
      call MPI_INIT(mpi_stat)

      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_stat)
C..
C.. Set up Cluster Sparse Solver control parameter
C..
      do i = 1, 64
         iparm(i) = 0
      enddo
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(8) = 2 ! numbers of iterative refinement steps
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(13) = 1 ! maximum weighted matching algorithm is switched-off
      iparm(40) = 0 ! Input: matrix/rhs/solution stored on master
      error = 0 ! initialize error flag
      msglvl = 1 ! print statistical information
      mtype = 13 ! symmetric, indefinite
C.. Initiliaze the internal solver memory pointer. This is only
C necessary for the FIRST call of the Cluster Sparse Solver.
      do i = 1, 64
         pt(i)%DUMMY = 0
      enddo
C.. Reordering and Symbolic Factorization, This step also allocates
C all memory that is necessary for the factorization
      phase = 11 ! only reordering and symbolic factorization
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) write(*,*) 'ERROR during symbolic factorization: ', error
         goto 999
      endif
      if (rank.eq.0) write(*,*) 'Reordering completed ... '

C.. Factorization.
      phase = 22 ! only factorization
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) write(*,*) 'ERROR during numerical factorization: ', error
         goto 999
      endif
      if (rank.EQ.0) write(*,*) 'Factorization completed ... '

C.. Back substitution and iterative refinement
      phase = 33 ! only solution
      if (rank.eq.0) then
        do i = 1, n
          b(i) = (1.d0,1.d0)
        enddo
      endif
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, b, x, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) write(*,*) 'The following ERROR was detected: ', error
         goto 999
      endif
      
      if (rank.eq.0) then
          write(*,*) 'Solve completed ... '
          write(*,*) 'The solution of the system is '
          do i = 1, n
              write(*,'("  x( ", I1, " ) = (", F19.16, ", ", F19.16,")")') i, x(i)
          enddo
          call mkl_zcsrgemv ('N', n, a, ia, ja, x, bs)
          res  = (0.d0,0.d0)
          res0 = (0.d0,0.d0)
          do i=1,n
              res = res + (bs(i)-b(i))*conjg((bs(i)-b(i)))
              res0 = res0 + b(i)*conjg(b(i))
          enddo
          print *, 'Relative residual = ', sqrt(abs(res))/sqrt(abs(res0))
      endif

C.. Termination and release of memory
      phase = -1 ! release internal memory
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error)
      if (error .NE. 0) then
         if (rank.EQ.0) write(*,*) 'The following ERROR was detected: ', error
         call MPI_FINALIZE(mpi_stat)
         stop 1
      endif

      if (rank .eq. 0) then
          if (sqrt(abs(res))/sqrt(abs(res0)).gt.1.d-10) then
             write(*,*) 'Error: residual is too high!'
             stop 1
          endif
          if (error.ne.0) then
              write(*,*) char(10), 'TEST FAILED'
          else
              write(*,*) char(10), 'TEST PASSED'
          endif
      endif

999   continue
      call MPI_FINALIZE(mpi_stat)
      if (error.ne.0.and.rank.eq.0) then
          stop 1
      endif
      end
