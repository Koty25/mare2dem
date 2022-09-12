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
*             Fortran example real, double precision, symmetric matrix
*
********************************************************************************
C----------------------------------------------------------------------
C Example program to show the use of the "CLUSTER_SPARSE_SOLVER" routine
C for symmetric linear systems with distributed CSR format of input data
C---------------------------------------------------------------------

      program cluster_sparse_solver_sym_distr_f
      implicit none
      include 'mkl_cluster_sparse_solver.fi'
      include 'mpif.h'
C.. Internal solver memory pointer for 64-bit architectures
      TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE) pt(64)
C.. All other variables
      INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
      INTEGER*4 rank, mpi_stat
      INTEGER  iparm(64)
      INTEGER, allocatable, dimension(:) :: ia
      INTEGER, allocatable, dimension(:) :: ja
      REAL*8, allocatable, dimension(:) :: a
      REAL*8, allocatable, dimension(:) :: b
      REAL*8, allocatable, dimension(:) :: x
      INTEGER i, idum(1)
      REAL*8 ddum(1)
C.. Fill all arrays containing matrix data.
      DATA n /5/, nrhs /1/, maxfct /1/, mnum /1/
C..
C.. Initialize MPI.
      call MPI_INIT(mpi_stat)
      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_stat)
C..
C.. Set up Cluster Sparse Solver control parameter
C..
        
      do i = 1, 64
         iparm(i) = 0
      end do
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(8) = 2 ! numbers of iterative refinement steps
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
      iparm(13) = 1 ! maximum weighted matching algorithm is switched-off
      iparm(27) = 1 ! Input: Matrix checker of initial matrix      
      iparm(40) = 1 ! Input: matrix/rhs distributed between processes, solution stored on master
      error = 0 ! initialize error flag
      msglvl = 1 ! print statistical information
      mtype = -2 ! symmetric, indefinite
C..
C.. Initialize matrix and rhs components on each process. 
C..      
      if (rank.eq.0) then
          ALLOCATE(ia(4),ja(6),a(6),b(3),x(5), stat=error)
          if (error.ne.0) then
            write(*,*) 'The following ERROR on allocation was detected on 0 rank: ', error
            write(*,*) char(10), 'TEST FAILED'
            stop 1
          endif

          ia(1) = 1
          ia(2) = 4
          ia(3) = 5
          ia(4) = 7

          ja(1) = 1
          ja(2) = 2
          ja(3) = 4
          ja(4) = 2
          ja(5) = 3
          ja(6) = 5

          a(1) = 6
          a(2) = -1
          a(3) = -3
          a(4) = 5
          a(5) = 3
          a(6) = 2

          b(1) = 1
          b(2) = 1
          b(3) = 0.25

          iparm(41) = 1
          iparm(42) = 3
      else if (rank.eq.1) then
          ALLOCATE(ia(4),ja(5),a(5),b(3), stat=error)
          if (error.ne.0) then
              write(*,*) 'The following ERROR on allocation was detected on 1 rank: ', error
              write(*,*) char(10), 'TEST FAILED'
              stop 2
          endif
          ia(1) = 1
          ia(2) = 4
          ia(3) = 5
          ia(4) = 6

          ja(1) = 3
          ja(2) = 4
          ja(3) = 5
          ja(4) = 4
          ja(5) = 5

          a(1) = 8
          a(2) = 5
          a(3) = 2
          a(4) = 10
          a(5) = 5

          b(1) = 0.75
          b(2) = 1
          b(3) = 1

          iparm(41) = 3
          iparm(42) = 5
      else
          iparm(41) = 2
          iparm(42) = 1
      endif
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
         if (rank.eq.0) then
             write(*,*) 'ERROR during symbolic factorization: ', error
         endif
         goto 999
      endif
      if (rank.eq.0) write(*,*) 'Reordering completed ... '

C.. Factorization.
      phase = 22 ! only factorization
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) then
             write(*,*) 'ERROR during symbolic factorization: ', error
         endif
         goto 999
      endif
      if (rank.eq.0) write(*,*) 'Factorization completed ... '

C.. Back substitution and iterative refinement
      phase = 33 ! only solution
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,
     1 idum, nrhs, iparm, msglvl, b, x, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) then
             write(*,*) 'ERROR during solution: ', error
         endif
         goto 999
      endif

      if (rank.eq.0) then
          write(*,*) 'Solve completed ... '
          write(*,*) 'The solution of the system is '
          do i = 1, n
              write(*,'("  x( ", I1, " ) = ", F19.16)') i, x(i)
          enddo
      endif

C.. Termination and release of memory
      phase = -1 ! release internal memory
      call cluster_sparse_solver (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,
     1 idum, nrhs, iparm, msglvl, ddum, ddum, MPI_COMM_WORLD, error)
      if (error.ne.0) then
         if (rank.eq.0) then 
             write(*,*) 'ERROR during release memory: ', error
         endif
      endif

999   continue
      if (rank.le.1) then
        deallocate(ia,ja,a,b)
      endif
      if (rank.eq.0) then
          deallocate(x)
          if (error.ne.0) then
              write(*,*) char(10), 'TEST FAILED'
          else
              write(*,*) char(10), 'TEST PASSED'
          endif
      endif
      CALL MPI_FINALIZE(mpi_stat)
      END
