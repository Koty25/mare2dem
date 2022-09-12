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

!   Content : Intel(R) Math Kernel Library (Intel(R) MKL) Cluster Sparse Solver
!             Fortran-90 example
!
!*******************************************************************************
!----------------------------------------------------------------------
! Example program to show the use of the "CLUSTER_SPARSE_SOLVER" routine 
! for solving symmetric linear systems
!---------------------------------------------------------------------
program cluster_sparse_solver_sym
use mkl_cluster_sparse_solver
implicit none
include 'mpif.h'
integer, parameter :: dp = kind(1.0D0)
!.. Internal solver memory pointer for 64-bit architectures
TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), allocatable  :: pt(:)
!.. All other variables
integer maxfct, mnum, mtype, phase, nrhs, error, msglvl
integer i, idum(1), n, nnz
integer*4 mpi_stat, rank
real(kind=dp)  ddum(1), norm_b, res

integer, allocatable :: iparm( : )
integer, allocatable :: ia( : )
integer, allocatable :: ja( : )
real(kind=dp), allocatable :: a( : )
real(kind=dp), allocatable :: b( : )
real(kind=dp), allocatable :: bs( : )
real(kind=dp), allocatable :: x( : )
real(kind=dp), allocatable :: r( : )

integer(4) MKL_COMM

data n /8/, nnz /18/, nrhs /1/, maxfct /1/, mnum /1/

!..
!.. Initialize task data on master MPI process
!..
MKL_COMM=MPI_COMM_WORLD
call mpi_init(mpi_stat)
call mpi_comm_rank(MKL_COMM, rank, mpi_stat)
!..
!.. Set task sizes
!..
allocate( pt ( 64 ) )
allocate( iparm ( 64 ) )
!..
!.. Set up Cluster Sparse Solver control parameters
!..
do i = 1, 64
  iparm(i) = 0
end do
!..
!.. Initiliaze internal Cluster Sparse Solver memory pointers.
!..
do i = 1, 64
  pt(i)%dummy = 0
end do
!..
!.. Initiliaze task data on master MPI process.
!..
if (rank.eq.0) then
    allocate( ia ( n + 1 ) )
    allocate( ja ( nnz ) )
    allocate( a ( nnz ) )
    allocate( b ( n ) )
    allocate( bs( n ) )
    allocate( r ( n ) )
    allocate( x ( n ) )
!
!.. Set the problem to be solved.
!
ia = (/1, 5, 8, 10, 12, 15, 17, 18, 19/)
ja = (/1,    3,       6,  7,    &
&         2, 3,    5,           &
&            3,              8, &
&               4,       7,     &
&                  5, 6, 7,     &
&                     6,     8, &
&                        7,     &
&                            8/)
a  = (/7.d0,     1.d0,          2.d0,7.d0,      &
&          -4.d0,8.d0,     2.d0,                &
&                1.d0,                     5.d0,&
&                     7.d0,     9.d0,           &
&                          5.d0,1.d0,5.d0,      &
&                              -1.d0,      5.d0,&
&                                    11.d0,     &
&                                          5.d0/)

!
!.. Fill in matrix data arrays
!
    iparm(1) = 1 ! no solver default
    iparm(2) = 3 ! fill-in reordering from METIS
    iparm(6) = 0 ! =0 solution on the first n compoments of x
    iparm(8) = 2 ! numbers of iterative refinement steps
    iparm(10) = 13 ! perturbe the pivot elements with 1E-13
    iparm(11) = 0 ! use nonsymmetric permutation and scaling MPS
    iparm(13) = 1 ! maximum weighted matching algorithm is switched-off 
    iparm(40) = 0 ! Input: matrix/rhs/solution stored on master
    error = 0  ! initialize error flag
    msglvl = 1 ! print statistical information
    mtype = -2  ! symmetric, indefinite
endif
!..
!.. Reordering and symbolic factorization
!..
phase = 11
call cluster_sparse_solver ( &
    pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
    idum, nrhs, iparm, msglvl, b, x, MKL_COMM, error )
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during symbolic factorization: ', error
    goto 999
endif
if (rank.eq.0) write(*,*) 'Reordering completed ... '
!..
!.. Numeric factorization
!..
phase = 22
call cluster_sparse_solver ( &
    pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
    idum, nrhs, iparm, msglvl, b, x, MKL_COMM, error )
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during numerical factorization: ', error
    goto 999
endif
if (rank.eq.0) write(*,*) 'Factorization completed ... '
!..
!.. Forward & backward substitution
!..

if (rank.eq.0) then
do i=1,n
  b(i)  = 1.
  bs(i) = 0.
enddo
endif

phase = 33
call cluster_sparse_solver ( &
    pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
    idum, nrhs, iparm, msglvl, b, x, MKL_COMM, error )
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during during solution: ', error
    goto 999
endif
!..
!.. Computing relative residual
!..
if (rank.eq.0) then
    call mkl_dcsrsymv('U', n, a, ia, ja, x, bs)
    res  = 0.d0
    norm_b = 0.d0
    do i = 1, n
      res  = res  + (bs(i)-b(i))**2
      norm_b = norm_b + b(i)**2
    enddo
    res = sqrt(res/norm_b)
    print *, 'Relative residual = ', res
endif
!..
!.. Termination and release of memory
!..
phase = -1
call cluster_sparse_solver ( &
     pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
     idum, nrhs, iparm, msglvl, ddum, ddum, MKL_COMM, error )
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during release memory: ', error
endif

999   continue
if (rank.eq.0) then
    if ( allocated( ia ) )      deallocate( ia )
    if ( allocated( ja ) )      deallocate( ja )
    if ( allocated( a ) )       deallocate( a )
    if ( allocated( b ) )       deallocate( b )
    if ( allocated( bs ) )      deallocate( bs )
    if ( allocated( r ) )       deallocate( r )
    if ( allocated( x ) )       deallocate( x )
endif
if ( allocated( pt ) )      deallocate( pt )
if ( allocated( iparm ) )   deallocate( iparm )

if (rank.eq.0) then
    if (res.gt.1.d-10) then
        print *, 'Error: residual is too high!'
        error = 1
    endif
    if (error.ne.0) then
        print *, char(10), 'TEST FAILED'
    else
        print *, char(10), 'TEST PASSED'
    endif
endif
call mpi_finalize(mpi_stat)
if (error.ne.0.and.rank.eq.0) then
    stop 1
endif
end
