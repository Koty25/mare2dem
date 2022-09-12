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
! Example program to show the use of the "CLUSTER_SPARSE_SOLVER" routine and
! the export functionality which provides the factors L, U and permutations P
! and Q such that P * A * Q = L * U.
!---------------------------------------------------------------------
program cluster_sparse_solver_sym
use mkl_cluster_sparse_solver
use ISO_C_BINDING
implicit none
include 'mpif.h'
integer, parameter :: dp = kind(1.0D0)
!.. Internal solver memory pointer for 64-bit architectures
TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), allocatable  :: pt(:)
!.. All other variables
integer maxfct, mnum, mtype, phase, nrhs, error, msglvl
! Variables and pointers for the export functionality specifically
integer l_nrows, l_nnz, u_nrows, u_nnz
integer, allocatable :: p_ptr( : )
integer, allocatable :: q_ptr( : )
integer, allocatable :: l_ia( : )
integer, allocatable :: l_ja( : )
real(kind=dp), allocatable :: l_vals( : )
integer, allocatable :: u_ia( : )
integer, allocatable :: u_ja( : )
real(kind=dp), allocatable :: u_vals( : )

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
    iparm(13) = 0 ! maximum weighted matching algorithm is switched-off 
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
!.. Getting the local number of rows and local number of nonzeros in
!.. factor L.
!..
call cluster_sparse_solver_get_csr_size(pt, SPARSE_PTLUQT_L, l_nrows, l_nnz, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during getting csr size for factor L: ', error
    goto 999
endif
print *, 'Local number of rows in the L factor: ', l_nrows
print *, 'Local number of nnz  in the L factor: ', l_nnz
!..
!.. Allocating memory for the csr arrays of the L factor
!..
allocate( l_ia   ( l_nrows + 1 ) )
allocate( l_ja   ( l_nnz ) )
allocate( l_vals ( l_nnz ) )
!..
!.. Getting the local number of rows and local number of nonzeros
!.. in factor U.
!..
call cluster_sparse_solver_get_csr_size(pt, SPARSE_PTLUQT_U, u_nrows, u_nnz, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during getting csr size for factor U: ', error
    goto 999
endif
print *, 'Local number of rows in the U factor: ', u_nrows
print *, 'Local number of nnz  in the U factor: ', u_nnz
if (rank.eq.0) write(*,*) 'Sizes of factors L and U have been received ... '
!..
!.. Allocating memory for the csr arrays of the L factor
!..
allocate( u_ia   ( u_nrows + 1 ) )
allocate( u_ja   ( u_nnz ) )
allocate( u_vals ( u_nnz ) )
!..
!.. Allocating memory for the permutations P and Q (represented as
!.. flat arrays since they are diagonal)
!..
allocate( p_ptr ( l_nrows ) )
allocate( q_ptr ( l_nrows ) )
!..
!.. Saving allocated pointers internally in the cluster sparse solver
!..
call cluster_sparse_solver_set_csr_ptrs(pt, SPARSE_PTLUQT_L, l_ia, l_ja, l_vals, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during setting csr ptrs for L: ', error
    goto 999
endif
call cluster_sparse_solver_set_csr_ptrs(pt, SPARSE_PTLUQT_U, u_ia, u_ja, u_vals, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during setting csr ptrs for U: ', error
    goto 999
endif
call cluster_sparse_solver_set_ptr(pt, SPARSE_PTLUQT_P, p_ptr, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during setting csr ptrs for P: ', error
    goto 999
endif
call cluster_sparse_solver_set_ptr(pt, SPARSE_PTLUQT_Q, q_ptr, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during setting csr ptrs for Q: ', error
    goto 999
endif
if (rank.eq.0) write(*,*) 'Pointers for L, U, P, Q have been set ... '

call cluster_sparse_solver_export(pt, SPARSE_PTLUQT, MKL_COMM, error);
if (error.ne.0) then
    if (rank.eq.0) print *, 'ERROR during exporting data: ', error
    goto 999
endif
if (rank.eq.0) write(*,*) 'Factors L, U and permutations P, Q have been exported ... '
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
if ( allocated( pt ) )          deallocate( pt )
if ( allocated( iparm ) )       deallocate( iparm )

if ( allocated( p_ptr  ) )      deallocate( p_ptr )
if ( allocated( q_ptr  ) )      deallocate( q_ptr )
if ( allocated( l_ia   ) )      deallocate( l_ia )
if ( allocated( l_ja   ) )      deallocate( l_ja )
if ( allocated( l_vals ) )      deallocate( l_vals )
if ( allocated( u_ia   ) )      deallocate( u_ia )
if ( allocated( u_ja   ) )      deallocate( u_ja )
if ( allocated( u_vals ) )      deallocate( u_vals )


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
