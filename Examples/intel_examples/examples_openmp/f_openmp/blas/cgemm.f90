!===============================================================================
! Copyright 2020 Intel Corporation.
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
!      Intel(R) oneAPI Math Kernel Library (oneMKL)
!      FORTRAN OpenMP offload examples for CGEMM
!*******************************************************************************

include "mkl_omp_offload.f90"
include "common_blas.f90"

program cgemm_example
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
use common_blas  

character*1 :: transa = 'N', transb = 'N'
integer :: i, j, m = 5, n = 3, k = 10
integer :: lda, ldb, ldc, passed
complex :: alpha = 1.0, beta = 1.0
complex,allocatable :: a(:,:), b(:,:), c(:,:), c_ref(:,:)

lda = m
ldb = k
ldc = m

allocate(a(lda,k))
allocate(b(ldb,n))
allocate(c(ldc,n))
allocate(c_ref(ldc,n))

if (.not. allocated(a)) goto 998
if (.not. allocated(b)) then
   deallocate(a)
   goto 998
end if
if (.not. allocated(c)) then
   deallocate(a)
   deallocate(b)
   goto 998
end if
if (.not. allocated(c_ref)) then
   deallocate(a)
   deallocate(b)
   deallocate(c)
   goto 998
end if

! initialize matrices
call cinit_matrix(transa, m, k, lda, a)
call cinit_matrix(transb, k, n, ldb, b)
call cinit_matrix('N', m, n, ldc, c)
call ccopy_matrix('N', m, n, ldc, c, c_ref)

! Calling cgemm on the CPU
call cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c_ref, ldc)

! Calling cgemm on the GPU
!$omp target data map(a,b,c)
!$omp target variant dispatch use_device_ptr(a,b,c)
call cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
!$omp end target variant dispatch
!$omp end target data

! Compare result of CPU and GPU implementation

passed = ccheck_matrix(m, n, ldc, c, c_ref)

if (passed.ne.0) then
   deallocate(a)
   deallocate(b)
   deallocate(c)
   deallocate(c_ref)
   goto 999
end if

print *, "Matrix computed on GPU:"
call cprint_matrix(m, n, ldc, c)

print *, "Matrix computed on CPU:"
call cprint_matrix(m, n, ldc, c_ref)

deallocate(a)
deallocate(b)
deallocate(c)
deallocate(c_ref)

stop

998 print *, 'Error: cannot allocate matrices' 
999 stop 1
end program
