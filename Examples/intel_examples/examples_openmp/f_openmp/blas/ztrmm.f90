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
!      FORTRAN OpenMP offload examples for ZTRMM
!*******************************************************************************

include "mkl_omp_offload.f90"
include "common_blas.f90"

program ztrmm_example
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
use common_blas  

character*1 :: side = 'L', uplo = 'U', trans = 'N', diag = 'U'
integer :: i, j, m = 5, n = 3
integer :: lda, ldb, passed
complex*16 :: alpha = 1.0
complex*16,allocatable :: a(:,:), b(:,:), b_ref(:,:)

lda = m
ldb = m

allocate(a(lda,m))
allocate(b(ldb,n))
allocate(b_ref(ldb,n))

if (.not. allocated(a)) goto 998
if (.not. allocated(b)) then
   deallocate(a)
   goto 998
end if
if (.not. allocated(b_ref)) then
   deallocate(a)
   deallocate(b)
   goto 998
end if

! initialize matrices
call zinit_matrix(trans, m, m, lda, a)
call zinit_matrix('N', m, n, ldb, b)
call zcopy_matrix('N', m, n, ldb, b, b_ref)

! Calling ztrmm on the CPU
call ztrmm(side, uplo, trans, diag, m, n, alpha, a, lda, b_ref, ldb)

! Calling ztrmm on the GPU
!$omp target data map(a,b)
!$omp target variant dispatch use_device_ptr(a,b)
call ztrmm(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
!$omp end target variant dispatch
!$omp end target data

! Compare result of CPU and GPU implementation

passed = zcheck_matrix(m, n, ldb, b, b_ref)

if (passed.ne.0) then
   deallocate(a)
   deallocate(b)
   deallocate(b_ref)
   goto 999
end if

print *, "Matrix computed on GPU:"
call zprint_matrix(m, n, ldb, b)

print *, "Matrix computed on CPU:"
call zprint_matrix(m, n, ldb, b_ref)

deallocate(a)
deallocate(b)
deallocate(b_ref)

stop

998 print *, 'Error: cannot allocate matrices' 
999 stop 1
end program
