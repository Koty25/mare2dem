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
!      FORTRAN OpenMP offload examples for ZGERU
!*******************************************************************************

include "mkl_omp_offload.f90"
include "common_blas.f90"

program zgeru_example
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
use common_blas  

integer :: m = 5, n = 3
integer :: lda, incx = 1, incy = 1, passed
complex*16 :: alpha = 1.0
complex*16,allocatable :: a(:,:), a_ref(:,:), x(:), y(:)

lda = m

allocate(a(lda,n))
allocate(a_ref(lda,n))
allocate(x(1 + (m - 1) * abs(incx)))
allocate(y(1 + (n - 1) * abs(incy)))

if (.not. allocated(a)) goto 998
if (.not. allocated(x)) then
   deallocate(a)
   goto 998
end if
if (.not. allocated(y)) then
   deallocate(a)
   deallocate(x)
   goto 998
end if
if (.not. allocated(a_ref)) then
   deallocate(a)
   deallocate(x)
   deallocate(y)
   goto 998
end if

! initialize matrices
call zinit_matrix('N', m, n, lda, a)
call zinit_vector(m, incx, x)
call zinit_vector(n, incy, y)
call zcopy_matrix('N', m, n, lda, a, a_ref)

! Calling zgeru on the CPU
call zgeru(m, n, alpha, x, incx, y, incy, a_ref, lda)

! Calling zgeru on the GPU
!$omp target data map(a,x,y)
!$omp target variant dispatch use_device_ptr(a,x,y)
call zgeru(m, n, alpha, x, incx, y, incy, a, lda)
!$omp end target variant dispatch
!$omp end target data

! Compare result of CPU and GPU implementation

passed = zcheck_matrix(m, n, lda, a, a_ref)

if (passed.ne.0) then
   deallocate(a)
   deallocate(x)
   deallocate(y)
   deallocate(a_ref)
   goto 999
end if

print *, "Matrix computed on GPU:"
call zprint_matrix(m, n, lda, a)

print *, "Matrix computed on CPU:"
call zprint_matrix(m, n, lda, a_ref)

deallocate(a)
deallocate(x)
deallocate(y)
deallocate(a_ref)

stop

998 print *, 'Error: cannot allocate matrices' 
999 stop 1
end program
