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
!      FORTRAN OpenMP offload examples for ZDOTC
!*******************************************************************************

include "mkl_omp_offload.f90"
include "common_blas.f90"

program zdotc_example
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
use common_blas  

integer :: n = 5, i
integer :: incx = 1, incy = 1, passed
complex*16 :: dot_cpu = 0, dot_gpu = 0
complex*16,allocatable :: x(:), y(:)
double precision :: diff

allocate(x(1 + (n - 1) * abs(incx)))
allocate(y(1 + (n - 1) * abs(incy)))

if (.not. allocated(x)) goto 998
if (.not. allocated(y)) then
   deallocate(x)
   goto 998
end if

! initialize matrices
call zinit_vector(n, incx, x)
call zinit_vector(n, incy, y)

! Calling zdotc on the CPU
dot_cpu = zdotc(n, x, incx, y, incy)

! Calling zdotc on the GPU
!$omp target data map(x,y)
!$omp target variant dispatch use_device_ptr(x,y)
dot_gpu = zdotc(n, x, incx, y, incy)
!$omp end target variant dispatch
!$omp end target data

! Compare result of CPU and GPU implementation
passed = zcheck_scalar(dot_gpu, dot_cpu)

if (passed.ne.0) then
   deallocate(x)
   deallocate(y)
   goto 999
end if

print *, "Result computed on GPU:", dot_gpu

print *, "Result computed on CPU:", dot_cpu

deallocate(x)
deallocate(y)

stop

998 print *, 'Error: cannot allocate matrices' 
999 stop 1
end program
