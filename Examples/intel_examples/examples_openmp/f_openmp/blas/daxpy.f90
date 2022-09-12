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
!      FORTRAN OpenMP offload examples for DAXPY
!*******************************************************************************

include "mkl_omp_offload.f90"
include "common_blas.f90"

program daxpy_example
#if defined(MKL_ILP64)
use onemkl_blas_omp_offload_ilp64
#else
use onemkl_blas_omp_offload_lp64
#endif
use common_blas  

integer :: n = 5
integer :: incx = 1, incy = 1, passed
double precision :: alpha = 1.0
double precision,allocatable :: x(:), y(:), y_ref(:)

allocate(x(1 + (n - 1) * abs(incx)))
allocate(y(1 + (n - 1) * abs(incy)))
allocate(y_ref(1 + (n - 1) * abs(incy)))

if (.not. allocated(x)) goto 998
if (.not. allocated(y)) then
   deallocate(x)
   goto 998
end if
if (.not. allocated(y_ref)) then
   deallocate(y)
   deallocate(x)
   goto 998
end if

! initialize matrices
call dinit_vector(n, incx, x)
call dinit_vector(n, incy, y)
call dcopy_vector(n, incy, y, y_ref)

! Calling daxpy on the CPU
call daxpy(n, alpha, x, incx, y_ref, incy)

! Calling daxpy on the GPU
!$omp target data map(x,y)
!$omp target variant dispatch use_device_ptr(x,y)
call daxpy(n, alpha, x, incx, y, incy)
!$omp end target variant dispatch
!$omp end target data

! Compare result of CPU and GPU implementation

passed = dcheck_vector(n, incy, y, y_ref)

if (passed.ne.0) then
   deallocate(x)
   deallocate(y)
   deallocate(y_ref)
   goto 999
end if

print *, "Vector computed on GPU:"
write( * , '(10G14.6)') y(:)

print *, "Vector computed on CPU:"
write( * , '(10G14.6)') y_ref(:)

deallocate(x)
deallocate(y)
deallocate(y_ref)

stop

998 print *, 'Error: cannot allocate matrices' 
999 stop 1
end program
