!===============================================================================
! Copyright 2015-2020 Intel Corporation.
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

! Content:
!       Example of using Intel(R) Math Kernel Library (Intel(R) MKL)
!       wrappers for FFTW3 MPI
!
!*******************************************************************************

program dp_c2c_1d
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'
  include 'fftw3-mpi.f03'

  ! Sizes of 3D transform
  integer(C_INTPTR_T), parameter :: N = 3424

  type(C_PTR) :: plan_fwd, plan_bwd, x_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: x(:)

  real(C_DOUBLE) :: r, l, loc_err_l2, loc_norm_l2, err_l2, norm_l2, rel_err
  real(C_DOUBLE), parameter :: eps = 1.d-14
  logical :: test_passed

  ! FFT local sizes
  integer(C_INTPTR_T) :: local_ni, local_i_start, local_no, local_o_start
  integer(C_INTPTR_T) :: i, j, local_size

  ! MPI stuff
  integer*4, parameter :: comm = MPI_COMM_WORLD
  integer*4 :: stat, irank, nrank

  ! MPI initialization
  call mpi_init(stat)
  call mpi_comm_size(comm, nrank, stat)
  call mpi_comm_rank(comm, irank, stat)

  if (0 == irank) then
      print '(" Number of MPI processes "I0)', nrank
      print '(" FFT data size: "I0)', N
  endif

  ! FFTW3 MPI initialization
  call fftw_mpi_init()

  ! allocating array
  local_size = fftw_mpi_local_size_1d(N, comm, FFTW_FORWARD, FFTW_ESTIMATE, &
      local_ni, local_i_start, local_no, local_o_start)
  x_ptr = fftw_alloc_complex(local_size)
  call c_f_pointer(x_ptr, x, [local_ni])

  print '(" FFT local data size on rank "I2": in "I4".."I4", out: "I4".."I4)', &
      irank, local_i_start, local_i_start + local_ni, &
      local_o_start, local_o_start + local_no

  plan_fwd = fftw_mpi_plan_dft_1d(N, x, x, comm, FFTW_FORWARD, FFTW_ESTIMATE)
  plan_bwd = fftw_mpi_plan_dft_1d(N, x, x, comm, FFTW_BACKWARD, FFTW_ESTIMATE)

  do i = 1, local_ni
      x(i) = foo(local_i_start + i, N)
  end do

  call fftw_mpi_execute_dft(plan_fwd, x, x)
  call fftw_mpi_execute_dft(plan_bwd, x, x)

  loc_err_l2 = 0
  loc_norm_l2 = 0
  do i = 1, local_ni
      l = foo(local_i_start + i, N)
      r = x(i)/N - l

      loc_err_l2 = loc_err_l2 + r**2
      loc_norm_l2 = loc_norm_l2 + l**2
  end do

  err_l2 = 0
  norm_l2 = 0
  call mpi_allreduce(loc_err_l2, err_l2, 1, MPI_DOUBLE, MPI_SUM, comm, stat)
  call mpi_allreduce(loc_norm_l2, norm_l2, 1, MPI_DOUBLE, MPI_SUM, comm, stat)
  err_l2 = sqrt(err_l2)
  norm_l2 = sqrt(norm_l2)
  rel_err = err_l2 / norm_l2

  call fftw_destroy_plan(plan_fwd)
  call fftw_destroy_plan(plan_bwd)
  call fftw_free(x_ptr)

  call fftw_mpi_cleanup()
  call mpi_finalize(stat)

  test_passed = rel_err < eps
  if (0 == irank) then
      print *, "=================================="
      print *, "computational relative error in L2 = ", rel_err
      if (test_passed) then
          print *," TEST PASSED"
      else
          print *," TEST FAILED"
      endif
  endif

  if (.not. test_passed) call exit(1)
  call exit(0)

contains

  pure complex(C_DOUBLE_COMPLEX) function foo(i, N)
      integer(C_INTPTR_T), intent(in) :: i, N

      foo = CMPLX( sin(1.1d0 * i / N), cos(1.2d0 * i / N) )
  end function foo

end program dp_c2c_1d
