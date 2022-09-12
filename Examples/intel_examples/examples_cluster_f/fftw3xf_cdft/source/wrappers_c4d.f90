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

program dp_c2c_4d
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'
  include 'fftw3-mpi.f03'

  ! Sizes of 3D transform
  integer(C_INTPTR_T) :: N(4) = [65, 128, 132, 16]

  type(C_PTR) :: plan_fwd, plan_bwd, x_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: x(:,:,:,:)

  real(C_DOUBLE) :: r, l, loc_err_l2, loc_norm_l2, err_l2, norm_l2, rel_err
  real(C_DOUBLE), parameter :: eps = 1.d-14
  logical :: test_passed

  ! FFT local sizes
  integer(C_INTPTR_T) :: i, j, k, m, local_size, local_n1, local_1_start

  ! MPI stuff
  integer*4, parameter :: comm = MPI_COMM_WORLD
  integer*4 :: stat, irank, nrank

  ! MPI initialization
  call mpi_init(stat)
  call mpi_comm_size(comm, nrank, stat)
  call mpi_comm_rank(comm, irank, stat)

  if (0 == irank) then
      print '(" Number of MPI processes "I0)', nrank
      print '(" FFT data size: "I0" x "I0" x "I0" x "I0)', &
          N(1), N(2), N(3), N(4)
      print *,"Data distribution across first dimension"
  endif

  ! FFTW3 MPI initialization
  call fftw_mpi_init()

  ! allocating array
  local_size = fftw_mpi_local_size(4, N, comm, local_n1, local_1_start)
  x_ptr = fftw_alloc_complex(local_size)
  call c_f_pointer(x_ptr, x, [N(4), N(3), N(2), local_n1])

  print '(" FFT local data size on rank "I2": "I3" x "I3" x "I3" x "I3)', &
      irank, local_n1, N(2), N(3), N(4)

  plan_fwd = fftw_mpi_plan_dft(4, N, x, x, comm, FFTW_FORWARD, FFTW_ESTIMATE)
  plan_bwd = fftw_mpi_plan_dft(4, N, x, x, comm, FFTW_BACKWARD, FFTW_ESTIMATE)

  ! Note reversed indexing: x(m,k,j,i) <- foo(i,j,k,m)
  forall (i=1:local_n1, j=1:N(2), k=1:N(3), m=1:N(4))
      x(m, k, j, i) = foo(local_1_start + i, j, k, m, N(1), N(2), N(3), N(4))
  end forall

  call fftw_mpi_execute_dft(plan_fwd, x, x)
  call fftw_mpi_execute_dft(plan_bwd, x, x)

  loc_err_l2 = 0
  loc_norm_l2 = 0
  do i = 1, local_n1
      do j = 1, N(2)
          do k = 1, N(3)
              do m = 1, N(4)
                  l = foo(local_1_start + i, j, k, m, N(1), N(2), N(3), N(4))
                  r = x(m, k, j, i)/(N(1) * N(2) * N(3) * N(4)) - l

                  loc_err_l2 = loc_err_l2 + r**2
                  loc_norm_l2 = loc_norm_l2 + l**2
              end do
          end do
      end do
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

  pure complex(C_DOUBLE_COMPLEX) function foo(i, j, k, m, N1, N2, N3, N4)
      integer(C_INTPTR_T), intent(in) :: i, j, k, m, N1, N2, N3, N4

      foo = CMPLX( sin(1.1d0 * i / N1 + 1.2d0 * j / N2), &
          cos(1.3d0 * k / N3 + 1.4d0 * m / N4) )
  end function foo

end program dp_c2c_4d
