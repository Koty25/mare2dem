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

!   Content : TR Solver F90 example
!
!********************************************************************************


module u_data
    type, public :: my_data
          integer a
          integer sum
    end type my_data
end module

!   nonlinear least square problem without boundary constraints
    include 'mkl_rci.f90'
program EXAMPLE_EX_NLSQP_BC_F90_X
    use MKL_RCI
    use MKL_RCI_TYPE
    use u_data
    implicit none
!   user's objective function
    external extended_powell 
!   n - number of function variables
!   m - dimension of function value
    integer n, m
    parameter (n = 4)
    parameter (m = 4)
!   precisions for stop-criteria (see manual for more details)
    real*8 eps(6)
!   solution vector. contains values x for f(x)
    real*8 x(n)
!   iter1 - maximum number of iterations
!   iter2 - maximum number of iterations of calculation of trial-step
    integer iter1, iter2
    parameter (iter1 = 1000)
    parameter (iter2 = 100)
!   initial step bound
    real*8 rs
!   reverse communication interface parameter
    integer RCI_Request
!   controls of rci cycle
    integer successful
!   function (f(x)) value vector
    real*8 fvec(m)
!   jacobi matrix
    real*8 fjac(m*n)
!   lower and upper bounds
    real*8 LW(n)
    real*8 UP(n)
!   number of iterations
    integer iter
!   number of stop-criterion
    integer st_cr
!   initial and final residuals
    real*8 r1, r2
!   TR solver handle
    type(HANDLE_TR) :: handle
!   cycle's counter
    integer i
!   results of input parameter checking
    integer info(6)

!   Additional users data
    type(my_data) :: m_data
    m_data%a = 1
    m_data%sum = 0
    rs = 0.0

!   set precisions for stop-criteria
    do i = 1, 6
        eps(i) = 0.00001
    end do
!   set the initial guess
    do i = 1, n/4
        x(4 * (i-1) + 1) = 3.0
        x(4 * (i-1) + 2) = -1.0
        x(4 * (i-1) + 3) = 0.0
        x(4 * (i-1) + 4) = 1.0
    end do
!   set initial values
    do i = 1, m
        fvec(i) = 0.0
    end do
    do i = 1, m * n
        fjac(i) = 0.0
    end do

!   set bounds
    do i = 1, n/4
        LW(4 * (i-1) + 1) = 0.1;
        LW(4 * (i-1) + 2) = -20.0;
        LW(4 * (i-1) + 3) = -1.0;
        LW(4 * (i-1) + 4) = -1.0;
        UP(4 * (i-1) + 1) = 100.0;
        UP(4 * (i-1) + 2) = 20.0;
        UP(4 * (i-1) + 3) = 1.0;
        UP(4 * (i-1) + 4) = 50.0;
    end do
!   initialize solver (allocate memory, set initial values)
!       handle       in/out: TR solver handle
!       n       in:     number of function variables
!       m       in:     dimension of function value
!       x       in:     solution vector. contains values x for f(x)
!       LW           in:             lower bound
!       UP           in:             upper bound
!       eps     in:     precisions for stop-criteria
!       iter1   in:     maximum number of iterations
!       iter2   in:     maximum number of iterations of calculation of trial-step
!       rs      in:     initial step bound
    if (dtrnlspbc_init (handle,n,m,x,LW,UP,eps,iter1,iter2,rs) /= TR_SUCCESS) then
!       if function does not complete successfully then print error message
        print*,'| error in dtrnlspbc_init'
!       Release internal Intel(R) Math Kernel Library (Intel(R) MKL) memory that might be used for computations.
!       NOTE: It is important to call the routine below to avoid memory leaks
!       unless you disable Intel MKL Memory Manager
        call MKL_FREE_BUFFERS
!       and exit
        stop 1
    end if
!       Checks the correctness of handle and arrays containing Jacobian matrix, 
!       objective function, lower and upper bounds, and stopping criteria.
    if (dtrnlspbc_check (handle, n, m, fjac, fvec, LW, UP, eps, info) /= TR_SUCCESS) then
!       if function does not complete successfully then print error message
        print*,'| error in dtrnlspbc_init'
!       Release internal Intel MKL memory that might be used for computations.
!       NOTE: It is important to call the routine below to avoid memory leaks
!       unless you disable Intel MKL Memory Manager
        call MKL_FREE_BUFFERS
!       and exit
        stop 1
    else
        if ( &
!           The handle is not valid. 
            info(1) /= 0 .or. &
!           The fjac array is not valid.
            info(2) /= 0 .or. &
!           The fvec array is not valid.
            info(3) /= 0 .or. &
!           The LW array is not valid.
            info(4) /= 0 .or. &
!           The UP array is not valid.
            info(5) /= 0 .or. &
!           The eps array is not valid.
            info(6) /= 0 ) then
            print*,'| input parameters for dtrnlspbc_solve are not valid'
!           Release internal Intel MKL memory that might be used for computations.
!           NOTE: It is important to call the routine below to avoid memory leaks
!           unless you disable Intel MKL Memory Manager
            call MKL_FREE_BUFFERS
!           and exit
            stop 1
        end if
    end if
!   set initial rci cycle variables
    RCI_Request = 0
    successful = 0
!   rci cycle
    do while (successful == 0)
!       call tr solver
!       handle               in/out: tr solver handle
!       fvec         in:     vector
!       fjac         in:     jacobi matrix
!       RCI_request in/out:  return number which denote next step for performing
        if (dtrnlspbc_solve (handle, fvec, fjac, RCI_Request) /= TR_SUCCESS) then
!           if function does not complete successfully then print error message
            print*, '| error in dtrnlspbc_solve'
!           Release internal Intel MKL memory that might be used for computations.
!           NOTE: It is important to call the routine below to avoid memory leaks
!           unless you disable Intel MKL Memory Manager
            call MKL_FREE_BUFFERS
!           and exit
            stop 1;
        end if
!       according with rci_request value we do next step
        if (-6 <= RCI_Request  .and. RCI_Request <= -1) then
!           exit rci cycle
            successful = 1
        end if
        if (RCI_Request == 1) then
!           recalculate function value
!               m            in:     dimension of function value
!               n            in:     number of function variables
!               x            in:     solution vector
!               fvec    out:    function value f(x)
            call extended_powell (m, n, x, fvec, m_data)
        end if
        if (RCI_Request == 2) then
!           compute jacobi matrix
!           extended_powell    in:     external objective function
!               n              in:     number of function variables
!               m              in:     dimension of function value
!               fjac           out:    jacobi matrix
!               x              in:     solution vector
!               jac_eps        in:     jacobi calculation precision !/
            if (djacobix (extended_powell,n,m,fjac,x,eps(1),%VAL(LOC(m_data))) &
                /= TR_SUCCESS) then
!               if function does not complete successfully then print error message
                print*, '| error in djacobix'
!               Release internal Intel MKL memory that might be used for computations.
!               NOTE: It is important to call the routine below to avoid memory leaks
!               unless you disable Intel MKL Memory Manager
                call MKL_FREE_BUFFERS
!               and exit
                stop 1;
            end if
        end if
    end do
!   get solution statuses
!   handle            in:        TR solver handle
!   iter              out:       number of iterations
!   st_cr             out:       number of stop criterion
!   r1                out:       initial residuals
!   r2                out:       final residuals
    if (dtrnlspbc_get (handle, iter, st_cr, r1, r2) /= TR_SUCCESS) then
!       if function does not complete successfully then print error message
        print*, '| error in dtrnlspbc_get'
!       Release internal Intel MKL memory that might be used for computations.
!       NOTE: It is important to call the routine below to avoid memory leaks
!       unless you disable Intel MKL Memory Manager
        call MKL_FREE_BUFFERS
!       and exit
        stop 1
    end if
    print*, 'Iterations : ',iter
    print*, 'Final residual : ',r2
    print*, 'Stop-criteria : ',st_cr
!   free handle memory
    if (dtrnlspbc_delete (handle) /= TR_SUCCESS) then
!       if function does not complete successfully then print error message
        print*, '| error in dtrnlspbc_delete'
!       Release internal Intel MKL memory that might be used for computations.
!       NOTE: It is important to call the routine below to avoid memory leaks
!       unless you disable Intel MKL Memory Manager
        call MKL_FREE_BUFFERS
!       and exit
        stop 1
    end if

!   Release internal Intel MKL memory that might be used for computations.
!   NOTE: It is important to call the routine below to avoid memory leaks
!   unless you disable Intel MKL Memory Manager
    call MKL_FREE_BUFFERS
!   if final residual less then required precision then print pass

    print*, 'User data ', m_data%sum

    if (r2 < 0.1) then
        print*, '|         dtrnlspbc Powell............PASS'
        stop 0
!   else print failed
    else
        print*, '|         dtrnlspbc Powell............FAILED'
        stop 1;
    end if
end program EXAMPLE_EX_NLSQP_BC_F90_X

!   nonlinear system equations without constraints
!   routine for extended Powell function calculation
!       m     in:     dimension of function value
!       n     in:     number of function variables
!       x     in:     vector for function calculating
!       f     out:    function value f(x)
!   user_data in: additional users data

subroutine extended_powell (m, n, x, f, user_data)
    use u_data
    implicit none
    integer m, n
    real*8  x(n), f(m)
    type(my_data) :: user_data
    integer i

    user_data%sum = user_data%sum + user_data%a

    do i = 1, n/4
        f(4*(i-1)+1) = x(4*(i-1)+1) + 10.0 * x(4*(i-1)+2)
        f(4*(i-1)+2) = 2.2360679774998 * (x(4*(i-1)+3) - x(4*(i-1)+4))
        f(4*(i-1)+3) = ( x(4*(i-1)+2) - 2.0 * x(4*(i-1)+3) )**2
        f(4*(i-1)+4) = 3.1622776601684 * (x(4*(i-1)+1) - x(4*(i-1)+4))**2
    end do
end subroutine extended_powell
