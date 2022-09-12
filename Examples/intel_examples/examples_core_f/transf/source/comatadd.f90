!===============================================================================
! Copyright 2010-2020 Intel Corporation.
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
!      C O M A T A D D  Example Program Text
!*******************************************************************************
    program  COMATADD_MAIN
    include 'mkl_trans.fi' 

    integer :: i, j
    complex :: alpha = (1.0,0.0), beta = (1.0,0.0)
!   source matrix A
    complex,dimension(3,5),parameter :: A=RESHAPE(&
                                         &(/(1.0,1.0), (5.0,5.0), ( 9.0, 9.0),& 
                                           &(2.0,2.0), (6.0,6.0), (10.0,10.0),& 
                                           &(3.0,3.0), (7.0,7.0), (11.0,11.0),&
                                           &(4.0,4.0), (8.0,8.0), (12.0,12.0),&
                                           &(0.0,0.0), (0.0,0.0), ( 0.0, 0.0)/),SHAPE(A))
!   source matrix B
    complex,dimension(4,2),parameter :: B=RESHAPE(&
                                         &(/(1.0,0.0), (4.0,0.0), (7.0,0.0), (10.0,0.0),&  
                                           &(2.0,0.0), (5.0,0.0), (8.0,0.0), (11.0,0.0)/),SHAPE(B))
    complex,dimension(3,2) :: C ! destination matrix 1 
    complex,dimension(2,3) :: D ! destination matrix 2

!   Executable Statements
    print*, ' Example of using mkl_comatadd transposition '
    print*, ' INPUT DATA: '
    print*, ' Source matrix A: '    
!   print source matrix           
    do i= 1, 3
          print 102, ( A(i, j), j=1, 5)
    enddo
    print*, ' Source matrix B: '    
!   print source matrix           
    do i= 1, 4
          print 102, ( B(i, j), j=1, 2)
    enddo
! 
!   Addition of transposed submatrix(3, 3) A and 
!   unchanged submatrix(3, 3) B     
!
    call MKL_COMATADD( 'C',    & ! column-major ordering
                     & 'C',    & ! A is transposed
                     & 'R',    & ! B is unchanged
                     &  3,     & ! rows of destination matrix
                     &  2,     & ! cols of destination matrix
                     &  alpha, & ! alpha
                     &  A,     & ! source matrix
                     &  3,     & ! lda
                     &  beta,  & ! beta
                     &  B,     & ! source matrix
                     &  4,     & ! ldb
                     &  C,     & ! destination matrix
                     &  3 )      ! ldc
!   New matrix: c = ( 
!     (2.,-1.),   (7., -5.),  
!     (6.,-2.),  (11., -6.),
!    (10.,-3.),  (15., -7.)
!   )  
    print*, 'OUTPUT DATA: '
    print*, 'Destination matrix: '      
    print*, 'Addition of transposed submatrix(3,3) of A &
            & and submatrix(3,3) of B:'
!   print destination matrix      
    do i= 1, 3
          print 102,  ( C(i, j), j = 1, 2)
    enddo
! 
!   Addition of transposed submatrix(3, 3) A and 
!   transposed submatrix(3, 3) B     
!
    call MKL_COMATADD( 'C',    & ! column-major ordering
                     & 'T',    & ! A is transposed
                     & 'T',    & ! B is transposed
                     &  2,     & ! rows of destination matrix
                     &  3,     & ! cols of destination matrix
                     &  alpha, & ! alpha
                     &  A,     & ! source matrix
                     &  3,     & ! lda
                     &  beta,  & ! beta
                     &  B,     & ! source matrix
                     &  4,     & ! ldb
                     &  D,     & ! destination matrix
                     &  2 )      ! ldd
!   New matrix: D = { 
!     (2., 1.),  (9.,  5.), (16.,  9.),
!     (4., 2.), (11.,  6.), (18., 10.)
!   )         
    print*, 'Destination matrix: '      
    print*, 'Addition of transposed submatrices(3,3) of A and B: '
!   print destination matrix      
    do i= 1, 2
          print 102, (D(i, j), j = 1,3)
    enddo
     
102 format(9x,10(f6.2,',',f6.2,'i',3x))
    
    stop
    end