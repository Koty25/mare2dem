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
!      D O M A T C O P Y 2  Example Program Text
!*******************************************************************************
    program  DOMATCOPY2_MAIN
    include 'mkl_trans.fi' 

    integer :: i, j
!   source matrix 
    double precision,dimension(3,5),parameter :: A = RESHAPE(&
                                                    &(/1, 4, 7,&  
                                                      &0, 0, 0,&
                                                      &2, 5, 8,& 
                                                      &0, 0, 0,&
                                                      &3, 6, 9/),SHAPE(A))

    double precision, dimension(4,2) :: B   ! destination matrix 
    double precision :: alpha = 1.0

!   Executable Statements
    print*, ' Example of using mkl_domatcopy2 transposition '
    print*, ' INPUT DATA: '
    print*, ' Source matrix: '    
!   print source matrix           
    do i= 1, 3
          print 102, ( A(i, j), j=1, 5)
    enddo
! 
!   Submatrix(2, 4) is transposed     
!
    call MKL_DOMATCOPY2('C',    & ! column-major ordering
                       &'T',    & ! A is transposed
                       & 2,     & ! rows of destination matrix
                       & 4,     & ! cols of destination matrix
                       & alpha, & ! alpha
                       & A,     & ! source matrix
                       & 3,     & ! lda
                       & 2,     & ! stridea
                       & B,     & ! destination matrix
                       & 4,     & ! ldb
                       & 1      ) ! strideb
!   New matrix: B =  
!   ( 1.0     7.0
!     0.0     0.0
!     2.0     8.0
!     0.0     0.0
!   )          
    print*, 'OUTPUT DATA: '
    print*, 'Destination matrix: '      
    print*, 'Submatrix(2, 4) is transposed'
!   print destination matrix      
    do i= 1, 4
          print 102,  ( B(i, j), j = 1, 2)
    enddo
    
102 format(9x,10(f8.3,2x))
    
    stop
    end
