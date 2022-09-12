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
!      Z I M A T C O P Y  Example Program Text
!*******************************************************************************
    program  ZIMATCOPY_MAIN
    include 'mkl_trans.fi' 

    integer :: i, j
!   source matrix 
    double complex, dimension(3,4) :: A =  RESHAPE(&
                                                     &(/1, 5, 9, &
                                                       &2, 6, 10,&
                                                       &3, 7, 11,&
                                                       &4, 8, 12/),SHAPE(A))
    double complex :: alpha = 1
!   Executable Statements
    print*, ' Example of using mkl_zimatcopy transposition '
    print*, ' INPUT DATA: '
    print*, ' Source matrix: '    
!   print source matrix           
    do i= 1, 3
          print 102, ( A(i, j), j=1, 4)
    enddo
! 
!   Submatrix(3, 2) is transposed 
!
    call MKL_ZIMATCOPY('C',    & ! column-major ordering
                      &'T',    & ! A is transposed
                      & 3,     & ! rows of destination matrix
                      & 2,     & ! cols of destination matrix
                      & alpha, & ! alpha
                      & A,     & ! source matrix
                      & 3,     & ! src_lda
                      & 3 )      ! dst_lda
!   New matrix: A =  
!  ( (1.0,0), (5.0,0), (9.0,0), (4.0,0),
!    (2.0,0), (6.0,0), (10.0,0),(8.0,0),
!    (9.0,0), (10.0,0),(11.0,0),(12.0,0)
!   )          
    print*, 'OUTPUT DATA: '
    print*, 'Destination matrix: '
    print*, 'Submatrix(3,2) is transposed'      
!   print destination matrix      
    do i= 1, 3
          print 102,  ( A(i, j), j = 1, 4)
    enddo

102 format(9x,10(f6.2,',',f6.2,'i',3x))
        
    stop
    end