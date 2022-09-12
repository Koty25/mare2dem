!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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
!      D R O T G  Example Program Text
!*******************************************************************************

      program  DROTG_MAIN

      use f95_precision, only: wp => dp
      use blas95, only: rotg

      implicit none

      real(wp) :: c, s
      real(wp) :: sa, sb

!     Executable Statements

      print*
      print*,'   D R O T G  EXAMPLE PROGRAM'

!     Read input data from input file
      read*
      read*, sa, sb

!     Print input data
      print*
      print*, '     INPUT DATA'
      print 100, sa, sb

!     Call DROTG subroutine
      call ROTG(sa, sb, c, s)

      print*
      print*, '     OUTPUT DATA'
      print 100, sa, sb
      print 101, c,  s

 100  format(7x,'SA=',f6.3,'  SB=',f6.3)
 101  format(7x,' C=',f6.3,'   S=',f6.3)
      stop
      end
