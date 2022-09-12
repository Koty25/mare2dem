!===============================================================================
! Copyright 1999-2020 Intel Corporation.
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
!      Z R O T G  Example Program Text
!*******************************************************************************

      program  ZROTG_MAIN
*
      double precision  c
      complex*16        s, ca, cb
*       External Subroutines
      external          ZROTG
*
*       Executable Statements
*
      print*
      print*,'   Z R O T G  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, ca, cb
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, ca, cb
*
*       Call ZROTG subroutine
      call ZROTG(ca,cb,c,s)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, ca
      print 102, c,  s
*
 100  format(7x,'CA=(',f6.3,',',f6.3,')  CB=(',f6.3,',',f6.3,')')
 101  format(7x,'CA=(',f6.3,',',f6.3,')')
 102  format(7x,' C= ',f6.3,11x,'S=(',f6.3,',',f6.3,')')
      stop
      end
