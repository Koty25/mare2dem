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
!     D R O T M G   Example Program Text
!*******************************************************************************

      program   DROTMG_MAIN

      double precision    param(5)
      double precision    dd1, dd2, dx1, dy1
      integer             i
*       External Subroutines
      external            DROTMG
*
*       Executable Statements
*
      print*
      print*,'   D R O T M G  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, dd1, dd2, dx1, dy1
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, dd1, dd2, dx1, dy1
*
*       Call DROTMG subroutine
      call DROTMG(dd1, dd2, dx1, dy1, param)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, dd1, dd2, dx1
      print 102, (param(i),i=1,5)
*
 100  format(7x,'DD1=',f5.2,'  DD2=',f5.2,'  DX1=',f5.2,'  DY1=',f5.2)
 101  format(7x,'DD1=',f6.3,'  DD2=',f6.3,'  DX1=',f6.3)
 102  format(7x,'VECTOR PARAM (some elements may be omitted)',/9x,5(f6.3
     $,1x))
      stop
      end
