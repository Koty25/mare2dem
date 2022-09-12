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
!      S R O T M G  Example Program Text
!*******************************************************************************

      program  SROTMG_MAIN
*
      real       param(5)
      real       sd1, sd2, sx1, sy1
*       External Subroutines
      external   SROTMG
*
*       Executable Statements
*
      print*
      print*,'   S R O T M G  EXAMPLE PROGRAM'
*       Read input data from input file
      read*
      read*, sd1, sd2, sx1, sy1
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, sd1, sd2, sx1, sy1
*
*       Call SROTMG subroutine
      call SROTMG(sd1,sd2,sx1,sy1,param)
*
      print*
      print*, '     OUTPUT DATA'
      print 101, sd1, sd2, sx1
      print 102, (param(i), i=1,5)
*
 100  format(7x,'SD1=',f4.1,'  SD2=',f4.1,'  SX1=',f4.1,
     $       '  SY1=',f4.1)
 101  format(7x,'SD1=',f5.2,'  SD2=',f5.2,'  SX1=',f5.2)
 102  format(7x,'VECTOR PARAM (some elements may be omitted)',/9x,5(f6.3
     $,1x))
      stop
      end
