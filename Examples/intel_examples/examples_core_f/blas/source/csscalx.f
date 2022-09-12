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
!      C S S C A L  Example Program Text
!*******************************************************************************

      program CSSCAL_MAIN
*
      integer      n, incx
      integer      nmax
      parameter    (nmax=10)
      complex      cx(nmax)
      real         sa
      integer      i
*       External Subroutines
      external     CSSCAL, PrintVectorC
*
*   Executable Statementcs
*
      print*
      print*,'   C S S C A L  EXAMPLE PROGRAM'
*
*       Read input data from input file
      read*
      read*, n, incx
      read*, sa
      if ( (1+(n-1)*abs(incx)).gt.nmax) then
        print*, ' Insufficient memory for arrays'
        goto 999
      end if
      read*, (cx(i),i=1,1+(n-1)*abs(incx))
*
*       Print input data
      print*
      print*, '     INPUT DATA'
      print 100, n
      print 101, sa
      call PrintVectorC(0,n,cx,incx,'X ')
*
*      Call CSSCAL subroutine
      call CSSCAL(n,sa,cx,incx)
*
      print*
      print*, '     OUTPUT DATA'
      call PrintVectorC(1,n,cx,incx,'X ')

      stop
 100  format(7x,'N=',i2)
 101  format(7x,'SA=',f5.2)
 999  stop 1
      end
