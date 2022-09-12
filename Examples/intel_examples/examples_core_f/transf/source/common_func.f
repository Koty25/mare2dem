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
!
!*******************************************************************************
!   Is'nt work yet :)
      subroutine printArray(a, n, m, data_type)

      integer       m, n
      character*1   data_type
      real          a(n,*)
      integer       i, j

      if (lge(data_type,'s')) then
        do i=1,n
           print 100, (a(i,j),j=1,m)
        end do
      else if (lge(data_type,'d')) then
             do i=1,n
                print 100, (a(i,j),j=1,m)
             end do
           else if (lge(data_type,'c')) then
                    do i=1,n
                       print 110, (a(i,j),j=1,m)              
                    end do
                else if(lge(data_type,'z')) then
                      do i=1,n
                        print 110, (a(i,j),j=1,m)
                      end do
                    end if

 100  format(9x, 10(f6.2,',',f6.2,3x)) ! for d, s
 110  format(9x, 10(f8.3,2x))          ! for c, z

      return
      end

