!===============================================================================
! Copyright 2001-2020 Intel Corporation.
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
!    vmlSet/Get/ClearErrorCallBack  Example Program Text
!*******************************************************************************

      program MKL_VML_TEST
      include "mkl_vml.f90"

      INTEGER errcb,terrcb

!     User callback function should be external in the program
!DEC$ ATTRIBUTES C, REFERENCE :: UserCallBack
      INTEGER,EXTERNAL :: UserCallBack

!     Arguments for vdLn function
      real(kind=8) :: a(1), r(1)
      integer      :: n

      print *,"Set/Get/Clear CallBack example program"
      print *,""

!     Testing vmlGetErrorCallback
      errcb=0
      errcb=VMLGETERRORCALLBACK()
      print 10,"Get CallBack returned address:",errcb,"h"

!     Registering user error callback function
      terrcb=VMLSETERRORCALLBACK(UserCallBack)
      errcb=VMLGETERRORCALLBACK()
      print 10,"User callback function address:",errcb,"h"
      print *,""
      print *,"Calling vdLn function..."

!     Call vdLn with invalid argument to test user callback function
      a(1) = 0.0
      n = 1
      call VDLN(n,a,r)

!     Testing vmlClearErrorCallback
      terrcb=VMLCLEARERRORCALLBACK()
      errcb=VMLGETERRORCALLBACK()
      print *,""
      print 10,"Callback address after calling Clear CallBack:",        &
     &          errcb,"h"

10    format (A,Z8.8,A)

      end program MKL_VML_TEST

      integer function UserCallBack( es )
!DEC$ ATTRIBUTES C, REFERENCE :: UserCallBack
      include "mkl_vml.f90"

      type(ERROR_STRUCTURE) :: es

      print 20,"In function vdLn argument a[",es%iindex,"]=",es%dba1,   &
     &         " is wrong."
      UserCallBack = 0

20    format (A,I1,A,F4.1,A)

      end function
