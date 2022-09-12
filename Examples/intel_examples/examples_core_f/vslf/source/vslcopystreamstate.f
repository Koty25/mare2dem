!===============================================================================
! Copyright 2003-2020 Intel Corporation.
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
!    vslCopyStreamState  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE (VSL_STREAM_STATE) :: stream
      TYPE (VSL_STREAM_STATE) :: streamCpy
      integer(kind=4) seed(6)
      integer(kind=4) r(1000)
      integer(kind=4) rCpy(1000)
      integer(kind=4) i,err,nn
      integer n,nnn
      integer brng,method,seedCpy
      integer(kind=4) errcode

      err=0
      n=1000
      nn=10
      brng=VSL_BRNG_MRG32K3A
      method=VSL_RNG_METHOD_UNIFORMBITS_STD
      nnn=6

!     ***** Initialize seeds *****
      seedCpy=1
      seed(1)=7777777
      do i=2,6
        seed(i)=seed(i-1)+11
      end do

!     ***** Initialize streams *****
      errcode=vslnewstreamex( stream,    brng,  nnn, seed )
      call CheckVslError(errcode)
      errcode=vslnewstream  ( streamCpy, brng,       seedCpy )
      call CheckVslError(errcode)
      errcode=vslcopystreamstate( streamCpy, stream )
      call CheckVslError(errcode)

!     ***** Call RNGs *****
      errcode=virnguniformbits( method, stream,   n, r )
      call CheckVslError(errcode)
      errcode=virnguniformbits( method, streamCpy, n, rCpy )
      call CheckVslError(errcode)

!     ***** Compare results *****
      do i=1,1000
        if (r(i) .NE. rCpy(i)) then
          err=err+1
        end if
      end do

!     ***** Printing results *****
      print *,"Sample of vslCopyStreamState"
      print *,"----------------------------"
      print *,""
      print *,"Parameters:"
      print 12,"    seed    = { ",seed(1)," ",                          &
     &  seed(2)," ",seed(3)," ",                                        &
     &  seed(4)," ",seed(5)," ",seed(6), " }"
      print 11,"    seedCpy =   ",seedCpy

      print *,""
      print *,"Results (first 10 of 1000):"
      print *,"---------------------------"
      do i=1,nn
        print 10, "r[",i-1,"]=0x",r(i)," rCpy[",i-1,"]=0x",rCpy(i)
      end do

      print *,""
      if (err>0) then
        print 13,"Error: ", err," values are incorrect!"
        stop 1
      else
        print *,"Results of original stream and its copy are identical."
      end if

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )
      call CheckVslError(errcode)
      errcode=vsldeletestream( streamCpy )
      call CheckVslError(errcode)

10    format(A,I1,A,Z8.8,A,I1,A,Z8.8)
11    format(A,I1)
12    format(A,I7,A,I7,A,I7,A,I7,A,I7,A,I7,A)
13    format(A,I8,A)

      end
