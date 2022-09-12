!===============================================================================
! Copyright 2019-2020 Intel Corporation.
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
!    vslSkipAheadStreamEx functions  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      TYPE NEW_TYPE
          TYPE (VSL_STREAM_STATE):: s1
      END TYPE NEW_TYPE

      TYPE (NEW_TYPE) :: stream
      TYPE (NEW_TYPE) :: streamS
      integer(kind=4) r(1000)
      integer(kind=4) rS(1000)
      integer(kind=4) i,j,k,err,s
      integer(kind=8) nskip
      integer(kind=8) skipTimes
      integer nParams
      integer(kind=8) params(2)

      integer n
      integer brng,method,seed
      integer(kind=4) errcode

!     ***** To skip 2^76 elements in the random stream SkipAheadStream(nskip) function should be called 2^14 times
!           with nskip equal to 2^62.

      nskip = scale(1.0,62)
      skipTimes = scale(1.0,14)

!     ***** To skip 2^76 elements in the random stream SkipAheadStreamEx(nskip) function should be called
!           with nskip represented as
!               nskip = 2^76 = 0 + 2^12 * 2^64
!           In general case
!               nskip = params[0] + params[1] * 2^64 + params[2] * 2^128 + ... *****
      nParams = 2
      params(1) = 0
      params(2) = scale(1.0,12)

      err=0
      n=1000
      brng=VSL_BRNG_MRG32K3A
      method=VSL_RNG_METHOD_UNIFORMBITS_STD
      seed=7777777

!     ***** Create main stream *****
      errcode=vslnewstream  ( stream%s1,   brng,  seed )
      call CheckVslError(errcode)

!     ***** Create skipahead streams as copies of the main one *****
      errcode=vslcopystream( streamS%s1,stream%s1 )
      call CheckVslError(errcode)

!     ***** Apply vslSkipAheadStreamEx to streamS with the nskip = 0 + 2^12 * 2^64 = 2^76 *****
      errcode=vslskipaheadstreamex( streamS%s1, nParams, params)
      call CheckVslError(errcode)

!     ***** Apply vslSkipAheadStream to stream 2^14 times with the nskip = 2^62 *****
      do i=1,skipTimes
        errcode=vslskipaheadstream( stream%s1, nskip )
        call CheckVslError(errcode)
      end do

!     ***** Generate random numbers for SkipAhead stream *****
      errcode=virnguniformbits( method, stream%s1, n, r )
      call CheckVslError(errcode)

!     ***** Generate random numbers for SkipAheadEx stream *****
      errcode=virnguniformbits( method, streamS%s1, n, rS)
      call CheckVslError(errcode)

!     ***** Compare results *****
      do i=1,1000
        if (r(i) .NE. rS(i)) then
          err=err+1
        end if
      end do

!     ***** Printing results *****
      print *,"Sample of vslSkipAheadStreamEx"
      print *,"----------------------------"
      print *,""
      print *,"Parameters:"
      print 11,"    seed   =  ",seed
      print *,""
      print *,"Results (first 10 of 1000):"
      print *,"---------------------------"
      do i=1,10
        print 10, "r[",i-1,"]=0x",r(i)," rS[",i-1,"]=0x",rS(i)
      end do

      print *,""
      if (err>0) then
        print 13,"Error: ", err," values are incorrect!"
        stop 1
      else
        print *,"Results of SkipAhead and SkipaheadEx streams are are
     & identical"
      end if

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream%s1 )
      call CheckVslError(errcode)
      errcode=vsldeletestream( streamS%s1 )
      call CheckVslError(errcode)

10    format(A,I1,A,Z8.8,A,I1,A,Z8.8)
11    format(A,I8)
12    format(A,I7,A,I1,A,I1,A,I1,A,I1,A,I1,A)
13    format(A,I8,A)

      end
