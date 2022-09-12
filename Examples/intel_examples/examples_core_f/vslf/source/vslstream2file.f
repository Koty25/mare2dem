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
!    stream2file functions  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      use MKL_VSL_TYPE
      USE MKL_VSL

      type (VSL_STREAM_STATE) :: stream
      integer N
      parameter (N=10)
      integer brng,method,seed

      real r_orig(N)
      real r_load(N)

      integer(kind=4) i, nerrors
      integer(kind=4) errcode

      brng=VSL_BRNG_R250
      method=VSL_RNG_METHOD_UNIFORM_STD
      seed=777
      nerrors=0

!     **** Create the original stream to be saved in a file ****
      errcode = vslnewstream( stream, brng, seed )
      call CheckVslError(errcode)

!     **** Save original stream to a file ****
      errcode = vslsavestreamf(stream, "vslstream2file.dat");
      call CheckVslError(errcode)


!     **** Generate random numbers using original stream ****
      errcode = vsrnguniform(method, stream, N, r_orig, 0.0, 1.0)
      call CheckVslError(errcode)

!     **** Delete original stream ****
      errcode = vsldeletestream(stream)
      call CheckVslError(errcode)

!     **** Load stream that is saved in a file ****
      errcode = vslLoadStreamF(stream, "vslstream2file.dat")
      call CheckVslError(errcode)

!     **** Generate random numbers using the stream loaded from file ****
      errcode = vsrnguniform(method, stream, N, r_load, 0.0, 1.0)
      call CheckVslError(errcode)

!     **** Delete stream loaded from file ****
      errcode = vsldeletestream(stream)
      call CheckVslError(errcode)

!     **** Compare random numbers from original and loaded stream. ****
!     **** Must be identical                                       ****
      do i=1,N
          print '(A,I2,A,F7.5,A,I2,A,F7.5)'," r_orig[",i,"]=",          &
     &          r_orig(i),"     r_load[",i,"]=",r_load(i)
          if ( r_orig(i) /= r_load(i) ) then
!             **** Here if results are not identical ****
              nerrors = nerrors + 1
          end if
      end do

      if (nerrors > 0) then
!         **** Here if results are not identical ****
          print *, "ERROR: Loaded stream differs from original stream."
          stop 1
      else
!         **** Here if results are identical ****
          print *, "PASS: Loaded stream identical with original stream."
      end if

      end program
