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
!    stream2memory functions  Example Program Text
!*******************************************************************************

      include 'mkl_vsl.f90'
      include "errcheck.inc"

      program MKL_VSL_TEST

      use MKL_VSL_TYPE
      USE MKL_VSL

      type (VSL_STREAM_STATE) :: stream_orig
      type (VSL_STREAM_STATE) :: stream_load
      integer DN, N
      parameter (DN=101)
      parameter ( N= 10)

      integer brng,method,seed

      real(kind=8) ad, bd
      real(kind=8) rd_orig(DN),rd_load(DN)

      integer membufsz
      parameter (membufsz=10000);
      integer(kind=1) membuf(membufsz);

      integer(kind=4) i, nerrors
      integer(kind=4) errcode

      ad = 0.0
      bd = 1.0
      seed=7777777
      brng=VSL_BRNG_SFMT19937
      method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2

!     **** Create the original stream to be saved in a memory ****
      errcode = vslnewstream( stream_orig, brng, seed )
      call CheckVslError(errcode)

!     **** Generate random numbers using original stream ****
      errcode = vdrnggaussian(method, stream_orig, DN, rd_orig, ad, bd)
      call CheckVslError( errcode )

!     **** Save stream to memory ****
!     **** Call function getstreamsize to compute memory size ****
!     **** to hold the stream ****
      errcode = vslsavestreamm( stream_orig, membuf )
      call CheckVslError( errcode )

!     **** Generate random numbers using original stream ****
      errcode = vdrnggaussian(method, stream_orig, DN, rd_orig, ad, bd)
      call CheckVslError( errcode )

!     **** Load stream from the memeory ****
      errcode = vslloadstreamm( stream_load, membuf )
      call CheckVslError( errcode )

!     **** Generate random numbers using loaded stream *****
      errcode = vdrnggaussian(method, stream_load, DN, rd_load, ad, bd)
      call CheckVslError( errcode )

!     **** Compare random numbers from original and loaded stream. ****
!     **** Must be identical                                       ****
       print *,"Gaussian numbers:"
      do i=1,N
          print '(A,I2,A,F8.5,A,I2,A,F8.5)'," rd_orig[",i,"]=",         &
     &          rd_orig(i),"     rd_load[",i,"]=",rd_load(i)
      end do

      nerrors = 0
      do i=1,DN
          if ( rd_orig(i) /= rd_load(i) ) then
!             **** Here if results are not identical ****
              nerrors = nerrors + 1
          end if
      end do

      if (nerrors > 0) then
!         **** Here if results are not identical ****
          print *, "ERROR: Loaded stream differs from original stream."
          stop 1
      end if

!     **** Delete original stream ****
      errcode = vsldeletestream(stream_orig)
      call CheckVslError( errcode )

!     **** Delete stream loaded from memory ****
      errcode = vsldeletestream(stream_load)
      call CheckVslError( errcode )

      print *, "PASS: Loaded stream identical with original stream."

      end program
