!===============================================================================
! Copyright 2013-2020 Intel Corporation.
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
!    Calculation of mean absolute deviation  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 3        ! Task dimension
      integer,parameter :: NN  = 1000     ! Number of observations

      TYPE(VSL_SS_TASK) task
      integer p, n, x_storage
      real(kind=4) x(NN,DIM)
      real(kind=4) mnad(DIM)
      real(kind=4) a, sigma
      integer i, j
      integer(kind=4) errcode
      integer errnums
      integer method, ed_estimate
      integer(kind=8) estimate
      real(kind=4) tD, tQ, tD2, s, sig, sD, deltaD, dn


!     ***** Initialize parameters of Summary Statistics task *****
      p           = DIM
      n           = NN
      x_storage   = VSL_SS_MATRIX_STORAGE_ROWS
      a           = 0.0
      sigma       = 1.0
      estimate    = VSL_SS_MNAD
      ed_estimate = VSL_SS_ED_MNAD
      method      = VSL_SS_METHOD_FAST
      errcode     = 0
      errnums     = 0
      dn          = NN

!     ***** Generate data set using VSL Gaussian RNG
!           with mean a = 0 and stdev = 1 *****
      errcode = sGenerateGaussianData( p, n, x, a, sigma )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vslsssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )


!     ***** Provide array of mean absolute deviation *****
      errcode = vslsssedittask( task, ed_estimate, mnad )
      call CheckVslError( errcode )

!    ***** Compute mean absolute deviation using FAST method *****
      errcode = vslssscompute( task, estimate, method )
      call CheckVslError( errcode )

!    ***** Check the correctness of computed MnAD *****
!    ***** Testing relies on property that for Gaussian distribution
!          standard deviation estimate ~= 1.2533 * mnad ******
      tD = sigma*sigma
      tQ = 720.0*sigma*sigma*sigma*sigma
      tD2=tD*tD
      s=((tQ-tD2)/dn)-(2.0*(tQ-2*tD2)/(dn*dn))+((tQ-3*tD2)/(dn*dn*dn))

      errnums = 0
      do i = 1, p
         sig = 1.2533 * mnad(i)
         sD = sig * sig
         deltaD = abs((tD-sD) / sqrt(s))
         if  ( deltaD > 3.0 ) then
             errnums = errnums + 1
         end if
      end do

!     ***** Printing results *****
      print *, 'Task dimension :         ', p
      print *, 'Number of observations : ', n
      print *, ''

!     ***** Printing computed mean absolute deviations *****
!     ***** Printing part of the initial matrix of observations *****

      print *,"Mean absolute deviation for all variables"

      print *, 'MnAD1 MnAD2 MnAD3'
      print 14, mnad(1), mnad(2), mnad(3)

!     ***** Printing summary of the test *****
      print *, ''
      if( errnums == 0 ) then
        print *, 'Computed mean absolute deviation',                    &
     &           ' agrees with theory.'
      else
        print *, 'Error: Computed mean absolute deviation',             &
     &           ' disagrees with theory.'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

14    format(F6.3     F6.3     F6.3)

      end
