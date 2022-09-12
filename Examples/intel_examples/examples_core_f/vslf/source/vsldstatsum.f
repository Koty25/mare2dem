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
!    Calculation of raw/central sums Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"
      include "statchars.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer,parameter :: DIM = 4        ! Task dimension
      integer,parameter :: NN  = 1000     ! Number of observations

      real(kind=8),parameter :: P_THRESHOLD = 0.01

      real(kind=8) C(DIM*DIM)
      data C / 1.0, 0.0, 0.0, 0.0,                                      &
     &         0.0, 1.0, 0.0, 0.0,                                      &
     &         0.0, 0.0, 1.0, 0.0,                                      &
     &         0.0, 0.0, 0.0, 1.0 /

      real(kind=8) a(DIM)
      data a / 5.0, 5.0, 5.0, 5.0 /

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      real(kind=8) x(DIM,NN)

      real(kind=8) mean(DIM), sum(DIM)

      real(kind=8) rs2(DIM), rs3(DIM), rs4(DIM)
      real(kind=8) rm2(DIM), rm3(DIM), rm4(DIM)

      real(kind=8) cs2(DIM), cs3(DIM), cs4(DIM)
      real(kind=8) cm2(DIM), cm3(DIM), cm4(DIM)

      integer i, j, order
      integer(kind=4) errcode
      integer errnums, method
      integer(kind=8) estimate

      real(kind=8) pval_mean(DIM)
      real(kind=8) pval_rm2(DIM), pval_rm3(DIM), pval_rm4(DIM)
      real(kind=8) pval_cm2(DIM), pval_cm3(DIM), pval_cm4(DIM)

!     ***** Initializing parameters for Summary Statistics task *****
      p               = DIM
      n               = NN
      x_storage       = VSL_SS_MATRIX_STORAGE_COLS

!     ***** Generate data set using VSL GaussianMV RNG *****
      errcode = dGenerateGaussianMVData( p, n, x, a, C )
      call CheckVslError( errcode )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Edit task parameters to be used in computations *****
      errcode = vsldssedittask( task, VSL_SS_ED_MEAN, mean )
      call CheckVslError( errcode )

!    ***** Edit task parameters for computation of sum, 2nd, 3rd
!          and 4th raw and central sums *****
      errcode = vsldsseditsums( task, sum, rs2, rs3, rs4,               &
     &                                     cs2, cs3, cs4 )
      call CheckVslError( errcode )


!     ***** Sum, 2nd, 3rd and 4th raw and central sums are included
!            in the list of estimates to compute *****
      estimate = VSL_SS_SUM
      estimate = ior( estimate, VSL_SS_2R_SUM )
      estimate = ior( estimate, VSL_SS_3R_SUM )
      estimate = ior( estimate, VSL_SS_4R_SUM )
      estimate = ior( estimate, VSL_SS_2C_SUM )
      estimate = ior( estimate, VSL_SS_3C_SUM )
      estimate = ior( estimate, VSL_SS_4C_SUM )

      method = VSL_SS_METHOD_1PASS

!     ***** Compute the estimates using 1PASS method *****
      errcode = vsldsscompute( task, estimate, method )
      call CheckVslError( errcode )
!    ***** Edit task parameters for computation of mean, 2nd, 3rd
!           and 4th raw and central moments *****
      errcode = vsldsseditmoments( task, mean, rm2, rm3, rm4,           &
     &                             cm2, cm3, cm4 )
      call CheckVslError( errcode )

!    ***** Convert sums into moments *****
      estimate = VSL_SS_MEAN
      estimate = ior( estimate, VSL_SS_2R_MOM )
      estimate = ior( estimate, VSL_SS_3R_MOM )
      estimate = ior( estimate, VSL_SS_4R_MOM )
      estimate = ior( estimate, VSL_SS_2C_MOM )
      estimate = ior( estimate, VSL_SS_3C_MOM )
      estimate = ior( estimate, VSL_SS_4C_MOM )

      method = VSL_SS_METHOD_SUM_TO_MOM
      errcode = vsldsscompute( task, estimate, method )
      call CheckVslError( errcode )



!     ***** Testing stat characteristics of computed estimates *****
      errnums = 0

!     Compute p-values for mean estimates
      call dComputePvalsMean( p, n, mean, a, C, pval_mean )
      order = 2
      call dComputePvalsRawMoments( p, n, rm2, order, a, C, pval_rm2 )
      order = 3
      call dComputePvalsRawMoments( p, n, rm3, order, a, C, pval_rm3 )
      order = 4
      call dComputePvalsRawMoments( p, n, rm4, order, a, C, pval_rm4 )
!     Compute p-values for central moments estimates
      order = 2
      call dComputePvalsCentralMoments( p, n, cm2, order, a, C,         &
     &                                  pval_cm2 )
      order = 3
      call dComputePvalsCentralMoments( p, n, cm3, order, a, C,         &
     &                                  pval_cm3 )
      order = 4
      call dComputePvalsCentralMoments( p, n, cm4, order, a, C,         &
     &                                  pval_cm4 )

!     ***** Checking the validity of p-values for all estimates *****
      errnums = 0
      do i = 1, p
        if (pval_mean(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_rm2(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_rm3(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_rm4(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cm2(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cm3(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
        if (pval_cm4(i) < P_THRESHOLD) then
          errnums = errnums + 1
        end if
      end do

!     ***** Printing results *****
      print 9, 'Task dimension :         ', p
      print 9, 'Number of observations : ', n
      print *, ''

!     ***** Printing computed minimum, maximum, mean and moments estimates *****
      print 10, '                    '
      print 10, 'SUM           2RS            3RS            '
      print 10, '4RS           2CS            3CS            4CS'
      print *, ''

      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 12, sum(i), ' '
        print 13, rs2(i), ' ', rs3(i), ' ', rs4(i), ' '
        print 13, cs2(i), ' ', cs3(i), ' ', cs4(i), ''
        print *, ''
      end do
      print *, ''

!     ***** Printing p-values for moments obtained from sums *****
      print *, 'P-values of the moments obtained from sums'
      print *, ''
      print *, ''
      print 10, '                    '
      print 10, 'MEAN          2RM           3RM           '
      print 10, '4RM           2CM           3CM           4CM'
      print *, ''

      do i = 1, p
        print 11, 'Variable #', i, ' '
        print 12, pval_mean(i), ''
        print 13, pval_rm2(i), '', pval_rm3(i), '', pval_rm4(i), ''
        print 13, pval_cm2(i), '', pval_cm3(i), '', pval_cm4(i), ''
        print *, ''
      end do
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
        print *, ' All the computed estimates agree with theory'
      else
        print *, ' Error: At least one of the computed estimates',      &
     &           ' disagrees with theory'
        stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

9     format(A,I4)
10    format(A,$)
11    format(A,I1,A,$)
12    format(F14.6,A,$)
13    format(F14.6,A,F14.6,A,F14.6,A,$)
14    format(4F10.6,$)

      end
