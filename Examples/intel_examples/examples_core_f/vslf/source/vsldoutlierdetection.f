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
!   Outlier detection function (without users weight)  Example Program Text
!******************************************************************************/

      include 'mkl_vsl.f90'
      include "errcheck.inc"
      include "generatedata.inc"

      program MKL_VSL_TEST

      USE MKL_VSL_TYPE
      USE MKL_VSL

      integer(kind=8),parameter :: DIM      = 10    ! Task dimension
      integer(kind=8),parameter :: NN       = 10000 ! Number of observations

!     Persent of outliers in observation
      integer,parameter :: EPSILON = 2
!     Mean of outliers in observation
      real(kind=8),parameter :: OUTL_MEAN = 30.0

!     BACON algorithm parameters
      integer,parameter :: NNPARAMS    = 3
      real(kind=8),parameter :: ALPHA  = 0.001 / NN
      real(kind=8),parameter :: BETA   = 0.005
      integer(kind=8),parameter :: INIT_METHOD =                        &
     &                             VSL_SS_METHOD_BACON_MAHALANOBIS_INIT

      TYPE(VSL_SS_TASK) task
      integer p
      integer n
      integer x_storage
      real(kind=8) x(DIM,NN)
      real(kind=8) weights(NN)
      real(kind=8) C(DIM,DIM)
      real(kind=8) a(DIM)
      real(kind=8) a_outl(DIM)
      real(kind=8) ro, eps
      integer i, j, k
      integer(kind=4) errcode
      integer cntoutl
      integer ndetoutl, ncorrdetoutl, nwrongdetoutl, nlostoutl
      real(kind=8) corrdetratio, wrongdetratio, lostratio
      integer task_method
      integer(kind=8) task_params
      integer errnums

      integer nparams
      real(kind=8) BACONparams(NNPARAMS)
      real(kind=8) BACONweights(NN)

!     ***** Initializing parameters for Summary Statistics task *****
      p           = DIM
      n           = NN
      x_storage   = VSL_SS_MATRIX_STORAGE_COLS
      ro          = 1.0 / ( p * 30.0 )
      eps         = EPSILON / 100.0
      task_params = VSL_SS_OUTLIERS
      task_method = VSL_SS_METHOD_BACON
      errcode     = 0
      errnums     = 0

!     ***** Generate data set *****
!     Definition of covariance matrix and means for input data
      do i = 1, p
        a(i) = 0.0
        a_outl(i) = OUTL_MEAN

        do j = 1, i - 1
          C(i, j) = ro
          C(j, i) = ro
        end do

        C(i, i) = 1.0
      end do

      errcode = dGenerateOutliersInput( p, n, x, eps, a, a_outl, C,     &
     &                                  weights, cntoutl )

!     ***** Create Summary Statistics task *****
      errcode = vsldssnewtask( task, p, n, x_storage, x )
      call CheckVslError( errcode )

!     ***** Initializing BACON algorithm parameters *****
      nparams = NNPARAMS

      BACONparams(1) = INIT_METHOD
      BACONparams(2) = ALPHA
      BACONparams(3) = BETA

      do i = 1, n
        BACONweights(i) = 0.0
      end do

!     ***** Initialization of the task parameters for Bacon algorithm *****
      errcode = vsldsseditoutliersdetection( task, nparams,             &
     &                                       BACONparams, BACONweights )
      call CheckVslError( errcode )


!     ***** Detect outliers using BACON method *****
      errcode = vsldsscompute( task, task_params, task_method )
      call CheckVslError( errcode )

!     ***** Compute stat characteristics of the test *****
      ndetoutl      = 0
      ncorrdetoutl  = 0
      nwrongdetoutl = 0
      nlostoutl     = 0

      do i = 1, n
        if (BACONweights(i) == 0.0) then
          ndetoutl = ndetoutl + 1
          if (weights(i) == 0.0) then
            ncorrdetoutl = ncorrdetoutl + 1
          else
            nwrongdetoutl = nwrongdetoutl + 1
          end if
        else
          if (weights(i) == 0.0) then
            nlostoutl = nlostoutl + 1
          end if
        end if
      end do

      errnums = nwrongdetoutl + nlostoutl

      if (cntoutl > 0) then
        corrdetratio = (100.0 * ncorrdetoutl) / cntoutl
        lostratio    = (100.0 * nlostoutl) / cntoutl
      else
        corrdetratio = 0.0
        lostratio    = 0.0
      end if

      if (ndetoutl > 0) then
        wrongdetratio = (100.0 * nwrongdetoutl) / ndetoutl
      else
        wrongdetratio = 0.0
      end if

!     ***** Printing results *****
      print 9,  'Task dimension :         ', p
      print 9,  'Number of observations : ', n
      print 10, 'Tail of Chi2 distribution alpha : ', ALPHA
      print 10, 'Stopping criteria beta :          ', BETA
      print *, ''

!     ***** Printing summary results of BACON algorithm *****
      print *, 'Results of BACON algorithm:'
      print 9, '    number of outliers in dataset : ', cntoutl
      print 9, '    number of detected outliers :   ', ndetoutl
      print *, ''

      print *, 'Check of the output of BACON algorithm:'
      print 11, '    Ratio of correctly detected  outliers  = ',        &
     &         corrdetratio, '%'
      print 11, '    Ratio of incorrectly detected outliers = ',        &
     &         wrongdetratio, '%'
      print 11, '    Ratio of missed outliers               = ',        &
     &         lostratio, '%'
      print *, ''
      print *, ''

!     ***** Printing summary of the test *****
      if ( errnums == 0 ) then
         print *, ' All outliers are correctly detected.'
      else
         print *, ' Error: Some of outliers are incorrectly detected',  &
     &            ' or missed.'
         stop 1
      end if

!     ***** Delete Summary Statistics task *****
      errcode = vslssdeletetask( task )
      call CheckVslError( errcode )

      call MKL_FREE_BUFFERS()

9     format (A,I5)
10    format (A,F10.7)
11    format (A,F6.2,A)

      end
