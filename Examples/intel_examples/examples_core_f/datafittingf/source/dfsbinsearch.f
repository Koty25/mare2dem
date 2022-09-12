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
!    Cells search in blocks Example Program Text
!*******************************************************************************

      include 'mkl_df.f90'
      include "errcheck.inc"
      include "rescheck.inc"
      include "generatedata.inc"

      PROGRAM MKL_DF_TEST

      USE MKL_DF_TYPE
      USE MKL_DF
      USE DF_GENERATE_DATA
      USE DF_EXAMPLE_RESCHECK


      ! number of breakpoints
      INTEGER,PARAMETER :: N           = 10
      ! number of sites for cell search in block
      INTEGER,PARAMETER :: NSITE_BLOCK = 15
      ! number of blocks
      INTEGER,PARAMETER :: NNBLOCKS    = 3
      ! length of the datahint array
      INTEGER,PARAMETER :: LDATA       = 5
      ! total number of sites for cell search
      INTEGER,PARAMETER :: NNSITE      = (NNBLOCKS*NSITE_BLOCK)
      ! left  limit of search interval
      REAL(4),PARAMETER :: LEFT_LIMIT  = -10.0
      ! right limit of search interval
      REAL(4),PARAMETER :: RIGHT_LIMIT = 10.0

      ! mean of the sites distribution
      REAL(4),PARAMETER :: MEAN_  = 0.0
      ! sandard deviation of the sites distribution
      REAL(4),PARAMETER :: SIGMA_ = 1.0

      ! Data Fitting task descriptor
      TYPE (DF_TASK) task

      ! number of break points
      INTEGER :: nx
      ! additional info about break points
      INTEGER :: xhint
      ! number of search sites in one block
      INTEGER :: nsite_bl
      ! total number of sites
      INTEGER :: nsite
      ! number of blocks of sites
      INTEGER :: nblocks
      ! search sites storage format
      INTEGER :: sitehint
      ! indices of cells containing search sites
      INTEGER,TARGET  :: cell(NNSITE)
      INTEGER,POINTER :: cell_ptr(:)
      ! array of break points
      REAL(4) :: x(N)
      ! array of search sites
      REAL(4),TARGET :: site(NNSITE)
      REAL(4),POINTER :: site_ptr(:)
      ! additional info about the structure of search sites
      REAL(4) :: datahint(LDATA)
      ! parameters of the search sites distribution
      REAL(4) :: a,sigma
      ! limits of the search interval
      REAL(4) :: left, right
      ! method that is used to perform calculations
      INTEGER :: method

      INTEGER :: test_cell(NNSITE)
      INTEGER :: i,j
      INTEGER :: errnums
      INTEGER(4) :: errcode
      REAL(4) :: r

      !***** Initializing parameters for Data Fitting task *****

      errcode = 0
      errnums = 0

      !***** Parameters describing search interval *****
      nx = N
      xhint = DF_NON_UNIFORM_PARTITION

      !***** Parameter describing number of search sites *****
      nsite_bl  = NSITE_BLOCK
      nsite     = NNSITE
      !  No additional info is provided in case of cell search
      sitehint  = DF_NO_HINT

      !***** Generate uniformly distributed search sites *****
      left = LEFT_LIMIT
      right = RIGHT_LIMIT
      errcode = sUniformRandSortedData( x, left, right, nx )
      CALL CheckDfError(errcode)

      !***** Parameter describing additional info about the search sites *****
      !  Presumably most likely cell is cell number (N/2 + 1)
      datahint(1) = 1.0
      datahint(2) = FLOAT(DF_APRIORI_MOST_LIKELY_CELL)
      datahint(3) = 0.0
      datahint(4) = 1.0
      datahint(5) = FLOAT(((nx/2) + 1))

      !***** Generate normally distributed break points *****
      a         = MEAN_;
      sigma     = SIGMA_;
      errcode = sPeakData( site, a, sigma, nsite )
      CALL CheckDfError(errcode)

      !***** Create Data Fitting task *****
      errcode = dfsNewTask1D( task, nx, x, xhint )
      CALL CheckDfError(errcode)

      !***** Perform cells search *****
      nblocks = NNBLOCKS

      method = DF_METHOD_STD

      DO i = 1, nblocks
        site_ptr => site((i-1)*nsite_bl+1:)
        cell_ptr => cell((i-1)*nsite_bl+1:)
        errcode = dfsSearchCells1D( task, method, nsite_bl,              &
     &     site_ptr, sitehint, datahint, cell_ptr )
        CALL CheckDfError(errcode)
      END DO

      !***** Check search results *****
      errcode = sFindCells( nx, x, nsite, site, test_cell )
      DO i = 1, nsite
        test_cell(i) = test_cell(i) - 1
        IF ( test_cell(i) /= cell(i) ) errnums = errnums+1
      END DO

      !***** Print given function *****
      WRITE (*,901) "Number of break points : ",nx
      WRITE (*,901) "Number of search sites : ",nsite
      WRITE (*,902) "  X:"
      DO j = 1, nx
        WRITE (*,903) " ",x(j)
      END DO
      WRITE (*,902) "Results of cell search:"
      WRITE (*,904) "    Site          Computed idx      Expected idx"
      DO i = 1, nsite
        WRITE (*,905) " ",site(i),"   ",cell(i),"   ",test_cell(i)
      END DO

      !***** Delete Data Fitting task *****
      errcode = dfDeleteTask( task )
      CALL CheckDfError(errcode)

      !***** Print summary of the test *****
      IF (errnums /= 0) THEN
        WRITE (*,906) "","Error: Computed cell indices are incorrect"
        STOP 1
      ELSE
        WRITE (*,906) "","Computed cell indices are correct"
      END IF
      STOP 0

901   FORMAT (A,I0)
902   FORMAT (/A)
903   FORMAT (A,SP,F11.6)
904   FORMAT (99A)
905   FORMAT (A,SP,F11.6,A,I15,A,I15)
906   FORMAT (//A)
      END PROGRAM
