!===============================================================================
! Copyright 2005-2020 Intel Corporation.
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

! Content:
!       Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT interface
!       example program (Fortran-interface)
!
!       Forward-Backward 2D complex transform for double precision data inplace.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX                  (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE                   (obligatory)
!           DFTI_DIMENSION      = 2                             (obligatory)
!           DFTI_LENGTHS        = {M,N}                         (obligatory)
!           DFTI_FORWARD_SCALE  = 1.0                           (default)
!           DFTI_BACKWARD_SCALE = 1.0/(M*N)                     (default=1.0)
!
!*****************************************************************************

      PROGRAM DM_COMPLEX_2D_DOUBLE_EX1

      USE MKL_CDFT
      USE CDFT_EXAMPLE_SUPPORT

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      COMPLEX(8), DIMENSION(:,:), ALLOCATABLE :: X_EXP
      COMPLEX(8), ALLOCATABLE :: X_IN(:,:), LOCAL(:), WORK(:)
      INTEGER   ELEMENTSIZE
      PARAMETER (ELEMENTSIZE = 16)

      TYPE(DFTI_DESCRIPTOR_DM), POINTER :: DESC

      INTEGER   NX,NX_OUT,START_X,START_X_OUT,SIZE

      INTEGER   STATUS
      REAL(8)   SCALE
      INTEGER   LENGTHS(2)

      REAL(8)   MAXERR
      REAL(8), PARAMETER :: EPS = DOUBLE_EPS

      INTEGER(MPI_KIND) M, N
      INTEGER(MPI_KIND) ROOTRANK
      INTEGER(MPI_KIND) MPI_ERR
      INTEGER(MPI_KIND) MPI_NPROC
      INTEGER(MPI_KIND) MPI_RANK
      INTEGER(MPI_KIND) COMM
      INTEGER(4)        MKL_COMM

      CHARACTER*1024    INFILE

      LOGICAL FAILURE

      FAILURE = .FALSE.

!
!     1. Initiate MPI by calling MPI_Init (Perform MPI initialization)
!
      CALL MPI_INIT(MPI_ERR)

      IF (MPI_ERR .NE. MPI_SUCCESS) THEN
         PRINT *, 'MPI initialization error'
         PRINT *, 'TEST FAILED'
         STOP 1
      END IF

      COMM     = MPI_COMM_WORLD
      MKL_COMM = MPI_COMM_WORLD

      CALL MPI_COMM_SIZE(COMM, MPI_NPROC, MPI_ERR)
      CALL MPI_COMM_RANK(COMM, MPI_RANK, MPI_ERR)
      IF (MPI_RANK .EQ. 0) THEN
         PRINT '(" Program is running on ",I2," processes"/)',MPI_NPROC
       END IF

!
!     Read input parameters from input file
!     m - size of transform  along first dimension
!     n - size of transform  along second dimension
!
      IF (MPI_RANK .EQ. 0) THEN
         CALL GETARG(1, INFILE)
         OPEN ( 10, FILE = INFILE, STATUS = 'OLD', ACTION = 'READ')
         READ ( 10, FMT = * )
         READ ( 10, FMT = * ) M
         READ ( 10, FMT = * ) N
         CLOSE (10)
      END IF

      CALL MPI_BCAST(M,1_MPI_KIND,MPI_INTEGER,0_MPI_KIND,COMM,MPI_ERR)
      CALL MPI_BCAST(N,1_MPI_KIND,MPI_INTEGER,0_MPI_KIND,COMM,MPI_ERR)

      IF (LEGEND_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'DM_COMPLEX_2D_DOUBLE_EX1'
         PRINT *
         PRINT *, 'Forward-Backward 2D complex transform for double '// &
     &                          'precision data inplace'
         PRINT *
         PRINT *, 'Configuration parameters:'
         PRINT *
         PRINT *, 'DFTI_FORWARD_DOMAIN       = DFTI_COMPLEX'
         PRINT *, 'DFTI_PRECISION            = DFTI_DOUBLE '
         PRINT *, 'DFTI_DIMENSION            = 2'
         PRINT '(" DFTI_LENGTHS              = {",I3,",",I3,"}")', M, N
         PRINT *, 'DFTI_FORWARD_SCALE        = 1.0 '
         PRINT *, 'DFTI_BACKWARD_SCALE       = 1.0/(M*N)'
         PRINT *
      END IF

      LENGTHS(1) = M
      LENGTHS(2) = N

!
!     Allocate dynamic arrays and put input data
!
      ALLOCATE(X_IN(M,N), X_EXP(M,N), STAT=STATUS)
      IF(STATUS/=0) GOTO 102

      X_IN = 0.0_8
      X_IN(1,1) = 1.0_8
      X_EXP = X_IN

!
!     Put input data and expected result
!
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *,'INPUT Global vector X, N columns)'
         CALL PRINT_DATA_2D_Z(X_IN, M, N, N)
      END IF

!
!     2. Allocate memory for the descriptor by calling DftiCreateDescriptorDM
!
      STATUS = DftiCreateDescriptorDM(MKL_COMM,DESC,DFTI_DOUBLE,        &
     &                                DFTI_COMPLEX,2,LENGTHS)
      IF (MPI_RANK .EQ. 0) PRINT *, 'CREATE=', STATUS
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 102
      END IF

!
!     3. Obtain some values of configuration parameters by calls to
!        DftiGetValueDM
!
      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_SIZE,SIZE)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GET=',STATUS,',SIZE=',SIZE
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_NX,NX)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GET=',STATUS,',NX=',NX
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_X_START,START_X)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GET=',STATUS,',START_X=',START_X
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_OUT_NX,NX_OUT)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GET=',STATUS,',NX_OUT=',NX_OUT
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_OUT_X_START,START_X_OUT)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GET=',STATUS,',START_X_OUT=',START_X_OUT
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

      ALLOCATE(LOCAL(SIZE), WORK(SIZE), STAT=STATUS)
      IF(STATUS/=0) GOTO 102

!
!     4. Specify a value(s) of configuration parameters by a call(s) to
!        DftiSetValueDM
!
      STATUS = DftiSetValueDM(DESC,CDFT_WORKSPACE,WORK)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'SET=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     5. Perform initialization that facilitates DFT computation by a call to
!        DftiCommitDescriptorDM
!
      STATUS = DftiCommitDescriptorDM(DESC)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'COMMIT=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     6. Create arrays for local parts of input and output data
!        (if it is needed) and fill the local part of input data with
!        values (for more information, see Distributing Data among Processes)
!
      ROOTRANK=0
      STATUS = MKL_CDFT_SCATTERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,2,      &
     &                                LENGTHS,X_IN,NX,START_X,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'Scatter=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     7. Compute the transform by calling
!        DftiComputeForwardDM or DftiComputeBackward
!        (Compute Forward transform)
!
      STATUS = DftiComputeForwardDM(DESC,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'ComputeForward=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Gather data among processors
!
      STATUS = MKL_CDFT_GATHERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,         &
     &                               2,LENGTHS,X_IN,NX,START_X,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GATHER=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Set Scale number for Backward transform
!
      SCALE = 1.0_8/(N*M)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0))THEN
         PRINT *,'DFTI_BACKWARD_SCALE  = 1/(M*N)'
      END IF
      STATUS = DftiSetValueDM(DESC,DFTI_BACKWARD_SCALE,SCALE)
      IF (MPI_RANK .EQ. 0) PRINT *, 'Set=',STATUS
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Commit DftiDM descriptor
!
      STATUS = DftiCommitDescriptorDM(DESC)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'COMMIT=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Spread data among processors
!
      STATUS = MKL_CDFT_SCATTERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,        &
     &                                2,LENGTHS,X_IN,NX,START_X,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'Scatter=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Compute Backward transform
!
      STATUS = DftiComputeBackwardDM(DESC,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'ComputeBackward=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     8. Gather local output data into the global array using MPI functions
!        or use them otherwise
!
      STATUS = MKL_CDFT_GATHERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,         &
     &                               2,LENGTHS,X_IN,NX,START_X,LOCAL)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'GATHER=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
         GOTO 101
      END IF

!
!     Print data after DftiComputeBackwardDM; data assembled together
!
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *,'Backward OUTPUT vector X, N columns)'
         CALL PRINT_DATA_2D_Z(X_IN, M, N, N)
      END IF

!
!     Check result
!
      MAXERR = MAXVAL(ABS(X_IN - X_EXP))
      IF (ACCURACY_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *
         PRINT '(" Accuracy = ",G15.6,/)', MAXERR
      END IF

      IF (MAXERR .LT. EPS) THEN
         IF (MPI_RANK .EQ. 0) THEN
            PRINT *, 'TEST PASSED'
         END IF
      ELSE
         IF (MPI_RANK .EQ. 0) THEN
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
      END IF

!
!     9. Release memory allocated for a descriptor by a call to
!        DftiFreeDescriptorDM
!        (Free DftiDM descriptor)
!
101   STATUS = DftiFreeDescriptorDM(DESC)
      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
         PRINT *, 'FreeDescriptor=',STATUS
      END IF
      IF (STATUS .NE. DFTI_NO_ERROR) THEN
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(STATUS)
            PRINT *, 'TEST FAILED'
         END IF
         FAILURE = .TRUE.
      END IF

!
!     Free memory for dynamic arrays
!
102   DEALLOCATE(X_IN, X_EXP, LOCAL, WORK)

!
!     Finalize MPI
!
      CALL MPI_FINALIZE(MPI_ERR)

      IF (FAILURE) THEN
         STOP 1
      END IF

      IF (MPI_RANK .EQ. 0) THEN
         PRINT *, 'END OF TEST'
      END IF

      END PROGRAM
