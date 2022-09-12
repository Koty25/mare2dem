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
!       Forward-Backward 2D real transform for double precision data inplace.
!
!*****************************************************************************
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_REAL                     (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE                   (obligatory)
!           DFTI_DIMENSION      = 2                             (obligatory)
!           DFTI_LENGTHS        = {M,N}                         (obligatory)
!           DFTI_FORWARD_SCALE  = 1.0                           (default)
!           DFTI_BACKWARD_SCALE = 1.0/(M*N)                     (default=1.0)
!
!*****************************************************************************

      PROGRAM DM_REAL_2D_DOUBLE_EX1

      USE MKL_CDFT
      USE CDFT_EXAMPLE_SUPPORT

      USE, INTRINSIC :: ISO_C_BINDING

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      COMPLEX(8), ALLOCATABLE, TARGET :: LOCAL(:), WORK(:)
      REAL(8), POINTER :: X_IN(:,:)
      COMPLEX(8), POINTER :: X_OUT(:,:)
      INTEGER   ELEMENTSIZE
      PARAMETER (ELEMENTSIZE = 16)

      TYPE(DFTI_DESCRIPTOR_DM), POINTER :: DESC

      INTEGER   NX,NX_OUT,START_X,START_X_OUT,SIZE,I,J

      INTEGER   STATUS
      REAL(8)   SCALE
      INTEGER   LENGTHS(2)

      REAL(8)   MAXERR, G_MAXERR
      REAL(8), PARAMETER :: EPS = DOUBLE_EPS

      INTEGER(MPI_KIND) M, N
      INTEGER M_PADDED
      INTEGER(MPI_KIND) MPI_ERR
      INTEGER(MPI_KIND) MPI_NPROC
      INTEGER(MPI_KIND) MPI_RANK
      INTEGER(MPI_KIND) COMM
      INTEGER(4)        MKL_COMM

      CHARACTER*1024    INFILE

      LOGICAL FAILURE

      FAILURE = .FALSE.

!
!     Initiate MPI by calling MPI_Init (Perform MPI initialization)
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
         PRINT *, 'DM_REAL_2D_DOUBLE_EX1'
         PRINT *
         PRINT *, 'Forward-Backward 2D real-to-complex transform '//    &
     &                          'for double precision data inplace'
         PRINT *
         PRINT *, 'Configuration parameters:'
         PRINT *
         PRINT *, 'DFTI_FORWARD_DOMAIN       = DFTI_REAL'
         PRINT *, 'DFTI_PRECISION            = DFTI_DOUBLE '
         PRINT *, 'DFTI_DIMENSION            = 2'
         PRINT '(" DFTI_LENGTHS              = {",I3,",",I3,"}")', M, N
         PRINT *, 'DFTI_FORWARD_SCALE        = 1.0 '
         PRINT *, 'DFTI_BACKWARD_SCALE       = 1.0/(M*N)'
         PRINT *
      END IF

      LENGTHS(1) = M
      LENGTHS(2) = N

      M_PADDED = M/2 + 1
!
!     Allocate memory for the descriptor by calling DftiCreateDescriptorDM
!
      STATUS = DftiCreateDescriptorDM(MKL_COMM,DESC,DFTI_DOUBLE,        &
     &                                DFTI_REAL,2,LENGTHS)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Descriptor for forward created'

!
!     Obtain some values of configuration parameters by calls to
!        DftiGetValueDM
!
      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_SIZE,SIZE)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      PRINT '(" Process ",I3,": size = ",I6,"")', MPI_RANK,SIZE

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_NX,NX)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      PRINT '(" Process ",I3,": NX = ",I6,"")', MPI_RANK,NX

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_X_START,START_X)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      PRINT '(" Process ",I3,": X start = ",I6,"")', MPI_RANK,START_X

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_OUT_NX,NX_OUT)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      PRINT '(" Process ",I3,": out NX = ",I6,"")', MPI_RANK,NX_OUT

      STATUS = DftiGetValueDM(DESC,CDFT_LOCAL_OUT_X_START,START_X_OUT)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      PRINT '(" Process ",I3,": out X start= ",I6,"")', MPI_RANK,       &
     &                                                  START_X_OUT

!
!     Allocate dynamic arrays and put input data
!
      ALLOCATE(LOCAL(SIZE), WORK(SIZE), STAT=STATUS)
      IF(STATUS/=0) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Memory for buffers allocated'

      CALL C_F_POINTER (C_LOC(LOCAL), X_IN, [2*M_PADDED,NX])
      CALL C_F_POINTER (C_LOC(LOCAL), X_OUT, [M_PADDED,NX_OUT])

!
!     Specify a value(s) of configuration parameters by a call(s) to
!        DftiSetValueDM
!
      STATUS = DftiSetValueDM(DESC,CDFT_WORKSPACE,WORK)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Work buffer set'

!
!     Perform initialization that facilitates DFT computation by a call to
!        DftiCommitDescriptorDM
!
      STATUS = DftiCommitDescriptorDM(DESC)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Descriptor for forward committed'

!
!     fill the local part of input data with values
!
      DO I = 1, M
         DO J = 1, NX
            X_IN(I,J) = INPUT_DATA_D(I,START_X+J-1)
         END DO
      END DO

!
!     Compute the transform by calling
!        DftiComputeForwardDM or DftiComputeBackward
!        (Compute Forward transform)
!
      STATUS = DftiComputeForwardDM(DESC,LOCAL)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Forward transform computed'

!
!     Set Scale number for Backward transform
!
      SCALE = 1.0_8/(N*M)
      STATUS = DftiSetValueDM(DESC,DFTI_BACKWARD_SCALE,SCALE)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Backward scale set'

!
!     Commit DftiDM descriptor
!
      STATUS = DftiCommitDescriptorDM(DESC)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Descriptor for backward committed'

!
!     Compute Backward transform
!
      STATUS = DftiComputeBackwardDM(DESC,LOCAL)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Backward transform computed'

!
!     Check result
!
      MAXERR = 0.0
      DO I = 1, M
         DO J = 1, NX
            MAXERR = MAX(MAXERR,                                        &
     &                   ABS(X_IN(I,J)-INPUT_DATA_D(I,START_X+J-1)))
         END DO
      END DO
      CALL MPI_ALLREDUCE(MAXERR, G_MAXERR, 1, MPI_DOUBLE_PRECISION,     &
     &                   MPI_MAX, COMM, MPI_ERR)
      MAXERR = G_MAXERR
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
!     Release memory allocated for a descriptor by a call to
!        DftiFreeDescriptorDM
!        (Free DftiDM descriptor)
!
      STATUS = DftiFreeDescriptorDM(DESC)

      FAILURE = GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)
      IF (FAILURE) GOTO 102

      IF (MPI_RANK .EQ. 0) PRINT *, 'Descriptor freed'

!
!     Free memory for dynamic arrays
!
102   DEALLOCATE(LOCAL, WORK)

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
