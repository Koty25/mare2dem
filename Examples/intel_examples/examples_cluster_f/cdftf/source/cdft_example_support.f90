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
!       Intel(R) Math Kernel Library (Intel(R) MKL) Cluster DFT example
!       support functions (Fortran-interface)
!
!*****************************************************************************

      MODULE CDFT_EXAMPLE_SUPPORT

      IMPLICIT NONE

      REAL      SINGLE_EPS
      PARAMETER (SINGLE_EPS = 1.0E-5)

      REAL(8)   DOUBLE_EPS
      PARAMETER (DOUBLE_EPS = 1.0E-11)

      LOGICAL    LEGEND_PRINT
      PARAMETER (LEGEND_PRINT = .TRUE.)

      LOGICAL    ADVANCED_DATA_PRINT
      PARAMETER (ADVANCED_DATA_PRINT = .TRUE.)

      LOGICAL    ACCURACY_PRINT
      PARAMETER (ACCURACY_PRINT = .TRUE.)

      INTEGER   MPI_KIND
      PARAMETER (MPI_KIND = MPI_KIND_)

      CONTAINS

      SUBROUTINE Dfti_Example_Status_Print(STATUS)

      USE MKL_DFTI

      INTEGER STATUS
      CHARACTER(DFTI_MAX_MESSAGE_LENGTH) Error_Message

      Error_Message = DftiErrorMessage(STATUS)
      PRINT *, ' Error message: ', Error_Message
      PRINT *, ' Error status = ', STATUS

      END SUBROUTINE

!---------------
      SUBROUTINE PRINT_DATA_2D_Z(X, M, N, C)

      COMPLEX(8) X(:,:)
      INTEGER(MPI_KIND) M, N, C

      INTEGER I, J

      DO J = 1, N
         PRINT '(" Row ",I3,":"/" ")', J
         DO I = 1, M
            PRINT '("(",F8.3,",",F8.3,")")',REAL(X(I,J)),AIMAG(X(I,J))
            IF ((MOD(I,C) .EQ. 0) .AND. (I .NE. M)) PRINT '(" ")'
         END DO
      END DO
      PRINT *

      END SUBROUTINE

!---------------
      SUBROUTINE PRINT_DATA_2D_C(X, M, N, C)

      COMPLEX(4) X(:,:)
      INTEGER(MPI_KIND) M, N, C

      INTEGER I, J

      DO J = 1, N
         PRINT '(" Row ",I3,":"/" ")', J
         DO I = 1, M
            PRINT '("(",F8.3,",",F8.3,")")',REAL(X(I,J)),AIMAG(X(I,J))
            IF ((MOD(I,C) .EQ. 0) .AND. (I .NE. M)) PRINT '(" ")'
         END DO
      END DO
      PRINT *

      END SUBROUTINE

!---------------
      INTEGER FUNCTION MKL_CDFT_DATA_S(COMM,ROOTRANK,ELEMENTSIZE,DIM,   &
     &                                 LENGTHS,GLOBAL,NX,START_X,LOCAL, &
     &                                 FLAG)

      INCLUDE 'mpif.h'

!     MPI related integer should have the kind MPI expect
      INTEGER(MPI_KIND), ALLOCATABLE :: COUNTS(:),DISPLS(:),BUF(:)
      INTEGER(MPI_KIND) I,TMP(2),REQ,STAT(MPI_STATUS_SIZE)
      INTEGER(MPI_KIND) COMM,ROOTRANK,NPROC,NRANK,MPI_ERR

      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X,FLAG,STATUS
      COMPLEX(4) GLOBAL(*),LOCAL(*)
      INTENT(IN) COMM,ROOTRANK,ELEMENTSIZE,DIM,LENGTHS,NX,START_X,FLAG

      INTEGER FD

      CALL MPI_COMM_RANK(COMM,NRANK,MPI_ERR)
      IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
      IF (NRANK==ROOTRANK) THEN
         CALL MPI_COMM_SIZE(COMM,NPROC,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
         ALLOCATE(COUNTS(NPROC),DISPLS(NPROC),BUF(2*NPROC), STAT=STATUS)
         IF(STATUS/=0) GOTO 100
      END IF

      FD=1
      DO I=1,DIM-1
         FD=FD*LENGTHS(I)
      END DO

      TMP(1)=NX*FD*ELEMENTSIZE
      TMP(2)=(START_X-1)*FD

      CALL MPI_GATHER(TMP,2_MPI_KIND,MPI_INTEGER,BUF,2_MPI_KIND,        &
     &                MPI_INTEGER,ROOTRANK,COMM,MPI_ERR)
      IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

      IF (NRANK==ROOTRANK) THEN
         COUNTS=BUF(1:2*NPROC-1:2)
         DISPLS=BUF(2:2*NPROC:2)
      END IF

      IF (FLAG==0) THEN

         CALL MPI_IRECV(LOCAL,TMP(1),MPI_BYTE,ROOTRANK,123_MPI_KIND,    &
     &                  COMM,REQ,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

         IF (NRANK==ROOTRANK) THEN
            DO I=0,NPROC-1
               CALL MPI_SEND(GLOBAL(DISPLS(I+1)+1),COUNTS(I+1),MPI_BYTE,&
     &                       I,123_MPI_KIND,COMM,MPI_ERR)
               IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
            END DO
         END IF

         CALL MPI_WAIT(REQ,STAT,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
      ENDIF

      IF (FLAG==1) THEN
         CALL MPI_ISEND(LOCAL,TMP(1),MPI_BYTE,ROOTRANK,222,COMM,REQ,    &
     &                  MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

         IF (NRANK==ROOTRANK) THEN
            DO I=0,NPROC-1
               CALL MPI_RECV(GLOBAL(DISPLS(I+1)+1),COUNTS(I+1),MPI_BYTE,&
     &                       I,222_MPI_KIND,COMM,STAT,MPI_ERR)
               IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
            END DO
         END IF
         CALL MPI_WAIT(REQ,STAT,MPI_ERR)
      END IF

      IF (NRANK==ROOTRANK) DEALLOCATE(COUNTS,DISPLS,BUF)

100   MKL_CDFT_DATA_S=MPI_ERR

      END FUNCTION

!---------------
      INTEGER FUNCTION MKL_CDFT_DATA_D(COMM,ROOTRANK,ELEMENTSIZE,DIM,   &
     &                                 LENGTHS,GLOBAL,NX,START_X,LOCAL, &
     &                                 FLAG)

      INCLUDE 'mpif.h'

!     MPI related integer should have the kind MPI expect
      INTEGER(MPI_KIND), ALLOCATABLE :: COUNTS(:),DISPLS(:),BUF(:)
      INTEGER(MPI_KIND) I,TMP(2),REQ,STAT(MPI_STATUS_SIZE)
      INTEGER(MPI_KIND) COMM,ROOTRANK,NPROC,NRANK,MPI_ERR


      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X,FLAG,STATUS
      COMPLEX(8) GLOBAL(*),LOCAL(*)
      INTENT(IN) COMM,ROOTRANK,ELEMENTSIZE,DIM,LENGTHS,NX,START_X,FLAG

      INTEGER FD

      CALL MPI_COMM_RANK(COMM,NRANK,MPI_ERR)
      IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

      IF (NRANK==ROOTRANK) THEN
         CALL MPI_COMM_SIZE(COMM,NPROC,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
         ALLOCATE(COUNTS(NPROC),DISPLS(NPROC),BUF(2*NPROC), STAT=STATUS)
         IF(STATUS/=0) GOTO 100
      END IF

      FD=1
      DO I=1,DIM-1
         FD=FD*LENGTHS(I)
      END DO

      TMP(1)=NX*FD*ELEMENTSIZE
      TMP(2)=(START_X-1)*FD

      CALL MPI_GATHER(TMP,2_MPI_KIND,MPI_INTEGER,BUF,2_MPI_KIND,        &
     &                MPI_INTEGER,ROOTRANK,COMM,MPI_ERR)
      IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

      IF (NRANK==ROOTRANK) THEN
         COUNTS=BUF(1:2*NPROC-1:2)
         DISPLS=BUF(2:2*NPROC:2)
      END IF

      IF (FLAG==0) THEN

         CALL MPI_IRECV(LOCAL,TMP(1),MPI_BYTE,ROOTRANK,123_MPI_KIND,    &
     &                  COMM,REQ,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

         IF (NRANK==ROOTRANK) THEN
            DO I=0,NPROC-1
               CALL MPI_SEND(GLOBAL(DISPLS(I+1)+1),COUNTS(I+1),MPI_BYTE,&
     &                       I,123_MPI_KIND,COMM,MPI_ERR)
               IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
            END DO
         END IF

         CALL MPI_WAIT(REQ,STAT,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
      ENDIF

      IF (FLAG==1) THEN

         CALL MPI_ISEND(LOCAL,TMP(1),MPI_BYTE,ROOTRANK,222_MPI_KIND,    &
     &                  COMM,REQ,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

         IF (NRANK==ROOTRANK) THEN
            DO I=0,NPROC-1
               CALL MPI_RECV(GLOBAL(DISPLS(I+1)+1),COUNTS(I+1),MPI_BYTE,&
     &                       I,222_MPI_KIND,COMM,STAT,MPI_ERR)
               IF (MPI_ERR/=MPI_SUCCESS) GOTO 100
            END DO
         END IF
         CALL MPI_WAIT(REQ,STAT,MPI_ERR)
         IF (MPI_ERR/=MPI_SUCCESS) GOTO 100

      END IF

      IF (NRANK==ROOTRANK) DEALLOCATE(COUNTS,DISPLS,BUF)

100   MKL_CDFT_DATA_D=MPI_ERR

      END FUNCTION

!---------------
      INTEGER FUNCTION MKL_CDFT_SCATTERDATA_S(COMM,ROOTRANK,ELEMENTSIZE,&
     &                                        DIM,LENGTHS,GLOBAL_IN,NX, &
     &                                        START_X,LOCAL_IN)

      INTEGER (MPI_KIND) ROOTRANK, COMM
      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X
      COMPLEX(4) GLOBAL_IN(*),LOCAL_IN(*)

      MKL_CDFT_SCATTERDATA_S=MKL_CDFT_DATA_S(COMM,ROOTRANK,ELEMENTSIZE, &
     &                                       DIM,LENGTHS,GLOBAL_IN,NX,  &
     &                                       START_X,LOCAL_IN,0)

      END FUNCTION

!---------------
      INTEGER FUNCTION MKL_CDFT_SCATTERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,&
     &                                        DIM,LENGTHS,GLOBAL_IN,NX, &
     &                                        START_X,LOCAL_IN)

      INTEGER (MPI_KIND) ROOTRANK, COMM
      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X
      COMPLEX(8) GLOBAL_IN(*),LOCAL_IN(*)

      MKL_CDFT_SCATTERDATA_D=MKL_CDFT_DATA_D(COMM,ROOTRANK,ELEMENTSIZE, &
     &                                       DIM,LENGTHS,GLOBAL_IN,NX,  &
     &                                       START_X,LOCAL_IN,0)

      END FUNCTION

!---------------
      INTEGER FUNCTION MKL_CDFT_GATHERDATA_S(COMM,ROOTRANK,ELEMENTSIZE, &
     &                                       DIM,LENGTHS,GLOBAL_OUT,NX, &
     &                                       START_X,LOCAL_OUT)

      INTEGER (MPI_KIND) ROOTRANK, COMM
      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X
      COMPLEX(4) GLOBAL_OUT(*),LOCAL_OUT(*)

      MKL_CDFT_GATHERDATA_S=MKL_CDFT_DATA_S(COMM,ROOTRANK,ELEMENTSIZE,  &
     &                                      DIM,LENGTHS,GLOBAL_OUT,NX,  &
     &                                      START_X,LOCAL_OUT,1)

      END FUNCTION

!---------------
      INTEGER FUNCTION MKL_CDFT_GATHERDATA_D(COMM,ROOTRANK,ELEMENTSIZE, &
     &                                       DIM,LENGTHS,GLOBAL_OUT,NX, &
     &                                       START_X,LOCAL_OUT)

      INTEGER (MPI_KIND) ROOTRANK, COMM
      INTEGER ELEMENTSIZE,DIM,LENGTHS(*),NX,START_X
      COMPLEX(8) GLOBAL_OUT(*),LOCAL_OUT(*)

      MKL_CDFT_GATHERDATA_D=MKL_CDFT_DATA_D(COMM,ROOTRANK,ELEMENTSIZE,  &
     &                                      DIM,LENGTHS,GLOBAL_OUT,NX,  &
     &                                      START_X,LOCAL_OUT,1)

      END FUNCTION

!---------------
      LOGICAL FUNCTION GLOBAL_CDFT_STATUS(COMM, MPI_RANK, STATUS)

      USE MKL_CDFT
      INCLUDE 'mpif.h'

      INTEGER (MPI_KIND) COMM
      INTEGER (MPI_KIND) MPI_RANK
      INTEGER STATUS, GLOBAL_STATUS
      INTEGER(MPI_KIND) MY_STATUS, MIN_STATUS, MAX_STATUS
      INTEGER(MPI_KIND) MPI_ERR
      LOGICAL FAILURE

      FAILURE = .FALSE.
      MY_STATUS = STATUS

!     Make sure all processes returned DFTI_NO_ERROR
      CALL MPI_ALLREDUCE(MY_STATUS, MAX_STATUS, 1_MPI_KIND, MPI_INTEGER,&
     &                   MPI_MAX, COMM, MPI_ERR)
      CALL MPI_ALLREDUCE(MY_STATUS, MIN_STATUS, 1_MPI_KIND, MPI_INTEGER,&
     &                   MPI_MIN, COMM, MPI_ERR)
      GLOBAL_STATUS = DFTI_NO_ERROR
      IF (MIN_STATUS .NE. DFTI_NO_ERROR) GLOBAL_STATUS = MIN_STATUS
      IF (MAX_STATUS .NE. DFTI_NO_ERROR) GLOBAL_STATUS = MAX_STATUS
      IF (GLOBAL_STATUS .NE. DFTI_NO_ERROR) THEN
         FAILURE = .TRUE.
         IF (MPI_RANK .EQ. 0) THEN
            CALL Dfti_Example_Status_Print(GLOBAL_STATUS)
            PRINT *, 'TEST FAILED'
         END IF
      END IF

      GLOBAL_CDFT_STATUS = FAILURE

      END FUNCTION

!---------------
      REAL(4) FUNCTION INPUT_DATA_S(I, J)

      REAL(4) DATA_VAL

      INTEGER I, J

      DATA_VAL = 0.0;
      IF ((I .EQ. 1) .AND. (J .EQ. 1)) DATA_VAL = 1.0

      INPUT_DATA_S = DATA_VAL

      END FUNCTION

!---------------
      REAL(8) FUNCTION INPUT_DATA_D(I, J)

      REAL(8) DATA_VAL

      INTEGER I, J

      DATA_VAL = 0.0;
      IF ((I .EQ. 1) .AND. (J .EQ. 1)) DATA_VAL = 1.0

      INPUT_DATA_D = DATA_VAL

      END FUNCTION

      END MODULE
