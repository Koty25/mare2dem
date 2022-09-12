!===============================================================================
! Copyright 2009-2020 Intel Corporation.
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

!
!  Content:
!      Intel(R) Math Kernel Library ScaLAPACK Fortran example source file
!
!      Performance and validation of PDGEQRF function
!
!*******************************************************************************
!
!===============================================================================
!==== ScaLAPACK PDGEQRF ========================================================
!===============================================================================
!
! Example performs QR factorization of general matrix with validation and
! performance measurement
!
! Example also demonstrates BLACS routines usage
!
! It works with block-cyclic distribution on 2D process grid
!
! List of routines demonstrated in the example:
!
! PDGEQRF
! PDLANGE
! PDLAMCH
! PDLASET
! PDLACPY
! PDGEMM
! PDUNGQR
! PDMATGEN
! BLACS_GET
! BLACS_PINFO
! BLACS_GRIDINIT
! BLACS_GRIDINFO
! BLACS_GRIDEXIT
! BLACS_EXIT
! BLACS_BARRIER
! DESCINIT
! IGEBS2D
! IGEBR2D
! SGEBS2D
! SGEBR2D
! NUMROC
! INDXG2P
!
! In this example the matrix generated by PDMATGEN is factorized by PDGEQRF.
! For performance measurement only PDGEQRF execution time is taken into
! account.
!
!*******************************************************************************

      PROGRAM PDGEQRF_DRIVER
*     ==== Parameters ==================================================
      INTEGER            CSRC_, CTXT_, DLEN_, DTYPE_, LLD_, MB_, M_,
     $                   NB_, N_, RSRC_, NIN
      DOUBLE PRECISION   ZERO, ONE
      CHARACTER*8        VERSION
      PARAMETER          ( DLEN_ = 9, DTYPE_ = 1, CTXT_ = 2, M_ = 3,
     $                     N_ = 4, MB_ = 5, NB_ = 6, RSRC_ = 7,
     $                     CSRC_ = 8, LLD_ = 9, NIN = 1,
     $                     ZERO = 0.0D+0, VERSION = '1.0',
     $                      ONE = 1.0D+0 )
*
*     ==== Local scalars ===============================================
      INTEGER            IAM, NPROCS, NPROW, NPCOL, ICTXT, M, N, MB, NB,
     $                   INFO, IASEED, MP, NQ, MYROW, MYCOL, PASSED,
     $                   NOUT, NGRIDS, NN, NNB, I, J, K, TESTS_FAILED,
     $                   TESTS_SKIPPED, TESTS_PASSED, TESTS_WOCHK, LNCH,
     $                   ITERATIONS, IAROW, IACOL, MP0, NQ0, LWORK
      CHARACTER*1024     INFILE, OUTFILE, PROGNAM
      DOUBLE PRECISION   GFLOPS, EPS, FRESID, ANORM, DIFNORM, TIME
      REAL               THRESH
      LOGICAL            CHECK
*
*     ==== Local arrays ================================================
      INTEGER            IWORK( 5 ), DESCA( DLEN_ ), DESCQ( DLEN_ ),
     $                   DESCR( DLEN_ )
      DOUBLE PRECISION   ASK_LWORK( 1 )
*
*     ==== Dynamic arrays ==============================================
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: A(:)
      DOUBLE PRECISION, ALLOCATABLE :: Q(:)
      DOUBLE PRECISION, ALLOCATABLE :: R(:)
      DOUBLE PRECISION, ALLOCATABLE :: TAU(:)
      DOUBLE PRECISION, ALLOCATABLE :: TIME_ITER(:)
      INTEGER, ALLOCATABLE :: PVAL(:)
      INTEGER, ALLOCATABLE :: QVAL(:)
      INTEGER, ALLOCATABLE :: MVAL(:)
      INTEGER, ALLOCATABLE :: NVAL(:)
      INTEGER, ALLOCATABLE :: MBVAL(:)
      INTEGER, ALLOCATABLE :: NBVAL(:)
*
*     ==== External subroutines ========================================
      EXTERNAL BLACS_GET, BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,
     $         BLACS_GRIDEXIT, BLACS_EXIT, BLACS_BARRIER, PDMATGEN,
     $         PDGEQRF, DESCINIT, IGEBS2D, IGEBR2D, SGEBS2D, SGEBR2D, 
     $         PDLASET, PDLACPY, PDGEMM, PDORGQR, CALC_AVERAGE_TIME
*
*     ==== External functions ==========================================
      DOUBLE PRECISION   DSECND, PDLANGE, PDLAMCH
      INTEGER            NUMROC, INDXG2P
      EXTERNAL           DSECND, PDLANGE, PDLAMCH, NUMROC, INDXG2P
*
*     ==== Intrinsic functions =========================================
      INTRINSIC          MIN, MAX, DBLE
*
*     ==== Executable statements =======================================
*
*
      CALL BLACS_PINFO( IAM, NPROCS )
*
*     Get program name and in/out file names
      CALL GETARG(0, PROGNAM)
      CALL GETARG(1, INFILE)
      CALL GETARG(2, OUTFILE)
*     Check whether in/out file names are valid
      IF( INFILE.EQ.'' .OR. INFILE.EQ.' ' .OR.
     $    OUTFILE.EQ.'' .OR. OUTFILE.EQ.' ' )
     $   GOTO 50
*
*     Open input file, skip header and read NOUT
      IF( IAM.EQ.0 ) THEN
         OPEN( NIN, FILE = INFILE, STATUS = 'OLD', ACTION = 'READ' )
         READ( NIN, FMT = * )
         READ( NIN, FMT = * ) NOUT
      ENDIF
*
*     Open output file
      IF( NOUT.NE.0 .AND. NOUT.NE.6 .AND. IAM.EQ.0 )
     $   OPEN( NOUT, FILE = OUTFILE, STATUS = 'UNKNOWN',
     $         ACTION = 'WRITE' )
*
*     Initialize temporary process grid
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS)
*
*     ==== Read and distribute data ====================================
      IF( IAM.EQ.0 ) THEN
         READ( NIN, FMT = * ) THRESH
         READ( NIN, FMT = * ) ITERATIONS
*        Read task dimensions array
         READ( NIN, FMT = * ) NN
         ALLOCATE( MVAL( NN ) )
         ALLOCATE( NVAL( NN ) )
         READ( NIN, FMT = * ) ( MVAL( I ), I = 1, NN )
         READ( NIN, FMT = * ) ( NVAL( I ), I = 1, NN )
*        Read block dimensions array
         READ( NIN, FMT = * ) NNB
         ALLOCATE( MBVAL( NNB ) )
         ALLOCATE( NBVAL( NNB ) )
         READ( NIN, FMT = * ) ( MBVAL( I ), I = 1, NNB )
         READ( NIN, FMT = * ) ( NBVAL( I ), I = 1, NNB )
*        Read grid arrays
         READ( NIN, FMT = * ) NGRIDS
         ALLOCATE( PVAL( NGRIDS ) )
         ALLOCATE( QVAL( NGRIDS ) )
         READ( NIN, FMT = * ) ( PVAL( I ), I = 1, NGRIDS )
         READ( NIN, FMT = * ) ( QVAL( I ), I = 1, NGRIDS )
*        Pack data into array
         IWORK( 1 ) = NN
         IWORK( 2 ) = NNB
         IWORK( 4 ) = NGRIDS
         IWORK( 5 ) = ITERATIONS
*        Send data to all processes
         CALL SGEBS2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1 )
         CALL IGEBS2D( ICTXT, 'All', ' ', 5, 1, IWORK, 5 )
         CALL IGEBS2D( ICTXT, 'All', ' ', NGRIDS, 1, PVAL, NGRIDS )
         CALL IGEBS2D( ICTXT, 'All', ' ', NGRIDS, 1, QVAL, NGRIDS )
         CALL IGEBS2D( ICTXT, 'All', ' ', NN, 1, MVAL, NN )
         CALL IGEBS2D( ICTXT, 'All', ' ', NN, 1, NVAL, NN )
         CALL IGEBS2D( ICTXT, 'All', ' ', NNB, 1, MBVAL, NNB )
         CALL IGEBS2D( ICTXT, 'All', ' ', NNB, 1, NBVAL, NNB )
      ELSE
*        Receive and unpack data
         CALL SGEBR2D( ICTXT, 'All', ' ', 1, 1, THRESH, 1, 0, 0 )
         CALL IGEBR2D( ICTXT, 'All', ' ', 5, 1, IWORK, 5, 0, 0 )
         NN = IWORK( 1 )
         NNB = IWORK( 2 )
         NGRIDS = IWORK( 4 )
         ITERATIONS = IWORK( 5 )
*        Allocate memory and receive arrays
         ALLOCATE( PVAL( NGRIDS ) )
         ALLOCATE( QVAL( NGRIDS ) )
         ALLOCATE( MVAL( NN ) )
         ALLOCATE( NVAL( NN ) )
         ALLOCATE( MBVAL( NNB ) )
         ALLOCATE( NBVAL( NNB ) )
         CALL IGEBR2D( ICTXT, 'All', ' ', NGRIDS, 1, PVAL, NGRIDS, 0, 0)
         CALL IGEBR2D( ICTXT, 'All', ' ', NGRIDS, 1, QVAL, NGRIDS, 0, 0)
         CALL IGEBR2D( ICTXT, 'All', ' ', NN, 1, MVAL, NN, 0, 0)
         CALL IGEBR2D( ICTXT, 'All', ' ', NN, 1, NVAL, NN, 0, 0)
         CALL IGEBR2D( ICTXT, 'All', ' ', NNB, 1, MBVAL, NNB, 0, 0)
         CALL IGEBR2D( ICTXT, 'All', ' ', NNB, 1, NBVAL, NNB, 0, 0)
      ENDIF
*
      IF( IAM.EQ.0 )
     $   CLOSE( NIN )
*
*     Destroy temporary grid
      CALL BLACS_GRIDEXIT( ICTXT )
*
*     Write table header to output
      IF( IAM.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9999 ) 'PDGEQRF performance results'
         WRITE( NOUT, FMT = 9996 )
         WRITE( NOUT, FMT = 9995 ) THRESH, ITERATIONS
         WRITE( NOUT, FMT = * )
         WRITE( NOUT, FMT = 9997 )
      ENDIF
*
*     Set variable values
      IASEED = 100
      TESTS_WOCHK = 0
      TESTS_SKIPPED = 0
      TESTS_FAILED = 0
      TESTS_PASSED = 0
      CHECK = ( THRESH.GT.0.0E+0 )
*
*     ==== Loop over grid dimensions ===================================
      DO 10 I = 1, NGRIDS
*
         NPROW = PVAL( I )
         NPCOL = QVAL( I )
*
*        Make sure grid is correct
         IF( NPROW.LT.0 ) THEN
            TESTS_SKIPPED = TESTS_SKIPPED + 1
            GOTO 10
         ELSE IF( NPCOL.LT.0 ) THEN
            TESTS_SKIPPED = TESTS_SKIPPED + 1
            GOTO 10
         ELSE IF( NPROW*NPCOL.GT.NPROCS ) THEN
            TESTS_SKIPPED = TESTS_SKIPPED + 1
            GOTO 10
         ENDIF
*
*        Initialize working process grid
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*        If my process out of grid, go to bottom
         IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $      GOTO 10
*
*        ==== Loop over task dimensions ================================
         DO 20 J = 1, NN
*
            M = MVAL( J )
            N = NVAL( J )
*
*           Make sure M and N are correct
            IF( M.LT.1 ) THEN
               TESTS_SKIPPED = TESTS_SKIPPED + 1
               GOTO 10
            ELSE IF( N.LT.1) THEN
               TESTS_SKIPPED = TESTS_SKIPPED + 1
               GOTO 10
            ENDIF
*
*           ==== Loop over block dimensions ============================
            DO 30 K = 1, NNB
*
               MB = MBVAL( K )
               NB = NBVAL( K )
*
*              Make sure MB and NB are correct
               IF( MB.LT.1 ) THEN
                  TESTS_SKIPPED = TESTS_SKIPPED + 1
                  GOTO 10
               ELSE IF( NB.LT.1 ) THEN
                  TESTS_SKIPPED = TESTS_SKIPPED + 1
                  GOTO 10
               ENDIF
*
*              ==== Factorization part =================================
*
*              Calculate precise block dimensions
               MP = NUMROC( M, MB, MYROW, 0, NPROW )
               NQ = NUMROC( N, NB, MYCOL, 0, NPCOL )
*
*              Initialize array descriptors
               CALL DESCINIT( DESCA, M, N, MB, NB, 0, 0, ICTXT,
     $                        MAX(1, MP), INFO )
               CALL DESCINIT( DESCQ, M, N, MB, NB, 0, 0, ICTXT,
     $                        MAX(1, MP), INFO )
               CALL DESCINIT( DESCR, M, N, MB, NB, 0, 0, ICTXT,
     $                        MAX(1, MP), INFO )
*
*
               ALLOCATE( TIME_ITER( ITERATIONS ) )
*
*              Loop for measurement stability, outer iterations
               DO 40 LNCH = 1, ITERATIONS
*
*                 Allocate memory for matrix A
                  ALLOCATE( A( MP*NQ ) )
*
*                 Generate matrix A
                  CALL PDMATGEN( ICTXT, 'No', 'No', DESCA( M_ ), 
     $                           DESCA( N_ ), DESCA( MB_ ),
     $                           DESCA( NB_ ), A, DESCA( LLD_ ),
     $                           DESCA( RSRC_ ), DESCA( CSRC_ ), IASEED,
     $                           0, MP, 0, NQ, MYROW, MYCOL, NPROW,
     $                           NPCOL )
*
                  CALL BLACS_BARRIER( ICTXT, 'All' )
*
*                 Allocate array for reflection multipliers
*                  ALLOCATE( TAU( MIN( M, N ) ) )
                  ALLOCATE( TAU( NUMROC( MIN( M, N ), NB, MYCOL, 0,
     $                      NPCOL ) ) )
*
*                 Ask necessary size of working array for factorization
                  CALL PDGEQRF( M, N, A, 1, 1, DESCA, TAU, ASK_LWORK,-1,
     $                          INFO )
                  LWORK = INT(ASK_LWORK(1))
                  ALLOCATE( WORK( LWORK ) )
*                 Factor matrix A
                  TIME = DSECND()
                  CALL PDGEQRF( M, N, A, 1, 1, DESCA, TAU, WORK,
     $                          LWORK, INFO )
                  CALL BLACS_BARRIER( ICTXT, 'All' )
                  TIME = DSECND() - TIME
                  DEALLOCATE( WORK )
*
                  TIME_ITER( LNCH ) = TIME
*
*                 ==== Checking part ===================================
*
                  IF( CHECK ) THEN
*
*                    Allocate memory for temporary matrices
                     ALLOCATE( Q( MP*NQ ) )
                     ALLOCATE( R( MP*NQ ) )
*
*                    Copy R and Q into different matrices
                     CALL PDLASET( 'L', DESCR( M_ ), DESCR( N_ ), ZERO,
     $                             ZERO, R, 1, 1, DESCR )
                     CALL PDLACPY( 'U', DESCA( M_ ), DESCA( N_ ), A, 1,
     $                             1, DESCA, R, 1, 1, DESCR )
                     CALL PDLACPY( 'L', DESCA( M_ ), DESCA( N_ ), A, 1,
     $                             1, DESCA, Q, 1, 1, DESCQ )
*
*                    Generate matrix Q from reflections                    
                     CALL PDORGQR( M, MIN( M, N ), MIN( M, N ), Q, 1, 1,
     $                             DESCQ, TAU, ASK_LWORK, -1, INFO )
                     LWORK = INT(ASK_LWORK( 1 ))
                     ALLOCATE( WORK( LWORK ) )
                     CALL PDORGQR( M, MIN( M, N ), MIN( M, N ), Q, 1, 1,
     $                             DESCQ, TAU, WORK, LWORK, INFO )
                     DEALLOCATE( WORK )
*
*                    Regenerate matrix A
                     CALL PDMATGEN( ICTXT, 'No', 'No', DESCA( M_ ),
     $                              DESCA( N_ ), DESCA( MB_ ),
     $                              DESCA( NB_ ), A, DESCA( LLD_ ),
     $                              DESCA( RSRC_ ), DESCA( CSRC_ ),
     $                              IASEED, 0, MP, 0, NQ, MYROW, MYCOL,
     $                              NPROW, NPCOL )
*
*                    Calculate necessary size of working array
                     IAROW = INDXG2P( 1, DESCA( MB_ ), MYROW, DESCA(
     $                                RSRC_ ), NPROW )
                     IACOL = INDXG2P( 1, DESCA( NB_ ), MYCOL, DESCA(
     $                                CSRC_ ), NPCOL )
                     MP0 = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYROW,
     $                             IAROW, NPROW )
                     NQ0 = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYCOL,
     $                             IACOL, NPCOL )
*                    For I-norm LWORK = MP0, for 1-norm LWORK = NQ0
                     LWORK = MP0
*
*                    Calculate norm of the matrix A ( 1 or infinity)
                     ALLOCATE( WORK( LWORK ) )
                     ANORM = PDLANGE( 'I', M, N, A, 1, 1, DESCA, WORK )
*
*                    Compute matrix Q*R-A -> A
                     CALL PDGEMM( 'N', 'N', M, N, MIN( M, N ), ONE,
     $                            Q, 1, 1, DESCQ, R, 1, 1, DESCR,
     $                            -ONE, A, 1, 1, DESCA )
*                    Deallocate memory for temporary matrices
                     DEALLOCATE( Q )
                     DEALLOCATE( R )
*
*                    Calculate norm of the matrix Q*R-A ( 1 or infinity)
                     DIFNORM = PDLANGE( 'I', M, N, A, 1, 1, DESCA, WORK)
                     DEALLOCATE( WORK )      
*
*                    Calculate machine epsilon
                     EPS = PDLAMCH( ICTXT, 'e' )
*                    Calculate FRESID = ||Q*R-A||/(N*||A||*eps)
                     FRESID = DIFNORM / ( DBLE( N )*EPS*ANORM )
*
*                    Check, whether FRESID passed the THRESHold
                     PASSED = 1
                     IF( FRESID.GT.THRESH )
     $                  PASSED = 0
                     TESTS_PASSED = TESTS_PASSED + PASSED
                     TESTS_FAILED = TESTS_FAILED + ( 1 - PASSED )
                  ELSE
                     FRESID = 0.0D+0
                     PASSED = 1
                     TESTS_WOCHK = TESTS_WOCHK + 1
                  ENDIF
*
*                 Deallocate memory for matrix A, reflection multipliers
                  DEALLOCATE( A )
                  DEALLOCATE( TAU )
*
   40          CONTINUE
*
*              Calculate average time and free memory
               CALL CALC_AVERAGE_TIME( TIME_ITER, ITERATIONS, TIME )
               DEALLOCATE( TIME_ITER )
*
*              Floating point operations per second
               IF( TIME.GT.0 ) THEN
                  GFLOPS = 2.0D+0 * ( ( DBLE( MIN( M, N ) )**2 *
     $                     DBLE( 3*MAX( M, N ) - MIN( M, N ) ) )
     $                     / 3.0D+0 ) / TIME*1.0D-9
               ELSE
                  PRINT *,'Invalid timer'
                  IF( IAM.EQ.0 ) THEN
                     STOP 3
                  ELSE
                     STOP
                  ENDIF
               ENDIF
*
*              Write results into output file
               IF( IAM.EQ.0 )
     $            WRITE( NOUT, FMT = 9998 )
     $               M, N, MB, NB, NPROW, NPCOL, TIME, GFLOPS,
     $               PASSED, FRESID
   30       CONTINUE
   20    CONTINUE
*
*        Close process grid
         CALL BLACS_GRIDEXIT( ICTXT )
   10 CONTINUE
*
*     Deallocate grid, task dimensions and block dimensions arrays
      DEALLOCATE( PVAL )
      DEALLOCATE( QVAL )
      DEALLOCATE( MVAL )
      DEALLOCATE( NVAL )
      DEALLOCATE( MBVAL )
      DEALLOCATE( NBVAL )
*
*     Close output file
      IF( NOUT.NE.6 .AND. NOUT.NE.0 .AND. IAM.EQ.0 )
     $   CLOSE ( NOUT )
*
*     ..Destroy grid..
      CALL BLACS_EXIT(0)
*
 9999 FORMAT( 1X, 80A )
 9998 FORMAT( 1X, I7,', ', I7,', ', I5,', ', I5,', ', I3,', ', I3,', ',
     $        1pe11.3, ', ', 1pe11.3,', ', I6,', ', 1pe11.3 )
 9997 FORMAT
     $( 1X,'      M,       N,    MB,    NB,   P,   Q,   FACT_TIME,
     $      GFLOPS, PASSED,     RESIDUE' )
 9996 FORMAT( 1X, '  Threshold, Time_Iterations')
 9995 FORMAT( 1X, 1pe11.3, ', ', I15 )
*
      GOTO 600
   50 CONTINUE
      PRINT *, 'Usage: ', TRIM( PROGNAM ), ' <input file> <output file>'
      IF( IAM.EQ.0 ) THEN
         STOP 2
      ELSE
         STOP
      ENDIF
  600 CONTINUE
*
      IF( IAM.EQ.0 ) THEN
         IF( TESTS_FAILED.EQ.0 ) THEN
            STOP
         ELSE
            STOP 1
         ENDIF
      ELSE
         STOP
      ENDIF
*     ==================================================================
*     ==== END OF PDGEQRF_DRIVER =======================================
*     ==================================================================
      END