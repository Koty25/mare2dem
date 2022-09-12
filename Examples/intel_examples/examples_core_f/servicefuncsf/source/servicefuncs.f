!===============================================================================
! Copyright 2008-2020 Intel Corporation.
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
!      mklservicefunctions example program demonstrate a lot of
!      Intel(R) MKL service functions.
!*******************************************************************************

      program mklservicefunctions

      include "mkl_service.fi"  !! Interfaces for Intel(R) MKL service functions
      include "mkl_lapack.fi"   !! Here is the interface for DSECND

      real*8     freq, SECONDS_S, SECONDS_E

      REAL*8     A,B,C  !! DGEMM's arrays
      POINTER    (A_PTR,A(1)), (B_PTR,B(1)), (C_PTR,C(1))

      INTEGER    N, I

#ifdef _IA32
      integer*4 ALLOC_SIZE, NUM_ELEM, SIZE_ELEM
#else
      integer*8 ALLOC_SIZE, NUM_ELEM, SIZE_ELEM
#endif

      integer*4 NTHRS
      integer*8 ALLOCATED_BYTES
      integer*4 ALLOCATED_BUFFERS
      integer*4 ALIGNMENT

      character*198    buf

      N = 1000
      ALIGNMENT = 128

! Information and defaults

      SECONDS_S = DSECND()

      write(*,'(/A)') 'Intel(R) MKLservice functions example started'
      write(*,'(/A)') 'Intel(R) MKLrelease version:'

      CALL mkl_get_version_string(buf)
      write(*,'(a)') buf

      write(*,'(/A)') 'Information and defaults'

      freq = mkl_get_cpu_frequency()
      write(*,'(A,F8.4,A)') '    Current CPU frequency:', freq,'GHz'

      freq = mkl_get_max_cpu_frequency()
      write(*,'(A,F8.4,A)') '    Maximum CPU frequency:', freq,'GHz'

      freq = mkl_get_clocks_frequency()
      write(*,'(A,F8.4,A)') '    Frequency:', freq,'GHz'

      CALL aux_print_n_threads()

! Memory functions
      ALLOCATED_BYTES = MKL_PEAK_MEM_USAGE(MKL_PEAK_MEM_ENABLE)
      write(*,'(/A)') 'Memory functions'
      write(*,'(A)')  '    Allocate DGEMM''s arrays'
      NUM_ELEM = N*N
      SIZE_ELEM = 8
      ALLOC_SIZE = NUM_ELEM*SIZE_ELEM
      A_PTR = MKL_MALLOC(ALLOC_SIZE,ALIGNMENT)
      if(A_PTR.eq.0) then
         write(*,*) 
     $     'Cannot allocate pointer array of the length ',
     $     ALLOC_SIZE
         stop 1
      endif
      B_PTR = MKL_MALLOC(NUM_ELEM,ALIGNMENT) ! Allocates a memory buffer of smaller size to realloc it to appropriate size later
      if(B_PTR.eq.0) then
         write(*,*)
     $     'Cannot allocate pointer array of the length ',
     $     NUM_ELEM
         stop 1
      endif
      C_PTR = MKL_CALLOC(NUM_ELEM,SIZE_ELEM,ALIGNMENT)
      if(C_PTR.eq.0) then
         write(*,*)
     $     'Cannot allocate pointer array of the length ',
     $     NUM_ELEM
         stop 1
      endif

      B_PTR = MKL_REALLOC(B_PTR,ALLOC_SIZE)
      if(B_PTR.eq.0) then
         write(*,*)
     $     'Cannot reallocate pointer array of the length ',
     $     ALLOC_SIZE
         stop 1
      endif

      write(*,'(A)') '    CALL DGEMM'
      CALL aux_call_dgemm(N,A,B,C,0)
      write(*,'(A)') '    ...Done'

      ALLOCATED_BYTES = MKL_PEAK_MEM_USAGE(MKL_PEAK_MEM)
      write(*,'(A,I8,A)')
     $  '    Peak memory allocated by Intel(R) MKL allocator :',
     $  ALLOCATED_BYTES, ' bytes.'

      ALLOCATED_BYTES = MKL_MEM_STAT(ALLOCATED_BUFFERS)
      write(*,'(A,I8,A,I3,A)')
     $  '    Currently allocated by Intel(R) MKL allocator :',
     $  ALLOCATED_BYTES, ' bytes in',
     $  ALLOCATED_BUFFERS,' buffers.'

      CALL MKL_FREE_BUFFERS()
      ALLOCATED_BYTES = mkl_mem_stat(ALLOCATED_BUFFERS)
      write(*,'(A)') '    After MKL_FREE_BUFFERS was called:'
      write(*,'(A,I8,A,I3,A)')
     $  '    Currenlt allocated by Intel(R) MKL allocator :',
     $  ALLOCATED_BYTES, ' bytes in',
     $  ALLOCATED_BUFFERS,' buffers.'
      if (MKL_PEAK_MEM_USAGE(MKL_PEAK_MEM_DISABLE)<0) then
        write(*,'(/A)') 'Peak memory statistics is not disabled'
      end if

! Threading functions

      write(*,'(/A)') 'DGEMM & Intel(R) MKL threading'

! DGEMM on N thread (default)

      CALL aux_call_dgemm(N,A,B,C,1)

! DGEMM on 1 thread

      NTHRS = 1
      i = MKL_DOMAIN_SET_NUM_THREADS(NTHRS,MKL_DOMAIN_BLAS)
      CALL aux_call_dgemm(N,A,B,C,1)

! DGEMM on XXX threads. MKL_DYNAMIC

      write(*,'(/A)') 'MKL_DYNAMIC experiment'

      write(*,'(A)')  '    Force MKL_DYNAMIC=TRUE'
      CALL MKL_SET_DYNAMIC(MKL_DYNAMIC_TRUE)

      write(*,'(A)')  '        Set Intel(R) MKL BLAS-N-Threads to 64'
      NTHRS = 64
      i = MKL_DOMAIN_SET_NUM_THREADS(NTHRS,MKL_DOMAIN_BLAS)
      CALL aux_print_n_threads()
      CALL aux_call_dgemm(N,A,B,C,1)

      write(*,'(/A)') '    Switch off MKL_DYNAMIC facility'
      CALL MKL_SET_DYNAMIC(MKL_DYNAMIC_FALSE)

      write(*,'(A)')  '        Set Intel(R) MKL BLAS-N-Threads to 64'
      i = MKL_DOMAIN_SET_NUM_THREADS(NTHRS,MKL_DOMAIN_BLAS)
      CALL aux_print_n_threads()
      CALL aux_call_dgemm(N,A,B,C,1)

! Free all work arrays

      write(*,'(/A)') '    Free DGEMM''s arrays'
      CALL MKL_FREE(A_PTR)
      CALL MKL_FREE(B_PTR)
      CALL MKL_FREE(C_PTR)

      SECONDS_E = DSECND()
      SECONDS_E = (SECONDS_E-SECONDS_S)

      write(*,'(/A,F8.4,A)')
     $  'Intel(R) MKLservice functions example finished at',
     $   SECONDS_E,' seconds'

      STOP
      END

      SUBROUTINE aux_call_dgemm(N,A,B,C,PRINT_CLOCKS)

      include 'mkl_service.fi'
      include 'mkl_blas.fi'

      INTEGER    N, PRINT_CLOCKS, I, J
      REAL*8     A(N,N),B(N,N),C(N,N)
      REAL*8     ALPHA, BETA
      integer*8  DGEMM_S, DGEMM_E
      integer*8  DGEMM_CLOCKS, DGEMM_CB, DGEMM_CT, DGEMM_CU

      ALPHA = 1.1; BETA = -1.2
      DO I=1,N
      DO J=1,N
        A(I,J) = I
        B(I,J) = -I
      END DO
      END DO

      CALL MKL_GET_CPU_CLOCKS(DGEMM_S)
      CALL DGEMM('N','N',N,N,N,ALPHA,A,N,B,N,BETA,C,N);
      CALL mkl_get_cpu_clocks(DGEMM_E)
      DGEMM_CLOCKS = DGEMM_E-DGEMM_S

      if (PRINT_CLOCKS == 1) THEN
          DGEMM_CB = DGEMM_CLOCKS/1000000;
          DGEMM_CT = (DGEMM_CLOCKS-(DGEMM_CB*1000000))/1000
          DGEMM_CU = (DGEMM_CLOCKS-(DGEMM_CB*1000000)-(DGEMM_CT*1000))
          write(*,'(A,I4,A,I2,A,I8.3,A,I3.3,A,I3.3,A)')
     $      '    DGEMM (',N,') on ',
     $      MKL_Domain_Get_Max_Threads(MKL_DOMAIN_BLAS), ' thread(s):',
     $      DGEMM_CB,'.',DGEMM_CT,'.',DGEMM_CU,' clocks'
      end if

      RETURN
      END

      SUBROUTINE aux_print_n_threads()

      include 'mkl_service.fi'

      if (MKL_GET_DYNAMIC() == 0) then
          write(*,'(/A,I1)') '        MKL_DYNAMIC : FALSE'
      else
          write(*,'(/A,I1)') '        MKL_DYNAMIC : TRUE'
      end if

      write(*,'(A)')
     $  '   Intel(R) MKL Number of THREADS: ALL BLAS FFT PARDISO VML'
      write(*,'(A,5I5)') '                             ',
     $      MKL_Get_Max_Threads(),
     $      MKL_Domain_Get_Max_Threads(MKL_DOMAIN_BLAS),
     $      MKL_Domain_Get_Max_Threads(MKL_DOMAIN_FFT),
     $      MKL_Domain_Get_Max_Threads(MKL_DOMAIN_PARDISO),
     $      MKL_Domain_Get_Max_Threads(MKL_DOMAIN_VML)

      RETURN
      END
