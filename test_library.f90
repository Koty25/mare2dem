 program test_library 
!
! utility for testing the mare2dem forward library that can be linked from external codes
!
! build notes:  in mare2dem_source:  "make CLUSTER=macpro" to make MARE2DEM and mare2dem_lib.o
!                                      then "make CLUSTER=macpro test_library" to make this test program
!                                      then run in mare2dem_examples/Demo/inversion_MT with "mpirun -n 4 ../../../mare2dem_source/test_library"
!
! Test the following subroutines in mare2dem_lib_interface.f90:
!
! 1) mare2dem_load_files()         - Loads a resistivity file into mare2dem's memory,
!                                    returns number of free regions and number of parameters per region (e.g., for isotropic = 2, tiz anisotropy = 2, etc)
!                                    Should only be called on rank 0 processor (for a given local communicator)
!
! 2) mare2dem_get_params()         - Returns the log10 resistivity parameters in each free parameter region (length = #regions x #params_per_region)
! 
! 3) mare2dem_forward()            - Computes forward model response for model loaded with mare2dem_load_files. You must pass in log10resitivity
!                                     array to use for the free parameters, and the number of processors per MPI team.
!
!
! For Julia interface, you should be albe to use 1 to load the files and get the number of free parameter regions. 
! Then you can allocate memory for them in julia an initialize the arrays, then use 2 to return the actually arrays for the 
! free parameter log10 resistivity and the free parameter region centroids. So 1 returns the number of array elements, 2 returns
! the actual arrays. Then call 3 whenever you need a forward call to get the chi^2.
!
!
 
	 
	use mare2dem_lib_interface
	use mpi
	
	use, intrinsic :: iso_fortran_env, only: int64, real64
	
	implicit none
	 
    integer         :: rank,ierr,nproc,i
    integer(int64)  :: nproc_per_team,nFreeRegions,nParamsPerRegion
    character(256)  :: testResistivityFile, testSaveFile
    real(real64)    :: chiSquaredMisfit
    
    real(real64),dimension(:),allocatable    :: nLog10RhoFree  
    real(real64),dimension(:,:),allocatable  :: freeCentroids   ! nFreeREgions x 2 free region centroids (y,z)
	
	!real(real64),dimension(:), allocatable :: freeResistivity
	
	
! KWK debug: hard wired for simple testing of a model in the Demo folder (MT or CSEM):	
    testResistivityFile ='Demo.2' ! note mare2dem expects this without the .resistivity extension
    
!
! Initialize MPI:   
!
    call mpi_init( ierr )
    
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )	 ! here mcomm defined to be mpi_world in mare2dem_mpi, but could be passed in instead
    call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )  
    
 
    
!
! Test 1: load MARE2DEM input files into memory only on rank 0 processor (mare2dem_forward() will dole out data to worker processors)
!	 
    if (rank == 0)	then
        call mare2dem_load_files(testResistivityFile,nFreeRegions,nParamsPerRegion)    ! only need to load files on manager processor
	    
	    !write(*,*) 'nFreeRegions: ',nFreeRegions
	    !write(*,*) 'nParamsPerRegion: ',nParamsPerRegion 
	
    endif

!
! Test 2: get free parameters and parameter centroids:
!
    if (rank == 0) then
        call mare2dem_get_params(nLog10RhoFree,freeCentroids)   
         write(*,*)    '10**nLog10RhoFree(1:10):',    10**nLog10RhoFree(1:10) 	
         
        
             
	    ! also test routine that will save a .resistivity file:
	    testSaveFile = 'testSave'
	    call mare2dem_save_resistivity(testSaveFile,nLog10RhoFree,0d0) 
	    
	    
    endif
!     do i=1,size(freeCentroids,1)
!         write(*,*) freeCentroids(i,1:2)
!     enddo
   
  
!
! Test 3: call the forward kernel 
!	
    nproc_per_team = nproc  ! for testing, use only one team consisting of all processors. For Dan's PT code, each PT chain will have 
                            ! some number of processors set asside for MARE2DEM forward calls (like nproc_per_team = 4 for 6 PT chains on a 24-core habanero node)
         
    call mare2dem_forward(nproc_per_team,nLog10RhoFree,chiSquaredMisfit) 
 
!
! Notes: Steps 1 and 2 above could be called by the local rank 0 processor on startup of the Julia code.
!        Then Test 3 will be called for each forward evaluation and is called by all processors. nproc_per_team is used to determine 
!        the local MPI communicator team. nLog10RhoFree only needs to be defined on the team's rank=0 process since in MARE2DEM it 
!        passed the resistivity array to the workers.
!    
    
 
!
! Finalize MPI
!    
 
    call mpi_finalize(ierr)
 
    
    !call mare2dem_deallocate ! only need to this at very end of code, not require afer each fwd call
     
end program test_library

 
!    
	