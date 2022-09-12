!-----------------------------------------------------------------------------------------------------------------------------------
    module mare2dem_lib_interface
!
! Interface for calling MARE2DEM from external codes for forward computations
!
! General sequence of calling from external routine:
!    (1) mare2dem_load_files                -  loads mare2dem format data and model files
!    (2) mare2dem_forward(freeResistivity)  - computes forward response and misfit for freeResistivity array of free parameter vals
!
! Only need to do step (1) once at the start. Then step (2) can be done any number of times for varying freeResistivity arrays.
!    
! Optional: 
!    mare2dem_get_params()        - can be called after (1) to get y,z of free parameter centroids. 
!    mare2dem_save_resistivity()  - save free parameters into a mare2dem .resistivity file 
!
!-----------------------------------------------------------------------------------------------------------------------------------

    use occam  
    use mare2d_mpi_definitions    ! for MPI communicator tags etc
    use mare2dem_global           ! for setting defaults and storing worker status array
    use mare2dem_io
 

    use, intrinsic :: iso_fortran_env, only: int64, real64
    
    implicit none
 
    public :: mare2dem_load_files_Test, mare2dem_load_files, mare2dem_get_params, mare2dem_forward, mare2dem_save_resistivity

    contains 
    
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------mare2dem_load_files_Test   
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_load_files_Test(nFreeRegions,nParamsPerRegion) bind(C, name="mare2dem_load_files_Test")

    ! test version of routine below
    
    integer(int64), intent(out)    :: nFreeRegions,nParamsPerRegion
 
!    
! Copy to model shared variable:
!
    resistivityFile = '../mare2dem_examples/Demo/inversion_MT/Demo.2'
 
!
! Load MARE2DEM .resistivity, .poly, .settings and data files:
!
    call displayBanner    
    call readModel 
    call readData
    if ( nFree < 1000 ) lUseInversionMeshCoarsening = .false. ! don't use coarsening if small model

    nParamsPerRegion = nRhoPerRegion
    nFreeRegions     = nFree/nRhoPerRegion                             

    end subroutine mare2dem_load_files_Test

!-----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------mare2dem_load_files   
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_load_files(inputResistivityFile,nFreeRegions,nParamsPerRegion) bind(C, name="mare2dem_load_files")
!
! Loads the .resistivity, .poly, data and settings files into mare2dem memory
! Returns number of free regions and the number of resistivity params per region (so total free params = nFreeRegions*nParamsPerRegion)
!

    character(1), dimension(256)   :: inputResistivityFile            ! e.g. Demo.0 . note mare2dem wants this without the .resistivity extension
    integer(int64), intent(out)    :: nFreeRegions,nParamsPerRegion
 
!    
! Copy to model shared variable:
!
    resistivityFile = inputResistivityFile(1)(1:256) ! clunky way to deal with bind(C) only allowing character(1) string, 
    ! so we use array of 256 of them, and this inserts those into a character(256) string
    
    write(*,*) 'resistivityFile: ',trim(resistivityFile) !kwk debug
 
!
! Load MARE2DEM .resistivity, .poly, .settings and data files:
!
    call displayBanner    
    call readModel 
    call readData
    if ( nFree < 1000 ) lUseInversionMeshCoarsening = .false. ! don't use coarsening if small model

    nParamsPerRegion = nRhoPerRegion
    nFreeRegions     = nFree/nRhoPerRegion                             

    end subroutine mare2dem_load_files
    
!-----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------mare2dem_save_resistivity   
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_save_resistivity(sFileName,nLog10RhoFreeIn,misfit) bind(C, name="mare2dem_save_resistivity")
!
! Loads the .resistivity, .poly, data and settings files into mare2dem memory
! Returns number of free regions and the number of resistivity params per region (so total free params = nFreeRegions*nParamsPerRegion)
!
    character(1), dimension(256), intent(in)     :: sFileName 
    real(real64), dimension(nFree), intent(in)   :: nLog10RhoFreeIn 
    real(real64),                   intent(in)   :: misfit
 
!    
! Copy to model shared variable:
!
    outputFileRoot = sFileName(1)(1:256) ! clunky way to deal with bind(C) only allowing character(1) string, 
    ! so we use array of 256 of them, and this inserts those into a character(256) string
    
    write(*,*) 'mare2dem_save_resistivity will save to file: ',trim(outputFileRoot) ! kwk debug
    ! note outputFileRoot is global variable defined in maredem_common.f90 and used by subroutine writeResistivityFile() below
    
    call insertFreeParams(nLog10RhoFreeIn) ! in em2d.f90, inserts nLog10RhoFreeIn into rhoParams array used
                                           ! by writeResistivityFile
    !
    ! Save to file:
    !
    ! note that this routine writes to a file named: <outputFileRoot>.<nCurrentIter>.resistivity and we're setting nCurrentIter=0 below
    
    call writeResistivityFile(0,misfit,0d0,0d0,0d0)
    ! writeResistivityFile(nCurrentIter,misfit,misfitRequested,lagrange,roughness)
    
    end subroutine mare2dem_save_resistivity
       
!-----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------mare2dem_get_params
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_get_params(nLog10RhoFree,freeCentroids)  bind(C, name="mare2dem_get_params")

    use kx_io 
    
    real(real64), dimension(:),allocatable,intent(out)   :: nLog10RhoFree
    real(real64), dimension(:,:),allocatable,intent(out) :: freeCentroids
    
    integer :: i,j,ict
!
! Copy free parameters to output:
!
    allocate(nLog10RhoFree(nFree))
    ict = 0
    do i = 1,nRegions   
        do j = 1,nRhoPerRegion
            if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then ! note this is NOT setup for cole-cole ip 
                ict = ict + 1
                nLog10RhoFree(ict) = log10(rhoParams((i-1)*nRhoPerRegion+j)) 
            endif
        enddo
    enddo  
    
!
! Get free region centroids:
!    
    call getFreeRegionCentroids()
    allocate(freeCentroids(size(freeRegionCentroids,1),2))
    freeCentroids = freeRegionCentroids

    end subroutine mare2dem_get_params 
!-----------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------mare2dem_forward    
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine mare2dem_forward(nproc_per_team,nLog10RhoFreeIn,chiSquaredMisfit) bind(C, name="mare2dem_forward")
!
! Routine for computing MARE2DEM forward response that can be called from external libraries.
! You must first call mare2dem_load_files(inputResistivityFile) on the rank 0 processor 
!
! Inputs:   nproc_per_team   variable to use to create local mpi communicator in Fortran in same way as done in 
!                            external calling code, which may be running multiple forward models in parallel 
!                            on separate teams of processors. e.g.,Julia  PT chains running using nproc_per_team cores 
!                            for each PT chain. The calling routine needs to have teams set up identically 
!                            (we recreate the MPI communicator here so its available to mare2dem's fortran routines)
!           nLog10RhoFreeIn  1D array (size nFree) of free parameter log10 resistivities (log10 ohm-m)
!
! Outputs:  chiSquaredMisfit data fit to input model using nRhoFreeIn 
!

    integer(int64),                 intent(in)   :: nproc_per_team
    real(real64), dimension(nFree), intent(in)   :: nLog10RhoFreeIn
    real(real64),                   intent(out)  :: chiSquaredMisfit
    
    integer             :: rank,nproc,nteams,team, ierr
    
!
! Intel MKL specific: make sure one thread per MPI process so mkl solver calls don't oversubscribe
!
    call mkl_set_num_threads ( 1 )  

!
! Get rank and setup communicator team(s)
!
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
    call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )   
    if (nproc == 1) then
        write(*,*) ' '
        write(*,*) '!!!!! Error: MARE2DEM needs to be run using at least 2 processors !!!'
        write(*,*) '             You specified only one, sorry, stopping!'
        write(*,*) ' '
        call exitMARE2DEM()
    endif
    nworkers = nproc - 1   
    
    nteams = ceiling(dble(nworkers)/dble(nproc_per_team))
    team   = (rank-1)/nproc_per_team
    
    call mpi_comm_split(MPI_COMM_WORLD, team, rank, mcomm ,ierr) ! assigns this rank to be part of mcomm for team.
                                                                 ! so each team has different mcomm
    
    write(*,*) 'rank, team: ', rank,team
    
    call mpi_barrier(MPI_COMM_WORLD,ierr) 
       
    !
    ! Launch manager and worker controllers:
    !
    if (rank /= manager) then ! worker
        
        call mpi_worker_controller
        
    else ! manager
    
        !
        ! Initialize the worker status array:
        !
        allocate (lWorker_status(nworkers))      
        lWorker_status  = .true.      
        
        !
        ! Compute forward response:
        !
        linversion = .false. ! tell mare2dem to NOT compute Jacobians otherwise it'll run way more slowly 

        call computeFwd( linversion, nLog10RhoFreeIn )
        
        !  
        ! Compute chi-squared misfit
        !
        chiSquaredMisfit =  sum( ((d - dm)/sd)**2 ) 
        
        ! debug check:
        write(*,*) 'rms: ',sqrt(chiSquaredMisfit/size(d))
     
        !
        ! Tell workers to exit mpi_worker_controller:
        !
        call mpi_shutDownWorkers()
        
    endif
    
    !
    ! Broadcast chiSquaredMisfit to all processors on local communicator
    !
    call mpi_bcast(chiSquaredMisfit, 1, MPI_DOUBLE_PRECISION, manager, mcomm, ierr)
    
    !debug: write(*,*) 'rank, chiSquaredMisfit: ',rank, chiSquaredMisfit
    
    end subroutine mare2dem_forward    
    


    
end module mare2dem_lib_interface   
    