    program runMARE2DEM
 
    use Occam  
    use mare2d_mpi_definitions    ! for MPI communicator tags etc
    use mare2dem_global           ! for setting defaults and storing worker status array
    use mare2dem_io
    use mare2dem_input_data_params

    implicit none

    logical         :: lFwdOnly             = .false.   ! use command line argument -F to set this to true. 
    logical         :: lSaveJacobian        = .false.   ! use command line argument -FJ to set this to true
    logical         :: lSaveSensitivity     = .false.   ! use command line argument -FS to set this to true
    logical         :: lrunMARE2DEM, lCompJ = .false.
    real(8)         :: timeOffset, time0 ! Timer temp
    integer         :: ierr, nproc, myID,iWorker !  Local MPI variables
    character(12)   :: ctemp   
 
!
! Initialize MPI:
!
    call mpi_init( ierr )
    
!
! Start the timer:
! 
    call get_time_offset(0d0, time0)
    
!
! Get each processor's rank:
!
    call mpi_comm_rank( mcomm, myID, ierr )
!
! For Intel math kernel library, set #threads to 1. Occam matrix routines will override this since they are done 
! after parallel MARE2DEM calls have finished: 
#if defined(__INTEL_COMPILER)  
    call mkl_set_num_threads ( 1 )  ! Intel MKL specific  
! 25 Feb 2013: note that threading is only used on small machines, MARE2DEM now calls ScaLAPACK
! if a large number of processors are being used. This is handled in module Occam.
#endif
   

!
! Get total number of processes and number of workers:
!
    call mpi_comm_size( MPI_COMM_WORLD, nProc, ierr )   
    if (nProc == 1) then
        write(*,*) ' '
        write(*,*) '!!!!! Error: MARE2DEM needs to be run using at least 2 processors !!!'
        write(*,*) '             You specified only one, sorry, stopping!'
        write(*,*) ' '
        call exitMARE2DEM()
    endif
    nWorkers = nproc - 1   


!
! Launch the worker controllers:
!
    if (myID /= manager) call mpi_worker_controller

!
! Read in input files using the manager process:
!
    if (myID == manager) then 

        if (lPrintBanner) call displayBanner()    

        !
        ! Allocate the worker status array:
        !
        allocate (lWorker_status(nworkers))      
        lWorker_status  = .true.  
    
        !
        ! Tell the user about the cluster:
        !
        write(*,*) ' ' 
        write(ctemp,'(i8)') nworkers
        write(*,fmt='(a,a,a)') ' MARE2DEM is using one manager node and ',trim(adjustl(ctemp)), ' compute nodes'
        write(*,*) ' '     

    
        !
        ! Get any command line arguments:
        !
        call getMARE2DEMCommandLineArguments(lFwdOnly,lSaveJacobian,lSaveSensitivity,currentIteration)    
    
        !
        ! Read in the model files (.resistivity, .poly,  .settings and .penalty):
        !
        call readModel 
    
        if ( nFree < 1000 ) lUseInversionMeshCoarsening = .false. ! don't use coarsening if small model
    
        if (nParams == 0) then
            lFwdOnly = .true. ! No free parameters to invert for...
            write(*,*) ' ' 
            write(*,*) ' No free parameters found in the input model, therefore I will only '
            write(*,*) ' compute the forward response of the input model. '
            write(*,*) ' ' 
            lSaveSensitivity = .false. ! can't compute these if no free parameters. Further, if these 
            lSaveJacobian    = .false. ! are true you will get scratch file error messages, so just don't even try...
        endif
        write(*,fmt='(a32,a)') 'Scratch folder:  ',trim(adjustl(scratchFolder))    
        write(*,*) ' '    
        
        ! 
        ! Read in the data file:
        ! 
        call readData
 
        !
        ! Distribute a few of the input arrays for parallel processing with ScaLAPACK:
        !       
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
        do iWorker = 1,nworkers
            call mpi_send( iflag_TellWorkerToBcast, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
        enddo
        call mpi_all_Bcast
#endif
!
! The manager node does everything below here while the workers have been sent off to listen for 
! forward modeling task assignments.
!
 
            
    !
    ! Forward Call:
    !
    if (lFwdOnly) then
    
        !
        ! If Jacobian and/or sensitivity requested, initialize wj matrix:
        !
        if ((lSaveJacobian).or.(lSaveSensitivity)) lCompJ = .true.
        
        if (lCompJ) then
        
            allocate( wj(nd,nParams) , stat = ierr)
            if (ierr .ne. 0) then
                write(*,*) ' Out of memory.  Too many data and free parameters (', nd,nParams, ')'
                call exitMARE2DEM()
            endif
            wj = 0
            
            !
            ! Initialize the scratch folders
            ! 
            if ((nTxCSEM > 0).or.(nTxDC>0)) call mpi_manager_initializeScratch   
        
        
        endif
        
        !
        ! Compute forward response:
        ! 
        
        call computeFwd( lCompJ, transformToBound(pm,lowerBound,upperBound,lBoundMe)) 
        
        ! Show the misfit  
        write(*,*) ' '
        if (any(abs(d)>0)) then ! Only show misfit if input data is sensible:
            write(ctemp,'(f12.2)')   sqrt(sum( ((d-dm)/sd)**2 )/nd )  
            call print2col32('Forward Model Misfit:', ctemp,6)
        endif 
        call get_time_offset( time0,timeOffset)
        write(ctemp,'(f8.2)')  timeOffset
        call print2col32(' Total Time (s):', ctemp,6)
       
        ! Display any joint inversion misfits:
        call displayJointInvMisfits()
            
        ! Write out the model response:
        call writeResponse(currentIteration) 
        write(*,*) ' '
       
        ! Write out the model Jacobian matrix:
        if (lSaveJacobian) call writeJacobian(currentIteration+1) ! +1 since it subtracts one to make it for starting model 
       
        ! Write out the model Jacobian matrix:
        if (lSaveSensitivity) call writeSensitivity(currentIteration+1) ! +1 since it subtracts one 
                
        if (lCompJ) deallocate (wj)       
        
    !
    ! Inversion Call:
    !    
    else 
     
        !
        ! Initialize the scratch folders
        ! 
        if ((nTxCSEM > 0).or.(nTxDC>0)) call mpi_manager_initializeScratch   
        
        !
        ! Open a logfile:
        !
        call openOccamLogFile(outputFileRoot)
        
        !
        ! Initialize the inversion settings:
        !
        lrunMARE2DEM       = .true.
        
        !
        ! Run the main loop
        !
        do while (lrunMARE2DEM)
            
            currentIteration = currentIteration + 1 
            
            !
            ! Compute an Occam iteration:
            !
             call computeOccamIteration(lSaveJacobian,lSaveSensitivity) 
             
            !
            ! Write out the results:
            ! 
              if ( convergenceFlag /= 4 .and. convergenceFlag /= 5 ) then
                
                ! Convert pm from unbound to bound space and insert into nRhoParams:
                    
                call insertFreeParams(transformToBound(pm,lowerBound,upperBound,lBoundMe))
                ! in em2d.f90, inserts pm into rhoParams array used by writeResistivityFile
                ! we also transform from Occam's unbound space to real resistivity in bound space
                ! since pm is always in unbounded space
                                           
                ! Display any joint inversion misfits:
                call displayJointInvMisfits()
                 
                ! output an iteration file of current rhoParams array:
                call writeResistivityFile(currentIteration,modelRMS,targetRMS,modelMu,modelRoughness)
                
                ! output the model response:
                call writeResponse(currentIteration)            

            endif

            !
            ! Check the convergence flag:
            ! 
            if (convergenceFlag > 1) then
                exit
            end if

        enddo
        
        ! Close the Occam log file:
        call closeOccamLogFile


    endif
    
    !
    ! Deallocate memory
    !
        write(*,*) 'Deallocating memory...'
        call deallocateOccam
        call mare2dem_deallocate
            
    !
    ! Shut down nicely:
    !
        call get_time_offset( time0,timeOffset)
        
        call exitMARE2DEM()
        
    !
    ! End of Manager node section  
    !
    endif ! (myid == 0)
   
    call mpi_finalize(ierr)
    
end program runMARE2DEM


!==================================================================================================================================! 
!=================================================================================================== getMARE2DEMCommandLineArguments
!==================================================================================================================================!      
    subroutine getMARE2DEMCommandLineArguments(lFwdOnly,lSaveJacobian,lSaveSensitivity,nCurrentIter)
    !
    ! Subroutine to get command line arguments for MARE2DEM.
    ! 
    use mare2dem_io
    use mare2dem_global, only : resistivityFile, outputFileRoot, lFwdFields, scratchFolder
    use string_helpers
 
    implicit none
    
    logical, intent(out)   :: lFwdOnly, lSaveJacobian, lSaveSensitivity
    integer, intent(out)   :: nCurrentIter
    
    
    character(256) :: arg, sExt1, rfile='', ofile='', rbase
    integer        :: n, iExt1, idot, ict
    logical        :: lokay1, lBad=.false., lrfile = .false.
        
    n = command_argument_count()
    
    if ( (n == 0) ) then
        write(*,*) ' '
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) 'MARE2DEM error, no command line arguments given!' 
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) ' '
        call mare2dem_help
    endif
    
    write(*,*) ' Checking input arguments... '
    ! 
    ! n > 0:
    !
    ict = 1
    
    do while (ict <= n)
    
        call get_command_argument(ict, arg)
    
        select case (arg)
    
        case ('?','-help','/help',char(92)//'help') ! Help info requested:
            call mare2dem_help
            
        case ('-f','-F','/F','/f',char(92)//'F',char(92)//'f') ! Forward run only
            lFwdOnly = .true.
            
         !    if (n == 1) then
!                 write(*,*) ' '
!                 write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
!                 write(*,*) 'MARE2DEM error, no command line arguments given!' 
!                 write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
!                 write(*,*) ' '
!                 call mare2dem_help
!             endif
!  
            write(*,*) ' -F argument given. This is a forward only call'
 
            ict = ict + 1
        case ('-fj','-FJ','/FJ','/fj',char(92)//'FJ',char(92)//'fj') ! Forward run + Jacobian computation only
            lFwdOnly = .true.
            lSaveJacobian = .true.
 
            write(*,*) ' -FJ argument given. This is a forward + Jacobian matrix only call'
 
            ict = ict + 1

        case ('-j','-J','/J','/j',char(92)//'J',char(92)//'j') ! Save Jacobian from start of each iteration
            lSaveJacobian = .true.
 
 
            write(*,'(a)') '  -J argument given. Jacobian from the STARTING model of each iteration will be saved.'
 
            ict = ict + 1
            
        case ('-S','-s','/S','/s',char(92)//'S',char(92)//'s') ! Save sensitivity from start of each iteration
            lSaveSensitivity = .true.
 
 
            write(*,'(a)') '  -S argument given. Sensitivity from the STARTING model of each iteration will be saved.'
 
            ict = ict + 1
      
            
        case  ('-scratch','-SCRATCH',char(92)//'scratch',char(92)//'SCRATCH') ! Scratch folder
            ict = ict + 1    
            call get_command_argument(ict, scratchFolder)                    
            ict = ict + 1           
            write(*,*) '-scratch arg given. The 2.5D inversion scratch folder is: ', trim(scratchFolder)
            
        case default
            if (.not.lrfile) then
                rfile = arg
                lrfile = .true.   
            else
                ofile = arg
            endif
            ict = ict + 1  
            
        end select   
        
    enddo 

     if (.not.lrfile) then
        write(*,*) ' '
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) 'MARE2DEM error, no resistivity file given!' 
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ' 
        write(*,*) ' '
        call mare2dem_help
    endif              
    
    !
    ! Now parse the input names to remove any .resistivity extensions:
    !
    rfile = adjustl(rfile)
    ofile = adjustl(ofile)
    
    ! Remove .resistivity from filename:        
    resistivityFile = trim(rfile)
    idot = index(resistivityFile,'.',.true.)
    if (idot > 0) then 
        sExt1 = resistivityFile(idot+1:)
        if (trim(sExt1) =='resistivity')  resistivityFile = resistivityFile(1:idot-1)
    endif
    
    ! At this point there is no resistivity extension. Now get the iteration number:
    idot = index(resistivityFile,'.',.true.)
    if (idot > 0) then 
        sExt1 = resistivityFile(idot+1:)
        call convertStrToInteger(sExt1,iExt1,lokay1)
        if (lokay1) then ! extension is an integer, omit it
            nCurrentIter = iExt1
            rbase = resistivityFile(1:idot-1)
        else
            lBad = .true.
        endif
    else
        lBad = .true.
    endif                 
    if (lBad) then
        write(*,*) ' ' 
        write(*,*) 'Error with command line arguments, no iteration number in resistivity file: '
        write(*,*) trim(adjustl(rfile))
        write(*,*) ' ' 
        write(*,*) 'Try again!'
        call exitMARE2DEM()
    endif
 
    if (len_trim(ofile) == 0)  then
        outputFileRoot = rbase
    else
        outputFileRoot = ofile
    endif
    
    
        
 
    write(*,*) ' '
    write(*,*) ' Input resistivity file name root: ', trim(resistivityFile)
    write(*,*) 'Output resistivity file name root: ', trim(outputFileRoot)
    write(*,*) ' '

    end subroutine  getMARE2DEMCommandLineArguments


    