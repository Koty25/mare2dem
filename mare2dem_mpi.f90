!-----------------------------------------------------------------------
!
!    Copyright 2008-2017
!    Kerry Key
!    Scripps Institution of Oceanography
!    kkey@ucsd.edu
!
!    This file is part of MARE2DEM.
!
!    MARE2DEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    MARE2DEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------
!==================================================================================================================================! 
!============================================================================================================ mare2d_mpi_definitions
!==================================================================================================================================!
#include <scorep/SCOREP_User.inc>

    module mare2d_mpi_definitions

#if defined(__INTEL_COMPILER)
    use IFPORT          ! DGM 2/23/2012 for the sleepqq routine
#endif
    
#if defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
    include 'mpif.h'
#else
    use mpi
#endif

!
! MPI definitions common to manager and Workers:
!
    
    integer, parameter :: manager = 0        ! manager is 0, all others are Workers
    integer            :: mcomm = MPI_COMM_WORLD
    integer            :: status(MPI_STATUS_SIZE)   
              
    integer, parameter :: tag_M_SendingData              = 1 ! tag for manager sending data     
    integer, parameter :: tag_W_SendingData              = 2 ! tag for worker sending data   
            
    ! iflags are like tags, but i put them in the data field of the mpi command
    integer, parameter :: iflag_TellWorkerToQuit            = 0
    integer, parameter :: iflag_TellWorkerToRunSubset       = 1   
    integer, parameter :: iflag_TellWorkerToScaLapackMM     = 2      
    integer, parameter :: iflag_TellWorkerToScaLapackMS     = 3
    integer, parameter :: iflag_TellWorkerToScratch         = 4
    integer, parameter :: iflag_TellWorkerToBcast           = 5
           
    end module mare2d_mpi_definitions

 
!==================================================================================================================================! 
!====================================================================================================================== mpi_mare2dem
!==================================================================================================================================!
    subroutine mpi_mare2dem()
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu 
!
    use mare2d_mpi_definitions
    use mare2dem_global    
    use mare2dem_input_data_params
    
    implicit none
 
    logical :: lkeepgoing, lflag
    
    integer :: ierr,  iWorker, nProc 

    call mpi_comm_size( mcomm, nProc, ierr )  
         
!
! Set up the wavenumbers, starting mesh and the refinement group arrays:
!
    call getWavenumbers    
    call getRxTxGroups
    call getFqGroups   

    call getRefinementGroups
    
    if (lprintDecomposition) call printDecomposition 
    
    if (lPrintSetup) call display_MARE2DEM_Params
    
!
! Launch the manager controller:
!
    do    
                
        status = 0
        
        !
        ! Check for a message from any Worker:
        !
        call  mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mcomm, lflag, status, ierr)
        
        !
        ! If message from Worker, then decode it:
        !
        if (lflag) then 
        
            iWorker = status(MPI_SOURCE)

            if  (status(MPI_TAG) == tag_W_SendingData) then
                
                !
                ! Receive data from Worker (results for this Worker's subset)
                !
                call mpi_manager_receive_results(iWorker)  
        
                                     
            endif 
            
        endif ! lflag 
                
        !
        ! Try sending out a job to the first available worker in lWorker_status
        !
        call mpi_manager_send_job()
        
        !
        ! David Myer's tweak for idling the manager, which can be useful when running MARE2DEM on a small
        ! multicore desktop or laptop. This keeps the manager from consuming 100% cpu while waiting for
        ! the workers to finish:
        !
#if defined(__INTEL_COMPILER)  
        if (nProc <= 12) call sleepqq(1)
#endif
   
        !
        ! Check to see if we're done with all the models.  If so then terminate the worker controllers.
        !
        lkeepgoing = .true.
        if ( all( lWorker_status ) ) then
             lkeepgoing = .false.
             exit
        endif
     
        
    enddo ! manager controller loop
!        
! Done with the adaptive mesh refinements 
!

!
! Deallocate all internal arrays:
!
    deallocate( RxTxGroups, FqGroups, refinementGroups )
    ! this also deallocates the allocatable fields within the derived types
    
    if (allocated(wavenum))       deallocate (wavenum)

    if (allocated(inputmesh%attr) ) call deallocate_trimesh(inputmesh,.false.)
    
    if (lprintDebug) write(*,*) 'leaving mpi_mare2dem '
    
    
    end subroutine mpi_mare2dem

!==================================================================================================================================! 
!==================================================================================================== mpi_manager_initializeScratch
!==================================================================================================================================!
    subroutine mpi_manager_initializeScratch

    use mare2d_mpi_definitions
    use mare2dem_global    

    implicit none
    
    integer :: ierr, iWorker 

    
    do iWorker = 1,nWorkers
        call mpi_send( iflag_TellWorkerToScratch,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
        call mpi_send( scratchFolder,             256, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )        
        
    enddo
    
    end subroutine mpi_manager_initializeScratch

!==================================================================================================================================! 
!============================================================================================================== mpi_manager_send_job
!==================================================================================================================================!
    subroutine mpi_manager_send_job

 
    use mare2d_mpi_definitions
    use mare2dem_global    

    implicit none
    
    integer :: i, iWorker
    
    !
    ! Send out work until there are no more workers are available:
    !
    
    do 
        
        !
        ! Find first available worker:
        !
        iWorker = 0
        do i = 1,nWorkers
            if (lWorker_status(i)) then
                iWorker = i
                exit
            endif
        enddo

        if (iWorker == 0) return ! no workers available

        
        !
        ! Now find a job for that worker:
        !
        if ( iPtr_refGroups <= nRefinementGroups ) then

            call mpi_manager_send_subset(iWorker, iPtr_refGroups)            
            iPtr_refGroups = iPtr_refGroups + 1
            cycle
      
        endif    
     
        !
        ! There are no tasks to do, return to main loop
        !
        return
        
    enddo
    
    end subroutine mpi_manager_send_job

!==================================================================================================================================! 
!=========================================================================================================== mpi_manager_send_subset
!==================================================================================================================================!
   subroutine mpi_manager_send_subset(iWorker,iGrp)
   
    use mare2d_mpi_definitions
    use kx_io
    use mare2dem_global
    use mare2dem_input_data_params
    use Occam, only : d,numForwardCalls,currentIteration
    
    implicit none
    
    integer, intent(in) :: iWorker, iGrp
    
    !
    ! Local:
    !
    integer         :: nRR, nTT, ierr ,nnod, nele, nn, iRxTxG, iFqG, nFq, rx0,rx1,tx0,tx1,fq0,fq1, nsend, ntr
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
       
    call cpu_time(tStart) 
    
    lWorker_status(iWorker) = .false. 
    
    if (lprintDebug) write(*,'(a,2(i6,1x))') 'manager sending job, iGrp, iWorker: ', iGrp, iWorker
    
    !
    ! Tell the worker that your going to send it a model to run:
    !
    call mpi_send( iflag_TellWorkerToRunSubset, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
 
    !
    ! Tell the worker which refinement group is being sent:
    !
    call mpi_send( iGrp,     1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr ) 
 
    iRxTxG  = refinementGroups(iGrp)%iRxTxGroup 
    
    nRR     = RxTxGroups(iRxTxG)%nRx          
    rx0     = RxTxGroups(iRxTxG)%iRx(1)
    rx1     = RxTxGroups(iRxTxG)%iRx(nRR)
    
    nTT     = RxTxGroups(iRxTxG)%nTx  
    tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
    tx1     = RxTxGroups(iRxTxG)%iTx(nTT)

    iFqG    = refinementGroups(iGrp)%iFqGroup
        
    nFq     = FqGroups(iFqG)%nFq
    fq0     = FqGroups(iFqG)%iFq(1)
    fq1     = FqGroups(iFqG)%iFq(nFq)
   
    if (lprintDebug) write(*,'(a,6(i6,1x))') 'mpi_manager_send_subset: rx0,rx1,tx0,tx1,fq0,fq1:', rx0,rx1,tx0,tx1,fq0,fq1
    
    sType = refinementGroups(iGrp)%sType
 
    call mpi_send( sType,             2, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lMTscatteredField, 1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )   
                
    if (sType == 'mt') then 
    
        ! send MT receivers:
        call mpi_send( nRR,                                  1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( xRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( yRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( zRxMT(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( ThetaRxMT(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( AlphaRxMT(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( BetaRxMT(RxTxGroups(iRxTxG)%iRx),   nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( lengthRxMT(RxTxGroups(iRxTxG)%iRx), nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
        call mpi_send( nQuadRxMT,                            1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )    
 
        ! Send frequencies:    
        nFq = FqGroups(iFqG)%nFq
        call mpi_send(   nFq,                              1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )  
        call mpi_send( fTxMT(FqGroups(iFqG)%iFq),    nFq, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr ) 
    
    elseif  (sType == 'cs') then 
    
        ! Send CSEM receivers:
        call mpi_send( nRR,                                    1, MPI_INTEGER,          iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( xRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( yRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( zRxCSEM(RxTxGroups(iRxTxG)%iRx),      nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( ThetaRxCSEM(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( AlphaRxCSEM(RxTxGroups(iRxTxG)%iRx),  nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( BetaRxCSEM(RxTxGroups(iRxTxG)%iRx),   nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( lengthRxCSEM(RxTxGroups(iRxTxG)%iRx), nRR, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( nQuadRxCSEM,                            1, MPI_INTEGER,          iWorker, tag_M_SendingData,mcomm,ierr) 
        
        ! Send CSEM transmitters:
        call mpi_send( nTT,                                     1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( xTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( yTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( zTxCSEM(RxTxGroups(iRxTxG)%iTx),       nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( azimuthTxCSEM(RxTxGroups(iRxTxG)%iTx), nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( dipTxCSEM(RxTxGroups(iRxTxG)%iTx),     nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( lengthTxCSEM(RxTxGroups(iRxTxG)%iTx),  nTT, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( cSourceType(RxTxGroups(iRxTxG)%iTx), 8*nTT, MPI_CHARACTER,       iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( nQuadTxCSEM,                             1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )
 
       ! Send wavenumbers:
        call mpi_send(   nwave,                                 1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )  
        call mpi_send( wavenum,                             nwave, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( nKxPerGroup,                             1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )
        
        ! Send frequencies:
        nFq = FqGroups(iFqG)%nFq
        call mpi_send( nFq,                                     1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr ) 
        call mpi_send( fTxCSEM(FqGroups(iFqG)%iFq),           nFq, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
    
    elseif  (sType == 'dc') then 
        
        ! Send Rx and Tx arrays that point to specific electrode pairs:
        call mpi_send( nRR,                                     1, MPI_INTEGER,         iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( RxTxGroups(iRxTxG)%iRxDC,            nRR*2, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( nTT,                                     1, MPI_INTEGER,         iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send( RxTxGroups(iRxTxG)%iTxDC,            nTT*2, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )  
!         write(*,'(a,*(i1,x))') ' RxTxGroups(iRxTxG)%iRxDC: ' , RxTxGroups(iRxTxG)%iRxDC
!         write(*,'(a,*(i1,x))') ' RxTxGroups(iRxTxG)%iTxDC: ' , RxTxGroups(iRxTxG)%iTxDC
             
        ! Send Rx and Tx electrodes:
        ntr = size(RxTxGroups(iRxTxG)%iRxTrodes,1)
        call mpi_send( ntr,                                              1, MPI_INTEGER,         iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send(trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes,1:3), 3*ntr, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )  
        
        ntr = size(RxTxGroups(iRxTxG)%iTxTrodes,1)
        call mpi_send( ntr,                                              1, MPI_INTEGER,         iWorker, tag_M_SendingData,mcomm,ierr)
        call mpi_send(trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes,1:3), 3*ntr, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )  
              
        ! Send wavenumbers:
        call mpi_send(   nwave,                                 1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )  
        call mpi_send( wavenum,                             nwave, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )
        call mpi_send( nKxPerGroup,                             1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr )

        ! Send frequencies:
        nFq = 1
        call mpi_send( nFq,                                     1, MPI_INTEGER,         iWorker,tag_M_SendingData,mcomm,ierr ) 
        call mpi_send( 0.,                                    nFq, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr )

    endif
 
    !
    ! Send the data parameter array:
    !
    call mpi_send( refinementGroups(iGrp)%nd,      1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )      
    call mpi_send( 4,                              1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )  
    nsend = refinementGroups(iGrp)%nd*4
 
    call mpi_send( refinementGroups(iGrp)%dp(:,1:4),  nsend, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )    
    call mpi_send(  d(refinementGroups(iGrp)%dp(:,5)) ,  refinementGroups(iGrp)%nd, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )   
    
    !
    ! Send the boundary model: 
    !    
    nnod   = inputmodel%nnod
    nele   = inputmodel%nele
    call mpi_send( nnod, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( nele, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
          
    call mpi_send( inputmodel%y,              nnod, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%z,              nnod, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%attr,           nele, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%emap,         3*nele, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%neighborlist, 3*nele, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
        
    call mpi_send( inputmodel%numberofpointattributes,    1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%numberofcorners,            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%numberoftriangleattributes, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%numberofsegments,           1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%numberofholes,              1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( inputmodel%numberofregions,            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )


    nn =  inputmodel%numberofpointattributes*inputmodel%nnod
    if (nn > 0 ) call mpi_send( inputmodel%pointattributelist, nn, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

    call mpi_send( inputmodel%pointmarkerlist,   nnod,                         MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)
    call mpi_send( inputmodel%segmentlist,       2*inputmodel%numberofsegments, MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)
    call mpi_send( inputmodel%segmentmarkerlist,   inputmodel%numberofsegments, MPI_INTEGER, iWorker, tag_M_SendingData,mcomm,ierr)

    if ( inputmodel%numberofholes > 0 ) then
        call mpi_send( inputmodel%holelist, 2*inputmodel%numberofholes, MPI_DOUBLE_PRECISION,iWorker,tag_M_SendingData,mcomm,ierr)
    endif

    call mpi_send( inputmodel%regionlist, 4*inputmodel%numberofregions, MPI_DOUBLE_PRECISION, iWorker,  &
          & tag_M_SendingData, mcomm,ierr )

    ! send fixed and free resistivity values:
    call mpi_send( cAnisotropy,        24, MPI_CHARACTER,        iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( nRhoParams ,         1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( rhoParams , nRhoParams, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( iFreeParam, nRhoParams, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send( nFree ,              1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    
    ! send settings for mesh refinement:
    call mpi_send( maxnadapt,                            1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( errortolerance,  1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( idual_func,      1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( minQangle,       1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )

    call mpi_send( outputFileRoot,         256, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lSaveMeshFiles,            1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send( linversion,                1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )  
    call mpi_send( lDisplayRefinementStats ,  1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )
   
    call mpi_send( linearSolver,  32, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lprintDebug,    1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )    
	call mpi_send( lUseInversionMeshCoarsening, 1, MPI_LOGICAL,   iWorker, tag_M_SendingData, mcomm,ierr )      
	call mpi_send( phaseConvention,4, MPI_CHARACTER, iWorker, tag_M_SendingData, mcomm,ierr )
	
    call mpi_send( minRangeProtector, 1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( minArea,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( maxMeshNodes,      1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( ecutoff,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( hcutoff,           1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( max_nsubrefine,    1, MPI_INTEGER,          iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( pct_refine,        1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( reciprocityUsed,   4, MPI_CHARACTER,        iWorker, tag_M_SendingData, mcomm,ierr )  

    call mpi_send( numForwardCalls,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lReuseRefine,      1, MPI_LOGICAL, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( maxnadapt_default,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    !call mpi_send( t0_em2dkxTrace,  1, MPI_DOUBLE_PRECISION, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( lprintTrace_em2dkx,      1, MPI_LOGICAL, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( currentIteration,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
    call mpi_send( nbOccItToReuseRefine,   1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm,ierr )
       
    !
    ! Worker will then compute the 2D EM response and send its results back to the manager...
    !           
    call cpu_time(tEnd)     
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,    69x , a32,  2x,a19)'    
    if (lprintMPItimers) write(*,sFmt)  'Proc:',0,'Group:',iGrp,   'mpi_manager_send_subset:',   trim(stime)

    end subroutine mpi_manager_send_subset
!==================================================================================================================================!
!=========================================================================================================== mpi_manager_send_subset
!==================================================================================================================================!
    subroutine mpi_manager_receive_results(iWorker)
 
    use em2dkx_mod 
    use mare2d_mpi_definitions
    use mare2dem_global
 
    
    use Occam, only : wj, dm   
    
    implicit none
    
    integer, intent(in) :: iWorker
    
    integer :: nrcv, iGrp, ierr,i,  nFreeRecv
    
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
    
    integer, dimension(:),  allocatable :: iFreeInd
    real(8), dimension(:),  allocatable :: dbuff
    real(8), dimension(:,:),allocatable :: wjbuff
    
 
    call cpu_time(tStart) 
!
! First receive the refinement group index:
!
    call mpi_recv( iGrp,   1, MPI_INTEGER, iWorker, tag_W_SendingData , mcomm, status, ierr )
    
    nrcv = refinementGroups(iGrp)%nd
      
    allocate(dbuff(nrcv))
 

    call mpi_recv( dbuff,   nrcv, MPI_DOUBLE_PRECISION, iWorker, tag_W_SendingData , mcomm, status, ierr )

    dm(refinementGroups(iGrp)%dp(:,5)) = dbuff

    deallocate(dbuff)
   
    if (linversion) then
            
        if (( refinementGroups(iGrp)%sType /= 'mt').and.(lUseInversionMeshCoarsening)) then

            call mpi_recv( nFreeRecv,         1, MPI_INTEGER, iWorker, tag_W_SendingData , mcomm, status, ierr )
            allocate(iFreeInd(nFreeRecv))
            call mpi_recv( iFreeInd,   nFreeRecv, MPI_INTEGER, iWorker, tag_W_SendingData , mcomm, status, ierr )

        else
            allocate(iFreeInd(nFree))
            iFreeInd(1:nFree) = [1:nFree]
            nFreeRecv = nFree
        endif
       
!       write(*,*) 'manager recv iGrp: ',iGrp
!       do i = 1,nFree
!        write(*,*) i,iFreeInd(i)
!        
!        enddo

             
                
        nrcv = refinementGroups(iGrp)%nd*nFreeRecv
        
        allocate(wjbuff(refinementGroups(iGrp)%nd,nFreeRecv) )
        
        call mpi_recv( wjbuff, nrcv, MPI_DOUBLE_PRECISION, iWorker, tag_W_SendingData , mcomm, status, ierr )
        
        wj(refinementGroups(iGrp)%dp(1:refinementGroups(iGrp)%nd,5),iFreeInd) = wjbuff
        
        
!        write(*,*) 'manager recv index global data: ',iGrp         
!        do i = 1, refinementGroups(iGrp)%nd  
!         write(*,*) i,refinementGroups(iGrp)%dp(i,5), wjbuff(i,1)  
!        
!        enddo
        
        deallocate(wjbuff,iFreeInd) 
    endif

    call mpi_recv( refinementGroups(iGrp)%timer, 1, MPI_DOUBLE_PRECISION,  iWorker, tag_W_SendingData , mcomm, status,ierr )
 
      
    call cpu_time(tEnd)     
    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,    a24,2x,i7,  14x,  54x,   2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',0,'Group:',iGrp,   'Manager recv results:',0,   trim(stime )
!
! Set the availability status for this worker:
!
    lWorker_status(iWorker) = .true.   
        
    end subroutine mpi_manager_receive_results 
        
!==================================================================================================================================! 
!============================================================================================================= mpi_worker_controller
!==================================================================================================================================!
    subroutine mpi_worker_controller

    use mare2d_mpi_definitions    
    use occam 
 
    
    implicit none
    
    logical :: lflag
    
    integer :: myID, ierr, iflag
    real(8) :: muDummy
    
    !
    ! Determine which worker I am:
    !
    call mpi_comm_rank( mcomm, myID, ierr )
    
    !
    ! Create an infinite loop for this worker controller:
    ! 
    do 

        !
        ! David Myer's tweak for idling the worker, which can be useful when running MARE2DEM on a small
        ! multicore desktop or laptop. This keeps the worker from consuming 100% cpu while waiting for
        ! a new task from the manager:
        !
#if defined(__INTEL_COMPILER)            
        !
        ! Check for a message from the manager:
        !
        call  mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, mcomm, lflag, status, ierr)
        
        !
        ! If no message from the manager, take a valium and try later:
        !
        if (.not.lflag) then 
            call sleepqq(1) ! kwk debug: this is specific to the intel compiler
            cycle
        endif
#endif        
        !
        ! Get the command from the manager:
        !
        call mpi_recv( iflag, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
        
        !
        ! What is the manager telling me?
        !
        
        select case (iflag)
        
        case ( iflag_TellWorkerToQuit ) 
           
            write(*,'(a12,1x,i5,1x,a32)') ' Worker: ', myID, ' is done with mare2dem...'     
            
            exit    ! Worker now leaves this subroutine
           
        case ( iflag_TellWorkerToRunSubset )     
            
            call mpi_worker_run_subset()  
            
         case ( iflag_TellWorkerToScratch )     
            
            call mpi_worker_intializeScratch()                      
            
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
        case ( iflag_TellWorkerToScaLapackMM )     
            
            call mpi_worker_ScaLapackMM()

        case ( iflag_TellWorkerToScaLapackMS )     
            
            muDummy = 0d0
            call mpi_worker_ScaLapackMS(ierr,muDummy)

        case ( iflag_TellWorkerToBcast )     
            
            call mpi_all_Bcast 
#endif                           
        case default
            write(*,*) 'Error in Worker controller, received unknown command from the manager !'
            write(*,*) 'iflag is: ',iflag
            write(*,*) ' quit mpi and debug the code ...'
            stop
        end select

    enddo   ! the infinite while loop
    
    end subroutine mpi_worker_controller
    
!==================================================================================================================================! 
!======================================================================================================= mpi_worker_intializeScratch
!==================================================================================================================================!  
    subroutine mpi_worker_intializeScratch


    use mare2d_mpi_definitions    
    use mare2dem_global
   
    implicit none
 
   
    integer         :: myID, ierr
    
    character(256)  :: cFilename, strID
    integer         :: ios
    
    logical         :: lstatus
    !
    ! Determine which worker I am:
    !
    call mpi_comm_rank( mcomm, myID, ierr )

    !
    ! Get scratch folder base name:
    ! 
    call mpi_recv( scratchFolder,         256, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
             
             

    scratchFolder = trim(scratchFolder)//'/'  
 
    
    ! 
    ! Test for scratch folder existence by trying to write a dummy file:
    !
    write(strID,*) myID       
    cFilename  = trim(scratchFolder)//'test'//trim(adjustl(strID)) 
 
    open(unit=21,file=trim(cFilename),form="unformatted",status="replace",access="stream",iostat=ios)
              
    if (ios /= 0) then
    
        !
        ! No scratch folder, so create one:
        !

#if defined(__INTEL_COMPILER)  
 
        if( .not. makedirqq(trim(scratchFolder)) ) then
        
            write(*,'(a)') ' ------------------------------------------------------------------------------------------------------'
            write(*,'(a)') ' Error on processor: ',myID
            write(*,'(a)') ' MARE2DEM cannot create a new scratch folder. Stopping! '
            write(*,'(a,a)') ' The scratch folder MARE2DEM tried to use is: ',trim(scratchFolder)
            write(*,'(a)') ' '
            write(*,'(a)') ' What has gone wrong? '
            write(*,'(a)') '  (1) The default scratch folder does not exist on your cluster. '
            write(*,'(a)') '  or '
            write(*,'(a)') '  (2) You specified a scratch folder in the settings file, but it does not exist on your system.  '
            write(*,'(a)') '     Solution: create this folder before running MARE2DEM. '
            write(*,*)    ' '
            write(*,'(a)') ' Note that on large clusters, an existing local folder for each node should be used (for example /tmp) '
            write(*,'(a)') ' so that the scratch data is not transmitted across the network and instead only uses the local disk.   '
            write(*,'(a)') ' If you are running MARE2DEM on a single computer (laptop or desktop), then you can use a local '
            write(*,'(a)') ' directory for scratch space. '
            write(*,'(a)') ' Or if you are running MARE2DEM on cluster with a superfast parallel filesystem, you might be able to use'
            write(*,'(a)') ' a networked directory for scratch space, but a local directory is still probably more efficient. '            
            write(*,*)    ' '
            write(*,*)    ' So long and thanks for all the fish! '
            write(*,*)    ' '
            stop
        endif
 
!#elif defined(__GNUC__) 
!        
!        call system("mkdir scratch")
                  
#endif
    
    
    else
    
        !
        ! A-okay, folder exists, delete the dummy file:
        !
        close(21,status="delete")
        
#if defined(__INTEL_COMPILER)  
        
        ! Erase its contents:
        lstatus = delfilesqq(trim(scratchFolder)//'mare2dem_tmp*')
 
 
!#elif defined(__GNUC__) 
!        
!        call system("rmdir -fR scratch") ! this should work but i've not tested in on gfortran yet...
!        call system("mkdir scratch")
      
                  
#endif   
     
    endif 
   
 
    
    end subroutine mpi_worker_intializeScratch  

!==================================================================================================================================!
!============================================================================================================= mpi_worker_run_subset
!==================================================================================================================================!
    subroutine mpi_worker_run_subset

    use mare2d_mpi_definitions 
    use mare2dem_worker
    use mare2dem_global
    
    implicit none
 
    integer(8)  :: t0,t1

    !SCOREP REGIONS DEFINITION
#ifdef TRACE
    SCOREP_USER_REGION_DEFINE(fwd_compute)
#endif
  
    call system_clock(count_rate=clockRate)    
    call system_clock(t0)        
!
! Receive data from the manager:
!
    call mpi_worker_receive_subset 
     
!
! Coarsen input model if this is CSEM data:  
!
    if ( (sType /= 'mt') .and. (lUseInversionMeshCoarsening) ) call coarsenModel
!
! Generate the starting mesh:
!    
    if (lprintDebug) write(*,*) '...call getStartingMesh   ',myid
    call getStartingMesh   
    
!
! Compute the 2D EM response:
!
#ifdef TRACE
    SCOREP_USER_REGION_BEGIN(fwd_compute, "fwd_compute", SCOREP_USER_REGION_TYPE_COMMON)
#endif
    if (lprintDebug) write(*,*) '...call worker_EM2D   ',myid
    call worker_EM2D

#ifdef TRACE
    SCOREP_USER_REGION_END(fwd_compute)
#endif
!
! Send results back to the manager
!
    call system_clock(t1)
    
    call mpi_worker_send_results( (t1-t0)/clockRate )   
!
! Deallocate:
!  
    call mare2dem_deallocate  ! this routine deallocates everything from the worker
      
    end subroutine mpi_worker_run_subset
    
!==================================================================================================================================! 
!========================================================================================================= mpi_worker_receive_subset
!==================================================================================================================================!
    subroutine mpi_worker_receive_subset

    use mare2d_mpi_definitions    
    use kx_io
    use mare2dem_global
    use mare2dem_input_data_params, only : phaseConvention, cnDataParams, reciprocityUsed
    use occam, only : numForwardCalls,currentIteration
    use mare2dem_worker
    
    implicit none
 
    integer        :: ierr, nele, nnod, nn
 
    real(8)         :: tStart, tEnd
    character(256)  :: sFmt
    character(32)   :: stime
    
    call cpu_time(tStart) 
       
    call mpi_comm_rank( mcomm, myID, ierr )
 
    nRxIn        = 0
    nfrequencies = 0   
       
    !
    ! Receive data for em2dkx from the manager:
    !
    call mpi_recv( iRefinementGrp,    1, MPI_INTEGER,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( sType,             2, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( lMTscatteredField, 1, MPI_INTEGER,   manager, tag_M_SendingData, mcomm, status, ierr )  

    if ( (sType == 'mt') .or. (sType == 'cs')) then ! receive the EM/MT receivers:
    
        call mpi_recv( nRxIn,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  

        allocate( xRxIn(nRxIn),yRxIn(nRxIn),zRxIn(nRxIn),thetaRxIn(nRxIn),alphaRxIn(nRxIn),betaRxIn(nRxIn),lengthRxIn(nRxIn)  )
      
        call mpi_recv( xRxIn,      nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( yRxIn,      nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( zRxIn,      nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( thetaRxIn,  nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( alphaRxIn,  nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( betaRxIn,   nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )    
        call mpi_recv( lengthRxIn, nRxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )    
        call mpi_recv( nquad_rx,       1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
      
        nKxPerGroup = 1  
        nTxIn       = 1  
        
    endif
    
    if (sType == 'cs') then  ! Receive CSEM transmitters:
       
        call mpi_recv( nTxIn,  1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
 
        allocate ( azimuthTxIn(nTxIn),dipTxIn(nTxIn),xTxIn(nTxIn),yTxIn(nTxIn),zTxIn(nTxIn), &
                 & lengthTxIn(nTxIn), TxTypeIn(nTxIn) )    ! x2 for negative wavenumbers (filled in later)
                     
        call mpi_recv( xTxIn,        nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
        call mpi_recv( yTxIn,        nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( zTxIn,        nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( azimuthTxIn,  nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( dipTxIn,      nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( lengthTxIn,   nTxIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( TxTypeIn ,  8*nTxIn, MPI_CHARACTER,        manager, tag_M_SendingData, mcomm, status, ierr )      
        call mpi_recv( nquad_tx,         1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  

        ! Receive the wavenumbers:
        call mpi_recv( nwave,               1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate( wavenum(nwave) )
        call mpi_recv( wavenum,         nwave, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( nKxPerGroup,         1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
                 
    endif
    
    if (sType == 'dc') then ! receiver DC Rx,Tx and electrodes:
        
        ! Receive Rx and Tx arrays that point to specific electrode pairs:
        call mpi_recv( nRxDCIn,                1, MPI_INTEGER,  manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate(RxDCIn(nRxDCIn,2))
        call mpi_recv( RxDCIn,         2*nRxDCIn, MPI_INTEGER,  manager, tag_M_SendingData, mcomm, status, ierr ) 
        
        call mpi_recv( nTxDCIn,                1, MPI_INTEGER,  manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate(TxDCIn(nTxDCIn,2))
        call mpi_recv( TxDCIn,         2*nTxDCIn, MPI_INTEGER,  manager, tag_M_SendingData, mcomm, status, ierr ) 
        
!         write(*,*)' m2d recv: '
!         do i = 1,nTxDCIn
!             write(*,'(a,3(i,1x))') 'i,TxDCIn(i,1:2):',i,TxDCIn(i,1),TxDCIn(i,2)
!         enddo
!         do i = 1,nRxDCIn
!             write(*,'(a,3(i,1x))') 'i,RxDCIn(i,1:2):',i,RxDCIn(i,1),RxDCIn(i,2)
!         enddo   
    
        ! Receive Rx and Tx electrodes:
        
        call mpi_recv( nTrodesRxDCIn,                1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate(trodes_RxDCIn(nTrodesRxDCIn,3))
        call mpi_recv( trodes_RxDCIn,  3*nTrodesRxDCIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        
        call mpi_recv( nTrodesTxDCIn,                1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate(trodes_TxDCIn(nTrodesTxDCIn,3))
        call mpi_recv( trodes_TxDCIn,  3*nTrodesTxDCIn, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
   
!          write(*,*)' m2d recv: '
! 
!         do i = 1,nTrodesRxDCIn
!             write(*,'(a,1x,i,1x,2(g,1x))') 'i,trodes_RxDCIn(i,2:3):',i,trodes_RxDCIn(i,2),trodes_RxDCIn(i,3)
!         enddo  
!         
!         do i = 1,nTrodesTxDCIn
!             write(*,'(a,1x,i,1x,2(g,1x))') 'i,trodes_TxDCIn(i,2:3):',i,trodes_TxDCIn(i,2),trodes_TxDCIn(i,3)
!         enddo
                      
        ! Receive the wavenumbers:
        call mpi_recv( nwave,               1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
        allocate( wavenum(nwave) )
        call mpi_recv( wavenum,         nwave, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
        call mpi_recv( nKxPerGroup,         1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
    
    endif
    
    ! 
    ! Receive the frequencies:
    !
    call mpi_recv( nfrequencies,            1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
    allocate( frequencies(nfrequencies) )
    call mpi_recv( frequencies,  nfrequencies, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )                          
 
    !
    ! Receive the data parameter array:
    !    
    call mpi_recv( nd_local,                     1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
    call mpi_recv( cnDataParams,                 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )  
    allocate(dp_local(nd_local,cnDataParams),d_local(nd_local))
    call mpi_recv( dp_local, nd_local*cnDataParams, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  
    call mpi_recv( d_local,               nd_local, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )  
  
    
    !
    ! Receive the model:
    !
    call mpi_recv( inputmodel%nnod, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%nele, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )

    nnod = inputmodel%nnod
    nele = inputmodel%nele

    allocate ( inputmodel%y(nnod),inputmodel%z(nnod) )
    allocate ( inputmodel%emap(3,nele), inputmodel%attr(nele) )
    allocate ( inputmodel%neighborlist(3,nele) ) 

    call mpi_recv( inputmodel%y,              nnod, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%z,              nnod, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%attr,           nele, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%emap,         3*nele, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%neighborlist, 3*nele, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )

    call mpi_recv( inputmodel%numberofpointattributes,    1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%numberofcorners,            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%numberoftriangleattributes, 1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%numberofsegments,           1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%numberofholes,              1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( inputmodel%numberofregions,            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )

    nn = inputmodel%numberofpointattributes*inputmodel%nnod
    allocate (  inputmodel%pointattributelist(nn), inputmodel%pointmarkerlist(inputmodel%nnod) ) 
    inputmodel%pointattributelist = 0

    if (nn > 0 ) call mpi_recv( inputmodel%pointattributelist, nn, MPI_DOUBLE_PRECISION, manager, & 
                             & tag_M_SendingData, mcomm, status,ierr )

    call mpi_recv( inputmodel%pointmarkerlist, nnod, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )

    allocate (inputmodel%segmentlist(2*inputmodel%numberofsegments), inputmodel%segmentmarkerlist(inputmodel%numberofsegments)) 

    call mpi_recv( inputmodel%segmentlist,    2*inputmodel%numberofsegments,MPI_INTEGER,manager,tag_M_SendingData,mcomm,status,ierr)
    call mpi_recv( inputmodel%segmentmarkerlist,inputmodel%numberofsegments,MPI_INTEGER,manager,tag_M_SendingData,mcomm,status,ierr)

    allocate (inputmodel%holelist (2*inputmodel%numberofholes) )
    if (inputmodel%numberofholes > 0 ) then
        call mpi_recv( inputmodel%holelist, 2*inputmodel%numberofholes,MPI_DOUBLE_PRECISION, & 
                       & manager,tag_M_SendingData,mcomm,status,ierr)
    endif

    allocate (inputmodel%regionlist(4*inputmodel%numberofregions) ) 
    call mpi_recv( inputmodel%regionlist, 4*inputmodel%numberofregions, MPI_DOUBLE_PRECISION, manager,  &
          & tag_M_SendingData, mcomm, status, ierr )


    call mpi_recv( cAnisotropy, 24, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( nRhoParams,  1, MPI_INTEGER,   manager, tag_M_SendingData, mcomm, status, ierr ) 

    allocate ( rhoParams(nRhoParams) ,iFreeParam(nRhoParams))

    call mpi_recv( rhoParams ,  nRhoParams, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( iFreeParam , nRhoParams, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( nFree,                1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr )  

    call mpi_recv( maxnadapt,                            1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( errortolerance, 1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv( idual_func ,    1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv( minQangle,      1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 


    call mpi_recv( outputFileRoot,        256, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( lSaveMeshFiles,          1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( linversion,              1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv( lDisplayRefinementStats, 1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr ) 


    call mpi_recv( linearSolver,  32, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
    call mpi_recv( lprintDebug ,   1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
	call mpi_recv( lUseInversionMeshCoarsening, 1, MPI_LOGICAL,   manager, tag_M_SendingData, mcomm, status, ierr )
	call mpi_recv( phaseConvention,4, MPI_CHARACTER, manager, tag_M_SendingData, mcomm, status, ierr)
	
    if (lPrintDebug) then
        lprintDebug_em2dkx = .true.
    else
        lprintDebug_em2dkx = .false.
    endif        

    call mpi_recv(minRangeProtector, 1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(minArea,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(maxMeshNodes,      1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(ecutoff,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(hcutoff,           1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(max_nsubrefine,    1, MPI_INTEGER,          manager, tag_M_SendingData, mcomm, status, ierr ) 
    call mpi_recv(pct_refine,        1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr ) 
    
    call mpi_recv( reciprocityUsed,  4, MPI_CHARACTER,        manager, tag_M_SendingData, mcomm, status, ierr )
    lUseReciprocity = .false.
    select case (trim(reciprocityUsed))
    case ('yes')
           lUseReciprocity = .true.          
    end select

    call mpi_recv(numForwardCalls ,1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )    
    call mpi_recv(lReuseRefine    ,1, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv(maxnadapt_default ,1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    !call mpi_recv(t0_em2dkxTrace, 1, MPI_DOUBLE_PRECISION, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv(lprintTrace_em2dkx    ,1, MPI_LOGICAL, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv(currentIteration ,1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    call mpi_recv(nbOccItToReuseRefine ,1, MPI_INTEGER, manager, tag_M_SendingData, mcomm, status, ierr )
    
    call cpu_time(tEnd)     

    write(stime,'(a8,f9.4,a2)') ' Timer: ',tEnd - tStart, ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,   69x, a32,     2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,   'mpi_worker_receive_subset:',   trim(stime )
      
    end subroutine mpi_worker_receive_subset

!==================================================================================================================================! 
!=========================================================================================================== mpi_worker_send_results
!==================================================================================================================================!
    subroutine mpi_worker_send_results(taskTime)

    use kx_io
    use mare2d_mpi_definitions    
    use mare2dem_global
    use mare2dem_worker

    implicit none
    
    real(8), intent(in) :: taskTime
    
    integer(8)  :: t0,t1
        
    integer         :: ierr
    character(256)  :: sFmt
    character(32)   :: stime
    
    call system_clock(t0)        

    call mpi_comm_rank( mcomm, myID, ierr )
        
    call mpi_send( iRefinementGrp, 1, MPI_INTEGER,  manager, tag_W_SendingData, mcomm, ierr )  
 
    call mpi_send( dm, nd_local, MPI_DOUBLE_PRECISION,  manager, tag_W_SendingData, mcomm, ierr ) 
   
    if (linversion) then
            if ((sType /= 'mt').and.(lUseInversionMeshCoarsening)) then
                call mpi_send( nFree, 1, MPI_INTEGER,  manager, tag_W_SendingData, mcomm, ierr )  
                call mpi_send( iFreeNewToOld(1:nFree), nFree, MPI_INTEGER,  manager, tag_W_SendingData, mcomm, ierr )  
            endif

         call mpi_send( wj, nd_local*nFree, MPI_DOUBLE_PRECISION,  manager, tag_W_SendingData, mcomm, ierr )      
    endif

    call mpi_send(taskTime, 1, MPI_DOUBLE_PRECISION, manager, tag_W_SendingData, mcomm, ierr)

    call system_clock(t1)   
    
    write(stime,'(a8,f9.4,a2)') ' Timer: ',(t1-t0)/clockRate , ' s'
    sFmt = '(a5,2x,i6,2x,       a6,2x,i6,2x,   69x, a32,     2x,a19)'     
    if (lprintMPItimers) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,   'mpi_worker_send_results:',   trim(stime )
     
 ! these should be deallocated in mpi routines...
    if ( allocated( dp_local ) )            deallocate( dp_local ) 
    if ( allocated( d_local ) )             deallocate( d_local )     
    if ( allocated( dm ) )                  deallocate( dm )     
    if ( allocated( wj ) )                  deallocate( wj )        
      
           
    end subroutine mpi_worker_send_results

!==================================================================================================================================! 
!===================================================================================================================== mpi_all_Bcast
!==================================================================================================================================!
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))       
    subroutine mpi_all_Bcast

! Distribute a few of the input arrays for parallel processing with ScaLAPACK. 
 
    use mare2d_mpi_definitions      
    use Occam  
    
    integer :: ierr,myID
    
    call mpi_comm_rank( mcomm, myID, ierr )
    
        
    call mpi_bcast(nParams, 1, MPI_INTEGER, manager, mcomm, ierr)
 
    if (myID /= manager) allocate(prewts(nParams))
 
    call mpi_bcast(prewts, nParams, MPI_DOUBLE_PRECISION, manager, mcomm, ierr)
 
    call mpi_bcast(nDiff, 1, MPI_INTEGER, manager, mcomm, ierr)
    if (nDiff > 1) then
         if (myID /= manager) allocate(preDiffwts(nParams),ijDiff(nParams,2))     ! note array has nParams size, but only nDiff non zero entries
         call mpi_bcast(preDiffwts, nParams, MPI_DOUBLE_PRECISION, manager, mcomm, ierr)
         call mpi_bcast(ijDiff, 2*nParams, MPI_INTEGER, manager, mcomm, ierr)
    endif
    
    ! KWK debug: only coding CSR format here...
    call mpi_bcast(pen%nnz,   1, MPI_INTEGER, manager, mcomm, ierr)
    call mpi_bcast(pen%nrows, 1, MPI_INTEGER, manager, mcomm, ierr)
    if (myID /= manager)  allocate(pen%colind(pen%nnz), pen%val(pen%nnz), pen%rowptr(pen%nrows+1) )   
 
    call mpi_bcast(pen%colind, pen%nnz    , MPI_INTEGER,           manager, mcomm, ierr) 
    call mpi_bcast(pen%val   , pen%nnz    , MPI_DOUBLE_PRECISION , manager, mcomm, ierr) 
    call mpi_bcast(pen%rowptr, pen%nrows+1, MPI_INTEGER,           manager, mcomm, ierr) 
    
    call mpi_bcast( lMGS,     1, MPI_LOGICAL,          manager, mcomm,ierr )   
    call mpi_bcast( beta_mgs, 1, MPI_DOUBLE_PRECISION, manager, mcomm,ierr )   
     
    end subroutine mpi_all_Bcast
    
#endif  
     
!==================================================================================================================================! 
!=============================================================================================================== mpi_shutDownWorkers
!==================================================================================================================================!  
    subroutine mpi_shutDownWorkers()
    
    use mare2d_mpi_definitions
    use mare2dem_global
    
    implicit none
 
    integer :: ierr, iworker
    
    write(*,*) ' '
    write(*,*) ' Shutting down the worker processors:'
    write(*,*) ' '
    
    do iworker= 1,nworkers
        
        ! Tell worker to shut down:
        call mpi_send( iflag_TellWorkerToQuit, 1, MPI_INTEGER, iworker, tag_M_SendingData, mcomm, ierr )    
       
    enddo    

    write(*,*) ' '  

    end subroutine mpi_shutDownWorkers    

!==================================================================================================================================! 
!====================================================================================================================== exitMARE2DEM
!==================================================================================================================================!    
    subroutine exitMARE2DEM()
!    
! Handles a graceful and elegant exit from MARE2EM rather than letting 
! chaos ensue 
!    
    integer :: ierr
  
!
! Shut down all workers:
!
    call mpi_shutDownWorkers()
    
!
! Finalize MPI
!    
    call mpi_finalize(ierr)
   
!
! Stop...hammer time:
! 
    stop
  
! ..can't touch this...
    
    end subroutine exitMARE2DEM
    
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! 
! Wrapper routines for Occam matrix-matrix and linear solves using Scalapack:
!
! multiplyATA_scalapack()
! solveLinearSystem_scalapack()
!
!-----------------------------------------------------------------------------------------------------------------------------------
#if .not. (defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64))
    subroutine multiplyATA_scalapack()
!----------------------------------------------------------------------------------------------------------------------------------- 
    
    use mare2d_mpi_definitions
    use mare2dem_global
    use occam 
    
    implicit none 

!
! Local variables:
!
    integer :: iWorker, ierr
          
    !
    ! Launch Scalapack Matrix Multiply subroutine on Workers
    !
    do iWorker= 1,nWorkers
        call mpi_send( iflag_TellWorkerToScaLapackMM, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
    enddo 
    
    !
    ! Launch Scalapack subroutine on manager as well
    !    
    call mpi_Worker_ScaLapackMM( )

 
    end subroutine multiplyATA_scalapack

!-----------------------------------------------------------------------------------------------------------------------------------     
    subroutine solveLinearSystem_scalapack(istat,mu)
!----------------------------------------------------------------------------------------------------------------------------------- 
    
    use mare2d_mpi_definitions  
    use mare2dem_global   
    use occam  
    
    implicit none
    
    integer, intent(out) :: istat
    real(8), intent(in)  :: mu
    
    integer :: iWorker, ierr
       
    !
    ! Launch Scalapack Matrix Solver subroutine on Workers
    !
    do iWorker= 1,nWorkers
        call mpi_send( iflag_TellWorkerToScaLapackMS, 1, MPI_INTEGER, iWorker, tag_M_SendingData, mcomm, ierr )
    enddo 
    
    !
    ! Launch Scalapack subroutine on manager as well
    !    
    call mpi_Worker_ScaLapackMS(istat,mu)
        
          
    end subroutine solveLinearSystem_scalapack
    
#endif

