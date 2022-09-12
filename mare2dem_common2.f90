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
!=================================================================================================================== mare2dem_global
!==================================================================================================================================! 
    module mare2dem_global
 
    use EM_constants
 
    use triangle_mesh  ! for trimesh structure

    character(80)   :: m2d_version = 'Version: 4.2b October 25, 2019'
     
!
!  Variables read in from the .resistivity file that are stored here  
!
    character(256)  :: resistivityFile  = ''     
    character(256)  :: outputFileRoot   = ''     
    character(256)  :: modelFileName    = ''      
    character(256)  :: penaltyFileName  = ''          
    character(256)  :: dataFileName     = ''
    character(256)  :: settingsFileName = ''                 
    real(8)         :: lowerBoundGlobal, upperBoundGlobal  ! There are used to initialize the bound arrays from module Occam
    character(256)  :: cParamHeader     = ''    ! This is not used by MARE2D, just passes thru to the output files.


    integer                               :: nRhoPerRegion = 1  ! 1, 2 or 3 depending on anisotropy   
    integer                               :: nRegions               
    integer                               :: nFree = 0         ! Number of conductivities in rhoParams that are free parameters     
    real(8), dimension(:,:), allocatable  :: boundsTemp, PrejTemp, ratioTemp  ! input arrays, keep them around for easy output 
    
! 
! Variables that can be changed using optional inputs in the settings file:
!
    logical         :: lprintDecomposition     = .true.   ! set to true to display the parallel decomposition settings                                            
    logical         :: lprintMPItimers         = .true.  ! set to true to display timing info on MPI send and recv commands
    logical         :: lPrintDebug             = .true.   
    logical         :: lPrintSetup             = .true.   
    logical         :: lPrintGroups            = .true.   ! prints all Fq, RxTx, Refinement and KxFq groups for debugging
    logical         :: lSaveLoadBalanceTimers  = .true.   ! saves loadBalanceTimers.txt with load timers for each worker process    
    logical         :: lSaveTaskTimers         = .true.   ! saves TaskTimer.txt with timers for each task   
    logical         :: lPrintData              = .true.    ! display data file when read in
    logical         :: lPrintBanner            = .true.    ! display banner at start
    
    integer         :: maxnadapt  = 30         ! maximum number of refinements        
    integer         :: nwave      = 30         ! # of wavenumbers for 2.5D Fourier transformation. 
    real(8)         :: loglower   = -5         ! log10 lower limit of wavenumbers
    real(8)         :: logupper   = -1         ! log10 upper limit of wavenumbers
    
    character(256)  :: scratchFolder  = '/tmp'   
    
    logical         :: lUseInversionMeshCoarsening = .true.  ! Use a moving window to coarsen the inversion mesh to the left and 
                                                             ! right of the region defined by a parallel subset of CSEM receivers 
                                                             ! and transmitters on a given processor. ***Only free parameters are 
                                                             ! coarsened and this is only applied for CSEM modeling since it could
                                                             ! be detrimental for MT which has long wavelength lateral dependencies.
                                                             ! Set this to false if you are using larger blocky parameters
                                                             ! for inversion (e.g., a boundary model).
   
    real(8), parameter :: meshCoarsenPadding_CS = 5000. ! m of padding to add to left/right sides of Rx-Tx array for mesh subset window.
    real(8), parameter :: meshCoarsenPadding_DC = 200.  ! m of padding to add to left/right sides of Rx-Tx array for mesh subset window.
   
    !
    ! Group parameters for data decomposition parallelization:
    !
    integer         :: nTxPerGroupCSEM   = 10     ! number of transmitters to group together for mesh refinement 
    integer         :: nRxPerGroupCSEM   = 40     ! number of receivers to group together for mesh refinement
    integer         :: nRxPerGroupMT     = 40     ! number of receivers to group together for mesh refinement
    integer         :: nTxPerGroupDC     = 1      
    integer         :: nRxPerGroupDC     = 8                      
    integer         :: nKxPerGroup       = 5      ! number of wavenumbers to group together for mesh refinement
    integer         :: nFreqPerGroupCSEM = 1      ! number of frequencies to group together for mesh refinement
    integer         :: nFreqPerGroupMT   = 1      ! number of frequencies to group together for mesh refinement
    
    !
    ! For finite Tx and Rx dipoles we need to specify the quadrature order to use (can be overridden by mare2dem.settings file):
    !   
    ! * Use an odd number:
    integer         :: nQuadTxCSEM = 9, nQuadRxCSEM = 9,  nQuadRxMT = 9          

!
! Variables set from command line arguments:
!
    logical         :: lFwdFields  = .false.  ! output all possible field components to a file. use MARE2DEM -FF to make true.  
     
 
!-----------------------------------------------------------------------------------------------------------------------------------
! Everything below here is an internal variable only
!    

!
! MPI worker status array ( true means ready for next task)
!   
    integer                                 :: nworkers
    logical,dimension(:), allocatable       :: lworker_status 
     
    type(trimesh)                           :: inputmodel, inputmesh ! see call_triangle.f90 for the definition of derived type


    integer                                 :: nMT, nCSEM, nDC           ! number of MT and CSEM data    
   
    logical, dimension(:,:,:), allocatable  :: lDataMaskCSEM ! nRxCSEM x nFreqCSEM x nTx. true where data file has data    
    logical, dimension(:,:),   allocatable  :: lDataMaskMT   ! nRxMT   x nFreqMT. true where data file has data    
    logical, dimension(:,:),   allocatable  :: lDataMaskDC   ! nRxDC   x nTxDC. true where data file has data  
  
         
    logical                                 :: linversion  = .false.  
 
    real(8), dimension(:),allocatable       :: wavenum    ! array of wavenumbers 
                       
!
! Arrays for RxTx and Frequency groups:
!
    type  :: RxTxG   ! Derived type to store each checkboard subgrid of Rx and Tx combinations:
        integer                             :: nTx      ! number of transmitters in this group 
        integer, dimension(:), allocatable  :: iTx      ! index of global transmitter number, used to get x,y,z of tx
        character(2)                        :: sType    !'cs','mt','dc'
        integer                             :: nRx
        integer, dimension(:), allocatable  :: iRx      ! index of global receiver numbers           
        
        integer, dimension(:),   allocatable  :: iRxTrodes  ! DC: local to global trode index. Used for passing trode x to worker
        integer, dimension(:,:), allocatable  :: iRxDC      ! DC: [nRx x 2] index to local trode array. Pass this to the worker
        integer, dimension(:),   allocatable  :: iTxTrodes  ! DC: local to global trode index. Used for passing trode x to worker
        integer, dimension(:,:), allocatable  :: iTxDC      ! DC: [nRx x 2] index to local trode array. Pass this to the worker   
             
    end type
 
    integer                                 :: nRxTxGroups
    type(RxTxG), dimension(:), allocatable  :: RxTxGroups   
    
    type  :: FqG   ! Derived type to store each checkboard subgrid of Rx and Tx combinations:
        integer                             :: nFq      ! number of frequencies in this group   
        integer, dimension(:), allocatable  :: iFq      ! index of global frequencies for this group
        !integer                             :: iFqRefine  ! local index to iFq for refinement frequency for this group         
        character(2)                        :: sType    !'cs', 'mt' or 'dc' 
    end type
 
    integer                                 :: nFqGroups 
    type(FqG), dimension(:), allocatable    :: FqGroups   

    type  :: RefG   ! Derived type refinement groups storing the index of the RxTxGroup, Frequency group and subset of dp() 
        integer                              :: iRxTxGroup    ! index to the group of Rx and Tx
        integer                              :: iFqGroup      ! index to the group of Kx and Fq
        character(2)                         :: sType         !'cs','mt','dc'
        integer                              :: nd            ! number of input data in this RefG
        integer, dimension(:,:), allocatable :: dp            ! (nd, cnDataParams) array subset of Occam's dp. Mapped to local iRx and iTx
        logical                              :: lReturned = .false. ! true when returned from worker         
        real(8)                              :: timer                            
    end type
    
    integer                                 :: nRefinementGroups   
    integer                                 :: link1, link2
    integer                                 :: iPtr_refGroups   ! points to next refinement group to send out
    type(RefG), dimension(:), allocatable   :: refinementGroups ! nRefinementGroups rows of  iRxTxGroup , iFqGroup

   
!
!  A few variables used by the worker nodes:
!
    integer     :: iRefinementGrp
    
    integer, dimension(:), allocatable    :: iFreeNewToOld
    
    real(8), dimension(:,:), allocatable  :: freeRegionCentroids
       
    contains
    
    subroutine deallocate_mare2dem_global  

    if ( allocated( boundsTemp ) )          deallocate( boundsTemp )
    if ( allocated( PrejTemp ) )            deallocate( PrejTemp )
    if ( allocated( ratioTemp ) )           deallocate( ratioTemp ) 
    if ( allocated( lworker_status ) )      deallocate (lworker_status)
    if ( allocated( inputmodel%attr ) )     call deallocate_trimesh(inputmodel,.false.)
    if ( allocated( inputmesh%attr ) )      call deallocate_trimesh(inputmesh,.false.)
    if ( allocated( inputmesh%area ) )      deallocate ( inputmesh%area  )
    if ( allocated( inputmodel%area ) )     deallocate ( inputmodel%area  )
    if ( allocated( wavenum ) )             deallocate( wavenum ) 
    if ( allocated( RxTxGroups ) )          deallocate( RxTxGroups )
    if ( allocated( FqGroups ) )            deallocate( FqGroups ) 
    if ( allocated( refinementGroups ) )    deallocate( refinementGroups )  
    if ( allocated( iFreeNewToOld ) )       deallocate( iFreeNewToOld )  
    if ( allocated( freeRegionCentroids ) ) deallocate( freeRegionCentroids )  
      
            
    end subroutine deallocate_mare2dem_global    
        
    end module mare2dem_global   
 
!==================================================================================================================================! 
!======================================================================================================== mare2dem_input_data_params
!==================================================================================================================================! 
    module mare2dem_input_data_params
    
! While the data parameter array dp is passed in/out of Occam, it is only used
! in the forward routines.  So we set cnDataParams in the forward routines
! and then allocate the dp array in readData.
!        
    integer         :: cnDataParams
 
    character(256)  :: cUTMLine        = ''         ! the UTM line from the data file (not used by MARE2D, just passes thru)

    character(4)    :: phaseConvention = 'lag'      ! 'lag' or 'lead'
    character(4)    :: reciprocityUsed = ' '        ! 'yes' or 'no'. if yes, em2dkx uses (-i/omega) scaling for e_rx from m_tx i*omega for m_rx from j_tx
    
!
! Data dependent Tx-Rx variables read from data file in subroutine readData:
!
    ! CSEM:
    integer                                   :: nTxCSEM=0, nRxCSEM=0, nFreqCSEM=0
    real(8), dimension(:), allocatable        :: azimuthTxCSEM,dipTxCSEM,lengthTxCSEM
    real(8), dimension(:), allocatable        :: xTxCSEM,yTxCSEM,zTxCSEM,xRxCSEM,yRxCSEM,zRxCSEM        
    real(8), dimension(:), allocatable        :: ThetaRxCSEM,AlphaRxCSEM, BetaRxCSEM  
    real(8), dimension(:), allocatable        :: lengthRxCSEM      ! nRx  for x,y,z dipoles all same length. Use new Rx if otherwise 
    real(8), dimension(:), allocatable        :: fTxCSEM  

    character(8),   dimension(:), allocatable :: cSourceType  
    character(128), dimension(:), allocatable :: cRxNamesCSEM
    character(128), dimension(:), allocatable :: cTxNamesCSEM  
    
    ! MT:
    integer                                   :: nRxMT=0, nFreqMT=0
    real(8), dimension(:), allocatable        :: xRxMT,yRxMT,zRxMT   
    real(8), dimension(:), allocatable        :: ThetaRxMT,AlphaRxMT, BetaRxMT, lengthRxMT      
    real(8), dimension(:), allocatable        :: fTxMT  
    integer, dimension(:), allocatable        :: iSolveStatic  ! 1 = TE + TM, 2 = TE only, 3 = TM only. otherwise ignored
   
    character(128), dimension(:), allocatable :: cRxNamesMT
        
    ! DC Resistivity:
    integer                                   :: nTxDC=0, nRxDC=0,nTrodesRxDC=0,nTrodesTxDC=0
    real(8), dimension(:,:), allocatable      :: trodes_RxDC,trodes_TxDC  ! n x 3 positions of electrodes for Rx and Tx
    integer, dimension(:,:), allocatable      :: RxDC,TxDC  ! nRx x 2, nTx x 2 listings of electrodes for Rx and Tx    

    character(128), dimension(:), allocatable :: cRxNamesDC
    character(128), dimension(:), allocatable :: cTxNamesDC
    
    contains
    
    subroutine deallocate_mare2dem_input_data_params
    
    if (allocated (fTxCSEM))        deallocate( fTxCSEM ) 
    if (allocated (xTxCSEM) )       deallocate( xTxCSEM, yTxCSEM, zTxCSEM, azimuthTxCSEM, dipTxCSEM, cSourceType, lengthTxCSEM )  
    if (allocated (cRxNamesCSEM) )  deallocate( cRxNamesCSEM )
    if (allocated (cTxNamesCSEM) )  deallocate( cTxNamesCSEM )

    if (allocated (fTxMT))          deallocate( fTxMT )    
    if (allocated (xRxCSEM))        deallocate( xRxCSEM, yRxCSEM, zRxCSEM, ThetaRxCSEM, AlphaRxCSEM, BetaRxCSEM, lengthRxCSEM )   
    if (allocated (xRxMT))          deallocate( xRxMT, yRxMT, zRxMT, ThetaRxMT, AlphaRxMT, BetaRxMT, lengthRxMT )         
    if (allocated (cRxNamesMT) )    deallocate( cRxNamesMT )
    if (allocated (iSolveStatic) )  deallocate( iSolveStatic )
        
    if (allocated (trodes_RxDC) )   deallocate( trodes_RxDC, trodes_TxDC, RxDC, TxDC)  
    if (allocated (cRxNamesDC) )    deallocate( cRxNamesDC )
    if (allocated (cTxNamesDC) )    deallocate( cTxNamesDC )
        
    end subroutine deallocate_mare2dem_input_data_params
        
    end module mare2dem_input_data_params
     
!==================================================================================================================================! 
!=============================================================================================================== mare2dem_deallocate
!==================================================================================================================================!     
    subroutine mare2dem_deallocate
!
! Deallocates all I/O variables for MARE2DEM and a few others. This routine is called at the end of Occam.
!        
    use kx_io 
    use mare2dem_global
    use mare2dem_input_data_params
    
    call deallocate_kx_io 
 
    call deallocate_mare2dem_global
    
    call deallocate_mare2dem_input_data_params

      
    end subroutine mare2dem_deallocate

!==================================================================================================================================! 
!=========================================================================================================== display_MARE2DEM_Params
!==================================================================================================================================!     
    subroutine display_MARE2DEM_Params
 
    use kx_io
    use mare2dem_global
    
    implicit none
    
    character(80)   :: ctemp

    write(6,*) '========== Parameters settings for MARE2DEM =============================='
    write(*,*) ' '     
 
 
       
    write(*,*) ' '
    write(*,*) ' ----------------------------'
    write(*,*) '  Adaptive refinement params:'
    write(*,*) ' ----------------------------'    
    write(ctemp,'(i5)') maxnadapt     
    write(6,fmt='(a32,a)') 'Max # refinements:  ',trim(adjustl(ctemp))         
    
    write(ctemp,'(f8.3)') errortolerance
    write(*,fmt='(a32,a)') 'Tolerance (%):  ',trim(adjustl(ctemp))     

    
    write(ctemp,'(f5.1)') minQangle 
    write(6,fmt='(a32,a)') 'Mesh quality angle:  ',trim(adjustl(ctemp))      
   
    
    write(ctemp,'(i5)') idual_func      
    write(*,fmt='(a32,a)') 'Dual function:  ', trim(adjustl(ctemp))         

    if (lUseInversionMeshCoarsening) then
        write(6,fmt='(a32,a3)') 'Use mesh coarsening:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Use mesh coarsening:  ','no' 
    endif  

     
    if (lSaveMeshFiles) then
        write(6,fmt='(a32,a3)') 'Save meshes:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Save meshes:  ','no' 
    endif 
    if (lprintDebug) then
        write(6,fmt='(a32,a3)') 'Print debug:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Print debug:  ','no' 
    endif     
    if (lDisplayRefinementStats) then
        write(6,fmt='(a32,a3)') 'Print adaptive:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Print adaptive:  ','no' 
    endif  

    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '         Wavenumber settings:'
    write(*,*) ' ----------------------------'     
    write(*,fmt='(a32,f6.1,1a,f6.1,1a,i6,32a)') 'Wavenumbers:  ', loglower,',', logupper,',', nwave, '   ! log10 lower, upper, # '  
!    write(6,fmt='(a32,a)') ' (Co)Sine Transform Filters:  ', trim(FCTfilter)  
    write(6,fmt='(a32,a)') ' Linear Solver:  ', trim(linearSolver)      
 
    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '         MT settings:'
    write(*,*) ' ----------------------------'    
            
    if (lMTscatteredField) then
        write(6,fmt='(a32,a3)') 'Use mt scattered field:  ','yes'  
    else
        write(6,fmt='(a32,a2)') 'Use mt scattered field:  ','no' 
    endif   
    
    
    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '    Error Estimator settings:'
    write(*,*) ' ----------------------------'     
    write(ctemp,'(f8.2)') minRangeProtector
    write(6,fmt='(a32,a)') 'Minimum error range:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(i5)') max_nsubrefine
    write(6,fmt='(a32,a)') 'Max # subrefinements:  ', trim(adjustl(ctemp))

    write(ctemp,'(f8.2)') pct_refine
    write(6,fmt='(a32,a)') 'Percent refine:  ', trim(adjustl(ctemp))

    write(ctemp,'(e11.4)') minArea
    write(6,fmt='(a32,a)') 'Minimum area:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(i0)') maxMeshNodes
    write(6,fmt='(a32,a)') 'Max # mesh nodes:  ', trim(adjustl(ctemp)) 
    
    write(ctemp,'(e11.4)') ecutoff
    write(6,fmt='(a32,a)') 'E noise floor:  ', trim(adjustl(ctemp)) 

    write(ctemp,'(e11.4)') hcutoff
    write(6,fmt='(a32,a)') 'H noise floor:  ', trim(adjustl(ctemp))   


    write(*,*) ' '    
    write(*,*) ' ----------------------------'
    write(*,*) '      Finite dipole settings:'
    write(*,*) ' ----------------------------'   
    
    write(ctemp,'(i5)') nQuadTxCSEM
    write(6,fmt='(a32,a)') ' Transmitter quadrature order:  ', trim(adjustl(ctemp))
     
    write(ctemp,'(i5)') nQuadRxCSEM
    write(6,fmt='(a32,a)') 'CSEM receiver quadrature order: ', trim(adjustl(ctemp)) 
    
    write(ctemp,'(i5)') nQuadRxMT
    write(6,fmt='(a32,a)') ' MT receiver quadrature order:  ', trim(adjustl(ctemp))    
    
    write(*,*) ' ' 
    

     
    end subroutine display_MARE2DEM_Params
    
!==================================================================================================================================! 
!======================================================================================================================== checkModel
!==================================================================================================================================!       
    subroutine checkModel 
    
    use kx_io
    use mare2dem_global    
    use triangle_mesh
    
    implicit none
 
    logical :: lHasSlivers
         
    character(32) :: cend
    !
    write(*,*) '======== Inspecting Input Values ========================================='
    write(*,*) ' '     
    
    !
    ! First check to make sure that if model has negative conductivities (i.e., indices to params), then the rho file has
    ! been defined and read in:
    !
 
    if ( nRhoParams == 0) then
 
        write(*,*) ' '
        write(*,*) ' Error, no resistivity parameters defined! '
        write(*,*) ' Please include a resistivity file!'
        write(*,*) ' Stopping!'
        stop        
  
    endif
    
    !
    ! Check the model for slivers:
    !
 
    write(*,*) ' ... Checking the input model for slivers ...'

    call check_model(inputmodel,minQangle, lHasSlivers)
 
    if (lHasSlivers) then
    
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '        
        write(cend,fmt='(F4.0)') minQangle
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) '!!!     Warning Model has slivers... stopping        !!!'
        write(*,*) '!!! Please fix the mesh at the locations noted above !!!'
        write(*,*) '!!! so that segment intersections are > ',trim(adjustl(cend)),' degrees  !!!'  
        write(*,*) '!!!    >>> Better luck next time, bucko <<<          !!!' ! creds to CJW
        write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '
        write(*,*) ' '                
        stop
    else
        write(*,*) ' '
        write(*,*) ' A-okay buddy:'
        write(*,*) ' '
        write(*,*) ' Your model has been checked for slivers and  '
        write(*,*) ' none were found.  '
        write(*,*) ' '   
    endif

    end subroutine checkModel       
         
!==================================================================================================================================! 
!==================================================================================================================== getWavenumbers
!==================================================================================================================================!       
    subroutine getWavenumbers
 
    use mare2dem_global
    use mare2dem_input_data_params ! for nFreqCSEM
    
    implicit none
 
    integer :: i 
    real(8) ::dlg
    
    if ((nfreqCSEM < 1).and.(nTxDC<1)) return
    
    if (lprintSetup) then
        write(6,*) '========== Setting up wavenumbers ==========='
        write(*,*) ' '     
    endif

    allocate (wavenum(nwave))
    dlg = abs(logupper - loglower) / (nwave - 1)  
    
    do i = 1,nwave
        wavenum(i) = 10**(min(loglower,logupper) +(i-1)*dlg)
    enddo  
    
    if (lprintSetup) then  
        write(*,'(a32)') ' Wavenumbers:  '
        do i = 1,nwave            
             write(*,'(i32,1x,es12.3)') i,wavenum(i)
        enddo    
        write(*,*) ' '  
    endif    
    
    end subroutine getWavenumbers

!==================================================================================================================================! 
!===================================================================================================================== getRxTxGroups
!==================================================================================================================================!       
    subroutine getRxTxGroups 
!
! Sets up the Rx-Tx groups:    
!
    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none
 
 
    integer ::  nTxGroupsCSEM, nRxGroupsCSEM, nRxGroupsMT, nTxGroupsDC, nRxGroupsDC, iGroup
    integer ::  itxG, ntxg, nt0, irxG, nrx0, nrxg, i,j,icnt, iAn,iBn
 
    integer, dimension(:),allocatable :: renumTrodes ! for DC resistivity trode wrangling

!
! Get total number of RxTxGroups:
!
    nTxGroupsCSEM = 0
    if (nTxCSEM > 0 ) nTxGroupsCSEM = ceiling( dble(nTxCSEM) / dble(nTxPerGroupCSEM) )   
    
    nTxGroupsDC = 0
    if (nTxDC > 0 )  nTxGroupsDC = ceiling( dble(nTxDC) / dble(nTxPerGroupDC) )   
     
    nRxGroupsCSEM     = ceiling( dble(nRxCSEM) / dble(nRxPerGroupCSEM) )   
    nRxGroupsMT       = ceiling( dble(nRxMT)   / dble(nRxPerGroupMT) )   
    nRxGroupsDC       = ceiling( dble(nRxDC)   / dble(nRxPerGroupDC) )  
    
    nRxTxGroups       = nRxGroupsCSEM*nTxGroupsCSEM + nRxGroupsMT + nTxGroupsDC*nRxGroupsDC
 
!
! Allocate storage for RxTxGroups
!   
    allocate ( RxTxGroups(nRxTxGroups) )
 
!
! CSEM RxTxGroups:
!
    iGroup = 0  
    do itxG = 1,nTxGroupsCSEM
    
        nt0     = (itxG-1)*nTxPerGroupCSEM               ! last input transmitter number from previous group
        ntxg    = min(nTxPerGroupCSEM, nTxCSEM - nt0)    ! number of transmitters in this group
 
        do irxG = 1,nRxGroupsCSEM
    
            nrx0    =  (irxG-1)*nRxPerGroupCSEM                 ! last input receiver number from previous group or 0 if first group
            nrxg    =  min(nRxPerGroupCSEM, nRxCSEM - nrx0 )    ! number of receivers in this group
            
            !
            ! Insert arrays into RxTxGroups(iGroup):
            !
            
            iGroup = iGroup + 1
            
            RxTxGroups(iGroup)%sType = 'cs'
                        
            !
            ! Tx indices:
            !
            RxTxGroups(iGroup)%nTx  = ntxg
            
            allocate( RxTxGroups(iGroup)%iTx(ntxg) )
 
            RxTxGroups(iGroup)%iTx  = [1:ntxg] + nt0    
             
            
            !
            ! Rx indices:
            !            
            RxTxGroups(iGroup)%nRx  = nrxg
            
            allocate( RxTxGroups(iGroup)%iRx(nrxg) )
            
            RxTxGroups(iGroup)%iRx  = [1:nrxg] + nrx0    
 
        enddo
    enddo
    
!
! MT RxTxGroups:
!
    if (nFreqMT > 0 ) then
    
        do irxG = 1,nRxGroupsMT
            
            nrx0  =  (irxG-1)*nRxPerGroupMT             ! last input receiver number from previous group or 0 if first group
            nrxg =  min(nRxPerGroupMT, nRxMT - nrx0 )   ! number of receivers in this group    
        
            iGroup = iGroup + 1
    
            !
            ! Tx stuff:
            !                           
            RxTxGroups(iGroup)%nTx      = 1         
                
            allocate ( RxTxGroups(iGroup)%iTx(1)  )  !kwk debug: need to update this for hybrid MT where E and H at separate Rx
            RxTxGroups(iGroup)%iTx(1)   = 0   
            RxTxGroups(iGroup)%sType = 'mt'
         
            !
            ! Rx stuff:
            !
            RxTxGroups(iGroup)%nRx      = nrxg
            allocate( RxTxGroups(iGroup)%iRx(nrxg) )
                
            RxTxGroups(iGroup)%iRx      =  [1:nrxg] + nrx0  
                      
                        
        enddo
    endif    

!
! DC RxTxGroups:
!
    do itxG = 1,nTxGroupsDC
    
        nt0     = (itxG-1)*nTxPerGroupDC               ! last input transmitter number from previous group
        ntxg    = min(nTxPerGroupDC, nTxDC - nt0)      ! number of transmitters in this group
 
        do irxG = 1,nRxGroupsDC
    
            nrx0    =  (irxG-1)*nRxPerGroupDC               ! last input receiver number from previous group or 0 if first group
            nrxg    =  min(nRxPerGroupDC, nRxDC - nrx0 )    ! number of receivers in this group
            
            !
            ! Insert arrays into RxTxGroups(iGroup):
            !
            
            iGroup = iGroup + 1
            
            RxTxGroups(iGroup)%sType = 'dc'
            
            !
            ! Tx indices:
            !
            
            RxTxGroups(iGroup)%nTx  = ntxg
            
            allocate( RxTxGroups(iGroup)%iTx(ntxg) )
 
            RxTxGroups(iGroup)%iTx  = [1:ntxg] + nt0    
             
                                   
            !
            ! Rx indices:
            !            
            RxTxGroups(iGroup)%nRx  = nrxg
            
            allocate( RxTxGroups(iGroup)%iRx(nrxg) )
            
            RxTxGroups(iGroup)%iRx  = [1:nrxg] + nrx0    
            
            
            ! New:
             
              
            ! Map Tx electrode indices to local electrodes:
            
            ! loop through electrode pair listing:
            allocate(renumTrodes(nTrodesTxDC))
            renumTrodes = 0
            icnt = 0
            do i = 1,ntxg
                do j = 1,2
                    ! mark electrodes that are used:
                    if (renumTrodes(TxDC(i+nt0,j)) == 0 ) then
                        icnt = icnt+1
                        renumTrodes(TxDC(i+nt0,j)) = icnt   
                    endif
                enddo            
            enddo
            ! renumTrodes is global to local index
            ! now make local to global index
            allocate(RxTxGroups(iGroup)%iTxTrodes(icnt))
     
            do i = 1,nTrodesTxDC
                if (renumTrodes(i)/=0) then
                    RxTxGroups(iGroup)%iTxTrodes(renumTrodes(i)) = i  ! local to global index
                    !write(*,'(a,3(i4,1x))') 'i,icnt,renumTrodes(i):',i,icnt,renumTrodes(i)
                endif
            enddo
            !write(*,*) 'iGroup,icnt,ntxg: ',iGroup,icnt,ntxg
            ! Now make new Tx array for with indices to local trodes:
            allocate(RxTxGroups(iGroup)%iTxDC(ntxg,2))
            do i = 1,ntxg
                do j = 1,2                    
                    RxTxGroups(iGroup)%iTxDC(i,j) = renumTrodes(TxDC(i+nt0,j)) 
                enddo
!                write(*,'(i4,1x,a,2(i4,1x),a,2(i4,1x))') i,' old: ',TxDC(i+nt0,1:2),'new: ',RxTxGroups(iGroup)%iTxDC(i,1:2)
            enddo
            deallocate(renumTrodes)
            ! So now we have:
            !  RxTxGroups(iGroup)%iTxTrodes  - the local to global trode index. Will use this for passing trode positions to worker
            !  RxTxGroups(iGroup)%iTxDC      - [ntxg x 2] the index to local trode array. Will pass this to the worker. 
                      
            ! Map Rx electrode indices to local electrodes:
            
            ! loop through electrode pair listing:
            allocate(renumTrodes(nTrodesRxDC))
            renumTrodes = 0
            icnt = 0
            do i = 1,nrxg
                do j = 1,2
                    ! mark electrodes that are used:
                    if (renumTrodes(RxDC(i+nrx0,j)) == 0 ) then
                        icnt = icnt+1
                        renumTrodes(RxDC(i+nrx0,j)) = icnt   
                    endif
                enddo            
            enddo
            ! renumTrodes now has 1:ict where trodes are used in this RxTxGroup
            ! now make index array for local to global trode index
            allocate(RxTxGroups(iGroup)%iRxTrodes(icnt))
            do i = 1,nTrodesRxDC
                if (renumTrodes(i)/=0) then
                    RxTxGroups(iGroup)%iRxTrodes(renumTrodes(i)) = i  ! local to global index
                endif
            enddo
            ! Now make new Rx array for with indices to local trodes:
            allocate(RxTxGroups(iGroup)%iRxDC(nrxg,2))
            do i = 1,nrxg
                do j = 1,2                    
                    RxTxGroups(iGroup)%iRxDC(i,j) = renumTrodes(RxDC(i+nrx0,j))
                enddo
            enddo
            deallocate(renumTrodes)
            ! So now we have:
            !  RxTxGroups(iGroup)%iRxTrodes  - the local to global trode index. Will use this for passing trode positions to worker
            !  RxTxGroups(iGroup)%iRxDC      - [n x 2] the index to local trode array. Will pass this to the worker. 
            
            !
            ! Print Check:
            !    
!             write(*,*) 'Transmitters:'
!             do i = 1,ntxg
!          
!                 iAn = RxTxGroups(iGroup)%iTxDC(i,1)
!                 iBn = RxTxGroups(iGroup)%iTxDC(i,2)
!                 write(*,*) 'iAn,iBn: ',iAn,iBn
!                  write(*,'(a,i4,1x,a,2(i4,1x),a,2(i4,1x))') 'i+nt0: ',i+nt0, 'Old: ',TxDC(i+nt0,1:2),'New:', RxTxGroups(iGroup)%iTxTrodes(iAn),RxTxGroups(iGroup)%iTxTrodes(iBn)
!             enddo
            
        enddo
    enddo
    
    if (lPrintGroups) then
    
        write(*,*) ' ' 
        write(*,*) '============== Receiver and Transmitter Groups: =========================='
        write(*,*) ' ' 
        write(*,*) ' # Groups: ',nRxTxGroups
        do i = 1,nRxTxGroups
            write(*,*) ' ' 
            write(*,'(a15,1x,i4)') 'Group:',i
            write(*,'(a15,3x,a2)') ' sType:', trim(RxTxGroups(i)%sType)
            write(*,'(a15,1x,i4)') ' nRx:', RxTxGroups(i)%nRx
            write(*,'(a15,1x,i0,a,i0)') 'Receivers:',RxTxGroups(i)%iRx(1),' to ',RxTxGroups(i)%iRx(RxTxGroups(i)%nRx)
            write(*,'(a15,1x,i4)') ' nTx:', RxTxGroups(i)%nTx
            write(*,'(a15,1x,i0,a,i0)') 'Transmitters:',RxTxGroups(i)%iTx(1),' to ',RxTxGroups(i)%iTx(RxTxGroups(i)%nTx)
!             if (RxTxGroups(i)%sType =='dc') then ! print out electrodes...
!             
!             endif
!         integer, dimension(:),   allocatable  :: iRxTrodes  ! DC: local to global trode index. Used for passing trode x to worker
!         integer, dimension(:,:), allocatable  :: iRxDC      ! DC: [nRx x 2] index to local trode array. Pass this to the worker
!         integer, dimension(:),   allocatable  :: iTxTrodes  ! DC: local to global trode index. Used for passing trode x to worker
!         integer, dimension(:,:), allocatable  :: iTxDC      ! DC: [nRx x 2] index to local trode array. Pass this to the worker   

        enddo 
        
    endif
        
    end subroutine getRxTxGroups
    
!==================================================================================================================================! 
!======================================================================================================================= getFqGroups
!==================================================================================================================================!        
    subroutine getFqGroups 
!
! Sets up the Fq groups:    
!
    use mare2dem_global
    use mare2dem_input_data_params
       
    implicit none
 
    integer :: nFreqGroupsCSEM,nFreqGroupsMT,nFreqGroupsDC         ! Total number of frequency groups    

    integer :: iGroup, i 
    integer :: ifreqG, nf0, nfg
 
    
!
! Get total number of FqGroups:
!    
 
    nFreqGroupsCSEM = 0
    if (nTxCSEM > 0 )  nFreqGroupsCSEM = ceiling( dble(nfreqCSEM)  / dble(nFreqPerGroupCSEM)   )    

    nFreqGroupsDC = 0
    if (nTxDC > 0 )  nFreqGroupsDC = 1
     
    nFreqGroupsMT = 0         
    if (nfreqMT > 0 ) nFreqGroupsMT = ceiling( dble(nfreqMT)    / dble(nFreqPerGroupMT) ) 
 
    nFqGroups = nFreqGroupsCSEM + nFreqGroupsMT + nFreqGroupsDC
  
    allocate( FqGroups(nFqGroups) ) 
    
 !
 ! CSEM:
 !   
    iGroup = 0 
    
    do ifreqG = 1,nFreqGroupsCSEM       
            
        nf0 =  (ifreqG-1)*nFreqPerGroupCSEM          ! last input frequency number from previous group
        nfg = min(nFreqPerGroupCSEM, nfreqCSEM - nf0)  ! number of freqs in this group
     
        iGroup = iGroup + 1
         
        FqGroups(iGroup)%nFq = nfg
        
        allocate ( FqGroups(iGroup)%iFq(nfg) )
        
        do i = 1,nfg
            
            FqGroups(iGroup)%iFq(i) = nf0 + i    
            
            !if ( mod(i+ floor(dble(nfg)/dble(2) ),nfg) == 0 ) FqGroups(iGroup)%iFqRefine = i
              
        enddo
                        
        FqGroups(iGroup)%sType = 'cs'     
                            
    enddo
    
 !
 ! MT:
 !   
    do ifreqG = 1,nFreqGroupsMT   
    
        nf0 = (ifreqG-1)*nFreqPerGroupMT          ! last input frequency number from previous group
        nfg = min(nFreqPerGroupMT, nfreqMT - nf0) 

        iGroup = iGroup + 1

        FqGroups(iGroup)%nFq = nfg
                
        allocate ( FqGroups(iGroup)%iFq(nfg) )
 
        do i = 1,nfg
        
            FqGroups(iGroup)%iFq(i) = nf0 + i    
            !if ( mod(i+ floor(dble(nfg)/dble(2) ),nfg) == 0 ) FqGroups(iGroup)%iFqRefine = i
          
        enddo     
        
        FqGroups(iGroup)%sType = 'mt'     
             
    enddo
    
 !
 ! DC:
 !   
    if (nFreqGroupsDC == 1) then
 
        iGroup = iGroup + 1

        FqGroups(iGroup)%nFq = 1
        
        allocate ( FqGroups(iGroup)%iFq(1) )
        FqGroups(iGroup)%iFq(1) = 1 
         
        FqGroups(iGroup)%sType = 'dc'     
                 
    endif  
    
    if (lPrintGroups) then
    
        write(*,*) ' ' 
        write(*,*) '======================= Frequency Groups: ================================'
        write(*,*) ' ' 
        write(*,*) ' # Groups: ',nFqGroups
        do i = 1,nFqGroups
            write(*,*) ' ' 
            write(*,'(a15,1x,i4)') 'Group:',i
            write(*,'(a15,3x,a2)')  ' sType:', trim(FqGroups(i)%sType)    
            write(*,'(a15,1x,i4)') ' nFq:', FqGroups(i)%nFq
            write(*,'(a15,1x,*(i4,1x))') ' Frequencies:',FqGroups(i)%iFq(1:FqGroups(i)%nFq)
            !write(*,'(a15,1x,i4)') ' Refine index:',   FqGroups(i)%iFqRefine
        enddo 
        
    endif
        
    end subroutine getFqGroups     
!==================================================================================================================================! 
!========================================================================================================================== between! 
!==================================================================================================================================!     
    logical function between(iTest,iLow,iHigh) 
    
    integer, intent(in) :: iTest, iLow, iHigh
    if ( (iTest >= iLow) .and.  (iTest <= iHigh) ) then
        between = .true.
    else
        between = .false.
    endif
    
    end function between    
    
          
!==================================================================================================================================! 
!============================================================================================================ getFreeRegionCentroids
!==================================================================================================================================!       
    subroutine getFreeRegionCentroids()
!
! Gets free parameter region centroids. Assumes convex parameter shape, so don't make crazy meshes with non convex parameters 
! since the centroid can lie outside the parameter. 
!    
 
    use mare2dem_global
    use kx_io
    use fem2d_utilities
 
    implicit none
 
    character (24) :: cend,tricommand 
    
    integer :: e,iregion, iparamnum(4),i,j,iFree, n(3)
    
    integer, dimension(:), allocatable :: iTriCounter, iFreeRegion
    
    real(8) :: area,ye(3),ze(3),a(3),b(3),c(3)
    real(8), dimension(:),   allocatable  :: sumArea
    logical :: lFree
   
    
    !
    ! Create a triangulation of input model:
    !
        
    tricommand = 'q0panjA'//CHAR(0)
           
    call call_triangle(tricommand,inputmodel,mesh)  
    
    !
    ! Create lookup table for free region indices: (note polygon regions have one or more associated parameters (e.g. > 1 for anisotropy) )
    !
    allocate( iFreeRegion(nRegions))
    iFreeRegion = 0
    
    iFree = 0 ! counter for free regions
    do i=1,nRegions
        lFree = .false.
        do j = 1,nRhoPerRegion
            if (iFreeParam((i-1)*nRhoPerRegion+j)>0) lFree = .true.
        enddo    
        if (lFree) then
            iFree = iFree + 1
            iFreeRegion(i) = iFree
        endif
    enddo
    
    !
    ! Traverse mesh, average triangle centroids for each parameter to get its centroid:
    !   
    allocate( sumArea(nFree/nRhoPerRegion),freeRegionCentroids(nFree/nRhoPerRegion,2) )

    freeRegionCentroids = 0
    sumArea = 0
 
    do e = 1,mesh%nele
           
        ! get polygon region containing element e:
        iregion = nint(abs(mesh%attr(e)))
 
       
        iFree = iFreeRegion(iregion)
        if (iFree > 0)  then 
        
            ! get area of element e:
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)   
            call get_abc_coeffs(ye,ze,a,b,c,area)   
            
            sumArea(iFree) = sumArea(iFree) + area
            
            freeRegionCentroids(iFree,1) = freeRegionCentroids(iFree,1) + area*sum(ye)/3d0
            freeRegionCentroids(iFree,2) = freeRegionCentroids(iFree,2) + area*sum(ze)/3d0
             
        endif
    
        
    enddo
    
    do i = 1,nFree/nRhoPerRegion
        freeRegionCentroids(i,1:2) = freeRegionCentroids(i,1:2) / sumArea(i)
        !write(*,*) freeRegionCentroids(i,1:2)
    enddo
     
    deallocate(iFreeRegion,sumArea)
       

    
    !deallocate
    call deallocate_trimesh(mesh,.false.)
    
    end subroutine getFreeRegionCentroids   
    
!==================================================================================================================================! 
!====================================================================================================================== getDataMasks
!==================================================================================================================================!     
    subroutine getDataMasks
    
    use mare2dem_global
    use mare2dem_input_data_params
    use Occam, only: nd,dp

    implicit none
    
    integer :: i,iRx, iTx, iFq
            
    !
    ! Specify CSEM data mask array:
    !
    if (nTxCSEM > 0 ) then
        
        allocate(lDataMaskCSEM(nRxCSEM,nFreqCSEM,nTxCSEM))  ! kwk debug: this could be generated at bottom of readData
        
        lDataMaskCSEM = .false.
        
        do i = 1,nd
            iRx = dp(i,4)  
            iTx = dp(i,3) 
            iFq = dp(i,2)            
            if ( dp(i,1) < 100 ) lDataMaskCSEM(iRx,iFq,iTx) = .true. 
        enddo

        if (lFwdFields) lDataMaskCSEM = .true.    
        
    endif        

    !
    ! DC resistivity data masks:
    ! 

    if (nTxDC > 0 ) then
        
        allocate(lDataMaskDC(nRxDC,nTxDC))  ! kwk debug: this could be generated at bottom of readData
        
        lDataMaskDC = .false.
        
        do i = 1,nd
            iRx = dp(i,4)  
            iTx = dp(i,3) 
                    
            if ( (dp(i,1) > 200).and.(dp(i,1) < 300) ) lDataMaskDC(iRx,iTx) = .true. 
        enddo

        if (lFwdFields) lDataMaskDC = .true.    
        
    endif     
    
    !
    ! MT data masks:
    !
    if (nFreqMT > 0 ) then
    
        allocate(lDataMaskMT(nRxMT,nFreqMT))  
        
        !         Type        Freq#        Tx#           Rx#   
        lDataMaskMT = .false.
        do i = 1,nd
            iRx = dp(i,4)
            iTx = dp(i,3)
            iFq = dp(i,2)
            
            if ( (dp(i,1) > 100).and.(dp(i,1) < 200) ) then ! only consider MT data here:
                 lDataMaskMT(iRx,iFq) = .true. 
                 if (iTx > 0 ) lDataMaskMT(iTx,iFq) = .true.  ! for magnetics at receiver iTx
            endif
        enddo
        
        if (lFwdFields) lDataMaskMT = .true.
        
    endif
    
    end subroutine getDataMasks        
!==================================================================================================================================! 
!=============================================================================================================== getRefinementGroups
!==================================================================================================================================!      
    subroutine getRefinementGroups
!
! Sets up the refinement groups by looking at RxTx and Fq combos that
! have input data. Also creates the data mask arrays. 
!
    use mare2dem_global
    use mare2dem_input_data_params
    use Occam, only: nd,dp
    
    implicit none  
    
    integer :: iRxTxG, iFqG, rx0,rx1,tx0,tx1,fq0,fq1, iPass, iGrp,i,j,iA,iB,iM,iN
    integer :: nFq, nRx, nTx, iD, ndCount
    logical :: between
    
    integer, dimension(:), allocatable :: g2L_iFq, g2L_iTx, g2L_iRx 
 
    character(2) :: sType_RxTx, sType_Fq  
    
 !
 ! Generate the data mask arrays which indicate which Rx,Tx,Freq combos have input data:
 !
     write(*,*) 'antes nd: ',nd

     call getDataMasks

     write(*,*) 'depois nd: ',nd
     
 !
 ! Use a first pass to count how many groups to allocate, then we use a second pass to allocate and insert:
 !   
 
    do iPass = 1,2
    
        if (iPass == 2) then
            allocate( refinementGroups(nRefinementGroups) )
        endif    
        
        nRefinementGroups = 0 
        
        do iRxTxG = 1,nRxTxGroups
        
            nRx = RxTxGroups(iRxTxG)%nRx
            nTx = RxTxGroups(iRxTxG)%nTx
            
            rx0 = RxTxGroups(iRxTxG)%iRx(1)
            rx1 = RxTxGroups(iRxTxG)%iRx(nRx)
            tx0 = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
            tx1 = RxTxGroups(iRxTxG)%iTx(nTx)
            
            sType_RxTx = RxTxGroups(iRxTxG)%sType
             
            do iFqG = 1,nFqGroups
        
                nFq = FqGroups(iFqG)%nFq
                fq0 = FqGroups(iFqG)%iFq(1)
                fq1 = FqGroups(iFqG)%iFq(nFq)
               
                sType_Fq = FqGroups(iFqG)%sType
                
                !
                ! MT refinement groups
                !    
                
                if ( ( sType_RxTx == 'mt' ) .and. ( sType_Fq == 'mt' ) ) then
                
                    if ( any( lDataMaskMT(rx0:rx1,fq0:fq1)) )  then
                        
                        nRefinementGroups = nRefinementGroups + 1  
                        
                        if (iPass == 2) then
                        
                            refinementGroups(nRefinementGroups)%iRxTxGroup = iRxTxG
                            refinementGroups(nRefinementGroups)%iFqGroup   = iFqG
                            refinementGroups(nRefinementGroups)%sType      = 'mt'                        
                            
                        endif    
                
                    endif  
            
                !
                ! CSEM refinement groups
                !
                elseif  ( ( sType_RxTx == 'cs' ) .and. ( sType_Fq == 'cs' ) ) then 
               
                    if ( any( lDataMaskCSEM(rx0:rx1,fq0:fq1,tx0:tx1)) )  then

                        nRefinementGroups = nRefinementGroups + 1  
                        
                        if (iPass == 2) then
                        
                            refinementGroups(nRefinementGroups)%iRxTxGroup = iRxTxG
                            refinementGroups(nRefinementGroups)%iFqGroup   = iFqG
                            refinementGroups(nRefinementGroups)%sType      = 'cs' 
                                        
                        endif    
                     
                    endif
                    
                !
                ! DC refinement groups
                !
                elseif  ( ( sType_RxTx == 'dc' ) .and. ( sType_Fq == 'dc' ) ) then
               
                    if ( any( lDataMaskDC(rx0:rx1,tx0:tx1)) )  then

                        nRefinementGroups = nRefinementGroups + 1  
                        
                        if (iPass == 2) then
                        
                            refinementGroups(nRefinementGroups)%iRxTxGroup = iRxTxG
                            refinementGroups(nRefinementGroups)%iFqGroup   = iFqG
                            refinementGroups(nRefinementGroups)%sType      = 'dc' 
                                        
                        endif    
                     
                    endif         
                          
                endif
                
            enddo ! iFq
                       
        enddo ! iRxTxG
 
    enddo ! iPass   

!
! Point to first group:
!
    iPtr_refGroups = 1

    write(*,*) 'nd aqui1: ',nd
    
!
! Lastly, create the local dp arrays for each refinement group:
!    
    allocate( g2L_iFq(max(nFreqCSEM,nFreqMT,1)), g2L_iTx(max(nTxCSEM,nTxDC)), g2L_iRx(max(nRxCSEM,nRxMT,nRxDC)) )

    write(*,*) 'nd aqui2: ',nd

    do iGrp = 1,nRefinementGroups
    
        iRxTxG  = refinementGroups(iGrp)%iRxTxGroup 

        nRx     = RxTxGroups(iRxTxG)%nRx          
        rx0     = RxTxGroups(iRxTxG)%iRx(1)
        rx1     = RxTxGroups(iRxTxG)%iRx(nRx)

        nTx     = RxTxGroups(iRxTxG)%nTx  
        tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
        tx1     = RxTxGroups(iRxTxG)%iTx(nTx)

        iFqG    = refinementGroups(iGrp)%iFqGroup

        nFq     = FqGroups(iFqG)%nFq
        fq0     = FqGroups(iFqG)%iFq(1)
        fq1     = FqGroups(iFqG)%iFq(nFq)
        
        write(*,'(a,*(i,1x))') ' nRx,rx0,rx1: ', nRx,rx0,rx1
        write(*,'(a,*(i,1x))') ' nTx,tx0,tx1: ', nTx,tx0,tx1
        write(*,'(a,*(i,1x))') ' nFq,fq0,fq1: ', nFq,fq0,fq1

        ! Get global to local mappings:
        g2L_iFq = 0
        g2L_iTx = 0
        g2L_iRx = 0 
        
        if (refinementGroups(iGrp)%sType == 'mt' ) then
            
            tx0 = rx0
            tx1 = rx1
            
            do i = 1,nFreqMT
                do j=1,nFq
                    if (i ==  FqGroups(iFqG)%iFq(j) ) g2L_iFq(i) = j 
                enddo 
            enddo
 
            do i = 1,nRxMT
                do j=1,nRx
                    if (i ==  RxTxGroups(iRxTxG)%iRx(j) ) g2L_iRx(i) = j 
                enddo             
            enddo
  
            ! Count number of input data for this group:
            ndCount = 0
            do iD = 1,nd
                if ( (between(dp(iD,2),fq0,fq1)) .and. &
                   & (between(dp(iD,3),tx0,tx1)) .and. &  !kwk debug: need to update this for hybrid MT responses
                   & (between(dp(iD,4),rx0,rx1)) ) then
                   if  ( (dp(iD,1)>100).and. (dp(iD,1)<200) )  ndCount = ndCount + 1
                endif
            enddo       

        elseif (refinementGroups(iGrp)%sType == 'cs' ) then
 
            do i = 1,nFreqCSEM
                do j=1,nFq
                    if (i ==  FqGroups(iFqG)%iFq(j) ) g2L_iFq(i) = j 
                enddo 
            enddo
            do i = 1,nRxCSEM
                do j=1,nRx
                    if (i ==  RxTxGroups(iRxTxG)%iRx(j) ) g2L_iRx(i) = j 
                enddo             
            enddo
            do i = 1,nTxCSEM
                do j=1,nTx
                    if (i ==  RxTxGroups(iRxTxG)%iTx(j) ) g2L_iTx(i) = j 
                enddo             
            enddo      

 
            ! Count number of input data for this group:
            ndCount = 0
            write(*,*) 'fq0: ',fq0
            write(*,*) 'fq1: ',fq1
            write(*,*) 'tx0: ',tx0
            write(*,*) 'tx1: ',tx1
            write(*,*) 'rx0: ',rx0
            write(*,*) 'rx1: ',rx1
            do iD = 1,nd
                if ( (between(dp(iD,2),fq0,fq1)) .and. &
                   & (between(dp(iD,3),tx0,tx1)) .and. & 
                   & (between(dp(iD,4),rx0,rx1)) ) then
                   if  (dp(iD,1)<100)  ndCount = ndCount + 1
                endif
            enddo
            write(*,*) 'Ver ndCount: ',ndCount
            write(*,*) 'nd: ',nd
                    
        elseif (refinementGroups(iGrp)%sType == 'dc' ) then
 
!             do i = 1,1
!                 do j=1,nFq
!                     if (i ==  FqGroups(iFqG)%iFq(j) ) g2L_iFq(i) = j 
!                 enddo 
!             enddo
            g2L_iFq(1) = 1 
             
            do i = 1,nRxDC
                do j=1,nRx
                    if (i ==  RxTxGroups(iRxTxG)%iRx(j) ) g2L_iRx(i) = j 
                enddo             
            enddo
            do i = 1,nTxDC
                do j=1,nTx
                    if (i ==  RxTxGroups(iRxTxG)%iTx(j) ) g2L_iTx(i) = j 
                enddo             
            enddo    
            ! Count number of input data for this group:
            ndCount = 0
            !write(*,*) 'fq0,fq1,tx0,tx1,rx0,rx1:',fq0,fq1,tx0,tx1,rx0,rx1
            do iD = 1,nd
                !write(*,*) 'iD,dp(iD,2:4): ' ,iD,dp(iD,2:4)
                if ( (between(dp(iD,2),fq0,fq1)) .and. &
                   & (between(dp(iD,3),tx0,tx1)) .and. &   
                   & (between(dp(iD,4),rx0,rx1)) ) then
                   if  ( (dp(iD,1)>200).and.(dp(iD,1)<300) ) ndCount = ndCount + 1
                endif
            enddo       
            !write(*,*) 'ndCount: ',ndCount                     
                                 
        endif

   
        ! Now allocate and fill in the dp values for the group using local indices:
        
        refinementGroups(iGrp)%nd = ndCount
        allocate(refinementGroups(iGrp)%dp(ndCount,cnDataParams+1))
        
        ndCount = 0
        do iD = 1,nd
            if ( (between(dp(iD,2),fq0,fq1)) .and. &
               & (between(dp(iD,3),tx0,tx1)) .and. & 
               & (between(dp(iD,4),rx0,rx1)) ) then
               
               !write(*,*) ' i am here: ',refinementGroups(iGrp)%sType,dp(iD,1)
               if ( ((refinementGroups(iGrp)%sType == 'mt') .and. (dp(iD,1) > 100).and.(dp(iD,1) < 200) ) .or. &
                  & ((refinementGroups(iGrp)%sType == 'cs') .and. (dp(iD,1) >   0).and.(dp(iD,1) < 100) ) .or. &
                  & ((refinementGroups(iGrp)%sType == 'dc') .and. (dp(iD,1) > 200).and.(dp(iD,1) < 300) ) ) then
    
                    ndCount = ndCount + 1
                    refinementGroups(iGrp)%dp(ndCount,1) = dp(iD,1)
                    refinementGroups(iGrp)%dp(ndCount,2) = g2L_iFq(dp(iD,2)) ! map from global to local index
                
                    if (refinementGroups(iGrp)%sType == 'mt')  then         !kwk debug: I think this only works when all MT Rx in same group...
                        refinementGroups(iGrp)%dp(ndCount,3) = g2L_iRx(dp(iD,3))  ! for different H location !kwk debug: need to fix stuff above to allow hybrid MT
                    else
                        refinementGroups(iGrp)%dp(ndCount,3) = g2L_iTx(dp(iD,3))
                    endif
                
                    !write(*,'(5(i5,1x))') iD, g2L_iFq(dp(iD,2)),g2L_iRx(dp(iD,3)),g2L_iRx(dp(iD,4))
                
                    refinementGroups(iGrp)%dp(ndCount,4) = g2L_iRx(dp(iD,4))
                    refinementGroups(iGrp)%dp(ndCount,5) = iD                ! index back to global data array      
                
                endif         
           endif
        enddo          
                              
    enddo

    deallocate( g2L_iFq, g2L_iTx, g2L_iRx )     
    
        
    if (lPrintGroups) then
        
        write(*,*) ' ' 
        write(*,*) '====================== Refinement Groups: ================================'
        write(*,*) ' ' 
        write(*,*) ' # Groups: ',nRefinementGroups
        do i = 1,nRefinementGroups
            write(*,*) ' ' 
            write(*,'(a15,1x,i4)') 'Group:',i
            write(*,'(a15,3x,a2)') ' sType:', refinementGroups(i)%sType
            iRxTxG = refinementGroups(i)%iRxTxGroup
            nRx     = RxTxGroups(iRxTxG)%nRx          
            rx0     = RxTxGroups(iRxTxG)%iRx(1)
            rx1     = RxTxGroups(iRxTxG)%iRx(nRx)   
            nTx     = RxTxGroups(iRxTxG)%nTx  
            tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
            tx1     = RxTxGroups(iRxTxG)%iTx(nTx)        
                  
            write(*,'(a15,1x,i4,2x,2(a,i0,a,i0))') ' iRxTxGroup:',iRxTxG,'  Rx: ',rx0,' to ',rx1,'   Tx: ',tx0,' to ',tx1
 
            iFqG    = refinementGroups(i)%iFqGroup
            nFq     = FqGroups(iFqG)%nFq
            fq0     = FqGroups(iFqG)%iFq(1)
            fq1     = FqGroups(iFqG)%iFq(nFq)        
        
            write(*,'(a15,1x,i4,2x,(a,i0,a,i0))') '   iFqGroup:',iFqG, '  Fq: ',fq0,' to ' ,fq1
            
            write(*,'(a15,1x,i4)') ' nd:',refinementGroups(i)%nd  
            
!             if (refinementGroups(i)%sType =='dc') then
!                 write(*,'(16x,10(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG', 'y_A','y_B','y_M','y_N'
!                 do j = 1,refinementGroups(i)%nd  
!                     iA = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),1) 
!                     iB = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),2)    
!                     iM = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),1) 
!                     iN = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),2)  
!                     write(*,*) iA,iB,iM,iN
!                     
!                     write(*,'(16x,6(i4,1x),4(f8.1,1x))') j,refinementGroups(i)%dp(j,:),&
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iA),2),  &
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iB),2), &
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iM),2), &  
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iN),2) 
!                 enddo
!             else 
                
                write(*,'(16x,6(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG'
                do j = 1,refinementGroups(i)%nd  
                    write(*,'(16x,*(i4,1x))') j,refinementGroups(i)%dp(j,:)
                enddo
            
!             endif
            
        enddo

        write(*,*) ' ' 
        write(*,*) '====================== modifing Refinement Groups: ================================'
        write(*,*) ' ' 
        write(*,*) ' # Groups: ',nRefinementGroups
        do i = 1,nRefinementGroups  
            write(*,*) ' ' 
            write(*,'(a15,1x,i4)') 'Group:',i
            write(*,'(a15,3x,a2)') ' sType:', refinementGroups(i)%sType
            iRxTxG = refinementGroups(i)%iRxTxGroup
            nRx     = RxTxGroups(iRxTxG)%nRx          
            rx0     = RxTxGroups(iRxTxG)%iRx(1)
            rx1     = RxTxGroups(iRxTxG)%iRx(nRx)   
            nTx     = RxTxGroups(iRxTxG)%nTx  
            tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
            tx1     = RxTxGroups(iRxTxG)%iTx(nTx)        
                  
            write(*,'(a15,1x,i4,2x,2(a,i0,a,i0))') ' iRxTxGroup:',iRxTxG,'  Rx: ',rx0,' to ',rx1,'   Tx: ',tx0,' to ',tx1
 
            iFqG    = refinementGroups(i)%iFqGroup
            nFq     = FqGroups(iFqG)%nFq
            fq0     = FqGroups(iFqG)%iFq(1)
            fq1     = FqGroups(iFqG)%iFq(nFq)        
        
            write(*,'(a15,1x,i4,2x,(a,i0,a,i0))') '   iFqGroup:',iFqG, '  Fq: ',fq0,' to ' ,fq1
            
            write(*,'(a15,1x,i4)') ' nd:',refinementGroups(i)%nd

            if (refinementGroups(i)%nd < 20) then
                ! right now the idea is to link the following small group to the previous none.


            
!             if (refinementGroups(i)%sType =='dc') then
!                 write(*,'(16x,10(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG', 'y_A','y_B','y_M','y_N'
!                 do j = 1,refinementGroups(i)%nd  
!                     iA = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),1) 
!                     iB = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),2)    
!                     iM = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),1) 
!                     iN = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),2)  
!                     write(*,*) iA,iB,iM,iN
!                     
!                     write(*,'(16x,6(i4,1x),4(f8.1,1x))') j,refinementGroups(i)%dp(j,:),&
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iA),2),  &
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iB),2), &
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iM),2), &  
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iN),2) 
!                 enddo
!             else 
                
                write(*,'(16x,6(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG'
                do j = 1,refinementGroups(i)%nd  
                    write(*,'(16x,*(i4,1x))') j,refinementGroups(i)%dp(j,:)
                enddo
            
!             endif
            
        enddo

        write(*,*) ' ' 
        write(*,*) '====================== Refinement Groups modified: ================================'
        write(*,*) ' ' 
        write(*,*) ' # Groupsm: ',nRefinementGroups
        do i = 1,nRefinementGroups
            write(*,*) ' ' 
            write(*,'(a15,1x,i4)') 'Group:',i
            write(*,'(a15,3x,a2)') ' sType:', refinementGroups(i)%sType
            iRxTxG = refinementGroups(i)%iRxTxGroup
            nRx     = RxTxGroups(iRxTxG)%nRx          
            rx0     = RxTxGroups(iRxTxG)%iRx(1)
            rx1     = RxTxGroups(iRxTxG)%iRx(nRx)   
            nTx     = RxTxGroups(iRxTxG)%nTx  
            tx0     = RxTxGroups(iRxTxG)%iTx(1)    ! note for MT, iTx(1) = 0
            tx1     = RxTxGroups(iRxTxG)%iTx(nTx)        
                  
            write(*,'(a15,1x,i4,2x,2(a,i0,a,i0))') ' iRxTxGroup:',iRxTxG,'  Rx: ',rx0,' to ',rx1,'   Tx: ',tx0,' to ',tx1
 
            iFqG    = refinementGroups(i)%iFqGroup
            nFq     = FqGroups(iFqG)%nFq
            fq0     = FqGroups(iFqG)%iFq(1)
            fq1     = FqGroups(iFqG)%iFq(nFq)        
        
            write(*,'(a15,1x,i4,2x,(a,i0,a,i0))') '   iFqGroup:',iFqG, '  Fq: ',fq0,' to ' ,fq1
            
            write(*,'(a15,1x,i4)') ' nd:',refinementGroups(i)%nd  
            
!             if (refinementGroups(i)%sType =='dc') then
!                 write(*,'(16x,10(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG', 'y_A','y_B','y_M','y_N'
!                 do j = 1,refinementGroups(i)%nd  
!                     iA = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),1) 
!                     iB = RxTxGroups(iRxTxG)%iTxDC(refinementGroups(i)%dp(j,3),2)    
!                     iM = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),1) 
!                     iN = RxTxGroups(iRxTxG)%iRxDC(refinementGroups(i)%dp(j,4),2)  
!                     write(*,*) iA,iB,iM,iN
!                     
!                     write(*,'(16x,6(i4,1x),4(f8.1,1x))') j,refinementGroups(i)%dp(j,:),&
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iA),2),  &
!                     trodes_TxDC(RxTxGroups(iRxTxG)%iTxTrodes(iB),2), &
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iM),2), &  
!                     trodes_RxDC(RxTxGroups(iRxTxG)%iRxTrodes(iN),2) 
!                 enddo
!             else 
                
                write(*,'(16x,6(a4,1x))') '#','iType','iFqL','iTxL','iRxL','idG'
                do j = 1,refinementGroups(i)%nd  
                    write(*,'(16x,*(i4,1x))') j,refinementGroups(i)%dp(j,:)
                enddo
            
!             endif
            
        enddo
        
    endif
    
    end subroutine getRefinementGroups
       
!==================================================================================================================================! 
!================================================================================================================ printDecomposition
!==================================================================================================================================!       
    subroutine printDecomposition

    use mare2dem_global
    use mare2dem_input_data_params
    
    implicit none  

    integer :: nTxGroupsCSEM, nTxGroupsDC, nRxGroupsCSEM, nRxGroupsDC, nRxGroupsMT
    integer :: nFreqGroupsCSEM,nFreqGroupsMT, nFreqGroupsDC          
    
    integer :: ict,i,j,k
    
    character (24) :: ctemp, str1, str2
         
    nTxGroupsCSEM = 0
    if (nTxCSEM > 0 ) nTxGroupsCSEM = ceiling( dble(nTxCSEM) / dble(nTxPerGroupCSEM) )   
    
    nTxGroupsDC = 0
    if (nTxDC > 0 )  nTxGroupsDC = ceiling( dble(nTxDC) / dble(nTxPerGroupDC) )   
    
    nRxGroupsCSEM     = ceiling( dble(nRxCSEM) / dble(nRxPerGroupCSEM) )   
    nRxGroupsMT       = ceiling( dble(nRxMT)   / dble(nRxPerGroupMT) )   
    nRxGroupsDC       = ceiling( dble(nRxDC)   / dble(nRxPerGroupDC) )  

    nRxTxGroups       = nRxGroupsCSEM*nTxGroupsCSEM + nRxGroupsMT + nTxGroupsDC*nRxGroupsDC

    nFreqGroupsCSEM = 0
    if (nTxCSEM > 0 )  nFreqGroupsCSEM   = ceiling( dble(nfreqCSEM)  / dble(nFreqPerGroupCSEM)   )    

    nFreqGroupsDC = 0
    if (nTxDC > 0 )  nFreqGroupsDC = 1
        
    nFreqGroupsMT = 0         
    if (nfreqMT > 0 )  nFreqGroupsMT = ceiling( dble(nfreqMT)    / dble(nFreqPerGroupMT) ) 
    
    nFqGroups = nFreqGroupsCSEM + nFreqGroupsMT  + nFreqGroupsDC
 
    if (lprintDecomposition) then
        write(*,*) ' '    
        write(6,*) '============== Parallel decomposition settings: =========================='
        write(*,*) ' '    
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nTxPerGroupCSEM
        if (nTxCSEM > 0 ) write(6,fmt='(a32,a)') 'CSEM Transmitters per group:  ', trim(adjustl(ctemp))
 
        if (nTxDC > 0 ) write(ctemp,'(i5)') nTxPerGroupDC
        if (nTxDC > 0 ) write(6,fmt='(a32,a)') 'DC Transmitters per group:  ', trim(adjustl(ctemp))
                
        write(ctemp,'(i5)') nRxPerGroupMT
        if (nfreqMT > 0 ) write(6,fmt='(a32,a)') 'MT Receivers per group:  ', trim(adjustl(ctemp))

        write(ctemp,'(i5)') nRxPerGroupCSEM
        if (nTxCSEM > 0 )  write(6,fmt='(a32,a)') 'CSEM Receivers per group:  ', trim(adjustl(ctemp))

        write(ctemp,'(i5)') nRxPerGroupDC
        if (nRxDC > 0 )  write(6,fmt='(a32,a)') 'DC Receivers per group:  ', trim(adjustl(ctemp))
                        
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nFreqPerGroupCSEM
        if (nTxCSEM > 0 ) write(6,fmt='(a32,a)') 'CSEM frequencies per group:  ', trim(adjustl(ctemp))
        
        
        if (nfreqMT > 0 )  write(ctemp,'(i5)') nFreqPerGroupMT  
        if (nfreqMT > 0 ) write(6,fmt='(a32,a)') 'MT frequencies per group:  ', trim(adjustl(ctemp))
    
        write(*,*) ' '
        write(*,*) ' Maximum possible groups if all Tx-Rx-Frequency '
        write(*,*) ' combinations are present in the input data: '
        write(*,*) ' '
        
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nTxGroupsCSEM
        if (nTxCSEM > 0 ) write(*,fmt='(a32,a)') '# CSEM transmitter groups:  ', trim(adjustl(ctemp))

        if (nTxDC > 0 ) write(ctemp,'(i5)') nTxGroupsDC
        if (nTxDC > 0 ) write(*,fmt='(a32,a)') '# DC transmitter groups:  ', trim(adjustl(ctemp))
        
        write(ctemp,'(i5)') nRxGroupsCSEM
        write(*,fmt='(a32,a)') '# CSEM receiver groups:  ', trim(adjustl(ctemp))
 
         write(ctemp,'(i5)') nRxGroupsDC
        write(*,fmt='(a32,a)') '# DC receiver groups:  ', trim(adjustl(ctemp))
        
        write(ctemp,'(i5)') nRxGroupsMT
        write(*,fmt='(a32,a)') '# MT receiver groups:  ', trim(adjustl(ctemp))
 
               
        if (nTxCSEM > 0 ) write(ctemp,'(i5)') nFreqGroupsCSEM
        if (nTxCSEM > 0 ) write(*,fmt='(a32,a)') '# CSEM frequency groups:  ', trim(adjustl(ctemp))
        
        if (nfreqMT > 0 ) write(ctemp,'(i5)') nFreqGroupsMT
        if (nfreqMT > 0 ) write(*,fmt='(a32,a)') '# MT frequency groups:  ', trim(adjustl(ctemp))
        
        write(*,fmt='(a32,a)') ' -','--'
        write(ctemp,'(i7)') nTxGroupsCSEM*nRxGroupsCSEM*nFreqGroupsCSEM + nRxGroupsMT*nFreqGroupsMT + nTxGroupsDC*nRxGroupsDC
        write(*,fmt='(a32,a)') ' # refinement groups possible:  ', trim(adjustl(ctemp))    
        write(*,*) ' '
        
        
        write(*,*) ' After inspecting the input data, I have found: '
        write(*,*) ' ' 
        write(ctemp,'(i7)') nRefinementGroups
        write(*,fmt='(a32,a)') ' # of refinement groups:  ', trim(adjustl(ctemp))      
        write(*,*) ' '  
        write(*,*) ' '  
        
        if (nTxCSEM > 0 ) then
    
            ! compute fill-in of CSEM data:
            ict = 0
            do i=1,nRxCSEM
                do j = 1,nFreqCSEM
                    do k = 1,nTxCSEM
                        if ( lDataMaskCSEM(i,j,k)) ict = ict + 1
                    enddo
                enddo
            enddo
 
            write(str1,*) ict
            write(str2,*) nRxCSEM*nFreqCSEM*nTxCSEM
            write(*,'(a32,a,1x,a6,1x,a,1x,a1,f6.1,a4)') ' CSEM Data fill in: ',trim(adjustl(str1)), 'out of', trim(adjustl(str2)),&
            &'(', real(ict)/real(nRxCSEM*nFreqCSEM*nTxCSEM)*100.,' % )'         
     
        endif
        
        if (nTxDC > 0 ) then
    
            ! compute fill-in of DC data:
            ict = 0
            do i=1,nRxDC
                do k = 1,nTxDC
                    if ( lDataMaskDC(i,k)) ict = ict + 1
                enddo
            enddo
 
            write(str1,*) ict
            write(str2,*) nRxDC*nTxDC
            write(*,'(a32,a,1x,a6,1x,a,1x,a1,f6.1,a4)') ' DC Data fill in: ',trim(adjustl(str1)), 'out of', trim(adjustl(str2)),&
            &'(', real(ict)/real(nRxDC*nTxDC)*100.,' % )'         
     
        endif
                
        if ( nFreqMT > 0 ) then
            
            ! compute fill-in of MT data:
            ict = 0
            do i=1,nRxMT
                do j = 1,nFreqMT
                    if ( lDataMaskMT(i,j)) ict = ict + 1 !kwk debug: can replace loops with count() command
                enddo
            enddo
 
            write(str1,*) ict
            write(str2,*) nRxMT*nFreqMT
            write(*,'(a32,a,1x,a6,1x,a,1x,a1,f6.1,a4)') ' MT Data fill in: ', trim(adjustl(str1)), 'out of', trim(adjustl(str2)), &
                 & '(', real(ict)/real(nRxMT*nFreqMT)*100.,' % )'
                            
        endif   
        write(*,*) '   '
    endif

    if (allocated(lDataMaskCSEM)) deallocate(lDataMaskCSEM) 
    if (allocated(lDataMaskMT))   deallocate(lDataMaskMT)  
    if (allocated(lDataMaskDC))   deallocate(lDataMaskDC)                
                  
    end subroutine printDecomposition

    
!==================================================================================================================================! 
!=================================================================================================================== get_time_offset
!==================================================================================================================================!  
    subroutine  get_time_offset(timein,timeout)
!    
! timein is the clock start time or 0.
! timeout is the time in seconds since the input time
!
! Uses date_and_time Fortran intrinsic
! ---------------------------------------------------------------------
   implicit none
   
   integer, dimension(8) :: values
   integer               :: i,j,k,mjd
   
   real(8) :: timein, timeout, fracday
    

 !
 ! New standard Fortran90 Time function:
 !
    call date_and_time(values=values) !this ouputs only values, ignoring other optional arguments
    
 ! Convert year, month day to modified julian day:
      
    i = values(1)
    j = values(2)
    k = values(3)
    mjd = -678927 + k + 1461*(i+(j-14)/12)/4 + 367*(j-2-12 * &
          & ((j-14)/12))/12 + (24002-12*i-j)/1200
    
               
  ! Add on fractional day:
                ! hour            ! minute          ! sec       ! millisec
    fracday = ((values(5)*60.d0 + values(6))*60.d0 + values(7) + values(8)/1000.d0 )/86400.d0
    
    timeout = mjd + fracday
    
  ! Finally, convert timeout to time difference between it and timein:  
    timeout =  timeout*86400.d0  - timein
                         
    end  subroutine  get_time_offset
          
!==================================================================================================================================! 
!============================================================================================================================ indexx
!==================================================================================================================================!
! From N.R.
       
    SUBROUTINE indexx(n,arr,indx) 
! Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that 
! arr(indx(j)) is in ascending order for j=1,2,...,N. 
! The input quantities n and arr are not changed.
      implicit none
      INTEGER n,indx(n),M,NSTACK
      REAL(8) arr(n)
      PARAMETER (M=7,NSTACK=1000)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL(8) a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

