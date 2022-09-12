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
!=============================================================================================================== dataTransformations
!==================================================================================================================================!
!
! Module for various data transformations
!
    module dataTransformations
    
    use EM_constants
    implicit none
    
    
    contains 

!
! MT Apparent resistivity
!
        pure real(8) function getAppRes(z,freq)
        
          implicit none
          complex(8), intent(in)   :: z    ! impedance Z, not C
          real(8), intent(in)      :: freq ! linear frequency
        
          
          getAppRes = abs(z)**2  / (2d0*pi*freq*mu0)
          
        end function

        elemental real(8) function getAppResDeriv(z,dz,freq)
        
          implicit none
          complex(8), intent(in)   :: z    ! impedance Z, not C
          complex(8), intent(in)   :: dz   ! derivative of impedance Z 
          real(8), intent(in)      :: freq ! linear frequency
          
          getAppResDeriv = 2d0 /  (2d0*pi*freq*mu0) * (dble(z)*dble(dz) + aimag(z)*aimag(dz))

        end function

!
! log10 MT Apparent Resistivity
!
        elemental real(8) function getlog10AppResDeriv(z,dz)
        
          implicit none
          complex(8), intent(in)   :: z    ! impedance Z, not C
          complex(8), intent(in)   :: dz   ! derivative of impedance Z 
 
          getlog10AppResDeriv = ( 2d0* log10(exp(1d0)) / (dble(z)*dble(z) + aimag(z)*aimag(z)) )* &
                                & (dble(z)*dble(dz) + aimag(z)*aimag(dz))
        
        end function
        
!
! Phase with unwrapping:
!  
        pure real(8) function getPhase(z,dataphi)
        
        ! JH hacky method of getting around lack of modulo subtraction in Occam.
        
          implicit none
          complex(8), intent(in)        :: z
          real(8), intent(in)           :: dataphi
          real(8)                       :: tmp 
    
          tmp = atan2( aimag(z), real(z) )*rad2deg  ! phase in degrees
          
          getPhase = tmp
          if( abs(tmp + 360d0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp + 360d0)
          if( abs(tmp - 360d0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp - 360d0)
          
    
        end function
        
        elemental real(8) function getPhaseDeriv(z,dz) 
        
          implicit none
          complex(8),intent(in) :: z, dz
          getPhaseDeriv = rad2deg * (dble(z)*aimag(dz) - aimag(z)*dble(dz) ) / (dble(z)*dble(z) + aimag(z)*aimag(z) )    
          
        end function          
         
!
! Absolute value:
!
        elemental real(8) function absDeriv(phi,dphi) 
        
          implicit none
          complex(8),intent(in) :: phi, dphi
    
          absDeriv = ( dble(phi)*dble(dphi) + aimag(phi)*aimag(dphi))  / sqrt( dble(phi)*dble(phi) + aimag(phi)*aimag(phi)) 
          
        end function

!
! Log10 Absolute value:
!
        elemental real(8) function log10absDeriv(phi,dphi) 
        
          implicit none
          complex(8),intent(in) :: phi, dphi
    
          log10absDeriv = (dble(phi)*dble(dphi) + aimag(phi)*aimag(dphi))  / &
                        & ( log(10d0)*(dble(phi)*dble(phi) + aimag(phi)*aimag(phi)) )
          ! log10absDeriv = absDeriv(phi,dphi)  / (log(10.)*abs(phi))
 
 
        end function
          
!
! Polarization ellipse parameters:
!

        pure real(8) function getPE(e1,e2,comp)
        
      ! Taken from Steve's old Occam routines - does not match  Smith and Ward 
      ! but probably works OK.
        implicit none
        complex(8), intent(in)   :: e1    ! E/B-field channel 1
        complex(8), intent(in)   :: e2    ! E/B-field channel 2
        character(4), intent(in)        :: comp  ! flag  'pmax' or 'pmin' for output 
        
      ! Local:
        real(8)                  :: x1, x2, y1, y2
        real(8)                  :: a, b, phi, s, c, p1, p2, pmin, pmax


      ! bust up the complex components of the electric field:
        x1 = real(e1)
        y1 = aimag(e1)
        x2 = real(e2)
        y2 = aimag(e2)
      
      ! find the critical angle
        a   = -(x1*y1 + x2*y2)
        b   = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe 
        s  = sin(phi)
        c  = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        
        s  = sin(phi+pi/2.d0)
        c  = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))

        pmax = max(p1,p2)
        pmin = min(p1,p2)
        
        select case (comp)
        
          case('pmin')
            getPE = pmin
            
          case('pmax')
            getPE = pmax
            
        end select    

        end function  
        
        elemental real(8) function getPEDeriv(e1,e2, de1, de2,comp)
      ! Taken from Steve's old Occam routines - does not match Smith and Ward 
      ! but probably works OK
        implicit none
        complex(8), intent(in)   :: e1    ! E/B-field channel 1
        complex(8), intent(in)   :: e2    ! E/B-field channel 2
        complex(8), intent(in)   :: de1   ! E/B-field deriv channel 1
        complex(8), intent(in)   :: de2   ! E/B-field deriv channel 2
        character(4), intent(in)        :: comp  ! flag  'pmax' or 'pmin' for output 
        
        ! Local:         
        real(8)                  :: x1, x2, y1, y2, dx1, dx2, dy1, dy2
        real(8)                  :: a, b, phi, s, c, p1, d1, p2, d2
        
      ! bust up the complex components of the electric field:
        x1  = real(e1)
        y1  = aimag(e1)
        x2  = real(e2)
        y2  = aimag(e2)
        dx1 = real(de1)
        dy1 = aimag(de1)
        dx2 = real(de2)
        dy2 = aimag(de2)
      
      ! find the critical angle
        a = -(x1*y1 + x2*y2)
        b = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe and derivatives
        s = sin(phi)
        c = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        d1 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d1 = d1/p1
        
        s = sin(phi+pi/2.d0)
        c = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))
        d2 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d2 = d2/p2


        if (p1 >= p2) then
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d2
                
                case('pmax')
                getPEDeriv = d1
            
            end select
        
        else ! ( p1 < p2)
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d1
                
                case('pmax')
                getPEDeriv = d2
            
            end select
        
        end if     
        
        end function  

        
    end module dataTransformations

!==================================================================================================================================! 
!=================================================================================================================== mare2dem_worker
!==================================================================================================================================! 
    module mare2dem_worker   
    
!
! Variables stored here are *ONLY* used at the local worker cpu level.
!
    use kx_io                           ! for passing data and settings
    use em2dkx_mod                      
    use dc2dkx_mod
    use mare2dem_global                 ! for settings 
    use spline_integ_kx
    
    implicit none   

    public
 
    integer                             :: nfrequencies
    real(8), dimension(:), allocatable  :: frequencies
     
    logical, dimension(:,:,:,:), allocatable  :: lDataMask  ! iFq,iTx,iRx,iComp(1:6)   true where data file has data    

!
! Receiver and transmitter parameters:
!
    integer, public :: nRxIn                                     ! Number of receivers       
    real(8), dimension(:), allocatable, public :: xRxIn, yRxIn, zRxIn, thetaRxIn, alphaRxIn, betaRxIn, lengthRxIn ! x,y,z location of the receivers 
 
    integer, public :: nTxIn                                         
    real(8), dimension(:), allocatable, public :: xTxIn, yTxIn, zTxIn, azimuthTxIn,dipTxIn, lengthTxIn ! x,y,z location of the transmitter 
    character(8), allocatable, public          :: TxTypeIn(:)  !'edipole' or 'bdipole'
                 
    ! DC Resistivity:
    integer                                   :: nTxDCIn=0, nRxDCIn=0,nTrodesRxDCIn=0,nTrodesTxDCIn=0
    real(8), dimension(:,:), allocatable      :: trodes_RxDCIn,trodes_TxDCIn  ! n x 3 positions of electrodes for Rx and Tx
    integer, dimension(:,:), allocatable      :: RxDCIn,TxDCIn  ! nRx x 2, nTx x 2 listings of electrodes for Rx and Tx    
       
       
    type  :: KxFqG   ! Derived type to store subgrid of Kx and Fq combinations:
        integer                             :: nKxFq    ! number of non-refinement wavenumber,frequency pairs in this group   
        integer, dimension(:), allocatable  :: iFq      ! index of non-refinement frequencies 
        integer, dimension(:), allocatable  :: iKx      ! index of wavenumbers for this group
        integer                             :: iFqRefine  ! local index to freq for refinement frequency for this group         
        integer                             :: iKxRefine  ! local index to wavenumber for refinement for this group  
    end type
 
    integer                                 :: nKxFqGroups
    type(KxFqG), dimension(:), allocatable  :: KxFqGroups      
        
    integer, parameter                      :: nComponentsPerScratch = 10
        
    logical     :: lRefine
 
    integer                                 :: nd_local
    integer, dimension(:,:), allocatable    :: dp_local  ! dp array local to worker processors
    real(8), dimension(:), allocatable      :: d_local, dm
    real(8), dimension(:,:), allocatable    :: wj
       
    complex(8), dimension(:,:,:), allocatable  ::  readBuffer   
    real(8),    dimension(:,:,:), allocatable  ::  readBufferDC   
    
    type :: sol   ! structure for storing Kx domain fields:
        integer                                 :: nComponents     
        integer, dimension(:,:), allocatable    :: Components      
        complex(8), dimension(:,:), allocatable :: Fields ! nComponents x nkx
        complex(8), dimension(:,:), allocatable :: dFieldsdRho ! nComponents x nFree for AFTER IFT
    end type
 
    
    type(sol), dimension(:), allocatable   :: solutions ! structure for each frequency
              
              
    type :: solT
        complex(8)                             :: field  
        complex(8), dimension(:), allocatable  :: dFieldsdRho   
    end type 
         
    type :: solG     
        type(solT), dimension(:,:,:), allocatable :: comps ! isite,itrans,icomp
    end type
    
    type(solG), dimension(:), allocatable   :: solGrid ! structure for each frequency
          
     character(14)                           :: cScratchPrefix ! date and time stamp to make scratch 
    ! files have unique names per MARE2DEM run, in case stupid user runs multiple jobs at same time 
    ! on same machine. eg cScratchPrefix='201503111015'.   
 
    public :: worker_EM2D
    
    real(8)     :: clockRate   
    
    contains
!==================================================================================================================================! 
!================================================================================================================= worker_deallocate
!==================================================================================================================================!       

    subroutine worker_deallocate
    
    call deallocate_kx_io  ! deallocate i/o arrays for em2dkx
  
    if ( allocated(frequencies) )           deallocate( frequencies )  

    if ( allocated(lDataMask) )             deallocate( lDataMask )  

    if ( allocated(KxFqGroups) )            deallocate( KxFqGroups )     
 
    if ( allocated(solutions) )             deallocate( solutions )    
    if ( allocated(solGrid) )               deallocate( solGrid )  
    
    if ( allocated(readBuffer) )            deallocate( readBuffer ) 
    if ( allocated(readBufferDC) )          deallocate( readBufferDC )                    
                               
    if (allocated(xRxIn))                   deallocate( xRxIn, yRxIn, zRxIn, thetaRxIn, alphaRxIn, betaRxIn, lengthRxIn)
    if (allocated(xTxIn))                   deallocate( xTxIn, yTxIn, zTxIn, azimuthTxIn,dipTxIn, lengthTxIn, TxTypeIn) 
    
    if (allocated(trodes_RxDCIn))           deallocate(trodes_RxDCIn,trodes_TxDCIn,RxDCIn,TxDCIn)
    
       
    end subroutine worker_deallocate    
            
!==================================================================================================================================! 
!======================================================================================================================= worker_EM2D
!==================================================================================================================================!       

    subroutine worker_EM2D()

   
    real(8) :: randNum
 
    
    call system_clock(count_rate=clockRate)  
    
    call random_seed
    call random_number(randNum)
    write(cScratchPrefix,'(i13)') floor(1d6*randNum)
    cScratchPrefix = adjustl(cScratchPrefix)    
    
    maxRefinements = maxnadapt 
   
!
! Create the KxGroups structure:
!
    if (lPrintDebug) write(*,*) 'call setup_KxFqGroups ', myID
    call setup_KxFqGroups  
          
!
! Create the lDataMask structure:
!
    if (lPrintDebug) write(*,*) 'call setup_lDataMask ', myID
    call setup_lDataMask    
  
! 
! Compute the 2.5D wavenumber domain responses for all wavenumbers and frequencies:
!
    if (lPrintDebug) write(*,*) 'call computeAllKxFq ', myID
    call computeAllKxFq
    
!
! Inverse Fourier transform from the (kx,y,z) domain to  the (x,y,z) domain
!    
    if (lPrintDebug) write(*,*) 'call inverseFourierTransform ', myID
    if (sType == 'cs') then
        call inverseFourierTransform
    elseif (sType == 'dc') then
        call inverseFourierTransform_DC
    endif
    
    if (linversion) then
        if (lPrintDebug) write(*,*) 'call inverseFourierTransform_Derivs ', myID
        if (sType == 'cs') then
            call inverseFourierTransform_Derivs
        elseif (sType == 'dc') then
            call inverseFourierTransform_Derivs_DC
        endif             
    endif
        
!
! Move solutions into nTx/2 by nRx arrays for easy lookup:
!
    if (lPrintDebug) write(*,*) 'call moveSolution ', myID
    call moveSolution
!
! Transform to requested model responses:
!
    if (lPrintDebug) write(*,*) 'call convert_to_data_format ', myID
    call convert_to_data_format
    
!
! Deallocate:
!
    if (lPrintDebug) write(*,*) 'call worker_deallocate ', myID
 
    call worker_deallocate
 
    if (lPrintDebug) write(*,*) 'leaving worker_EM2D ', myID
 
    end subroutine worker_EM2D
!==================================================================================================================================! 
!=================================================================================================================== getStartingMesh
!==================================================================================================================================!     
    subroutine getStartingMesh
!
 
    use triangle_mesh
    use mare2dem_global
    use em2dkx_mod
 
    implicit none
 
    character (24) :: cend,tricommand 
    
    integer :: i,j,n(3)
    
!
! Create the starting mesh:
!
    call copy_trimesh(inputmodel,inputmesh)

    if (inputmesh%nele == 0 ) then 
        if (lprintDebug) then
            write(*,*)  ' '
            write(*,*) ' Generating mesh from input polygon model by calling Triangle:'
            write(*,*) ' '
        endif
        
        write(cend,fmt='(F4.0)') minQangle
        tricommand = 'q'//trim(adjustl(cend))//'QpanjA'//CHAR(0)
           
        call call_triangle(tricommand,inputmesh)  
    else
        write(*,*) 'Error: input mesh already has been triangulated, stopping'
        write(*,*) ' '
        stop
    endif
 
    if (lprintDebug) then
        write(*,'(a38,1x,i6,a2,i6)') ' Starting mesh: # nodes, # elements: ',inputmesh%nnod,', ',inputmesh%nele
        write(*,'(a38,1x,g11.4,1x,g11.4)') ' Min/max resistivities: ',minval(rhoParams),maxval(rhoParams)
        write(*,*) ' '  
    endif
    
    !
    ! Check to make sure regions don't have 0 sig (usually regions that you've forgotten to assign during mesh design)
    !
    if ( any( inputmesh%attr == 0) )then
        write(*,*) ''
        write(*,*) ' !!! Error reading mesh !!!'
        Write(*,*) ' Some elements have 0 attributes!'
        do i = 1,inputmesh%nele
           if (inputmesh%attr(i) == 0 ) then
                n = inputmesh%emap(1:3,i)
                do j = 1,3
                    write(*,*) 'element: ', i,inputmesh%y(n(j)),inputmesh%z(n(j))   
                enddo
                write(*,*) ' '
           endif
        enddo
        write(*,*) ' Stopping. '
        stop    
    endif   
 
    if (lprintDebug) then
        write(*,*) '  '
        write(*,*) ' Set up done, starting 2.5D EM computations...'
        write(*,*) '  '
    endif
    
    
    end subroutine getStartingMesh

!==================================================================================================================================! 
!====================================================================================================================== coarsenModel
!==================================================================================================================================!       
    subroutine coarsenModel 
    
    use triangle_mesh
    use kdtree2_module    
    use fem2d_utilities  
       
    use em2dkx_mod
    use mare2dem_global    
    use mare2dem_input_data_params    
    implicit none
 
    integer     :: iseg, inod, n1, n2, icnt,j,i,nn,nFreeNew,ireg
    real(8)     :: ymin,ymax,meshCoarsenPadding
    
    integer, dimension(:), allocatable :: iOld2New, iNodeCnt
    logical, dimension(:), allocatable :: lNodeKeep, lSegKeep
    
    character(80)   :: ctemp
    
    type(trimesh)  :: tempMesh  
    character (24) :: cend,tricommand 

    type(kdtree2), pointer  :: tree 
 
    real(8), dimension(:,:),allocatable :: yz
    real(8), dimension(:),allocatable   :: sigTemp
    integer, dimension(:), allocatable  :: node2tri, nCnt,iInputToTemp 
    integer  ::  inele(1)
                    
    !
    ! Get data window defined by Tx and Rx y position limits plus some padding:
    !
    if (sType == 'cs') then
        ymin = min(minval(yRxIn),minval(yTxIn)) 
        ymax = max(maxval(yRxIn),maxval(yTxIn))  
        meshCoarsenPadding = meshCoarsenPadding_CS
    elseif (sType == 'dc') then
        ymin = min(minval(trodes_RxDCIn(:,2)),minval(trodes_TxDCIn(:,2)))  ! not allocated yet! 
        ymax = max(maxval(trodes_RxDCIn(:,2)),maxval(trodes_TxDCIn(:,2)))  
        meshCoarsenPadding = meshCoarsenPadding_DC
    else
        write(*,*) ' !!!!!!!!!!!' 
        write(*,*) ' Error in subroutine coarsenModel. Mesh coarsening is not supported for this data type: ',trim(sType)
        write(*,*) ' stopping!'
        stop
    endif
    ymin = ymin - meshCoarsenPadding
    ymax = ymax + meshCoarsenPadding
    
    if (lprintDebug) then    
        write(*,*) ' '
        write(*,*) ' ----------------------------'
        write(*,*) '  Inversion Model Coarsening:'
        write(*,*) ' ----------------------------'    
        write(*,'(a,f12.0,a,f12.0)') ' Footprint for this subset: ', ymin,' to ', ymax
        write(ctemp,'(i6,1x,i6)') inputmodel%nnod ,inputmodel%numberofsegments        
        write(6,fmt='(a32,a)') 'Before: # Nodes, #Segs:   ',ctemp
    endif
    
    !
    ! Mark nodes and segments for deletion:
    !
    allocate(lNodekeep(inputmodel%nnod), lSegKeep(inputmodel%numberofsegments) )
    lNodeKeep = .false.
    lSegKeep  = .false.
    
    do iseg = 1,inputmodel%numberofsegments

        n1 = inputmodel%segmentlist(2*iseg-1)
        n2 = inputmodel%segmentlist(2*iseg  )
        

        ! First case: fixed segments, always keep nodes and segs
        if (abs(inputmodel%segmentmarkerlist(iseg)) < 2 ) then ! 1 is for segments bounding fixed regions or outer boundary
            
            lSegKeep(iseg) = .true.
            
        else ! Free segment:
 
            ! Keep segment if it is inside or crosses the region

            if (    (min(inputmodel%y(n1),inputmodel%y(n2)) <= ymax).and. & 
                    (max(inputmodel%y(n1),inputmodel%y(n2)) >= ymin) ) then
                    lSegKeep(iseg) = .true.
            endif
          
        endif
        
    enddo
    
    ! second pass to set nodes, and to count number of segments incident on each node:
    allocate(iNodeCnt(inputmodel%nnod) )
    iNodeCnt = 0
    do iseg = 1,inputmodel%numberofsegments
        if (lSegKeep(iseg)) then
            n1 = inputmodel%segmentlist(2*iseg-1)
            n2 = inputmodel%segmentlist(2*iseg  )
            lNodeKeep(n1) = .true.
            lNodeKeep(n2) = .true.
            iNodeCnt(n1) = iNodeCnt(n1)  + 1
            iNodeCnt(n2) = iNodeCnt(n2)  + 1
        endif   
    enddo    
    
    ! Third pass to look for dangles (segments that don't form a closed loop, defined by nodes with one incidence count):
 
    do iseg = 1,inputmodel%numberofsegments
        if (lSegKeep(iseg)) then
            n1 = inputmodel%segmentlist(2*iseg-1)
            n2 = inputmodel%segmentlist(2*iseg  )    
            if (iNodeCnt(n1) == 1) then
                lNodeKeep(n1) = .false.
                lSegKeep(iseg) = .false.
            endif       
            if (iNodeCnt(n2) == 1) then
                lNodeKeep(n2) = .false.
                lSegKeep(iseg) = .false.
            endif   
        endif 
    enddo
    
    ! debug print check:
    !   do inod = 1,inputmodel%nnod
    !       write(*,*) inod, lNodeKeep(inod), inputmodel%y(inod)
    !   enddo
    !write(*,*) 'Nodes: ',inputmodel%nnod, ' -> ',count(lNodeKeep), 'Segments: ',inputmodel%numberofsegments, ' -> ',count(lSegKeep)

    !
    ! Now delete those nodes and segments:
    !
 
    allocate(iOld2New(inputmodel%nnod))
    iOld2New = -1
    
    icnt = 0
    
    do inod = 1,inputmodel%nnod
        
        if (lNodeKeep(inod)) then
        
            icnt = icnt + 1
            inputmodel%y(icnt) = inputmodel%y(inod)
            inputmodel%z(icnt) = inputmodel%z(inod) 
            inputmodel%pointmarkerlist(icnt) = inputmodel%pointmarkerlist(inod)
            
            !write(*,*) icnt, inputmodel%y(icnt),inputmodel%z(icnt),inputmodel%pointmarkerlist(icnt) 
            
            iOld2New(inod) = icnt
            
            do j = 1,inputmodel%numberofpointattributes
                inputmodel%pointattributelist((icnt-1)*inputmodel%numberofpointattributes+j) = inputmodel%pointattributelist((inod-1)*inputmodel%numberofpointattributes+j)
            enddo           
            
        endif
        
    enddo 
    inputmodel%nnod = icnt
 
    
    icnt = 0
    do iseg = 1,inputmodel%numberofsegments
    
        if (lSegKeep(iseg)) then
        
            icnt = icnt + 1
            n1 = inputmodel%segmentlist(2*iseg-1)
            n2 = inputmodel%segmentlist(2*iseg  )
            
            inputmodel%segmentlist(2*icnt-1) = iOld2New(n1)
            inputmodel%segmentlist(2*icnt  ) = iOld2New(n2)
            inputmodel%segmentmarkerlist(icnt) = inputmodel%segmentmarkerlist(iseg)
            
            ! debug check:
            if ( (iOld2New(n1) < 1) .or. (iOld2New(n2) < 1 ) ) then     
                write(*,*)' iOld2New negative, bug in code!'
                stop
            endif   
 
        endif
                
    enddo
    inputmodel%numberofsegments = icnt
    inputmodel%segmentlist(2*icnt+1:  ) = 0
    inputmodel%segmentmarkerlist(icnt+1:) = 0

    deallocate(lNodekeep, lSegKeep, iOld2New, iNodeCnt)
    
 
!        
! New Robust method: Generate CDT and record regions used in it.
!
    write(cend,fmt='(F4.0)') minQangle
    tricommand = 'q'//trim(adjustl(cend))//'QpanjA'//CHAR(0)
 
    call call_triangle(tricommand,inputmodel,tempMesh) 
 
    
    allocate(yz(2,size(tempMesh%y)))  !KWK Sept 25 2015: create temporary array to avoid memory segfaults when mesh gets large.
    yz(1,:) = tempMesh%y
    yz(2,:) = tempMesh%z    
    tree => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
    deallocate(yz)
    nullify(tree%the_data)
 
    allocate (node2tri(tempMesh%nnod))
    
    call getNode2Tri( tempMesh%nele, tempMesh%emap,tempMesh%nnod, node2tri)
   
    allocate (sigTemp(tempMesh%nele),iInputToTemp(inputmodel%numberofregions), nCnt(inputmodel%numberofregions))
    sigTemp = 1d0        
    
    nCnt = 0
    iInputToTemp = 0
 
    ! For each old region, find element in tempMesh containing y,z for the region. Get new region 
    ! number for that element. Increment counter for that region. Any region visited more than 
    ! once is then a merger of two or more old regions and hence should not be considered
    ! a free parameter in the coarsened model
    do i = 1,inputmodel%numberofregions
        ! Find element in tempMesh containing this location of this region from inputModel:
            
        call findElement(tree, tempMesh%nnod, tempMesh%y, tempMesh%z, tempMesh%nele, tempMesh%emap, tempMesh%neighborlist, &
                    & node2tri,  sigTemp, 1, inputmodel%regionlist(4*i - 3), inputmodel%regionlist(4*i - 2),inele)          

        
        ! Get tempMesh region label for that element:
        ireg = tempMesh%attr(inele(1))
        
        nCnt(ireg) =  nCnt(ireg) + 1                
        
        iInputToTemp(i) = ireg 
    
    enddo
 
    deallocate(node2tri,sigTemp)        
    call kdtree2_destroy(tree)  
      
    
    !
    ! Now set freeparameter = 0 for all regions outside the data window:
    !
   
    nFreeNew = 0
    allocate(iFreeNewToOld(nFree))
    iFreeNewToOld = 0
    
    do i = 1,inputmodel%numberofregions
 
        ireg = inputmodel%regionlist(4*i - 1)
        if (nCnt(iInputToTemp(ireg))>1) then  ! region has been merged, no long free parameter:
        
            select case (trim(cAnisotropy))

            case ('isotropic')
                iFreeParam(ireg) = 0
            case ('isotropic_ip')
                iFreeParam(4*(ireg-1) + 1) = 0 
                iFreeParam(4*(ireg-1) + 2) = 0 
                iFreeParam(4*(ireg-1) + 3) = 0  
                iFreeParam(4*(ireg-1) + 4) = 0  
            case ('isotropic_complex')
                iFreeParam(2*(ireg-1) + 1) = 0 
                iFreeParam(2*(ireg-1) + 2) = 0                               
            case ('triaxial')
                iFreeParam(3*(ireg-1) + 1) = 0 
                iFreeParam(3*(ireg-1) + 2) = 0 
                iFreeParam(3*(ireg-1) + 3) = 0 
            case ('tix','tiy','tiz')                               
                iFreeParam(2*(ireg-1) + 1) = 0 
                iFreeParam(2*(ireg-1) + 2) = 0 
            case default
                write(*,*) 'Error decoding anisotropy in mare2dem_worker.f90'
                write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
                stop

            end select
            
           ! write(*,*) 'region outside footprint, i, ireg: ',i,ireg
            
         else ! keep parameter, update new counter:
            select case (trim(cAnisotropy))
            
            case ('isotropic')
                nn = 1
            case ('isotropic_ip')
                nn = 4 
            case ('isotropic_complex')
                nn = 2                               
            case ('triaxial')
                nn = 3
            case ('tix','tiy','tiz')    
                nn = 2       
            case default
                write(*,*) 'Error decoding anisotropy in mare2dem_worker.f90'
                write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
                stop

            end select     
            do j = 1,nn
                if (iFreeParam( nn*(ireg-1)+j)> 0 )  then
                    nFreeNew = nFreeNew + 1
                    iFreeNewToOld(nFreeNew) = iFreeParam(nn*(ireg-1)+j)
                    iFreeParam(nn*(ireg-1)+j) = nFreeNew              
                endif
            enddo
                   
        endif
    enddo

    if (lprintDebug) then 
        write(ctemp,'(i6,1x,i6)') inputmodel%nnod ,inputmodel%numberofsegments        
        write(6,fmt='(a32,a)') 'After: # Nodes, #Segs:   ',ctemp   
        write(6,*) ' nFree before & after:', nFree,nFreeNew 
    endif      
    nFree = nFreeNew ! reset nFree, so we don't have to change other subroutines much

    deallocate(iInputToTemp,nCnt)
  
    if ( allocated( tempMesh%attr ) )      call deallocate_trimesh(tempMesh,.false.)
 
!         
  
         
    !if (lSaveMeshFiles) then
!        write(cadapt,'(i6)') myID  
!        cadapt = adjustl(cadapt)       
! 
!        filename  = trim(outputFileRoot)//'.'//trim(cadapt)//'.coarsened'
!        !write(*,*) 'writing file: ',trim(filename)
! 
!        call write_triangle(inputmodel,filename)
   ! endif
         
    end subroutine coarsenModel      
!==================================================================================================================================! 
!====================================================================================================================== moveSolution
!==================================================================================================================================!  
    subroutine moveSolution
    
    integer :: iFq, nc, ic, iRx, iTx, iComp 
   
    if (lprintDebug)  write(*,*) ' in moveSolution'     
    allocate( solGrid(nfrequencies))
    
    do iFq = 1,nfrequencies
    
        nc = solutions(iFq)%nComponents

        if (sType == 'cs') nc = nc/2 ! since nTx array is 2x larger for -ikx sources   
        
        if (nc == 0) continue
        
        if (sType == 'dc') then
            allocate(solGrid(iFq)%comps(nTrodesRxDCIn,nTrodesTxDCIn,1))
        else        
            allocate(solGrid(iFq)%comps(nRxIn,nTxIn,6))
        endif
        
        do ic = 1,nc
        
            iRx   = solutions(iFq)%Components(ic,4) ! 4 is index to RxIn or DC global Rx electrode
            iTx   = solutions(iFq)%Components(ic,5) ! 4 is index to TxIn or DC global Tx electrode
            iComp = solutions(iFq)%Components(ic,3)
            
           ! write(*,'(3(i4,1x),g8.1)') iRx,iTx,iComp, abs(solutions(iFq)%Fields(ic,1))
   
            solGrid(iFq)%comps(iRx,iTx,iComp)%field = solutions(iFq)%Fields(ic,1)
                       
            if (linversion) then
                allocate(solGrid(iFq)%comps(iRx,iTx,iComp)%dFieldsdRho(nFree) ) 
                solGrid(iFq)%comps(iRx,iTx,iComp)%dFieldsdRho = solutions(iFq)%dFieldsdRho(ic,:) 
                
                !write(*,'(i2,100(e8.1,1x))') ic, abs(solutions(iFq)%dFieldsdRho(ic,:)) 
            endif
       
            ! change phase convention for MT:
            if (sType == 'mt') then
                 solGrid(iFq)%comps(iRx,iTx,iComp)%field = conjg(solGrid(iFq)%comps(iRx,iTx,iComp)%field)
                if (linversion) solGrid(iFq)%comps(iRx,iTx,iComp)%dFieldsdRho = conjg(solGrid(iFq)%comps(iRx,iTx,iComp)%dFieldsdRho)
            endif
        enddo     

    enddo
    
    end subroutine moveSolution
    
!==================================================================================================================================! 
!================================================================================================================= getRotationMatrix
!==================================================================================================================================!       
    subroutine getRotationMatrix(theta,alpha,beta,RotR)
!
! Outputs the rotation matrix is for left multiplication: F_rotated = RotR * F_xyz
!    
    real(8), intent(in)                     :: theta, alpha, beta   ! site rotation angles
    real(8), dimension(3,3), intent(out)    :: RotR                 ! 3x3 rotation matrices
 
    real(8)                                 :: cct, sst, cca, ssa, ccb, ssb ! cos and sin values for rotation angles
    real(8) , dimension(3,3)                :: Rotx,Roty,Rotz          ! 3x3 rotation matrices
   
    
    cct = cos(deg2rad* theta)
    sst = sin(deg2rad* theta)
    cca = cos(deg2rad* alpha)
    ssa = sin(deg2rad* alpha)
    ccb = cos(deg2rad* beta)
    ssb = sin(deg2rad* beta)                

    Rotx = 0d0
    Rotx(1,1) = 1.d0
    Rotx(2,2) =  ccb
    Rotx(2,3) =  ssb
    Rotx(3,2) = -ssb
    Rotx(3,3) =  ccb

    Roty = 0d0
    Roty(1,1) =  cca  
    Roty(1,3) =  ssa
    Roty(2,2) =  1.d0 
    Roty(3,1) = -ssa
    Roty(3,3) =  cca

    Rotz = 0d0
    Rotz(1,1) =  cct
    Rotz(1,2) =  sst
    Rotz(2,1) = -sst
    Rotz(2,2) =  cct
    Rotz(3,3) =  1.d0    

    RotR = matmul(Roty,Rotz)  
    RotR = matmul(Rotx,RotR)       
    
    
    end subroutine getRotationMatrix 
    
!==================================================================================================================================! 
!================================================================================================================== setup_KxFqGroups
!==================================================================================================================================!       
    subroutine setup_KxFqGroups

    integer :: i, icnt, ikxG, nk0, nkg, ikx, iFq
    
!
! Count how many KxFq groups we need:
!
    if (sType == 'mt') then
        nKxFqGroups = 1
    else
        nKxFqGroups = ceiling( dble(nwave)      / dble(nKxPerGroup)  )   
    endif
!
! Allocate the structure:
!
    allocate ( KxFqGroups(nKxFqGroups) )
    
!
! Fill in the values:
!     
    if (sType == 'mt') then !MT
    
        KxFqGroups(1)%nKxFq = nfrequencies - 1
        
        allocate( KxFqGroups(1)%iFq(KxFqGroups(1)%nKxFq),KxFqGroups(1)%iKx(KxFqGroups(1)%nKxFq) )
        
        icnt = 0
        do i = 1,nfrequencies
            if ( mod(i+floor(dble(nfrequencies)/dble(2)),nfrequencies) ==0 ) then
                KxFqGroups(1)%iFqRefine = i
                KxFqGroups(1)%iKxRefine = 1
            else
                icnt = icnt + 1
                KxFqGroups(1)%iFq(icnt) = i
                KxFqGroups(1)%iKx(icnt) = 1
                if (lPrintDebug) write(*,'(a,2(i4,1x))') 'MT KxFqG refine icnt,iFq  : ',icnt,i
            endif
            
        enddo    
    
    else ! CSEM and DC:
    
        do ikxG = 1,nKxFqGroups   
        
            nk0 =  (ikxG-1)*nKxPerGroup
            nkg = min(nKxPerGroup, nwave - nk0)   ! number of wavenumbers in this group 
            
            KxFqGroups(ikxG)%nKxFq = nkg*nfrequencies - 1 ! -1 for the refinement pair
            
            allocate( KxFqGroups(ikxG)%iFq(KxFqGroups(ikxG)%nKxFq),KxFqGroups(ikxG)%iKx(KxFqGroups(ikxG)%nKxFq) )
            
            icnt = 0
            do iFq = 1,nfrequencies
            
                do i = 1,nkg
                    
                    ikx = nk0 + i
                    
                    if ( ( mod(iFq + floor(dble(nfrequencies)/dble(2)),nfrequencies) ==0 ) .and.   &
                    &    ( mod(i + floor(dble(nkg)/dble(2)),nkg) ==0 ) ) then
                        KxFqGroups(ikxG)%iFqRefine = iFq
                        KxFqGroups(ikxG)%iKxRefine = ikx
                        if (lPrintDebug) write(*,'(a,3(i4,1x))') 'CSEM/DC KxFqG refine ikxG,iFq,ikx  : ',ikxG,iFq,ikx  
                    else
                        icnt = icnt + 1
                        KxFqGroups(ikxG)%iFq(icnt) = iFq
                        KxFqGroups(ikxG)%iKx(icnt) = ikx
                        if (lPrintDebug) write(*,'(a,4(i4,1x))') 'CSEM/DC KxFqG compute ikxG,icnt,iFq,ikx  : ',ikxG,icnt,iFq,ikx  
                    endif         
                    
                   
                enddo
            enddo
        enddo
  
    endif
 
    end subroutine setup_KxFqGroups 

!==================================================================================================================================! 
!=================================================================================================================== computeAllKxFq
!==================================================================================================================================!       
    subroutine computeAllKxFq
 
    integer     :: ikx,ifq,ikxG,iKxFq, isubset, ifqLast
    
    logical :: lLocalRefineIn
 
    if (lPrintDebug) write(*,*) 'entering computeAllKxFq...'
    
    lLocalRefineIn = lLocalRefine
    isubset = 0
    ifqLast = 0
    do ikxG = 1,nKxFqGroups
        
        ikx = KxFqGroups(ikxG)%iKxRefine
        ifq = KxFqGroups(ikxG)%iFqRefine
        
        w   = 2*pi*frequencies(ifq)
        
        if (sType == 'mt') then
            kx = 0
        else
            kx = wavenum(ikx)        
        endif
        
        lRefine    = .true.
        lLocalRefine = lLocalRefineIn
        meshnumber = 1
        
        if (ifq /= ifqLast) then
            if (lPrintDebug) write(*,'(a,1x,3(i3,1x))') 'refine: ifq /= ifqLast: call make_iDataMask_Components: ',myid, ifq, ikx
            if (sType == 'dc') then
                call make_iDataMask_Components_DC(ifq)
            else
                call make_iDataMask_Components(ifq)
            endif
            
            if (lPrintDebug) write(*,*) 'call allocateSolutions(ifq): ',myid
            call allocateSolutions(ifq)
            ifqLast = ifq    
        endif
                     
        isubset = isubset + 1
        
        if (lPrintDebug) write(*,*) 'call refinementKernel(isubset): ',myid
        call refinementKernel(isubset,ikx,ifq)
        
        do iKxFq = 1,KxFqGroups(ikxG)%nKxFq
    
            ikx = KxFqGroups(ikxG)%iKx(iKxFq)
            ifq = KxFqGroups(ikxG)%iFq(iKxFq)
            
            w   = 2*pi*frequencies(ifq)
             
            if (sType == 'mt') then     
                kx = 0
            else
                kx = wavenum(ikx)        
            endif
            
            lRefine = .false.
            lLocalRefine = .false.
            
            if (lPrintDebug) write(*,*) 'call make_iDataMask_Components: ',myid, ifq, ikx
            if (ifq /= ifqLast) then
                if (lPrintDebug) write(*,'(a,1x,3(i3,1x))') 'nonrefine: ifq /= ifqLast: call make_iDataMask_Components: ',myid, ifq, ikx
                if (sType == 'dc') then
                    call make_iDataMask_Components_DC(ifq)
                else
                    call make_iDataMask_Components(ifq)
                endif
                call allocateSolutions(ifq)
                ifqLast = iFq
            endif
            
            isubset = isubset + 1
            
            if (lPrintDebug) write(*,*) 'call refinementKernel(isubset): ',myid
            call refinementKernel(isubset,ikx,ifq)

        enddo
         
    enddo
    
    lLocalRefine = lLocalRefineIn
 
    end subroutine computeAllKxFq
!==================================================================================================================================! 
!=================================================================================================================== setup_lDataMask
!==================================================================================================================================!       
    subroutine setup_lDataMask 

    implicit none
        
    include 'em_parameters.inc'
        
    integer :: i,isite,itrans,ifreq
    
   ! write(*,'(a,8(i6,1x))') 'setup_lDataMask: nfrequencies,nTxIn,nRxIn' ,nfrequencies,nTxIn,nRxIn
    !write(*,*) nd_local,nfrequencies, nTxDCIn,nRxDCIn
    
    if (sType =='dc') then
        allocate (lDataMask(nfrequencies,nTxDCIn,nRxDCIn,1), solutions(nfrequencies) )
    else
        allocate (lDataMask(nfrequencies,nTxIn,nRxIn,6), solutions(nfrequencies) )
    endif    
    lDataMask = .false.
 
    !
    ! Loop through dp_local array and fill in iDataMask as appropriate
    !
    do i = 1, nd_local
        
        isite  = dp_local(i,4)
        itrans = dp_local(i,3)
        ifreq  = dp_local(i,2)
        
        if (sType =='dc') ifreq = 1 ! override any input values...
        !write(*,*) isite,itrans,ifreq
        
        if (itrans <= 0) itrans = isite  ! if itrans not specified, make it the current receiver 
 
 
        select case (  dp_local(i,1) ) ! what data type is this?
                                     
            ! CSEM types:                                  
            case (indRealEx,indImagEx,indAmpEx,indPhsEx,indLog10AmpEx) 
                lDataMask(ifreq,itrans      ,isite,1) = .true. 
                
            case (indRealEy,indImagEy,indAmpEy,indPhsEy,indLog10AmpEy)  
                lDataMask(ifreq,itrans      ,isite,2) = .true.

            case (indRealEz,indImagEz,indAmpEz,indPhsEz,indLog10AmpEz)  
                lDataMask(ifreq,itrans      ,isite,3) = .true.
                
            case (indRealBx,indImagBx,indAmpBx,indPhsBx,indLog10AmpBx)  
                lDataMask(ifreq,itrans      ,isite,4) = .true.

            case (indRealBy,indImagBy,indAmpBy,indPhsBy,indLog10AmpBy)  
                lDataMask(ifreq,itrans      ,isite,5) = .true.

            case (indRealBz,indImagBz,indAmpBz,indPhsBz,indLog10AmpBz)                      
                lDataMask(ifreq,itrans      ,isite,6) = .true.

            case (iPEmax,iPEmin)
                lDataMask(ifreq,itrans      ,isite,1) = .true.
                lDataMask(ifreq,itrans      ,isite,2) = .true.

            case (iPBmax,iPBmin)
                lDataMask(ifreq,itrans      ,isite,4) = .true.
                lDataMask(ifreq,itrans      ,isite,5) = .true.
                
            ! MT types:   
            case (indRhoZXY,indPhsZXY,indRealZXY,indImagZXY,indlog10RhoZXY, indRealMZY, indImagMZY,indAmpMZY,indPhsMZY )    
                lDataMask(ifreq,1,isite, 1:3) = .true.
                lDataMask(ifreq,1,itrans,4:6) = .true.
 
            case (indRhoZYX,indPhsZYX,indRealZYX,indImagZYX,indlog10RhoZYX)    
                lDataMask(ifreq,1,isite, 1:3) = .true.
                lDataMask(ifreq,1,itrans,4:6) = .true.
                
!             case (indRealExMT,indImagExMT)    
!                 lDataMask(ifreq,1,isite,1) = .true.  
!                   
!             case (indRealEyMT,indImagEyMT)  
!                 lDataMask(ifreq,1,isite,2) = .true.  
! 
!             case (indRealEzMT,indImagEzMT)  
!                 lDataMask(ifreq,1,isite,3) = .true. 
! 
!             case (indRealHxMT,indImagHxMT)  
!                 lDataMask(ifreq,1,isite,4) = .true.  
!                   
!             case (indRealHyMT,indImagHyMT)  
!                 lDataMask(ifreq,1,isite,5) = .true.  
! 
!             case (indRealHzMT,indImagHzMT)  
!                 lDataMask(ifreq,1,isite,6) = .true. 

           ! KWK note: need to set all components for MT since convergence test in em2dkx assumes all MT components there: 
            case (indRealExMT,indImagExMT,indRealEyMT,indImagEyMT,indRealEzMT,indImagEzMT)    
                lDataMask(ifreq,1,isite, 1:3) = .true.
                lDataMask(ifreq,1,itrans,4:6) = .true.
                
            case (indRealHxMT,indImagHxMT,indRealHyMT,indImagHyMT,indRealHzMT,indImagHzMT)    
                lDataMask(ifreq,1,isite, 1:3) = .true.
                lDataMask(ifreq,1,itrans,4:6) = .true.                  
                                                                 
            ! DC resistivity
            case(indAppRes_DC,indlog10AppRes_DC)              
                lDataMask(ifreq,itrans,isite, 1 ) =  .true.
                lDataMask(ifreq,itrans,isite, 1 ) =  .true.
                
        end select  ! case dp(i,1)
        
    enddo
    
    !
    ! Also set TE and TM mode data logical flags if MT computations (used in em2dkx to skip mode if no data present)
    !
    do i = 1, nd_local
 
        select case (  dp_local(i,1) ) ! what data type is this?
        
        ! TE types:  
            case (indRhoZXY,indPhsZXY,indRealZXY,indImagZXY,indlog10RhoZXY, indRealMZY, indImagMZY,indAmpMZY,indPhsMZY )    
                l_has_TE_mode = .true.
               
            case (indRealExMT,indImagExMT,indRealHyMT,indImagHyMT,indRealHzMT,indImagHzMT )    
                l_has_TE_mode = .true.
                       
        ! TM types
            case (indRhoZYX,indPhsZYX,indRealZYX,indImagZYX,indlog10RhoZYX)    
                l_has_TM_mode = .true.
                
            case (indRealHxMT,indImagHxMT,indRealEyMT,indImagEyMT,indRealEzMT,indImagEzMT)    
                l_has_TM_mode = .true.
                                    
        end select  ! case dp(i,1)
        
    enddo
        
    end subroutine setup_lDataMask
!==================================================================================================================================! 
!========================================================================================================= make_iDataMask_Components
!==================================================================================================================================!       
    subroutine  make_iDataMask_Components(iFq)
    
    integer, intent(in) :: iFq

    
    integer :: isite,itrans,icomp,iPass, nc,nc0, ict
    
    real(8), dimension(3,6)     :: V =[1,0,0, 0,1,0, 0,0,1, 1,0,0, 0,1,0, 0,0,1]
    real(8)                     :: vec(3), theta, alpha, beta           ! site rotation angles
    real(8) , dimension(3,3)    :: RotR 
 
    integer, dimension(:), allocatable :: iG2LTx
  
    if (allocated(Components))  deallocate(Components)
 
!
! Array of transmitters for requested data:
!
    nTx = 0
    do itrans = 1,nTxIn
        if (any(lDataMask(iFq,itrans,:,:))) nTx  = nTx + 1
    enddo    
   ! write(*,*) 'nTx: ',nTx
    
    if (sType == 'mt') then 
        nTx = 1
        if (.not.(allocated(lNeg_ikx))) allocate(lNeg_ikx(nTx))
        lNeg_ikx = .false.    
        allocate(iG2LTx(nTx))    
        iG2LTx = 1
    else
        nTx = nTx*2 ! for -ikx
          
         ! write(*,*) 'nTx*2: ',nTx     
          
        if (allocated(xTx)) deallocate( xTx, yTx, zTx, azimuthTx,dipTx, lengthTx, TxType )
       
        allocate(xTx(nTx), yTx(nTx), zTx(nTx), azimuthTx(nTx), dipTx(nTx), lengthTx(nTx), TxType(nTx) )
         
        if (.not.(allocated(lNeg_ikx))) allocate(lNeg_ikx(nTx))  
        allocate(iG2LTx(nTxIn))
         
        iG2LTx = 0
        
        ict = 0
        
        do itrans = 1,nTxIn
            if (any(lDataMask(iFq,itrans,:,:))) then
                
                ict = ict + 1
                
                !write(*,*) 'itrans, ict: ',itrans,ict
                
                !+ikx:
                xTx(ict)            = xTxIn(itrans)
                yTx(ict)            = yTxIn(itrans)
                zTx(ict)            = zTxIn(itrans)
                azimuthTx(ict)      = azimuthTxIn(itrans)
                dipTx(ict)          = dipTxIn(itrans)
                lengthTx(ict)       = lengthTxIn(itrans)
                TxType(ict)         = TxTypeIn(itrans)
                lNeg_ikx(ict)       = .false.
                 
                iG2LTx(itrans)      = ict
               
                
                !-ikx:
                xTx(ict+nTx/2)            = xTxIn(itrans)
                yTx(ict+nTx/2)            = yTxIn(itrans)
                zTx(ict+nTx/2)            = zTxIn(itrans)
                azimuthTx(ict+nTx/2)      = azimuthTxIn(itrans)
                dipTx(ict+nTx/2)          = dipTxIn(itrans)
                lengthTx(ict+nTx/2)       = lengthTxIn(itrans)
                TxType(ict+nTx/2)         = TxTypeIn(itrans)     
                lNeg_ikx(ict+nTx/2)       = .true. 

            endif
        enddo
                
    endif        
            
!  
! now get receivers for each component:
!        
    nc0 = 0    
    do iPass =1,3
       
        if (iPass == 2) then
            
            nComponents = nc*2
            if (sType == 'mt') nComponents = nc
            
            allocate(Components(nComponents,5))                 
            
            if (allocated(xRx)) deallocate( xRx, yRx, zRx, azimuthRx,dipRx, lengthRx, RxType )
            if (allocated(iDataMask)) deallocate( iDataMask)
            
            allocate(xRx(nRx), yRx(nRx),zRx(nRx), azimuthRx(nRx),dipRx(nRx), lengthRx(nRx), RxType(nRx) )

            allocate(iDataMask(nRx,nTx))          
            iDataMask = 0

        
        endif
        
        nc = 0   
        nRx = 0
        do isite = 1,nRxIn
        
            do icomp = 1,6
            
                if (any(lDataMask(iFq,:,isite,icomp)))  then 
                
                    nRx = nRx + 1   ! add the receiver
            
                    if (iPass == 2) then

                   
                        xRx(nRx) = xRxIn(isite)
                        yRx(nRx) = yRxIn(isite)
                        zRx(nRx) = zRxIn(isite)
   
                        lengthRx(nRx) =  lengthRxIn(isite)
                        
                        if (icomp < 4) then
                            RxType(nRx) = 'edipole'
                        else
                            RxType(nRx) = 'bdipole'
                            lengthRx(nRx)  = 0
                        endif
                        
                        ! get sensor vector angles:
                        theta = ThetaRxIn(isite)
                        alpha = AlphaRxIn(isite)
                        beta  = BetaRxIn(isite)        
                        call getRotationMatrix(theta,alpha,beta, RotR)           
                        
                        !write(*,'(3(f6.1,1x))') -theta,-alpha,-beta
                        ! rotate unit vector along component:
                        vec = matmul(transpose(RotR),V(1:3,icomp))
                        azimuthRx(nRx) = rad2deg*atan2(vec(2),vec(1))
                        dipRx(nRx)     = rad2deg*asin(vec(3))
                        
                       ! write(*,'(i4,1x,i4,1x,2(f6.1,1x))') isite,icomp, azimuthRx(nRx),dipRx(nRx)                   
                    
                    endif
                    
                    ! Now loop through transmitters and add TF components:
                    if (iPass == 1) then
 
                        nc0 = 0
    
                    elseif (iPass==2) then               
                        nc0 = 0
                        
                    else
 
                        nc0 = nComponents/2             
                    endif
                    
                                        
                    do itrans = 1,nTxIn
                    
                        if (lDataMask(iFq,itrans,isite,icomp)) then 
                        
                            nc = nc + 1                                         
                            
                            !write(*,*) 'isite,itrans, nc:',isite,itrans, nc
                            if (iPass >1) then
                                
                                iDataMask(nRx,iG2LTx(itrans)+nTx/2*(iPass-2) ) =  2
                                
                                Components(nc+nc0,1) = nRx 
                                Components(nc+nc0,2) = iG2LTx(itrans)+nTx/2*(iPass-2)   
                                Components(nc+nc0,3) = icomp
                                Components(nc+nc0,4) = isite 
                                Components(nc+nc0,5) = itrans
                                
                               ! write(*,'(a,7(i3,1x))') 'nc,nc0,nRx,iG2LTx(itrans),icomp,isite,itrans: ',nc,nc0,nRx,iG2LTx(itrans)+nTx/2*(iPass-2),icomp,isite,itrans
  
                            endif
                            
                        endif
                    enddo ! itrans = 1,nTx
                    
                endif               
            
            enddo ! icomp = 1,6
           
        enddo ! isite = 1,nRxIn
        
        if ((sType=='mt').and.(iPass == 2)) exit
    enddo  ! iPass 
  
    if (allocated(iG2LTx)) deallocate(iG2LTx)     
    
    end subroutine make_iDataMask_Components    
!==================================================================================================================================! 
!===================================================================================================== make_iDataMask_Components_DC
!==================================================================================================================================!       
    subroutine  make_iDataMask_Components_DC(iFq)
    
    use kx_io
    
    implicit none
    
!
! Creates module kx_io arrays for dc2dkx call:
!  - iDataMask [Rx electrodes x Tx electrodes]
!  - Components
!
!   Components(nc+nc0,1) = local Rx electrode 
!   Components(nc+nc0,2) = local Tx electrode
!   Components(nc+nc0,3) = icomp - 1?
!   Components(nc+nc0,4) = isite - index to global Rx electrode
!   Components(nc+nc0,5) = itrans - index to global Tx electrode
!

  !
  ! Allocate iDataMask, then traverse dp_local to get iRx and iTx, get associate electrodes and
  ! set to 0 = no data, 1 = no data but include for refinement, 2=data present
  !
    
    integer, intent(in) :: iFq ! not used, but ...
    
    integer :: i,j
    integer, dimension(:), allocatable :: iG2LTx,iG2LRx 
    
 
   if (allocated(Components))  deallocate(Components)
    
!
! Create kx_io arrays electrode_tx and electrode_rx 
!       
    if (lPrintDebug) write(*,'(a,1x,3(i3,1x))') 'make_iDataMask_Components_DC:myid, nTrodesTxDCIn,nTrodesRxDCIn: ',myid,nTrodesTxDCIn,nTrodesRxDCIn
    
    allocate(iG2LTx(nTrodesTxDCIn),iG2LRx(nTrodesRxDCIn))
    iG2LTx = 0
    iG2LRx = 0
     
    ! count tx and rx electrodes actually used:
    ntrodes_tx = 0
    do i = 1,nTxDCIn
        if (any(lDataMask(iFq,i,:,1))) then ! this Tx is used.
            do j = 1,2
                if (iG2LTx(TxDCIn(i,j)) == 0) then
                    ntrodes_tx = ntrodes_tx + 1
                    iG2LTx(TxDCIn(i,j)) = ntrodes_tx
                endif
            enddo
        endif
    enddo
    ntrodes_rx = 0
    do i = 1,nRxDCIn
        if (any(lDataMask(iFq,:,i,1))) then ! this Rx is used.
            do j = 1,2
                if (iG2LRx(RxDCIn(i,j)) == 0) then
                    ntrodes_rx = ntrodes_rx + 1
                    iG2LRx(RxDCIn(i,j)) = ntrodes_rx
                endif
            enddo
        endif
    enddo
    
    if (lPrintDebug) write(*,'(a,1x,3(i3,1x))') 'make_iDataMask_Components_DC:myid, ntrodes_tx,ntrodes_rx: ',myid, ntrodes_tx,ntrodes_rx
    
    allocate(electrode_tx(ntrodes_tx,3),electrode_rx(ntrodes_rx,3),iDataMask(ntrodes_rx,ntrodes_tx))
    do i = 1,nTrodesTxDCIn
        if ( iG2LTx(i)>0 ) electrode_tx(iG2LTx(i),1:3) = trodes_TxDCIn(i,1:3)
    enddo
    do i = 1,nTrodesRxDCIn
        if ( iG2LRx(i)>0 ) electrode_rx(iG2LRx(i),1:3) = trodes_RxDCIn(i,1:3)
    enddo
    iDataMask = 0
    do i = 1,nTrodesTxDCIn
        do j = 1,nTrodesRxDCIn
            if ( (iG2LTx(i)>0).and.(iG2LRx(j)>0)) then ! this tx and rx used:
                iDataMask(iG2LRx(j),iG2LTx(i)) = 2
            endif
        enddo
    enddo    
        
!
! Create components array used by dc2dkx:
!        
    nComponents = count(iDataMask>0)
    allocate(Components(nComponents,5))      
    nComponents = 0
    do i = 1,nTrodesTxDCIn
        do j = 1,nTrodesRxDCIn
            if ( (iG2LTx(i)>0).and.(iG2LRx(j)>0)) then
                nComponents = nComponents + 1
                Components(nComponents,1) = iG2LRx(j)
                Components(nComponents,2) = iG2LTx(i)
                Components(nComponents,3) = 1
                Components(nComponents,4) = j 
                Components(nComponents,5) = i            
            endif
        enddo
    enddo
    
!   Components(nc+nc0,1) = local Rx electrode 
!   Components(nc+nc0,2) = local Tx electrode
!   Components(nc+nc0,3) = icomp - 1?
!   Components(nc+nc0,4) = isite - index to Rx electrode pair. or global Rx electrode?
!   Components(nc+nc0,5) = itrans - index to Tx electrode pair. or global Tx electrode?
!
 
    deallocate(iG2LTx,iG2LRx)
    
    
    end subroutine make_iDataMask_Components_DC    
!==================================================================================================================================! 
!================================================================================================================= allocateSolutions
!==================================================================================================================================!       
    subroutine allocateSolutions(ifq)
    
    integer, intent(in) :: ifq
    
    integer :: ierr
 
    ! out with the old and in with the new:
    if (allocated(Fields))           deallocate(Fields)
    if (allocated(dFieldsdRho))      deallocate(dFieldsdRho)
    
    if (lPrintDebug) write(*,'(a,3(i,1x))') 'in allocateSolutions(ifq): ifq, nComponents,nFree',ifq,nComponents,nFree
    
    !
    ! Allocate storage for fields and derivatives:
    !
    allocate(Fields(nComponents))
    if (linversion) then
        allocate(dFieldsdRho(nComponents,nFree))
    endif
        
    ! allocate solution arrays for this kxfq group:
    if (.not.allocated(solutions(ifq)%Components)) allocate(   solutions(ifq)%Components(nComponents,5), solutions(ifq)%fields(nComponents,nwave))
    solutions(ifq)%nComponents = nComponents
    solutions(ifq)%Components  = Components

    if (linversion) then
        if (.not.allocated(solutions(ifq)%dFieldsdRho)) allocate(solutions(ifq)%dFieldsdRho(nComponents,nFree), stat=ierr)
        if (ierr .ne. 0) then
            write(*,*)  'Out of memory in computeAllKxFq.  Too many free parameters: ', nComponents,nFree 
            stop 
        endif               
    endif
              
    if (lPrintDebug) write(*,*) 'leaving allocateSolutions'
              
    end subroutine allocateSolutions
!==================================================================================================================================! 
!================================================================================================================== refinementKernel
!==================================================================================================================================!       
    subroutine refinementKernel(isubset,ikx,ifq)

    integer, intent(in) :: isubset, ikx,ifq
    character(256)      ::  cGroup, cSubset     
    real(8)             ::  tstart, tend
 
!
! Compute the 2.5 response for this kx, and possibly do some mesh refinement:
!
    write(cGroup,'(i6)') iRefinementGrp
    write(cSubset,'(i6)') isubset      
    fileroot = trim(outputFileroot)//'.'//trim(adjustl(cGroup))//'.'//trim(adjustl(cSubset))   
    
!
! If first mesh, copy from input:
!
    !write(*,*) '***refinementKernel meshnumber:',meshnumber
    if (meshnumber == 1) then
        
        if ( allocated( mesh%attr ) ) call deallocate_trimesh(mesh,.false.)

        if (lPrintDebug) write(*,*) 'call copy_trimesh: ',myid
        call copy_trimesh(inputmesh,mesh)
    
    endif

!
! Compute EM response for current wavenumber using a priori local and goal-oriented adaptive refinement:
!    
    lCompDerivs  = linversion
    
    if (sType == 'dc') then
        call dc2dkx(iRefinementGrp,isubset,lrefine)  
    else
        call em2dkx(iRefinementGrp,isubset,lrefine)
    endif
    
!
! Move kx domain fields into holding arrays for temporary storage:
!
    solutions(ifq)%fields(:,ikx) = fields

!
! Save kx domain sensitivity derivatives to scratch space
!
    if (linversion) then
    
    !
    ! If MT, just copy the derivatives across dsig arrays:
    !
        if (sType=='mt') then
        
            solutions(iFq)%dFieldsdRho = dFieldsdRho ! move solution into structure
            
        else   
        
        !
        ! For CSEM and DC data, save to scratch space:
        !
            call cpu_time(tstart)

            call writeScratchFile(iFq,iKx)  

            call cpu_time(tend)

            if (lprintMPItimers)  write(*, '((a5,2x,i6,2x,a,f12.5))')  'Proc:',myID,  ' Worker writescratch:',  (tEnd - tStart)

        endif
             
     endif       
     
     end subroutine refinementKernel

!==================================================================================================================================! 
!================================================================================================================== writeScratchFile
!==================================================================================================================================! 
    subroutine writeScratchFile(iFq,iKx)
!
! Writes scratch files by peeling off nRxTxPerIFT Rx,Tx pairs from the RxTxGroup iGrp  
!    
    use em2dkx_mod
    
    implicit none 
    
    integer, intent(in)     :: iKx, iFq 

    integer                 :: iScratch, ios, posStart, iComp, i, icnt, nScratchGroups
    character(256)          :: cFilename     
 
    if (sType == 'cs') then
        nScratchGroups = ceiling(dble(nComponents/2)/dble(nComponentsPerScratch))  ! note the nComponents/2 since 2nd half is -ikx srcs
    elseif (sType == 'dc') then
         nScratchGroups = ceiling(dble(nComponents)/dble(nComponentsPerScratch)) 
    endif
    
    iComp = 0   
    do iScratch = 1,nScratchGroups
        
        !
        ! Get a name for the file:    
        !
        call getScratchFilename(iRefinementGrp,iScratch,iFq,iKx,cFilename)
        !write(*,*) ' ikx, nComponents,writeScratchFile: ',ikx,nComponents,trim(cFilename)
        
        !  
        ! Open the file:            
        !
        open(unit=21,file=trim(cFilename),form="unformatted",status="replace",access="stream",iostat=ios)
        if (ios /= 0) then
            write(*,*) 'w: error, opening scratch file: ',trim(cFilename),  ' ios: ', ios
            stop
        endif
        inquire(21,pos=posStart)
        write(21,iostat=ios) 0      ! dummy value for now     
        if (ios /= 0) then
            write(*,*) 'w: error, writing dummy value to scratch file: ',trim(cFilename), 'ios: ',ios
            stop
        endif
        
        !
        ! Write out the data:
        !
        icnt = 0
        do i = 1,nComponentsPerScratch
           
            iComp =  iComp  + 1
            
            !write(*,*) i,iComp,icnt
            if (sType == 'cs') then
                if (iComp > nComponents/2) exit
            else
                if (iComp > nComponents) exit
            endif 
            
            !
            ! Write data from this Rx,Tx to the file:
            !
            if (sType == 'dc') then ! save only real part since imag = 0
                write(21,iostat=ios)  dble(dFieldsdRho(iComp,1:nFree))   ! + ikx 
            else
                write(21,iostat=ios)  dFieldsdRho(iComp,1:nFree)   ! + ikx
            endif
            if (ios /= 0) then
                write(*,*) 'error, writing scratch file: ',trim(cFilename), ' ios: ',ios
                stop
            endif            
             
            if (sType == 'cs') then ! write negative wavenumbers too
                write(21,iostat=ios)  dFieldsdRho(iComp+nComponents/2,1:nFree)  ! -ikx
                if (ios /= 0) then
                    write(*,*) 'error, writing scratch file: ',trim(cFilename), ' ios: ',ios
                    stop
                endif
            endif
            icnt = icnt + 1
        enddo
        
        ! Rewind and write status indicator:
        write(21,iostat=ios, pos=posStart) icnt  ! insert this at the end of i/o so it can be used as a file i/o status indicator
        if (ios /= 0) then
            write(*,*) 'error, rewinding to start of scratch file: ',trim(cFilename), ' ios: ', ios
            stop
        endif   
 
        ! Close the file:
        close(21,iostat=ios)   
        if (ios /= 0) then
            write(*,*) 'error, closing scratch file: ',trim(cFilename),  ' ios: ', ios
            stop
        endif
     
    enddo
  
    end subroutine writeScratchFile
    
!==================================================================================================================================! 
!================================================================================================================ getScratchFilename
!==================================================================================================================================! 
    subroutine getScratchFilename(iGrp,iFileNum,iFq,iKx,cFilename)  
    
    implicit none 
    
    integer, intent(in)         :: iGrp,iKx,iFq,iFileNum
    character(256), intent(out) :: cFilename
    
    character(256) ::  cGrp, cFileNum, cKx, cFq, cID
 
     write(cID,'(i12)') myID  
    cID = adjustl(cID)   
    
    write(cGrp,'(i12)') iGrp  
    cGrp = adjustl(cGrp)   
     
    write(cFileNum,'(i12)') iFileNum  
    cFileNum = adjustl(cFileNum)    
    
     
    write(cFq,'(i12)') iFq  
    cFq = adjustl(cFq)     
    
    write(cKx,'(i12)') iKx  
    cKx = adjustl(cKx)    
    
    cFilename  = trim(scratchFolder)//'/mare2dem_tmp.'//trim(cID)//'.'//trim(cScratchPrefix)//'.'// &
                 &trim(cGrp)//'.'//trim(cFileNum)//'.'//trim(cFq)//'.'//trim(cKx) 
                   
    end subroutine getScratchFilename    

!==================================================================================================================================! 
!=========================================================================================================== inverseFourierTransform
!==================================================================================================================================! 
    subroutine inverseFourierTransform
 
    use SinCosFilters
    use mare2dem_input_data_params, only: phaseConvention
    
    integer                                 :: iFq,iComp,nift,nc, iRx, iTx
    real(8)                                 :: xr   
    complex(8), dimension(:), allocatable   :: fEven, fOdd  

    integer(8)      :: tStart,tEnd  
    character(256)  :: sFmt
    character(32)   :: stime
    
    if (sType=='mt') return
    
    call system_clock(tStart) 
    
 !
 ! Initialize digital filters
 !
    call spline_integ_kx_setup()
    
    allocate(fEven(nwave),fOdd(nwave))
 
    nift = 0
    
    do iFq = 1,nfrequencies
    
        nc = solutions(iFq)%nComponents
        
        do iComp = 1,nc/2 ! since nTx array is 2x larger for -ikx sources              

            fEven = (solutions(iFq)%Fields(iComp,:) + solutions(iFq)%Fields(iComp+nc/2,:))/2d0
            fOdd  = (solutions(iFq)%Fields(iComp,:) - solutions(iFq)%Fields(iComp+nc/2,:))/2d0

!           write(*,*) 'solutions:'
!            write(*,*) (wavenum(i0),dble(solutions(iFq)%Fields(iComp,i0)),aimag(solutions(iFq)%Fields(iComp,i0)),i0=1,nwave)
!            write(*,*) ' '
!                        
!            write(*,*) 'Even:'
!            write(*,*) (wavenum(i0),dble(fEven(i0)),aimag(fEven(i0)),i0=1,nwave)
!            write(*,*) ' '
!             write(*,*) 'Odd:'
!            write(*,*) (wavenum(i0),dble(fOdd(i0)),aimag(fOdd(i0)),i0=1,nwave)
!            write(*,*) ' ' 
                                  
            iRx = solutions(iFq)%Components(iComp,1)
            iTx = solutions(iFq)%Components(iComp,2)
            
            xr  = xRx(iRx) - xTx(iTx)
            
            call spline_integ_fkx(xr, fEven, 0, fEven(1))   
            call spline_integ_fkx(xr, fOdd,  1, fOdd(1))

            ! Store result back in Fields:
            solutions(iFq)%Fields(iComp,1) = fEven(1) + fOdd(1)
 
            nift = nift + 1
            
            select case (trim(phaseConvention))

            case ('lag')
                    ! do nothing since em2dkx uses lag by default
            case ('lead')
                solutions(iFq)%Fields(iComp,1) = conjg(solutions(iFq)%Fields(iComp,1))
            end select
            
        enddo
        
    enddo
    
    call spline_integ_kx_deallocate

    deallocate(fEven,fOdd)

    call system_clock(tEnd) 
 
    write(stime,'(a8,f9.3,a2)') ' Timer: ',(tEnd - tStart)/clockRate, ' s'

    sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x, 61x, 9x, a32,2x,i7,    2x,a19)'     
    if (lDisplayRefinementStats) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp, '# Field Transforms:',nift,   trim(stime)
                      
    end subroutine inverseFourierTransform

!==================================================================================================================================! 
!======================================================================================================== inverseFourierTransform_DC
!==================================================================================================================================! 
    subroutine inverseFourierTransform_DC
 
    use SinCosFilters

    integer                                 :: iFq,iComp,nift,nc, iRx, iTx
    real(8)                                 :: xr   
    complex(8), dimension(:), allocatable   :: fEven 

    integer(8)      :: tStart,tEnd  
    character(256)  :: sFmt
    character(32)   :: stime
 
    
    call system_clock(tStart) 
    
 !
 ! Initialize digital filters
 !
    call spline_integ_kx_setup()
    
    allocate(fEven(nwave))
 
    nift = 0
    

    iFq = 1

    nc = solutions(iFq)%nComponents

    do iComp = 1,nc            

        fEven = solutions(iFq)%Fields(iComp,:) 

        iRx = solutions(iFq)%Components(iComp,1) ! local Rx electrode
        iTx = solutions(iFq)%Components(iComp,2) ! local Tx electrode
    
        xr =  electrode_rx(iRx,1) - electrode_tx(iTx,1) 
    
        call spline_integ_fkx(xr, fEven, 0, fEven(1))   
 
        ! Store result back in Fields:
        solutions(iFq)%Fields(iComp,1) = fEven(1)  

        nift = nift + 1
    
    
    enddo
 
    call spline_integ_kx_deallocate

    deallocate(fEven)

    call system_clock(tEnd) 
 
    write(stime,'(a8,f9.3,a2)') ' Timer: ',(tEnd - tStart)/clockRate, ' s'

    sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x, 61x, 9x, a32,2x,i7,    2x,a19)'     
    if (lDisplayRefinementStats) write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp, '# Field Transforms:',nift,   trim(stime)
                      
    end subroutine inverseFourierTransform_DC


!==================================================================================================================================! 
!==================================================================================================== inverseFourierTransform_Derivs
!==================================================================================================================================! 
    subroutine inverseFourierTransform_Derivs
 
    use SinCosFilters
    use mare2dem_input_data_params, only: phaseConvention
    
    implicit none 
       
    integer :: iFq,nc,iScratch, nScratchGroups, nComps,iComp, iFree, iRx, iTx, nift, ic, ierr
    real(8) :: xr 
 
    integer(8)      :: tStart,tEnd 
     
    complex(8), dimension(:), allocatable   :: fEven, fOdd  
      
    character(256)  :: sFmt
    character(32)   :: stime
    
    if (sType == 'mt') return
    
    call system_clock(tStart) 
    
 !
 ! Initialize digital filters
 !
    call spline_integ_kx_setup()
    allocate(fEven(nwave),fOdd(nwave))
 
 !
 ! Loop over frequencies:
 !   
    nift = 0
        
    do iFq = 1,nfrequencies
    
        nc = solutions(iFq)%nComponents

        nScratchGroups = ceiling(dble(nc/2)/dble(nComponentsPerScratch)) 
               
        !
        ! Allocate Read Buffer:
        !        
        allocate( readBuffer(nFree,2*nComponentsPerScratch,nwave), stat=ierr)    
        if (ierr .ne. 0) then
            write(*,*)  'Error allocating readBuffer in inverseFourierTransform_Derivs. '
            write(*,*) ' Too many free parameters: ', nFree,  nComponentsPerScratch, nwave
            stop 
        endif       

        ic = 0
        do iScratch = 1,nScratchGroups
 
            !
            ! Read in buffer:
            !              
            call loadScratch(iScratch, iFq, nComps)  ! loods in all kx for this iFq, iScratch
 
            !
            ! Inverse Fourier transform:
            !
            do iComp = 1,nComps ! nComps is actual number read in
                ic = ic + 1
                
                do iFree = 1,nFree
                    !write(*,*) 'hh: ',iComp,ic,iFree
                    !readBuffer(1:nFree,2*nc,iKx)   2*iComp - 1 is +ikx, 2*iComp is -ikx
                    fEven = (readBuffer(iFree,2*iComp-1, :) + readBuffer(iFree,2*iComp  , :))/2d0
                    fOdd  = (readBuffer(iFree,2*iComp-1, :) - readBuffer(iFree,2*iComp  , :))/2d0
                                   
                    iRx = solutions(iFq)%Components(ic,1)
                    iTx = solutions(iFq)%Components(ic,2)
                    
                    xr  = xRx(iRx) - xTx(iTx)

                    call spline_integ_fkx(xr, fEven, 0, fEven(1))   
                    call spline_integ_fkx(xr, fOdd,  1, fOdd(1))

                    ! Store result back in Fields:
                    solutions(iFq)%dFieldsdRho(ic,iFree) = fEven(1) + fOdd(1)
                    
                    select case (trim(phaseConvention))

                    case ('lag')
                            ! do nothing since em2dkx uses lag by default
                    case ('lead')
                        solutions(iFq)%dFieldsdRho(ic,iFree) = conjg(solutions(iFq)%dFieldsdRho(ic,iFree))
                    end select

                    nift = nift + 1

                enddo
                !write(*,'(a,i2,100(e8.1,1x))') 'invft ',ic, abs(solutions(iFq)%dFieldsdRho(ic,:)) 
            enddo            
            
        enddo
        
        deallocate(readBuffer)
        
    enddo
   
    call spline_integ_kx_deallocate 
    deallocate(fEven,fOdd)
    
    call system_clock(tEnd) 
   
    write(stime,'(a8,f9.3,a2)') ' Timer: ',(tEnd - tStart)/clockRate, ' s'
 
    sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x, 61x, 9x, a32,2x,i7,    2x,a19)'   
    if (lDisplayRefinementStats) write(*,sFmt) 'Proc:',myID,'Group:',iRefinementGrp, '# Derivative Transforms:',nift,trim(stime)
    
                   
    end subroutine inverseFourierTransform_Derivs
    
!==================================================================================================================================! 
!================================================================================================= inverseFourierTransform_Derivs_DC
!==================================================================================================================================! 
    subroutine inverseFourierTransform_Derivs_DC
 
    use SinCosFilters
 
    implicit none 
       
    integer :: iFq,nc,iScratch, nScratchGroups, nComps,iComp, iFree, iRx, iTx, nift, ic, ierr
    real(8) :: xr 
 
    integer(8)      :: tStart,tEnd 
     
    complex(8), dimension(:), allocatable   :: fEven   
      
    character(256)  :: sFmt
    character(32)   :: stime
    
    if (sType == 'mt') return
    
    call system_clock(tStart) 
    
 !
 ! Initialize digital filters
 !
    call spline_integ_kx_setup()
    allocate(fEven(nwave))
 
 !
 ! Loop over frequencies:
 !   
    nift = 0
   
    iFq = 1   
 
    
    nc = solutions(iFq)%nComponents

    nScratchGroups = ceiling(dble(nc)/dble(nComponentsPerScratch)) 
           
    !
    ! Allocate Read Buffer:
    !        
    allocate( readBufferDC(nFree,nComponentsPerScratch,nwave), stat=ierr)    
    if (ierr .ne. 0) then
        write(*,*)  'Error allocating readBuffer in inverseFourierTransform_Derivs. '
        write(*,*) ' Too many free parameters: ', nFree,  nComponentsPerScratch, nwave
        stop 
    endif       

    ic = 0
    do iScratch = 1,nScratchGroups

        !
        ! Read in buffer:
        !              
        call loadScratch(iScratch, iFq, nComps)  ! loods in all kx for this iFq, iScratch

        !
        ! Inverse Fourier transform:
        !
        do iComp = 1,nComps ! nComps is actual number read in
            ic = ic + 1
            
            do iFree = 1,nFree
                !write(*,*) 'hh: ',iComp,ic,iFree
                 
                fEven = readBufferDC(iFree,iComp, :) ! dble to complex
                
                iRx = solutions(iFq)%Components(ic,1) ! local Rx electrode
                iTx = solutions(iFq)%Components(ic,2) ! local Tx electrode

                xr =  electrode_rx(iRx,1) - electrode_tx(iTx,1) 

                call spline_integ_fkx(xr, fEven, 0, fEven(1))   
                 
                ! Store result back in Fields:
                solutions(iFq)%dFieldsdRho(ic,iFree) = fEven(1)  
         
                nift = nift + 1

            enddo
            !write(*,'(a,i2,100(e8.1,1x))') 'invft ',ic, abs(solutions(iFq)%dFieldsdRho(ic,:)) 
        enddo            
        
    enddo
    
    deallocate(readBufferDC)
 
    call spline_integ_kx_deallocate 
    deallocate(fEven)
    
    call system_clock(tEnd) 
   
    write(stime,'(a8,f9.3,a2)') ' Timer: ',(tEnd - tStart)/clockRate, ' s'
 
    sFmt = '(a5,2x,i6,2x,        a6,2x,i6,2x,  5x, 61x, 9x, a32,2x,i7,    2x,a19)'   
    if (lDisplayRefinementStats) write(*,sFmt) 'Proc:',myID,'Group:',iRefinementGrp, '# Derivative Transforms:',nift,trim(stime)
    
                   
    end subroutine inverseFourierTransform_Derivs_DC    
    
!==================================================================================================================================! 
!======================================================================================================================= loadScratch
!==================================================================================================================================! 
    subroutine loadScratch(iScratch,iFq,nComps)
 
    integer, intent(in)  :: iScratch, iFq   
    integer, intent(out) :: nComps 
    
    integer              :: iKx 
    
    character(256)       :: cFilename  
     
 
    do iKx = 1,nwave
        
        ! Get a name for the file:     
        call getScratchFilename(iRefinementGrp,iScratch,iFq,iKx,cFilename)
 
        ! Read in the scratch file when it is ready to be read:
        call readScratchFileWhenReady(cFilename,iKx,nComps)
        
    enddo
       
    end subroutine loadScratch

!==================================================================================================================================! 
!========================================================================================================== readScratchFileWhenReady
!==================================================================================================================================! 
    subroutine readScratchFileWhenReady(cFilename, iKx, nComps)

#if defined(__INTEL_COMPILER)
    use IFPORT          ! DGM 8/23/2013 Windows Intel compiler requirement for the sleepqq routine
#endif
 
    integer, intent(in)         :: iKx
    character(256), intent(in)  :: cFilename  
    integer, intent(out)        :: nComps
 
    logical :: lReady
    integer :: nc, ios 
    
!
! First inquire to see if the file exists:
! 
    lReady = .false. 
    
    do while(.not.lReady)
    
        inquire(file=cFileName,exist=lReady)    

        if (lReady) exit 

        !
        ! Else take a chill pill:
        !   
#if defined(__INTEL_COMPILER)  
        call sleepqq(10)     ! sleepqq time is in milliseconds.
#endif    
    
    
    enddo

!
! Next open file and see if marker integer has been written, indication file i/o is done
!        
    do  
        !write(*,*) 'opening: ',trim(cFilename), nfree, nc,ikx
        open(unit=21,file=trim(cFilename),status='old',access="stream",form="unformatted",action="read",iostat=ios)
        if (ios /= 0) then
            write(*,*) 'r: error, opening scratch file: ',trim(cFilename),  ' iostat:' ,ios
            stop
        endif

        read(unit=21,iostat=ios) nc
         if (ios /= 0) then
            write(*,*) 'r: error, reading nc from scratch file: ',trim(cFilename),  ' iostat:' ,ios
            stop
        endif        
         
        if (nc > 0) exit

        !
        ! Else take a chill pill:
        !   
#if defined(__INTEL_COMPILER)  
        call sleepqq(10)     ! sleepqq time is in milliseconds.
#endif    
    
    
    enddo    
 
 
!
! Now read in the file's data:
!         
    if (sType == 'cs') then
        read(21,iostat=ios)  readBuffer(1:nFree,1:2*nc,iKx)   
    else ! 'dc'
        read(21,iostat=ios)  readBufferDC(1:nFree,1:nc,iKx)   
    endif
 
    if (ios /= 0) then
        write(*,*) 'error, reading from scratch file on unit: ',21
        write(*,*) 'cFilename: ',trim(cFilename)
        write(*,*) 'ikx: ',ikx
        write(*,*) 'nc: ',nc
        stop
    endif     
    
    nComps = nc
    
!
! Finally, close and delete the file:
!       

    close(21, status='delete')
    
    end subroutine readScratchFileWhenReady

            
!==================================================================================================================================! 
!============================================================================================================ convert_to_data_format
!==================================================================================================================================! 
    subroutine convert_to_data_format 
    
    use mare2dem_global
    use mare2dem_input_data_params    
    use dataTransformations
    use EM_constants

    implicit none
     
    include 'em_parameters.inc'
            
    integer    :: i,iFq,itrans,isite

    complex(8) :: ex,ey,ez,bx,by,bz,hx,hy,hz
    complex(8) :: zte, ztm, htipper
    
    real(8)     :: vMA,vMB,vNA,vNB,dV
    real(8)     :: appres, dMA,dMB,dNA,dNB,xA(3),xB(3),xM(3),xN(3)
    
    complex(8), dimension(:), allocatable :: dzte_dsig,dztm_dsig,dhtipper_dsig,dVdSig
 
    allocate(dm(nd_local))
 
!
! This is a giant block that handles many different data types:
! 
    do i = 1,nd_local
    
        iFq     =  dp_local(i,2)
        itrans  =  dp_local(i,3)
        isite   =  dp_local(i,4)         
        
!
! Special for MT data:
!        
        if (sType == 'mt') then
            
            if (itrans <= 0) itrans = isite  ! if itrans not specified, make it the current receiver 

            ex = solGrid(iFq)%comps(isite,1,1)%field 
            ey = solGrid(iFq)%comps(isite,1,2)%field 
            ez = solGrid(iFq)%comps(isite,1,3)%field 
            hx = solGrid(iFq)%comps(itrans,1,4)%field
            hy = solGrid(iFq)%comps(itrans,1,5)%field
            hz = solGrid(iFq)%comps(itrans,1,6)%field 
 
            zte     = ex / hy                                   ! Ex/Hy
            ztm     = ey / hx                                   ! Ey/Hx
            htipper = solGrid(iFq)%comps(isite,1,6)%field / hy  ! Hz/Hy    ! note Hz is isite local Hz, in case hybrid tipper
                                                                                                                        
        elseif  (sType == 'cs') then
        
            ex = solGrid(iFq)%comps(isite,itrans,1)%field 
            ey = solGrid(iFq)%comps(isite,itrans,2)%field 
            ez = solGrid(iFq)%comps(isite,itrans,3)%field 
            bx = solGrid(iFq)%comps(isite,itrans,4)%field*mu0
            by = solGrid(iFq)%comps(isite,itrans,5)%field*mu0
            bz = solGrid(iFq)%comps(isite,itrans,6)%field*mu0 
        
        elseif  (sType == 'dc') then 
            ! get potentials for the two receiver electrodes created from two source electrodes
            ! we did this to allow for handling arbitrary geometry while retaining 2.5D kx cosine tfm symmetry 
            ! (whereas in 3D could compute directly )
            ! injection electrodes A and B and receive electrodes M and N.
            
            
            vMA = solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,1),1)%field ! ie. potential at M due to current injection at A
            vMB = solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,2),1)%field
            vNA = solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,1),1)%field
            vNB = solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,2),1)%field 
            
            ! potential difference from current source AB measured across electrodes MN
            dV = (vMA - vMB) - (vNA - vNB)
         
            xM = trodes_RxDCIn(RxDCIn(isite,1) ,1:3)
            xN = trodes_RxDCIn(RxDCIn(isite,2) ,1:3)
            xA = trodes_TxDCIn(TxDCIn(itrans,1),1:3)
            xB = trodes_TxDCIn(TxDCIn(itrans,2),1:3)
            
            dMA = norm2(xM-xA)
            dMB = norm2(xM-xB)
            dNA = norm2(xN-xA)
            dNB = norm2(xN-xB)
            
            ! apparent resistivity:         
            appres = PI2*dble(dv)/(1d0/dMA - 1d0/dMB - 1d0/dNA + 1d0/dNB) 
            
            !write(*,'(i,1x,4(g12.5,1x),4(f5.1,1x),4(i,1x))') i,vMA,vMB,vNA,vNB,dMA,dMB,dNA,dNB,RxDCIn(isite,1:2),TxDCIn(itrans,1:2)
            !write(*,'(i,1x,g,1x,4(i,1x),4(f8.1,1x))') i,appres,RxDCIn(isite,1:2),TxDCIn(itrans,1:2),xM(2),xN(2),xA(2),xB(2)
            
        endif
                
        !write(*,*) i,itrans,isite,ey    
            
        select case (  dp_local(i,1) ) ! what data type is this? 
                
!
!  CSEM data:                                     
!                              
        ! Electric: Real and Imaginary Data:
                
                case (indRealEx)  
                  dm(i) = dble( ex )
                case (indImagEx)  
                  dm(i) = aimag( ex )
                case (indRealEy)  
                  dm(i) = dble( ey )
                case (indImagEy)  
                  dm(i) = aimag( ey )
                case (indRealEz)  
                  dm(i) = dble( ez )  
                case (indImagEz)  
                  dm(i) = aimag( ez ) 
                    
        ! Electric: Amplitude and Phase Data:
                           
                case (indAmpEx)   
                  dm(i) = abs( ex )    !     |Ex|
                case (indAmpEy)   
                  dm(i) = abs( ey )    !     |Ey|
                case (indAmpEz)   
                  dm(i) = abs( ez )    !     |Ez|
                case (indPhsEx) 
                  dm(i) = getPhase( ex, d_local(i) )
                case (indPhsEy) 
                  dm(i) = getPhase( ey, d_local(i))
                case (indPhsEz) 
                  dm(i) = getPhase( ez, d_local(i))

        ! Electric: Log10 Amplitude Data:
                           
                case (indLog10AmpEx)   
                  dm(i) = log10(abs(ex))    ! log10|Ex|
                case (indLog10AmpEy)   
                  dm(i) = log10(abs(ey))    ! log10|Ey|
                case (indLog10AmpEz)   
                  dm(i) = log10(abs(ez))    ! log10|Ez|                      
                  
        ! Magnetic: Real and Imaginary Data:

                case (indRealBx) 
                  dm(i) = dble(bx)   !        real(Bx)
                case (indImagBx) 
                  dm(i) = aimag(bx)  !        imag(Bx)
                case (indRealBy) 
                  dm(i) = dble(by)   !        real(By)
                case (indImagBy) 
                  dm(i) = aimag(by)  !        imag(By)
                case (indRealBz) 
                  dm(i) = dble(bz)   !        real(Bz)
                case (indImagBz) 
                  dm(i) = aimag(bz)  !        imag(Bz)
                    
        ! Magnetic: Amplitude and Phase Data:
                    
                case (indAmpBx) 
                  dm(i) = abs(bx)   !     |Bx|
                case (indAmpBy) 
                  dm(i) = abs(by)   !     |By|
                case (indAmpBz) 
                  dm(i) = abs(bz)   !     |Bz|
                case (indPhsBx) 
                  dm(i) = getPhase(bx,(d_local(i)))
                case (indPhsBy) 
                  dm(i) = getPhase(by,(d_local(i)))
                case (indPhsBz) 
                  dm(i) = getPhase(bz,(d_local(i)))
                  
        ! Magnetic: Log10 Amplitude Data:
                           
                case (indLog10AmpBx)   
                  dm(i) = log10(abs(bx))    ! log10|Bx|
                case (indLog10AmpBy)   
                  dm(i) = log10(abs(by))    ! log10|By|
                case (indLog10AmpBz)   
                  dm(i) = log10(abs(bz))    ! log10|Bz|                      
                  
                                 
         ! Polarization Ellipse Parameters:
         
                case (iPEmax)
                  dm(i) = getPE(ex,ey,'pmax') ! E-field max
                case (iPEmin)
                  dm(i) = getPE(ex,ey,'pmin') ! E-field min
                case (iPBmax)
                  dm(i) = getPE(bx,by,'pmax') ! B-field max
                case (iPBmin)
                  dm(i) = getPE(bx,by,'pmin') ! B-field min
                  
!
!  MT data:                                     
!                     
                  
              ! TE:            
                case (indRhoZXY)                                        ! TE Apparent resistivity
                    dm(i) = getAppRes(zte,frequencies(iFq))    
                case (indlog10RhoZXY)                                   ! log10 TE Apparent resistivity
                    dm(i) = log10(getAppRes(zte,frequencies(iFq)))                                                         
                case (indPhsZXY)                                        ! TE Phase
                    dm(i) = getPhase(zte,(d_local(i)))                  
                case (indRealZXY)                                       ! real(Zxy)
                    dm(i) = dble( zte ) 
                case (indImagZXY)                                       ! imag(Zxy)
                    dm(i) = aimag(zte )      
                                        
                !TM:    
                case (indRhoZYX)                                        ! TM Apparent resistivity
                    dm(i) = getAppRes( -ztm,frequencies(iFq))                         
                case (indlog10RhoZYX)                                   ! log10 TM Apparent resistivity
                    dm(i) = log10(getAppRes( -ztm,frequencies(iFq)))     
                case (indPhsZYX)                                        ! TM Phase (moved to first quadrant from third)   
                    dm(i) = getPhase( -ztm,(d_local(i)))                  
                case (indRealZYX)                                       ! real(Zyx)
                    dm(i) = dble( ztm )             
                case (indImagZYX)                                       ! imag(Zyx)
                    dm(i) = aimag( ztm ) 
                
                                             
                ! Magnetic Tipper:    
                case (indRealMZY)                                       ! real(Mzy) ! TE mode Hz = Mzy*Hy
                    dm(i) = dble(  htipper )        
                case (indImagMZY)                                       ! imag(Mzy)
                    dm(i) = aimag( htipper )         
                
                case (indAmpMZY)                                        ! |Mzy| ! TE mode Hz = Mzy*Hy
                    dm(i) = abs(htipper)     
                case (indPhsMZY)                                        ! phase(Mzy)
                    dm(i) = getPhase( htipper, (d_local(i)) )      
                    
                ! Raw MT fields (for model study purposes only!):
                
                case (indRealExMT)  
                  dm(i) = dble( ex )
                case (indImagExMT)  
                  dm(i) = aimag( ex )
                case (indRealEyMT)  
                  dm(i) = dble( ey )
                case (indImagEyMT)  
                  dm(i) = aimag( ey )
                case (indRealEzMT)  
                  dm(i) = dble( ez )  
                case (indImagEzMT)  
                  dm(i) = aimag( ez )  
                                                   
                case (indRealHxMT)  
                  dm(i) = dble( hx )
                case (indImagHxMT)  
                  dm(i) = aimag( hx )
                case (indRealHyMT)  
                  dm(i) = dble( hy )
                case (indImagHyMT)  
                  dm(i) = aimag( hy )
                case (indRealHzMT)  
                  dm(i) = dble( hz )  
                case (indImagHzMT)  
                  dm(i) = aimag( hz )                                                         
!
! DC resisivity:
!
                case(indAppRes_DC)
                    dm(i) = appres
                case(indlog10AppRes_DC)
                    dm(i) = log10(appres)
      
        end select  ! case dp(i,1)
 
    enddo  ! loop over nd
 
 if (linversion) then    
 
    allocate(wj(nd_local,nFree))
    if (sType=='dc') then
        allocate(dVdSig(nfree))
    elseif (sType=='mt') then
        allocate(dzte_dsig(nfree),dztm_dsig(nfree),dhtipper_dsig(nfree))
    endif 
            
!
! This is a giant block that handles many different data types:
! 
    do i = 1,nd_local
    
        iFq     =  dp_local(i,2)
        itrans  =  dp_local(i,3)
        isite   =  dp_local(i,4)         
        
!
! Special for MT data:
!        
        if (sType == 'mt') then

            if (itrans <= 0) itrans = isite  ! if itrans not specified, make it the current receiver 
  
            zte     = solGrid(iFq)%comps(isite,1,1)%field / solGrid(iFq)%comps(itrans,1,5)%field  ! Ex/Hy
            ztm     = solGrid(iFq)%comps(isite,1,2)%field / solGrid(iFq)%comps(itrans,1,4)%field  ! Ey/Hx
            htipper = solGrid(iFq)%comps(isite,1,6)%field / solGrid(iFq)%comps(itrans,1,5)%field  ! Hz/Hy    
        
                                          
            dzte_dsig  = (solGrid(iFq)%comps(isite,1,1)%dFieldsdRho - zte * solGrid(iFq)%comps(itrans,1,5)%dFieldsdRho ) / solGrid(iFq)%comps(itrans,1,5)%field 
            dztm_dsig  = (solGrid(iFq)%comps(isite,1,2)%dFieldsdRho - ztm * solGrid(iFq)%comps(itrans,1,4)%dFieldsdRho ) / solGrid(iFq)%comps(itrans,1,4)%field  
 
            dhtipper_dsig = (solGrid(iFq)%comps(isite,1,6)%dFieldsdRho - htipper * solGrid(iFq)%comps(itrans,1,5)%dFieldsdRho ) / solGrid(iFq)%comps(itrans,1,5)%field  
                  
        endif
 
        select case (  dp_local(i,1) )  ! what data type is this?
 
!
! CSEM data:
!              
           ! Electric: Real and Imaginary Data:
                      
              case (indRealEx)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho)    !        real(Ex)
              case (indImagEx)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho)   !        imag(Ex)
              case (indRealEy)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho)    !        real(Ey)
              case (indImagEy)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho)   !        imag(Ey)
              case (indRealEz)  
                  wj(i, :) = dble( solGrid(iFq)%comps(isite,itrans,3)%dFieldsdRho)   !        real(Ez)
              case (indImagEz)  
                  wj(i, :) = aimag( solGrid(iFq)%comps(isite,itrans,3)%dFieldsdRho)  !        imag(Ez)
                              
           ! Electric: Amplitude and Phase Data:
                                 
              case (indAmpEx)   
                  wj(i, :) = absDeriv( solGrid(iFq)%comps(isite,itrans,1)%field, solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho )   !     |Ex|
              case (indAmpEy)   
                  wj(i, :) = absDeriv( solGrid(iFq)%comps(isite,itrans,2)%field, solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho )   !     |Ey|
              case (indAmpEz)   
                  wj(i, :) = absDeriv( solGrid(iFq)%comps(isite,itrans,3)%field, solGrid(iFq)%comps(isite,itrans,3)%dFieldsdRho )   !     |Ez|
              case (indPhsEx) 
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,1)%field, solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho )   !     Ex phase 
              case (indPhsEy)
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,2)%field, solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho )   !     Ey phase 
              case (indPhsEz) 
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,3)%field, solGrid(iFq)%comps(isite,itrans,3)%dFieldsdRho )   !     Ez phase 

        ! Electric: Log10 Amplitude Data:
                        
              case (indLog10AmpEx)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,1)%field, solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho)  ! log10|Ex|
              case (indLog10AmpEy)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,2)%field, solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho)  ! log10|Ey|
              case (indLog10AmpEz)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,3)%field, solGrid(iFq)%comps(isite,itrans,3)%dFieldsdRho)  ! log10|Ez|    

           ! Magnetic: Real and Imaginary Data:

              case (indRealBx) 
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho)*mu0   !        real(Bx)
              case (indImagBx) 
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho)*mu0  !        imag(Bx)
              case (indRealBy) 
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho)*mu0   !        real(By)
              case (indImagBy) 
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho)*mu0  !        imag(By)
              case (indRealBz) 
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,itrans,6)%dFieldsdRho)*mu0   !        real(Bz)
              case (indImagBz) 
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,itrans,6)%dFieldsdRho)*mu0  !        imag(Bz)
                              
           ! Magnetic: Amplitude and Phase Data:
                          
              case (indAmpBx) 
                  wj(i, :) = absDeriv( solGrid(iFq)%comps(isite,itrans,4)%field*mu0, solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho*mu0 )   !     |Bx|
              case (indAmpBy) 
                  wj(i, :) = absDeriv( solGrid(iFq)%comps(isite,itrans,5)%field*mu0, solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho*mu0 )   !     |By|
              case (indAmpBz) 
                  wj(i, :) = absDeriv(solGrid(iFq)%comps(isite,itrans,6)%field*mu0 , solGrid(iFq)%comps(isite,itrans,6)%dFieldsdRho*mu0 )   !     |Bz|
              case (indPhsBx) 
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,4)%field, solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho )   !     Bx phase 
              case (indPhsBy) 
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,5)%field, solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho )   !     By phase 
              case (indPhsBz) 
                  wj(i, :) = getPhaseDeriv( solGrid(iFq)%comps(isite,itrans,6)%field, solGrid(iFq)%comps(isite,itrans,6)%dFieldsdRho )   !     Bz phase 

          ! Magnetic: Log10 Amplitude Data:
                        
              case (indLog10AmpBx)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,4)%field*mu0, solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho*mu0)  ! log10|Bx|
              case (indLog10AmpBy)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,5)%field*mu0, solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho*mu0)  ! log10|By|
              case (indLog10AmpBz)   
                  wj(i, :) = log10absDeriv(solGrid(iFq)%comps(isite,itrans,6)%field*mu0, solGrid(iFq)%comps(isite,itrans,6)%dFieldsdRho*mu0)  ! log10|Bz|    

       
           ! Polarization Ellipse Parameters:   
           
              case (iPEmax)
                  wj(i, :) = getPEDeriv(solGrid(iFq)%comps(isite,itrans,1)%field,solGrid(iFq)%comps(isite,itrans,2)%field, &
                                        solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho,solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho,'pmax') ! E-field max
              case (iPEmin)
                  wj(i, :) = getPEDeriv(solGrid(iFq)%comps(isite,itrans,1)%field,solGrid(iFq)%comps(isite,itrans,2)%field, &
                                        solGrid(iFq)%comps(isite,itrans,1)%dFieldsdRho,solGrid(iFq)%comps(isite,itrans,2)%dFieldsdRho,'pmin') ! E-field min
              case (iPBmax)
                  wj(i, :) = getPEDeriv(solGrid(iFq)%comps(isite,itrans,4)%field*mu0,solGrid(iFq)%comps(isite,itrans,5)%field*mu0, &
                                        solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho*mu0,solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho*mu0,'pmax') ! B-field max
              case (iPBmin)
                  wj(i, :) = getPEDeriv(solGrid(iFq)%comps(isite,itrans,4)%field*mu0,solGrid(iFq)%comps(isite,itrans,5)%field*mu0, &
                                        solGrid(iFq)%comps(isite,itrans,4)%dFieldsdRho*mu0,solGrid(iFq)%comps(isite,itrans,5)%dFieldsdRho*mu0,'pmin') ! B-field min

!
! MT data:
!
            ! TE mode:
                case (indRhoZXY)                                        ! TE Apparent resistivity
                    wj(i, :) = getAppResDeriv(zte,dzte_dsig(:),frequencies(iFq) )          
                case (indlog10RhoZXY)                                    
                    wj(i, :) = getlog10AppResDeriv(zte,dzte_dsig(:))  
                case (indPhsZXY)                                        ! TE Phase
                    wj(i, :) = getPhaseDeriv( zte,dzte_dsig(:) )        
                case (indRealZXY)                                       ! real(Zxy)
                   wj(i, :) = real( (dzte_dsig(:)) )
                case (indImagZXY)                                       ! imag(Zxy)
                    wj(i, :) = aimag( (dzte_dsig(:)) )                           

            ! TM mode:                                       
                case (indRhoZYX)                                        ! TM Apparent resistivity
                    wj(i, :) = getAppResDeriv( -ztm,-dztm_dsig(:),frequencies(iFq) )    
                case (indlog10RhoZYX)                                        ! TM Apparent resistivity
                    wj(i, :) = getlog10AppResDeriv( -ztm,-dztm_dsig(:))                      
                case (indPhsZYX)                                        ! TM Phase    
                    wj(i, :) = getPhaseDeriv( -ztm, -dztm_dsig(:) )    
                case (indRealZYX)                                       ! real(Zyx)
                    wj(i, :) = real( dztm_dsig(:) )     
                case (indImagZYX)                                       ! imag(Zyx)
                   wj(i, :) = aimag( dztm_dsig(:) )    
   
            ! Magnetic Tipper:    
                case (indRealMZY)                                       ! real(Mzy) ! TE mode Hz = Mzy Hy
                   wj(i, :) = real(  dhtipper_dsig(:) )         
                case (indImagMZY)                                       ! imag(Mzy)
                   wj(i, :) = aimag( dhtipper_dsig(:) )         

                case (indAmpMZY)                                        ! |Mzy| ! TE mode Hz = Mzy*Hy
                    wj(i, :) = absDeriv(htipper,dhtipper_dsig(:))      
                case (indPhsMZY)                                        ! phase(Mzy)
                    wj(i, :) = getPhaseDeriv(htipper,dhtipper_dsig(:))                 

           ! Raw MT fields (for model study purposes only!):
                      
              case (indRealExMT)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,1,1)%dFieldsdRho)    !        real(Ex)
              case (indImagExMT)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,1,1)%dFieldsdRho)   !        imag(Ex)
              case (indRealEyMT)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,1,2)%dFieldsdRho)    !        real(Ey)
              case (indImagEyMT)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,1,2)%dFieldsdRho)   !        imag(Ey)
              case (indRealEzMT)  
                  wj(i, :) = dble( solGrid(iFq)%comps(isite,1,3)%dFieldsdRho)   !        real(Ez)
              case (indImagEzMT)  
                  wj(i, :) = aimag( solGrid(iFq)%comps(isite,1,3)%dFieldsdRho)  !        imag(Ez)
                  
              case (indRealHxMT)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,1,4)%dFieldsdRho)    !        real(Hx)
              case (indImagHxMT)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,1,4)%dFieldsdRho)   !        imag(Hx)
              case (indRealHyMT)  
                  wj(i, :) = dble(solGrid(iFq)%comps(isite,1,5)%dFieldsdRho)    !        real(Hy)
              case (indImagHyMT)  
                  wj(i, :) = aimag(solGrid(iFq)%comps(isite,1,5)%dFieldsdRho)   !        imag(Hy)
              case (indRealHzMT)  
                  wj(i, :) = dble( solGrid(iFq)%comps(isite,1,6)%dFieldsdRho)   !        real(Hz)
              case (indImagHzMT)  
                  wj(i, :) = aimag( solGrid(iFq)%comps(isite,1,6)%dFieldsdRho)  !        imag(Hz)                  
                                      
!
! DC resisivity:
!
                case(indAppRes_DC,indlog10AppRes_DC)
            
                    ! potential difference from current source AB measured across electrodes MN
                    ! dVdSig = (vMA - vMB) - (vNA - vNB)
                    
                    dVdSig =   (   solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,1),1)%dFieldsdRho   &
                           &     - solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,2),1)%dFieldsdRho ) &
                           & - (   solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,1),1)%dFieldsdRho   &
                           &     - solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,2),1)%dFieldsdRho )

                    vMA = solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,1),1)%field ! ie. potential at M due to current injection at A
                    vMB = solGrid(iFq)%comps(RxDCIn(isite,1),TxDCIn(itrans,2),1)%field
                    vNA = solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,1),1)%field
                    vNB = solGrid(iFq)%comps(RxDCIn(isite,2),TxDCIn(itrans,2),1)%field 
            
                    ! potential difference from current source AB measured across electrodes MN
                    dV = (vMA - vMB) - (vNA - vNB)         
                                      
                    xM = trodes_RxDCIn(RxDCIn(isite,1) ,1:3)
                    xN = trodes_RxDCIn(RxDCIn(isite,2) ,1:3)
                    xA = trodes_TxDCIn(TxDCIn(itrans,1),1:3)
                    xB = trodes_TxDCIn(TxDCIn(itrans,2),1:3)
            
                    dMA = norm2(xM-xA)
                    dMB = norm2(xM-xB)
                    dNA = norm2(xN-xA)
                    dNB = norm2(xN-xB)
                   
                    appres = PI2*dble(dv)/(1d0/dMA - 1d0/dMB - 1d0/dNA + 1d0/dNB) 
            
                     
                    wj(i, :) = PI2*dble(dVdSig)/(1d0/dMA - 1d0/dMB - 1d0/dNA + 1d0/dNB)   
                    !write(*,*) 'b4:',i,sum(abs(wj(i,:))),abs(dVdSig)
                    if ( dp_local(i,1) == indlog10AppRes_DC ) wj(i, :) = wj(i, :)*log10(exp(1d0))/appres
                 
                    
           end select  ! case dp(i,1)

         end do
    
        if (allocated(dVdSig)) deallocate(dVdSig)
        if (allocated(dzte_dsig)) deallocate(dzte_dsig,dztm_dsig,dhtipper_dsig)
        
!        write(*,*) 'convert_to_data_format:'
!        do i = 1,size(wj,1)
!            write(*,'(100(e9.2,1x))') (wj(i,j),j = 1,size(wj,2))
!        enddo
        
      endif ! linversion
      
        
    end subroutine convert_to_data_format      


    end module mare2dem_worker   

