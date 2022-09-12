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

    module mare2dem_io
    
    use string_helpers
    use kx_io
    use mare2dem_global
    use mare2dem_input_data_params
    use Occam 

    implicit none
    
    private
    
    ! Contained subroutines:
    public :: mare2dem_help 
    public :: displayBanner
    public :: readModel
    public :: readData
    public :: readPoly
    ! public :: writeResponse
!     public :: writeJacobian
!     public :: writeSensitivity
    public :: writeResistivityFile
    public :: print2col32
    
    contains

!==================================================================================================================================! 
!===================================================================================================================== mare2dem_help
!==================================================================================================================================! 
    subroutine mare2dem_help
    
    write(*,'(a)') ' '
    write(*,'(a)') 'Usage:  MARE2DEM [-F] [-S] [-J] [-scratch <scratchfolder>] <resistivity file> <output root name>'
    write(*,'(a)') ' '
    write(*,'(a)') 'MARE2DEM has four optional flags (*normal inversion uses no flags):'
    write(*,'(a)') ' '
    write(*,'(a)') '      -F    Computes the forward response of the starting model only.'
    write(*,'(a)') '            The forward response is output to <resistivity file>.resp'
    write(*,'(a)') ' '
    write(*,'(a)') '      -S    Outputs the sensitivity vector of the STARTING model for each iteration' 
    write(*,'(a)') '            The sensitivity is output to <resistivity file>.<iter#>.sensitivity.'    
    write(*,'(a)') '            The sensitivity is computed by summing the absolute value of each '
    write(*,'(a)') '            row of the uncertainty weighted Jacobian matrix then dividing by  '   
    write(*,'(a)') '            the parameter area. See eq. 7 in Schwalenberg et al. (2002). '       
    write(*,'(a)') ' '      
    write(*,'(a)') '      -J    Outputs the full Jacobian matrix of the STARTING model for each iteration.'
    write(*,'(a)') '            The Jacaobian matrix is output to <resistivity file>.<iter#>.jacobianBin'
    write(*,'(a)') '            Note: J is unweighted (i.e. it has NOT been normalized by the data uncertainties).'             
    write(*,'(a)') ' '   
    write(*,'(a)') '      -scratch <scratchfolder>  Use the specified directory for the scratch files'   
    write(*,'(a)') '            required for 2.5D CSEM inversion (but not MT). Optimally '    
    write(*,'(a)') '            this should be a local directory on each compute node and not'
    write(*,'(a)') '            a networked directory.'      
    write(*,'(a)') ' '
    write(*,'(a)') 'MARE2DEM has one required parameter:'
    write(*,'(a)') ' '    
    write(*,'(a)') '      <resistivity file> - This is the name of the input resistivity  '
    write(*,'(a)') '      file. By convention, this file should have the extension .resistivity. '
    write(*,'(a)') '      For example inputModel.0.resistivity.  The model found by each'
    write(*,'(a)') '      inversion iteration is then output to a new resistivity file with the '
    write(*,'(a)') '      iteration number incremented.  For example: inputModel.1.resistivity, '
    write(*,'(a)') '      inputModel.2.resistivity, ... The corresponding model responses are '
    write(*,'(a)') '      written to inputModel.1.resp, inputModel.2.resp,... '
    write(*,'(a)') ' ' 
    write(*,'(a)') 'MARE2DEM has one optional parameter:'  
    write(*,'(a)') ' '              
    write(*,'(a)') '      <output file root> - With this option, the output files are '
    write(*,'(a)') '      named <outputfileroot>.1.resistivity, <outputfileroot>.1.resp, ' 
    write(*,'(a)') '      named <outputfileroot>.2.resistivity, <outputfileroot>.2.resp,... ' 
    write(*,'(a)') ' ' 
    
    call exitMARE2DEM()
     
    end subroutine mare2dem_help
    
!==================================================================================================================================! 
!===================================================================================================================== displayBanner
!==================================================================================================================================!       
    subroutine displayBanner()
 
    implicit none
 
    
    write(*,*) ' ' 
    write(*,*) '============================= MARE2DEM ==================================='
    write(*,*) ' '  
    write(*,*) ' MARE2DEM: Modeling with Adaptively Refined Elements for 2.5D EM'
    write(*,*) ''
    write(*,'(2x,a)') m2d_version
    write(*,*) ' '
    write(*,*) ' A parallel goal-oriented adaptive finite element forward and inverse'
    write(*,*) ' modeling code for electromagnetic fields from electric dipoles, magnetic'
    write(*,*) ' dipoles and magnetotelluric sources in triaxially anisotropic conducting'
    write(*,*) ' media. Iterative adaptive mesh refinement is accomplished using the'
    write(*,*) ' goal-oriented error estimation method described in Key and Ovall (2011) '
    write(*,*) ' Inversion is accomplished with Occam''s method (Constable et al., 1987).'
    write(*,*) ' Key (2016) describes most of the features in the current 2018 version '
    write(*,*) ' of the code.'
    write(*,*) ' '
    write(*,*) ' When citing the code, please use the most recent reference:'
    write(*,*) ' '
    write(*,*) ' Key, K. MARE2DEM: a 2-D inversion code for controlled-source electromagnetic '
    write(*,*) '     and magnetotelluric data. Geophysical Journal International 207, '
    write(*,*) '     571–588 (2016).  '  
    write(*,*) ''
    write(*,*) ' This work is currently supported by: '
    write(*,*) ''
    write(*,*) ' Electromagnetic Methods Research Consortium'
    write(*,*) ' Lamont-Doherty Earth Observatory'
    write(*,*) ' Columbia University'
    write(*,*) ' http://emrc.ldeo.columbia.edu'
    write(*,*) ' '
    write(*,*) ' Originally funded by:'
    write(*,*) ''
    write(*,*) ' Seafloor Electromagnetic Methods Consortium '
    write(*,*) ' Scripps Institution of Oceanography '
    write(*,*) ' University of California San Diego'
    write(*,*) ''
    write(*,*) ' Copyright (C) 2017-2019'
    write(*,*) ' Kerry Key'
    write(*,*) ' Lamont-Doherty Earth Observatory'
    write(*,*) ' Columbia University'
    write(*,*) ' http://emlab.ldeo.columbia.edu'
    write(*,*) ' '
    write(*,*) ' Copyright (C) 2008-2016'
    write(*,*) ' Kerry Key'
    write(*,*) ' Scripps Institution of Oceanography'
    write(*,*) ' University of California, San Diego'
    write(*,*) ''
    write(*,*) ''
    write(*,*) ' This file is part of MARE2DEM.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is free software: you can redistribute it and/or modify'
    write(*,*) ' it under the terms of the GNU General Public License as published by'
    write(*,*) ' the Free Software Foundation, either version 3 of the License, or'
    write(*,*) ' (at your option) any later version.'
    write(*,*) ''
    write(*,*) ' MARE2DEM is distributed in the hope that it will be useful,'
    write(*,*) ' but WITHOUT ANY WARRANTY; without even the implied warranty of'
    write(*,*) ' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    write(*,*) ' GNU General Public License for more details.'
    write(*,*) ''
    write(*,*) ' You should have received a copy of the GNU General Public License'
    write(*,*) ' along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>. '
    write(*,*) ' ' 
    write(*,*) '==========================================================================' 
    
    end subroutine displayBanner
        
!==================================================================================================================================! 
!=============================================================================================================== readResistivityFile
!==================================================================================================================================!       
    subroutine readResistivityFile
! 
! Reads in the resistivity file that specifies all resistivity regions
! for all fixed and free parameters.
! 
 
 
    integer             :: ierr, iLine, i, j, iAllocErr, iostat_err, iskip, ict, iFmt, iparam1, iparam2
    character(256)      :: filename, vers
    character(256)      :: sLine, sCode, sValue,sValueToComma 
    character(32)       :: cnum
    logical             :: bComment
    real(RealPrec)      :: pmtol, rtemp 
    real(RealPrec)      :: inputMisfit, inputRoughness ! only locally read in and displayed here since Occam will recompute them...

    !
    ! Open the file:
    !
    filename = trim(resistivityFile)//'.resistivity' 
    
    write(6,*) '=============== Reading in Resistivity File =============================='
    write(*,*) ' '            
    call print2col32('Opening file: ',filename,6)
    open(unit=10,file=trim(filename),status='old',iostat = ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening Resistivity File: ',trim(filename)
        call exitMARE2DEM
    end if  
    write(*,*) ' '  
    
    nRhoPerRegion  = 0
    lowerBoundGlobal = 0
    upperBoundGlobal = 0
    
    ! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
    
    iLine = 0
    do  ! infinite while loop
        
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( 10, '(A)', iostat = ierr ) sLine
        
        if (ierr /= 0) exit  ! end of file read, escape from while loop and 
                        ! proceed to checking that required inputs defined
        
        iLine = iLine + 1   ! advance line counter
        
        call ParseCode(  sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        
        !
        ! What do we have?
        !
        select case (trim(sCode))
        
        case ('version','format')
        
            vers = sValue(1:len(vers))    
            call Lower(vers)
            select case (trim(vers))
            
            case ('mare2dem_1.0')
                iFmt = 1
            case ('mare2dem_1.1')   
                iFmt = 2
            case default
           
                write(*,*) 'Error, unsupported resistivity file version: ',trim(vers)
                write(*,*) ' Try using mare2dem_1.1 format.   Stopping!'
                call exitMARE2DEM
            end select
                         
        case ('model file')   
        
            modelFileName = sValue(1:len(modelFileName))         
            call print2col32('Model File: ', modelFileName,6)
        
        case ('data file')   
        
            dataFileName = sValue(1:len(dataFileName))         
            call print2col32('Data File: ', dataFileName,6)
        
        case ('settings file')   
        
            settingsFileName = sValue(1:len(settingsFileName))         
            call print2col32('Settings File: ', settingsFileName,6)
                    
        case ('penalty file')   
        
            penaltyFileName = sValue(1:len(penaltyFileName))         
            call print2col32('Penalty File: ', penaltyFileName,6)
        
        case ('maximum iterations')   
            read(sValue,*) maxOccamIterations
        
        case ('bounds transform')
            cBoundsTransform = trim(sValue)              
            
            ! Check to make sure input is valid:
            !
            select case (cBoundsTransform)
            case ('exponential','bandpass')
                call print2col32('Model bounds transform: ', cBoundsTransform,6) 
            case default 
                write(*,*) '  Error, unrecognized model bounds transform:  ', trim(cBoundsTransform)
                call exitMARE2DEM
            end select
        
        case ('global bounds')    
            ! Limit value string must have two values: min, max
            ! Comma is required.
            i = index( sValue, ',' )
            
            if (i == 0) then
                write(*,*) 'Invalid "global Bounds" entry in resistivity file.'
                write(*,*) 'The format for this line is'
                write(*,*) '   Global Bounds: <min>,<max>'
                write(*,*) 'Do not include angle brackets.  Comma required.'
                write(*,*) 'Min must be less than max.'

                call exitMARE2DEM
            endif
             
            read(sValue(:i-1),*) lowerBoundGlobal
            
            ! Skip space padding
            do i=i+1,len_trim(sValue)
                if (sValue(i:i) .ne. ' ') exit
            enddo
            read(sValue(i:),*) upperBoundGlobal
            
            ! Check for mistake:
            if (lowerBoundGlobal > upperBoundGlobal) then
                rtemp = lowerBoundGlobal
                lowerBoundGlobal = upperBoundGlobal
                upperBoundGlobal = rtemp 
            endif
            
            write(cnum,'(es11.4,1x,es11.4)') lowerBoundGlobal, upperBoundGlobal
            call print2col32(' Global bounds: ',cnum,6)       
            
        case ('debug level', 'print level')
            read(sValue,*) occamPrintLevel
            write(cnum,'(i3)') occamPrintLevel
            call print2col32('Print Level: ',cnum,6)

        case ('roughness penalty')
  
 			i = index( sValue, ',' )
 
			if (i /= 0) then
				sValueToComma = sValue(1:i-1)
			else
				sValueToComma  = sValue
			endif
	
			call Lower(sValueToComma)
			select case (trim(adjustl(sValueToComma)))
        
            case ('mgs','MGS')
				lMGS = .true.  
				if ((i .ne. 0 ).and.(i > len_trim(sValueToComma)) ) then
					read(sValue(i+1:),*,iostat = ierr) beta_mgs
						if (ierr /= 0) then
							write(*,*) ' Error beta value for MGS roughness penalty!'
							write(*,*) ' Line read in = ',trim(sLine)
							write(*,*) ' Format should be: Roughness penalty: MGS,<betavalue> '
							write(*,*) ' Stopping'
							call exitMARE2DEM
						endif        
				endif            	
            	call print2col32('Roughness penalty: ',trim(sValue),6)   

            	    
            end select      	
        	
        case ('target misfit')
            read(sValue,*) targetRMS
            write(cnum,'(es11.4)') targetRMS
            call print2col32('Target Misfit: ',cnum,6)       
            
        case ('iteration')    
        ! kwk: iteration number is now in file name, ignore any input numbers
 
        case  ('misfit decrease threshold')
            read(sValue,*) rmsThreshold
            write(cnum,'(g11.4)') rmsThreshold     
            call print2col32('Misfit Decrease Threshold: ', cnum,6)           
 
        case  ('converge slowly')
            call Lower(sValue)
            select case (trim(sValue))
            case ('yes')
                  lConvergeSlowly = .true.
            end select  
            call print2col32('Converge Slowly: ', trim(sValue),6)
            
        case ('lagrange value')  
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) modelMu  
                write(cnum,'(es11.4)') modelMu
                call print2col32('Lagrange Value: ',cnum,6)
            endif
            
        case ('model roughness') 
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) inputRoughness 
                write(cnum,'(es11.4)') inputRoughness
                call print2col32('Model Roughness: ',cnum,6)
            endif
        
        case ('model misfit')
            if (len_trim(sValue) > 0 ) then
                read(sValue,*) inputMisfit
                write(cnum,'(es11.4)') inputMisfit
                call print2col32('Model Misfit: ',cnum,6)
            endif            
        case  ('date/time')
            ! skip
            
!        case  ('inversion method')
!            call Lower(sValue)
!            sInversionMethod = sValue(1:len(sInversionMethod))         
!            call print2col32('Inversion Method: ', sInversionMethod,6)
            
!        case  ('fixed mu cut')
!            read(sValue,*) rMuCut
!            write(cnum,'(es11.4)') rMuCut     
!            call print2col32('Fixed Mu Cut: ', cnum,6)                        
!        
        case  ('anisotropy') 
            cAnisotropy = sValue(1:len(cAnisotropy))  
            ! Check to make sure this makes sense:    
            select case (trim(cAnisotropy ))
                    
            case ('isotropic')
                nRhoPerRegion = 1
            case ('triaxial')
                nRhoPerRegion = 3
            case ('isotropic_ip') ! Cole-Cole Model
                nRhoPerRegion = 4        
            case ('isotropic_complex')
                nRhoPerRegion = 2     
            case ('tix','tiy','tiz')
                nRhoPerRegion = 2
            
            case default
                write(*,*) ' Error decoding anisotropy in resistivity file'
                write(*,*) ' Unknown or unsupported anisotropy code:', trim(cAnisotropy )
                write(*,*) ' Stopping!'
                call exitMARE2DEM
            end select    
            
            call print2col32('Anisotropy: ',cAnisotropy,6)
          
!    
        case  ('number of regions')     
            read(sValue,*) nRegions
            write(cnum,'(i8)') nRegions
            call print2col32('Number of regions: ',cnum,6)
                            
            nRhoParams =   nRegions*nRhoPerRegion             
            if (nRhoParams < 1)  then
                write(*,*) ' Error in resistivity file, nRegions, nRho per region: ',nRegions, nRhoPerRegion
                write(*,*) ' Stopping!'
                call exitMARE2DEM    
            endif
            
          !
          ! Now read in the regions:
          !
        
          !
          ! First allocate storage:
          !
            allocate(rhoParams(nRhoParams), iFreeParam(nRhoParams), stat=iAllocErr)
            if (iAllocErr .ne. 0) then
                write(*,*) ' Out of memory.  Too many resistivity parameters (', nRhoParams, ')'
                call exitMARE2DEM 
            endif      
            
            allocate (boundsTemp(nRegions,2*nRhoPerRegion), &
                      PrejTemp(nRegions,2*nRhoPerRegion), stat=iAllocErr)   
            if (iAllocErr .ne. 0) then
                write(*,*) ' Out of memory.  Too many resistivity parameters (', nRegions*nRhoPerRegion, ')'
                call exitMARE2DEM 
            endif
            
            if ( nRhoPerRegion > 1 ) then 
                allocate (ratioTemp(nRegions,nRhoPerRegion*(nRhoPerRegion-1)), stat=iAllocErr)  ! 2 or 6 columns
                if (iAllocErr .ne. 0) then
                    write(*,*) ' Out of memory.  Too many resistivity parameters (', nRegions*nRhoPerRegion, ')'
                    call exitMARE2DEM 
                endif            
                ratioTemp = 0 
            endif
                                   
            !
            ! Read in each line:
            !  
            ! Skip comment line:
            read( 10, '(A)', iostat = ierr ) cParamHeader
            
            if ( (iFmt == 1)  .or. (nRhoPerRegion == 1).or. (nRhoPerRegion == 4)  )then ! Format: ', ' MARE2DEM_1.0 or isotropic or isotropic_ip
                do i = 1,nRegions
                    read(10,*,iostat = iostat_err) iskip, ( rhoParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),   & ! rho is linear
                    &                                      (iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                    &                                      (boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &                                      (  prejTemp(i,j),j=1,2*nRhoPerRegion) 
 
                    if (iostat_err /= 0) then
                        write(*,*) ' Error reading parameters, stopping!'
                        call exitMARE2DEM
                    endif        
                enddo
            
            else  ! nRhoPerRegion > 1 and new format with anisotropy ratio support:
            
                do i = 1,nRegions
                    read(10,*,iostat = iostat_err) iskip, ( rhoParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),   & ! rho is linear
                    &                                      (iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                    &                                      (boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &                                      (  prejTemp(i,j),j=1,2*nRhoPerRegion), &
                    &                                      (  ratioTemp(i,j),j=1,nRhoPerRegion*(nRhoPerRegion-1))     
 
                    if (iostat_err /= 0) then
                        write(*,*) ' Error reading parameters, stopping!'
                        call exitMARE2DEM
                    endif        
                enddo
                            
            endif
                        
        case default
            write(*,*) 'Error reading resistivity file!'
            write(*,*) ' On line :', iLine            
            write(*,*) 'Unknown or unsupported code:'
            write(*,*) sCode
        
            call exitMARE2DEM
        
        end select 
        
    enddo ! read while loop
    
    write(*,*) ' ' 

!
! Close resistivity file
!
    close(10)
    
 !
 ! Perform a few odds and ends here:
 !    
    
    nFree = count(iFreeParam > 0)
    nParams = nFree
    
    write(*,*) 'Done reading resistivity file, here is what I found: '
    write(*,*) ' '
    write(cnum,'(i8)') nFree
    call print2col32('# Free parameters: ',cnum,6)
    write(cnum,'(i8)') nRhoParams-nFree
    call print2col32('# Fixed parameters: ',cnum,6)
    write(*,*) ' '
       
    !
    ! Now fill in the free parameter arrays
    !
    allocate( pm(nFree), lowerBound(nFree),upperBound(nFree),lBoundMe(nFree),premod(nFree),prewts(nFree), stat=iAllocErr)
    if (iAllocErr .ne. 0) then
        write(*,*) ' Out of memory.  Too many free parameters (', nFree, ')'
        call exitMARE2DEM 
    endif      
    
    ! Initialize these to the globals, then override if local value is nonzero:
    ! Also don't forget that everything needs to be converted to log10 here:
    
    lowerBound = 0 
    upperBound = 0 
    lBoundMe = .false.
    
    if ( (lowerBoundGlobal > 0) .and. (upperBoundGlobal > 0) ) then
        lowerBound = log10(lowerBoundGlobal)
        upperBound = log10(upperBoundGlobal)
        lBoundMe         = .true. 
    endif 
    
    premod     = 0
    prewts     = 0          
    
    select case (trim(cAnisotropy ))
 
        case ('isotropic', 'tix','tiy','tiz','triaxial')
                             
        ict = 0
        do i = 1,nRegions
            do j = 1,nRhoPerRegion
            
                if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
                
                    ict = ict + 1
                
                    if ( ( boundsTemp(i,2*j-1) > 0 ) .and. ( boundsTemp(i,2*j  ) > 0 ) ) then
                        lowerBound(ict) = log10(boundsTemp(i,2*j-1)) ! (i,1 3 5)
                        upperBound(ict) = log10(boundsTemp(i,2*j  )) ! (i,2 4 6)    
                        lBoundMe(ict)   = .true.
                    endif      
                
                    pm(ict)     = log10(rhoParams((i-1)*nRhoPerRegion+j)) ! inversion for log10(rho)
                
                    prewts(ict) =     abs(  prejTemp(i,2*j  ) ) ! (i,2 4 6) ! weights are linear      
                    if ( prewts(ict) > 0) then
                        premod(ict) = log10(prejTemp(i,2*j-1)) ! (i,1 3 5)
                    endif                
                    !write(*,*) ict, pm(ict),premod(ict),prewts(ict), lowerBound(ict),upperBound(ict),lBoundMe(ict)       
                endif
            enddo
        enddo
    
        if (ict /= nFree) write(*,*) 'Error in readResistivityFile, ict /=nFree: ', ict, nFree
    
        case ('isotropic_ip')
    
        ! Cole-Cole IP:    
        ict = 0
    
        do i = 1,nRegions
                   
            do j = 1,nRhoPerRegion
       
                if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
                
                    ict = ict + 1
                    
                    if (j == 1) then ! Rho 
                    
                        ! inversion parameter:
                        pm(ict)  = log10(rhoParams((i-1)*nRhoPerRegion+j)) ! inversion for log10(rho)
               
                        ! bounds:
                        if ( ( boundsTemp(i,2*j-1) > 0 ) .and. ( boundsTemp(i,2*j  ) > 0 ) ) then
                            lowerBound(ict) = log10(boundsTemp(i,2*j-1)) ! (i,1 3 5)
                            upperBound(ict) = log10(boundsTemp(i,2*j  )) ! (i,2 4 6)    
                            lBoundMe(ict)   = .true.      
                        endif
                        
                        ! prejudice:
                        prewts(ict) =     abs(  prejTemp(i,2*j  ) ) ! (i,2 4 6) ! weights are linear      
                        if ( prewts(ict) > 0) then
                            premod(ict) = log10(prejTemp(i,2*j-1)) ! (i,1 3 5)
                        endif                                
                    
                    else ! parameter is one if Eta, Tau, C     
                        
                        ! make sure values are positive:
                        !linear: 
                        !rhoParams((i-1)*nRhoPerRegion+j) = max(0d0,rhoParams((i-1)*nRhoPerRegion+j))
                        !log:
                        !rhoParams((i-1)*nRhoPerRegion+j) = max(1d-8,rhoParams((i-1)*nRhoPerRegion+j))
                        
                        ! inversion parameter:
                        !linear: 
                        pm(ict)  = (rhoParams((i-1)*nRhoPerRegion+j)) ! no conversion, parameters stay linear
                        ! log:
                        !pm(ict)  = log10(rhoParams((i-1)*nRhoPerRegion+j)) ! convert to log scaling
 
                        ! bounds:
                        ! if user has specified bounds:
                        if ( ( boundsTemp(i,2*j-1) > 0 ) .or. ( boundsTemp(i,2*j  ) > 0 ) ) then
                            lowerBound(ict) = (boundsTemp(i,2*j-1)) ! (i,1 3 5)
                            upperBound(ict) = (boundsTemp(i,2*j  )) ! (i,2 4 6)    
                            lBoundMe(ict)   = .true. 
                        else !
 
                            !linear:
                            if (j == 2) then ! eta
                                lowerBound(ict) = 0d0
                                upperBound(ict) = 1d0
                                lBoundMe(ict)   = .true.
                            elseif (j == 3) then! tau                               
                                lowerBound(ict) = 0d0
                                upperBound(ict) = 20d0 
                                lBoundMe(ict)   = .true.
                            elseif (j == 4) then! c
                                lowerBound(ict) = 0d0
                                upperBound(ict) = 1d0 
                                lBoundMe(ict)   = .true.
                            endif      
                                                                            
                        !log:
!                            if (j == 2) then ! eta
!                                lowerBound(ict) = -6d0
!                                upperBound(ict) = 0d0
!                            elseif (j == 3) then! tau
!                                lowerBound(ict) = -3d0
!                                upperBound(ict) = log10(3d0) 
!                            elseif (j == 4) then! c
!                                lowerBound(ict) = -3d0
!                                upperBound(ict) = 0d0 
!                            endif                          
                                                          
                          
                                            
                        endif
                                      
                        ! prejudice:
                        prewts(ict) =     abs(  prejTemp(i,2*j  ) ) ! (i,2 4 6) ! weights are linear     
                 
                        if ( prewts(ict) > 0) then
                            ! linear:
                            premod(ict) = (prejTemp(i,2*j-1)) ! (i,1 3 5)

                            ! log:
                            !premod(ict) = log10(max(1d-8,prejTemp(i,2*j-1))) ! (i,1 3 5) ! will fail on 0 input!
                        endif                 

                    endif      
                    ! write(*,'(i6,1x,i6,1x,6(g12.2,1x))') ict,j,rhoParams((i-1)*nRhoPerRegion+j), pm(ict), lowerBound(ict),upperBound(ict), premod(ict) ,prewts(ict)
                endif ! free param
            
                                                                 
            enddo ! j = 1,nRhoPerRegion                               
            
        enddo ! i = 1,nRegions
       
       
        case ('isotropic_complex')
    
        ! IP complex resistivity:
        
        ict = 0
    
        do i = 1,nRegions
                   
            do j = 1,nRhoPerRegion
       
                if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
                
                    ict = ict + 1
                  
                    ! inversion parameter:
                    pm(ict)  = log10(rhoParams((i-1)*nRhoPerRegion+j)) ! inversion for log10(rho)


                    if (j == 1) then ! Rho real
                    
                        ! bounds:
                        if ( ( boundsTemp(i,2*j-1) > 0 ) .and. ( boundsTemp(i,2*j  ) > 0 ) ) then
                            lowerBound(ict) = log10(boundsTemp(i,2*j-1)) ! (i,1 3 5)
                            upperBound(ict) = log10(boundsTemp(i,2*j  )) ! (i,2 4 6)    
                            lBoundMe(ict)   = .true.      
                        endif
   
                    else ! Rho imag  
                         
                      
                        ! bounds:
                        ! if user has specified bounds:
                        if ( ( boundsTemp(i,2*j-1) > 0 ) .and. ( boundsTemp(i,2*j  ) > 0 ) ) then
                            lowerBound(ict) = log10(boundsTemp(i,2*j-1)) ! (i,1 3 5)
                            upperBound(ict) = log10(boundsTemp(i,2*j  )) ! (i,2 4 6)    
                            lBoundMe(ict)   = .true. 
                        else ! do not use global bounds
                            lowerBound(ict) = 0d0
                            upperBound(ict) = 0d0
                            lBoundMe(ict)  = .false. 
                        endif  
                        
                    endif    
                      
                    ! prejudice:
                    prewts(ict) =     abs(  prejTemp(i,2*j  ) ) ! (i,2 4 6) ! weights are linear      
                    if ( prewts(ict) > 0) then
                        premod(ict) = log10(prejTemp(i,2*j-1)) ! (i,1 3 5)
                    endif  
              
                endif ! free param
                                                          
            enddo ! j = 1,nRhoPerRegion                      

        enddo ! i = 1,nRegions
                
    end select
    
    !
    ! If model is anisotropic, create parameter difference preference array (log of ratio), if any
    !
    nDiff = 0
    if ( ( nRhoPerRegion > 1 ).and. ( nRhoPerRegion < 4 ) ) then ! anisotropic or isotropic_complex and not Cole-Cole IP
        
        allocate(preDiff(nFree),preDiffwts(nFree),ijDiff(nFree,2)) ! over-allocating, only up to nDiff values will be set
        preDiffwts = 0
        preDiff    = 0
        
        ict = 0
        do i = 1,nRegions
            
            do j = 1,nRhoPerRegion

                if (iFreeParam((i-1)*nRhoPerRegion+j)>0) then
    
                    ict = ict + 1 ! free parameter counter
                    
                    if ( (iFreeParam((i-1)*nRhoPerRegion+j)>0).and.(iFreeParam((i-1)*nRhoPerRegion+eps(j+1))>0) ) then
                    ! if this parameter and next one in this region are free:
                        
                        if ( (nRhoPerRegion==2) .and. (j>1) )  then
                        
                         ! do nothing
                         
                         elseif  (  (ratioTemp(i,2*j) > 0) .and. (ratioTemp(i,2*j-1) > 0 ) ) then ! add a ratio penalty for two anisotropic components:
                                         
                            nDiff = nDiff + 1
                    
                            iparam1 = ict
                            iparam2 = ict+1
                            if (j==3) iparam2 = ict - 2  ! x/y, y/z, and z/x, so skip back to 1 if j==3
                    
                            preDiffwts(nDiff) = ratioTemp(i,2*j  ) 
                            preDiff(nDiff)    = log10(ratioTemp(i,2*j-1))
                            ijDiff(nDiff,1)   = iparam1
                            ijDiff(nDiff,2)   = iparam2 
                        
                           ! write(*,'(4(i8,1x),2(g12.4,1x))') j,nDiff,iparam1,iparam2,preDiff(nDiff),preDiffwts(nDiff) 
                        
                        endif                        
                        
                    endif    
                endif  
            enddo ! j
        enddo ! i       
        
    endif
    
!
! Double check the input pm to make sure its within the bounds, then transform it to the unbound parameters
!      
    do i = 1,nFree
       
        if (lBoundMe(i)) then
!
! August 2015: If an input parameter exceeds a bound, let the user know and stop execution:
!  
            if ( (pm(i) < lowerBound(i)) .or.  (pm(i) > upperBound(i) ) ) then
                
                write(*,*) ' '
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '              
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
                write(*,*) '  '
                write(*,*) ' Error in the input resistivity file !  '
                write(*,*) ' The resistivity of parameter: ',i, ', which has been designated a free parameter, exceeds '
                write(*,*) ' either the global or local bounds on resistivity! '
                write(*,*) '  '
                write(*,*) ' Please either adjust the resistivity to be within the bounds or make that parameter fixed. '
                write(*,*) '  ' 
                write(*,*) ' log10: lowerBound: ',lowerBound(i),' pm(i): ', pm(i), 'upperBound: ',upperBound(i)
                write(*,*) 
                call exitMARE2DEM               
                
                stop
            endif
        
            pmtol = (upperBound(i) - lowerBound(i)) / 232.       
            
            if (pm(i) <= lowerBound(i) + pmtol) pm(i) = lowerBound(i) + pmtol 
            if (pm(i) >= upperBound(i) - pmtol) pm(i) = upperBound(i) - pmtol
            
            ! now transform it:
            pm(i)  =  transformToUnbound(pm(i) ,lowerBound(i) ,upperBound(i),lBoundMe(i) ) 
            
            ! transform prejudice model too:
            if ( prewts(i) > 0) then
                ! only transform when prejudice is NOT at one of the bounds, otherwise
                ! it will go to infinity and implicitly be heavily weighted. This was a problem
                ! for IP Cole-Cole inversions with (0,1) bounds on eta where parameters were 
                ! prejudiced to 0 for stability.  The transform below would move the 0 value prejudice to some really large
                ! negative value that made it implicitly a large prejudice weight (due to P^T*P*prej term in eqn 18 in 
                ! MARE2DEM paper).
                if ( (premod(i) > lowerBound(i) + pmtol)  .and. (premod(i) < upperBound(i) - pmtol) ) then
                    premod(i)  =  transformToUnbound(premod(i) ,lowerBound(i) ,upperBound(i),lBoundMe(i) ) 
                endif
            endif  
            
        endif   
      
    enddo 
 
    end subroutine readResistivityFile           
!==================================================================================================================================! 
!============================================================================================================== writeResistivityFile
!==================================================================================================================================!       
    subroutine writeResistivityFile(nCurrentIter,misfit,misfitRequested,lagrange,roughness)
! 
! Saves resistivity file from array rhoParams that specifies all resistivity regions
! for all fixed and free parameters. 
! 
    integer,                      intent(in) :: nCurrentIter
    real(RealPrec),               intent(in) :: misfit,lagrange,roughness,misfitRequested 
    
    real(RealPrec) :: logrho 
    character(256) :: cNum, sFmt 
    integer        :: lerr, ioUnit = 21, i, ict, j
    character(80)  :: dateAndTime

!    
! Open resistivity file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.resistivity', iostat=lerr)
  
    call print2col32('Format: ', ' MARE2DEM_1.1',ioUnit)
    call print2col32('Model File: ', modelFileName,ioUnit)
    call print2col32('Data File: ', dataFileName,ioUnit)
    call print2col32('Settings File: ', settingsFileName,ioUnit)
    call print2col32('Penalty File: ', penaltyFileName,ioUnit)
    write (cNum,'(i6)') maxOccamIterations
    call print2col32('Maximum Iterations: ',cnum ,ioUnit)
    call print2col32('Bounds Transform: ', cBoundsTransform,ioUnit)
    write(cnum,'(es12.4,1a,1x,es12.4)') lowerBoundGlobal,',', upperBoundGlobal
    call print2col32('Global Bounds: ', cnum,ioUnit)
    write (cNum,'(i6)') occamPrintLevel
    call print2col32('Print Level: ', cnum,ioUnit)
    write(cnum,'(g12.4)') misfitRequested !
    call print2col32('Target Misfit: ',cnum,ioUnit) 
    
 
    write(cnum,'(g12.4)') rmsThreshold !
    call print2col32('Misfit Decrease Threshold: ',cnum,ioUnit)        
    if (lConvergeSlowly) then
        cnum = 'yes'
    else
        cnum = 'no'
    endif
    call print2col32('Converge Slowly: ',cnum,ioUnit)        
    if (lMGS) then
    	write(cnum,'(a4,e12.4)') 'MGS,',beta_mgs
    	call print2col32('Roughness Penalty: ',cnum,ioUnit)  
    endif
    write(cnum,'(e12.4)') lagrange
    call print2col32('Lagrange Value: ',cnum,ioUnit)
    write(cnum,'(g12.4)') roughness
    call print2col32('Model Roughness: ',cnum,ioUnit)    
    write(cnum,'(f11.4)') misfit
    call print2col32('Model Misfit: ',cnum,ioUnit)
    call datime(dateAndTime) 
    call print2col32('Date/Time: ',dateAndTime,ioUnit)
    
    ! kwk debug: this is currently a test feature and might not ever be a supported feature, so don't output 'fixed' flag unless
    !            it is being used:
!    if (trim(sInversionMethod) =='fixed') then
!        call print2col32('Inversion Method: ', sInversionMethod,ioUnit)
!        write(cnum,'(g12.4)')  rMuCut
!        call print2col32('Fixed Mu Cut: ', cnum,ioUnit)
!    endif
    
    call print2col32('Anisotropy: ',cAnisotropy,ioUnit) 
    write(cnum,'(i8)') nRegions  
    call print2col32('Number of regions: ',cnum,ioUnit)  

!       
    ! Now loop over the regions and write out a line for each:  
    sFmt = '(i8,1x'
    do j=1,nRhoPerRegion ! rho
        sFmt = trim(sFmt)//',es12.4,1x'
    enddo
    do j=1,nRhoPerRegion ! free
        sFmt = trim(sFmt)//',i8,1x'
    enddo   
    do j=1,2*2*nRhoPerRegion ! bounds, prej
        sFmt = trim(sFmt)//',es12.4,1x'
    enddo
    if ((nRhoPerRegion > 1).and. (nRhoPerRegion < 4) ) then
        do j=1,nRhoPerRegion*(nRhoPerRegion-1) ! ratio prejudice
            sFmt = trim(sFmt)//',es12.4,1x'
        enddo   
    endif
    sFmt = trim(sFmt)//')'   

!     write(*,*) sFmt
     write(ioUnit,'(a)') trim(adjustl(cParamHeader))
     do i = 1,nRegions
    
        if (nRhoPerRegion == 1) then
            write(ioUnit,sFmt)  i, (   rhoParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &  
                    &              (  iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                    &              (   boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                    &              (     prejTemp(i,j),j=1,2*nRhoPerRegion) 
        elseif  (nRhoPerRegion <= 3) then
            write(ioUnit,sFmt)  i, (   rhoParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &  
                &                  (  iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                &                  (   boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                &                  (     prejTemp(i,j),j=1,2*nRhoPerRegion),  &                
                &                  (    ratioTemp(i,j),j=1,nRhoPerRegion*(nRhoPerRegion-1) ) 
                
        elseif  (nRhoPerRegion == 4) then
            write(ioUnit,sFmt)  i, (   rhoParams((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &  
                &                  (  iFreeParam((i-1)*nRhoPerRegion+j),j=1,nRhoPerRegion),  &
                &                  (   boundsTemp(i,j),j=1,2*nRhoPerRegion),  &
                &                  (     prejTemp(i,j),j=1,2*nRhoPerRegion)             
        endif   
                  
     enddo  
      
! Close the file:    
    close(21)
    
    end subroutine writeResistivityFile    
    
!==================================================================================================================================! 
!============================================================================================================================ datime
!==================================================================================================================================!       
    subroutine datime(datetm)
    
    character(80) datetm
    character(10) cdate,ctime
    
    !
    ! Standard intrinsic Fortran90 data_and_time call:
    !
    call date_and_time(cdate,ctime) !this ouputs only values, ignoring other optional arguments
    
    datetm = cdate(5:6)//'/'//cdate(7:8)//'/'//cdate(1:4)//' '// ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:)               
    
    end subroutine datime  
    

  
!==================================================================================================================================! 
!========================================================================================================================== readPoly
!==================================================================================================================================!       
    subroutine readPoly
! 
! Reads in the .poly file used to describe the boundary model.
! October, 2011: retooled to use routines in call_triangle.f90 for most
! file reading, just does some minor adjustments to input model here
! 
 

    use triangle_mesh   
    use mare2dem_global, only: modelFileName,inputmodel
    
    implicit none
 
    character(256)      :: filename
       
    integer             :: idot
!
! Read in the .poly file (and optionally a .nod and .ele file if they exist):
!
 
    write(*,*) '======== Reading in Poly file ============================================'
    write(*,*) ' '     
 
    modelFileName = adjustl(modelFileName)
    idot = index(modelFileName,'.',.true.)
    filename = modelFileName(1:idot-1)

    if ( allocated( inputmodel%attr ) ) call deallocate_trimesh(inputmodel,.false.)
    call read_triangle(filename,inputmodel)
    write(*,*) 'Done reading Poly file'
    write(*,*) ' '
              
    end subroutine readPoly   
    
!==================================================================================================================================! 
!======================================================================================================================= readPenalty
!==================================================================================================================================!        
    subroutine readPenalty
 
    
    integer             :: ierr, i
    character(8)        :: cnum
    character(256)      :: sLine, sCode, sValue 
    logical             :: bComment 
    
    write(*,*) '======== Reading in Penalty file ========================================='
    write(*,*) ' '     

    call print2col32(' Opening Penalty file: ',penaltyFileName,6)
    open(unit=10,file=trim(penaltyFileName),status='old',iostat = ierr)
    
    read( 10, '(A)' ) sLine
    call ParseCode(  sLine, sCode, sValue, bComment )
    select case (trim(sCode))
    
    case ( 'format') ! only one 'format' for now, so no need to check the sValue
        
        call print2col32(' reading new CSR format... ','',6)
                
        read( 10, * ) pen%nnz, pen%nrows
       ! write( *, * ) pen%nnz, pen%nrows
        
        write(cnum,'(i8)') pen%nnz
        call print2col32('Number of penalties: ',cnum,6)
        
        allocate(pen%colind(pen%nnz), pen%val(pen%nnz), pen%rowptr(pen%nrows+1) ) 
         
        do i = 1,pen%nnz
            read( 10, * ) pen%colind(i), pen%val(i)
            !write(*,*) pen%colind(i), pen%val(i)
        enddo
        
        do i = 1,pen%nrows+1
            read( 10, * ) pen%rowptr(i)
            !write(*,*) pen%rowptr(i)
        enddo        
    
    case default ! old style penalty file:  npen \n [i,j,weight]
    
        read (sLine,*) npenalty
         
        write(cnum,'(i8)') npenalty
        call print2col32('Number of penalties: ',cnum,6)
        
        allocate(ijpenalty(npenalty,2),penaltywts(npenalty)) 
        do i = 1,npenalty
            read(10,*) ijpenalty(i,1:2),penaltywts(i)
            !write(*,*) i,ijpenalty(i,1:2),penaltywts(i)
        enddo
   
        
        if (any(ijpenalty > nFree)) then
            write(*,*) ' Error, penalty matrix points to parameter numbers larger than the input free parameters!'
            write(*,*) ' Largest index: ',maxval(ijpenalty),', Number of free paramters: ', nFree
            write(*,*) ' Stopping!'
            call exitMARE2DEM
        endif
        
    end select
    
    close(10)
   
    write(*,*) ' '
    write(*,*) 'Done reading penalty file '
    write(*,*) ' '
    
    end subroutine readPenalty
    
!==================================================================================================================================! 
!====================================================================================================================== readSettings
!==================================================================================================================================!      
    subroutine readSettings  
 
    integer             :: err, iLine = 0
    character(256)      :: sLine, sCode, sValue 
    logical             :: bComment
    
    integer             :: nRxPerGroup 
    
 
    write(6,*) '========== Reading in Parallel and Adaptive Refinement Settings =========='
    write(*,*) ' '     
    call print2col32('Reading settings file:  ',settingsFileName,6)
    write(*,*) ' '
        
    open (unit=10,file=settingsFileName,status='old',iostat=err)
    if (err .ne. 0) then
        write(*,*) ' Error opening settings file'
        call exitMARE2DEM
    end if
    
! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
! Certain sCodes are followed by lists such as model parameters, site locations, etc.
! Parsecode is D Myer's useful code tidebit for separating sCode and sValue pairs
 
    do  ! infinite while loop
      
      ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( 10, '(A)', iostat = err ) sLine
        
        if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
        iLine = iLine + 1   ! advance line counter
        
        call ParseCode(  sLine, sCode, sValue, bComment )
        if( bComment ) cycle
          
        !
        ! What do we have?
        !
        select case (trim(sCode))
       
        case ('version')
        
        case ('scratch folder','scratchfolder')    
            
            scratchFolder = trim(adjustl(sValue))
                                  
        case ('wavenumbers','csem wavenumbers')    
            
            read(sValue,*)  loglower, logupper, nwave 
   
        case ('max # refinement','max # refinements') 
        
            read(sValue,*) maxnadapt
            maxnadapt_default = maxnadapt
        
        case ('tolerance (%)','tolerance')
        
            read(sValue,*) errortolerance        

        case ('mesh quality angle')
        
            read(sValue,*) minQangle
            
            if (minQangle > 33) then
                write(6,fmt='(a32,f6.3)') 'Mesh quality angle:  ',minQangle  
                write(*,*) ' '
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write(*,*) ' !!!! Error: quality angle should be 10 < minQangle < 33 degrees !!!!'
                write(*,*) ' minimum angles > 33 degrees are nearly impossible to mesh ! '
                write(*,*) ' while angles < 10 degrees lead to numerical problems with the FE method '
                write(*,*) ' Sorry buddy, but try specifying a better angle!'
                write(*,*) ' ' 
                write(*,*) ' STOPPING!  '
                write(*,*) ' ' 
                write(*,*) ' ' 
                write(*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                call exitMARE2DEM
            endif           
 
        case ('save meshes')
            select case (trim(sValue))
            case ('yes')
              lSaveMeshFiles = .true.
            end select  
        ! lReuseRefine
        case ('reuse refine')
            select case (trim(sValue))
            case ('yes')
               lReuseRefine   = .true.
               lSaveMeshFiles = .true.
            end select

        case ('number of occam iterations to reuse refine')
            read(sValue,*) nbOccItToReuseRefine
            if (nbOccItToReuseRefine == 0) then
                nbOccItToReuseRefine = 1
            endif

        case ('dual function')

            read(sValue,*) idual_func        

        !
        ! Data decomposition parameters for parallel modeling:
        !
        case ('transmitters per group','csem transmitters per group')   
            read(sValue,*) nTxPerGroupCSEM 

        case ('receivers per group')   !  here for backwards compatibility       
            read(sValue,*) nRxPerGroup 
            nRxPerGroupMT   = nRxPerGroup
            nRxPerGroupCSEM = nRxPerGroup

        case ('mt receivers per group')   
            read(sValue,*) nRxPerGroupMT

        case ('csem receivers per group')    
            read(sValue,*) nRxPerGroupCSEM

        case ('dc transmitters per group')    
            read(sValue,*) nTxPerGroupDC   
        
        case ('dc receivers per group')    
            read(sValue,*) nRxPerGroupDC
            
        case ('csem frequencies per group')         
            read(sValue,*) nFreqPerGroupCSEM 

        case ('wavenumbers per group')         
            read(sValue,*) nKxPerGroup 

        case ('mt frequencies per group')         
            read(sValue,*) nFreqPerGroupMT 

        case ('csem custom groups file')
            read(sValue,*) refGrpsFileName

        case (' (Co)Sine Transform Filters','ct filters','ct filter','ctfilters','ctfilter')

        !           FCTfilter = sValue(1:len(FCTfilter))  ! sorry, no longer supporting this. Filters are fixed to be 601 pts...  

        case ('solver','linear solver','linearsolver')
            linearSolver =  sValue(1:len(linearSolver))  

        case ('use mt scattered field','mt scattered field')  
            select case (trim(sValue))
            case ('yes')
                lMTscatteredField = .true.
            case('no')
                lMTscatteredField = .false.                
            end select  
                
        case ('use bump fields')  
!         select case (trim(sValue))  May 2018: now deprecated
!         case ('yes')
!             lUseBumpFields = .true.
!         case('no')
!             lUseBumpFields = .false.                
!         end select  

        case ('use mesh coarsening', 'use inversion mesh coarsening')  
        
            select case (trim(sValue))
            case ('yes')
                lUseInversionMeshCoarsening = .true.
            case('no')
                lUseInversionMeshCoarsening = .false.                
            end select  


        !
        ! Printing options:
        !
        case ('print nothing')
            select case (trim(sValue))
            case ('yes')
                lprintDebug             = .false.
                lprintDebug_em2dkx      = .false.
                lprintDebug_dc2dkx      = .false.
                lprintSetup             = .false.
                lprintDecomposition     = .false.
                lprintData              = .false.
                lprintGroups            = .false.
                lDisplayRefinementStats = .false.
                lprintTrace_em2dkx      = .false.
                
            end select

        case ('print em2dkx trace')
            select case (trim(sValue))
            case ('yes')
                !lprintDebug             = .false.
                !lprintDebug_em2dkx      = .false.
                !lprintDebug_dc2dkx      = .false.
                !lprintSetup             = .false.
                !lprintDecomposition     = .false.
                !lprintData              = .false.
                !lprintGroups            = .false.
                !lDisplayRefinementStats = .false.
                lprintTrace_em2dkx      = .true.

            end select            

        case ('print data')
            select case (trim(sValue))
            case ('yes')
                lprintData = .true.  
            case ('no')
                lprintData = .false.   
            end select            
        
                
        case ('print debug','printdebug','debug')  
            select case (trim(sValue))
            case ('yes')
                lprintDebug = .true.
                lprintDebug_em2dkx = .true.
                lprintDebug_dc2dkx = .true.
            end select 
        
        case ('print setup','printsetup')  
            select case (trim(sValue))
            case ('yes')
                lprintSetup = .true.
            end select 
        
        case ('print decomposition')  
            select case (trim(sValue))
            case ('yes')
                lprintDecomposition = .true.
            end select    
                                    
        case ('print groups')  
            select case (trim(sValue))
            case ('yes')
                lprintGroups = .true.
            end select    
        
        case ('print adaptive')  
            select case (trim(sValue))
            case ('yes')
                lDisplayRefinementStats = .true.
            case('no')
                lDisplayRefinementStats = .false.
            end select             
        
         
        case ('save task timer file','save tasktimer file')   ! KWK March 2018: this is probably deprecated
            select case (trim(sValue))
            case ('yes')
                lSaveTaskTimers = .true.
            end select   
        
        case ('save loadbalancefile','save load balance file')  ! KWK March 2018: this is probably deprecated
            select case (trim(sValue))
            case ('yes')
                lSaveLoadBalanceTimers = .true.
            end select        
        
        !
        ! Items below are for fine tuning the error estimator and adaptive refinement: advanced users only!
        !    
        case ('minimum error range')  
            read(sValue,*) minRangeProtector  
 
        case ('minimum area')
            read(sValue,*) minArea
             
        case ('max # mesh nodes','maxMeshNodes')
            read(sValue,*) maxMeshNodes 
        
        case ('max # subrefinements')
            read(sValue,*) max_nsubrefine 
                        
        case ('percent refine')
            read(sValue,*) pct_refine 
            
        case ('e noise floor')  
            read(sValue,*) ecutoff 
               
        case ('h noise floor')  
            read(sValue,*) hcutoff
                          
        case ('transmitter quadrature order')
             read(sValue,*) nQuadTxCSEM  
             
        case ('mt quadrature order','mt receiver quadrature order') 
            read(sValue,*) nQuadRxMT  
            
        case ('csem quadrature order','csem receiver quadrature order') 
            read(sValue,*) nQuadRxCSEM  

 
                                 
        case default
            write(*,*) 'Error reading RUNFILE file!'
            write(*,*) ' On line :', iLine            
            write(*,*) 'Unknown or unsupported code:'
            write(*,*) sCode

            call exitMARE2DEM
        
        end select 
        
    enddo ! read while loop
!
! Close settings file:
!
    close(10)
    
    write(*,*) 'Done reading settings file '
    write(*,*) ' '
 
    end subroutine readSettings

!==================================================================================================================================!
!===================================================================================================================== readRefGroups
!==================================================================================================================================!
    subroutine readRefGroups

    integer             :: err, iAllocErr, rxStartIdx = 1, iLine = 0
    character(256)      :: sLine, sCode, sValue
    logical             :: bComment

    write(6,*) '========== Reading in Custom Refinement Group Settings =========='
    write(*,*) ' '
    call print2col32('Reading ref. Groups file:  ',refGrpsFileName,6)
    write(*,*) ' '

    open (unit=90,file=refGrpsFileName,status='old',iostat=err)
    if (err .ne. 0) then
        write(*,*) ' Error opening ref. groups settings file',refGrpsFileName, ')'
        call exitMARE2DEM
    end if

    ! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
    ! Certain sCodes are followed by lists such as model parameters, site locations, etc.
    ! Parsecode is D Myer's useful code tidebit for separating sCode and sValue pairs

    do  ! infinite while loop

      ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are
        !   stripped from the value portion.

        read( 90, '(A)', iostat = err ) sLine

        if (err /= 0) exit  ! end of file read, escape from while loop and
                            ! proceed to checking that required inputs defined

        iLine = iLine + 1   ! advance line counter

        call ParseCode(  sLine, sCode, sValue, bComment )
        if( bComment ) cycle

        !
        ! What do we have?
        !
        select case (trim(sCode))

          case ('custom number of receiver groups')
            read(sValue,*) customNRrxGroups
            allocate (customNbRxPerGroup(customNRrxGroups), stat=iAllocErr)  ! Allocate ref groups start array
            if (iAllocErr .ne. 0) then
                write(*,*) ' Out of memory.  Too many custom refinement groups (', customNRrxGroups, ')'
                call exitMARE2DEM
            endif

          case ('group number of recceivers')
            if (customNRrxGroups == 0) then
                write(*,*) ' Error: custom number of receiver groups not set before setting receiver group start indices'
                call exitMARE2DEM
            endif
            read(sValue,*) customNbRxPerGroup(rxStartIdx)
            if (lprintDebug) then
                write(6,fmt='(a32,a3)') ' customNbRxPerGroup (', customNbRxPerGroup(rxStartIdx), ')'
            endif
            rxStartIdx = rxStartIdx + 1

          case default
            write(*,*) 'Error reading RUNFILE file!'
            write(*,*) ' On line :', iLine
            write(*,*) 'Unknown or unsupported code:'
            write(*,*) sCode

            call exitMARE2DEM

        end select

    enddo ! read while loop

    write(*,*) 'Done reading custom ref. groups settings file '
    write(*,*) ' '

    end subroutine readRefGroups
    
!==================================================================================================================================! 
!========================================================================================================================= readModel
!==================================================================================================================================!       
    subroutine readModel
!
! Subroutine to read in the MARE2DEM model file and its dependencies:
!

! Step 1: read the .resistivity file: 
    call readResistivityFile
 
! Step 2: read in the Model Boundary .poly file:
    call readPoly
     
! Step 3: read in the Model Penalty .penalty file:
    call readPenalty
    
! Step 4: read the MARE2DEM settings file:
    call readSettings

! Step 5: read the custom refinement groups file:
    if (refGrpsFileName /= '') then
        call readRefGroups
    endif
! 
!  Check boundary model for slivers
!
    call checkModel
           
    end subroutine readModel 
    
!==================================================================================================================================! 
!========================================================================================================================== readData
!==================================================================================================================================!  
    subroutine readData
!
! Reads in the EMData_2.0 format data file used in the MARE2DEM inversion
!
 
    implicit none

    include 'em_parameters.inc'
          
!
! Local variables:
! 
    integer         :: i, err, iTxRead, iRxRead, iFreqRead, iDataRead, nAllocErr
    character(180)  :: sLine, sCode, sValue ! These are for reading the lines of the file
    character(180)  :: sFields(9),  sFmt
    logical         :: bComment, lerror 
    character(32)   :: cNum
    character(180)  :: sDataFormat 
    
    integer, parameter :: iof = 15   ! I/O File identifier  
    

    if (lprintData) then
        write(*,*) ' '    
        write(6,*) '=============== Reading in Data File ====================================='
        write(*,*) ' '    
    endif
!
! Open the data file:
!
    open (iof, file=trim(dataFileName), status='old', iostat=err)
    if (err /= 0) then
        write(*,*) ' Error opening data file:', trim(dataFileName)
        call exitMARE2DEM
    else
        if (lprintData) call print2col32('Reading data from file: ', dataFileName,6)
    end if  
    
!
! Loop through the model file looking for the Format and Number of Layers flags
! and skipping over any comments
!
    do while (.true.) 
    
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( iof, '(A180)', iostat = err ) sLine
       
        if (err /= 0) exit  ! end of file read, escape from while(.true.) loop and 
                            ! proceed to checking that required inputs are defined
        call ParseCode( sLine, sCode, sValue, bComment )

        if( bComment ) cycle
       
!
! What do we have?
!
        select case (trim(sCode))
       
        case ('format')         ! Set cnDataParams for the input data format
        
            call lower(sValue) 
            sDataFormat = trim(sValue)
            
            select case (trim(sDataFormat) )
            
            case ('emdata_2.0','emdata_2.1','emdata_2.2')  
                cnDataParams = 4    
                
 
            case default
            
                write(*,*) ' Error: data format unsupported: ', trim(sDataFormat)
                write(*,*) ' Try using EMDATA_2.2'
                close(iof)
                call exitMARE2DEM
            end select
            
            if (lprintData) call print2col32('Format: ',trim(sValue),6) 
     
            
        case ('utm of x,y origin (utm zone, n, e, 2d strike)','utm','origin')
            ! These aren't used in MARE2DEM, but save the line so we can include it in the .resp files:
            cUTMLine = trim(adjustl(sValue))
            if (lprintData) call print2col32('UTM of origin (zone,N,E,theta): ',trim(sValue),6) 
            
        case ('phase','phase convention')  
            call lower(sValue)
            
            select case (trim(sValue))
            
            case ('lag')
                    phaseConvention = 'lag'
            case ('lead')
                    phaseConvention = 'lead'
            case default
                write(*,*) ' Error: unrecognized phase convenction: ', trim(sValue)
                write(*,*) ' Try using lag or lead'
                close(iof)
                call exitMARE2DEM          
            end select
            if (lprintData) call print2col32('Phase convention: ', phaseConvention,6)

        case ('reciprocity','reciprocity used')  
            call lower(sValue)
          
            select case (trim(sValue))
            case ('yes')
                   reciprocityUsed = 'yes'
            case ('no')
                   reciprocityUsed = 'no'
            case default
            
                if (len_trim(sValue) > 0 ) then ! only error if a non-blank value is given
                    write(*,*) ' Error: unrecognized reciprocity used value: ', trim(sValue)
                    write(*,*) ' Try using yes or no'
                    close(iof)
                    call exitMARE2DEM     
                endif     
            end select
            if (lprintData) call print2col32('Reciprocity used: ', reciprocityUsed,6)      
                        
                               
        case ('# transmitters','# csem transmitters')
        
            read(sValue,*) nTxCSEM
            
            if (lprintData) write(*,*) ' '
            write(cnum,'(i3)') nTxCSEM
            if (lprintData) call print2col32('# Transmitters: ',cnum,6)
            !write(*,*) ''    
            !write(*,*) '  Transmitter Positions: '

            sFmt = '(a1,1x,7(a12,2x),a12)' 
            if (lprintData) write(*,sFmt) ' ','X','Y','Z','Azimuth','Dip','Length','Type','Name'

            
             
            !
            !  Allocate arrays dependent on nTx
            !
            allocate ( azimuthTxCSEM(nTxCSEM),dipTxCSEM(nTxCSEM), &
                     & xTxCSEM(nTxCSEM),yTxCSEM(nTxCSEM),zTxCSEM(nTxCSEM), lengthTxCSEM(nTxCSEM), &
                     & cSourceType(nTxCSEM), cTxNamesCSEM(nTxCSEM),   stat=nAllocErr )   
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many transmitters (', nTxCSEM, ')'
                call exitMARE2DEM 
            endif        
            cTxNamesCSEM = ' '
            lengthTxCSEM = 0
                
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            
            iTxRead = 0
            do while (iTxRead < nTxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iTxRead  = iTxRead + 1    
                
                select case (trim(sDataFormat) )
                
                case ('emdata_2.0')  
                    
                    read(sLine,*) xTxCSEM(iTxRead),yTxCSEM(iTxRead),zTxCSEM(iTxRead), &
                                & azimuthTxCSEM(iTxRead),dipTxCSEM(iTxRead),cSourceType(iTxRead)                     
 
                    
                case ('emdata_2.1') 
                
                    call parseFields( sLine, 7, sFields)
                     
                    read(sFields(1),*) xTxCSEM(iTxRead)
                    read(sFields(2),*) yTxCSEM(iTxRead) 
                    read(sFields(3),*) zTxCSEM(iTxRead) 
                    
                    read(sFields(4),*) azimuthTxCSEM(iTxRead)
                    read(sFields(5),*) dipTxCSEM(iTxRead)
                    read(sFields(6),*) cSourceType(iTxRead)

                    if (len_trim( sFields(7) ) > 0 ) then
                        read(sFields(7),*) cTxNamesCSEM(iTxRead)    
                    else
                        cTxNamesCSEM(iTxRead)     = ' '
                    endif
 
 
                case ('emdata_2.2') 
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xTxCSEM(iTxRead)
                    read(sFields(2),*) yTxCSEM(iTxRead) 
                    read(sFields(3),*) zTxCSEM(iTxRead) 
                    
                    read(sFields(4),*) azimuthTxCSEM(iTxRead)
                    read(sFields(5),*) dipTxCSEM(iTxRead)
                    read(sFields(6),*) lengthTxCSEM(iTxRead)
                    read(sFields(7),*) cSourceType(iTxRead)

                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cTxNamesCSEM(iTxRead)    
                    else
                        cTxNamesCSEM(iTxRead)     = ' '
                    endif
                

                end select                
                              
                 
                sFmt = '(2x,6(f12.1,2x),a12,2x,a12)'            

                if (lprintData) then
                    if ( iTxRead < 101 ) then
                        write(*,sFmt) xTxCSEM(iTxRead),yTxCSEM(iTxRead),zTxCSEM(iTxRead), &
                                & azimuthTxCSEM(iTxRead),dipTxCSEM(iTxRead), lengthTxCSEM(iTxRead), &
                                & trim(adjustl(cSourceType(iTxRead))), trim(adjustl(cTxNamesCSEM(iTxRead)))
                    elseif( iTxRead == 101 ) then
                        write(*,*) 'More than 101 CSEM transmitters. Not writing the rest.'  
                    endif                
                endif
                                                         
            enddo
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iTxRead/=nTxCSEM) then
                write(*,*) ' Error reading transmitters: iTxRead/=nTxCSEM:',iTxRead, nTxCSEM
                close(iof)
                call exitMARE2DEM            
            
            endif       
            
        
                
        case ('# csem frequencies')
            
            read(sValue,*) nFreqCSEM
            write(cnum,'(i3)') nFreqCSEM
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32('# CSEM Frequencies: ',cnum,6) 
            allocate (fTxCSEM(nFreqCSEM), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many CSEM frequencies (', nFreqCSEM, ')'
                call exitMARE2DEM 
            endif   
            
            !
            ! Now read in block of frequencies, skipping any comment lines:
            !
            iFreqRead = 0
            do while (iFreqRead < nFreqCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iFreqRead  = iFreqRead + 1                  
                read(sLine,*)  fTxCSEM(iFreqRead)
                if (lprintData) write(*,'(g12.3)') fTxCSEM(iFreqRead)
               
        
            enddo
            !
            ! Check and make sure all the frequencies were read in:
            !
            if (iFreqRead/=nFreqCSEM) then
                write(*,*) ' Error reading CSEM frequencies: iFreqRead/=nFreqCSEM:',iFreqRead, nFreqCSEM
                close(iof)
                call exitMARE2DEM                        
            endif          
             
        case ('# csem receivers')
         
            read(sValue,*) nRxCSEM
            write(cnum,'(i6)') nRxCSEM
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32(' # CSEM Receivers: ',cnum,6) 
    
            if (lprintData) write(*,'(a1,1x,7(a12,2x),a12)') ' ','X','Y', 'Z','Theta','Alpha','Beta','Length', 'Name' 
               
        
            !
            ! Allocate the field arrays: 
            !
            allocate(xRxCSEM(nRxCSEM), yRxCSEM(nRxCSEM), zRxCSEM(nRxCSEM), lengthRxCSEM(nRxCSEM),  &
                      ThetaRxCSEM(nRxCSEM),AlphaRxCSEM(nRxCSEM),BetaRxCSEM(nRxCSEM), cRxNamesCSEM(nRxCSEM),  &
                      stat=nAllocErr )   ! rxnames go here
                   
            cRxNamesCSEM = ' '
            lengthRxCSEM = 0
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many CSEM receivers (', nRxCSEM, ')'
                call exitMARE2DEM 
            endif         

            
 
         
            !
            ! Now read in block of receiver locations, skipping any comment lines:
            !
            iRxRead = 0
 
            do while (iRxRead < nRxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                select case (trim(sDataFormat) )
            
                case ('emdata_2.0')  
                
     
                    call parseFields( sLine, 6, sFields)
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)

  
                case ('emdata_2.1')  
                
                    call parseFields( sLine, 7, sFields)
                     
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)

                    if (len_trim( sFields(7) ) > 0 ) then
                        read(sFields(7),*) cRxNamesCSEM(iRxRead)
                    else
                        cRxNamesCSEM(iRxRead) = ' '
                    endif
  
                  case ('emdata_2.2')  
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xRxCSEM(iRxRead)
                    read(sFields(2),*) yRxCSEM(iRxRead)
                    read(sFields(3),*) zRxCSEM(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxCSEM(iRxRead)
                    read(sFields(5),*) AlphaRxCSEM(iRxRead)
                    read(sFields(6),*) BetaRxCSEM(iRxRead)
                    read(sFields(7),*) lengthRxCSEM(iRxRead) 
                    
                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cRxNamesCSEM(iRxRead)
                    else
                        cRxNamesCSEM(iRxRead) = ' '
                    endif
                      
                end select
                if (lprintData) then
                    if ( iRxRead < 101 ) then
                        write(*,'(2x,7(f12.1,2x),a12)')  xRxCSEM(iRxRead),yRxCSEM(iRxRead),zRxCSEM(iRxRead), & 
                             & ThetaRxCSEM(iRxRead),AlphaRxCSEM(iRxRead),BetaRxCSEM(iRxRead),lengthRxCSEM(iRxRead), & 
                             & trim(adjustl(cRxNamesCSEM(iRxRead)))
                    elseif( iRxRead == 101 ) then
                        write(*,*) 'More than 101 CSEM receivers. Not writing the rest.'  
                    endif
                endif
                                          
            enddo
            !
            ! Check and make sure all the rx were read in:
            !
            if (iRxRead/=nRxCSEM) then
                write(*,*) ' Error reading receivers: iRxRead/=nRxCSEM:',iRxRead,nRxCSEM
                close(iof)
                call exitMARE2DEM                        
            endif       
            
        case ('csem receiver names') 
           
            cRxNamesCSEM = ''   
            iRxRead = 0
            write(*,*) ' '
            if (lprintData) call print2col32('CSEM Receiver Names: ', ' ',6)
            do while (iRxRead < nRxCSEM) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                 cRxNamesCSEM(iRxRead) = trim(adjustl(sLine))
                 if (lprintData) write(*,'(1x,a)') trim(cRxNamesCSEM(iRxRead) )
  
                
            enddo  

!
! DC Resistivity inputs:
!   
   
        case ('# dc transmitters')
        
            read(sValue,*) nTxDC
            
            if (lprintData) write(*,*) ' '
            write(cnum,'(i6)') nTxDC
            if (lprintData) call print2col32('# DC Transmitters: ',cnum,6)

            sFmt = '(a1,1x,3(a12,2x))' 
            if (lprintData) write(*,sFmt) ' ','Electrode A','Electrode B','Name'

            allocate ( TxDC(nTxDC,2),cTxNamesDC(nTxDC),   stat=nAllocErr )   
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many DC transmitters (', nTxDC, ')'
                call exitMARE2DEM 
            endif        
            cTxNamesDC = ' '
                
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            
            iTxRead = 0
            do while (iTxRead < nTxDC) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iTxRead  = iTxRead + 1    
                
!                 select case (trim(sDataFormat) )
!  
!                 case ('emdata_2.2') 
                
                    call parseFields( sLine, 3, sFields)
                     
                    read(sFields(1),*) TxDC(iTxRead,1)
                    read(sFields(2),*) TxDC(iTxRead,2) 
 
                    if (len_trim( sFields(3) ) > 0 ) then
                        read(sFields(3),*) cTxNamesDC(iTxRead)    
                    else
                        cTxNamesDC(iTxRead)     = ' '
                    endif
                         
                ! end select                
                if (lprintData) then               
                sFmt = '(2x,2(i12,2x),a12)'            
            
                if ( iTxRead < 101 ) then
                    write(*,sFmt) TxDC(iTxRead,1:2), trim(adjustl(cTxNamesDC(iTxRead)))
                elseif( iTxRead == 101 ) then
                    write(*,*) 'More than 101 DC transmitters. Not writing the rest.'  
                endif
                
                endif
                                                         
            enddo
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iTxRead/=nTxDC) then
                write(*,*) ' Error reading DC transmitters: iTxRead/=nTxDC:',iTxRead, nTxDC
                close(iof)
                call exitMARE2DEM            
            
            endif       
            
        case ('# dc receivers')
         
            read(sValue,*) nRxDC
            write(cnum,'(i6)') nRxDC
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32(' # DC Receivers: ',cnum,6) 
            
            sFmt = '(a1,1x,3(a12,2x))' 
            if (lprintData) write(*,sFmt) ' ','Electrode M','Electrode N','Name'
            
            allocate ( RxDC(nRxDC,2),cRxNamesDC(nRxDC),   stat=nAllocErr )   
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many DC receivers (', nRxDC, ')'
                call exitMARE2DEM 
            endif        
            cRxNamesDC = ' '
                
            !
            ! Now read in block of DC receivers, skipping any comment lines:
            !
            
            iRxRead = 0
            do while (iRxRead < nRxDC) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1    
                
!                 select case (trim(sDataFormat) )
!  
!                 case ('emdata_2.2') 
                
                    call parseFields( sLine, 3, sFields)
                     
                    read(sFields(1),*) RxDC(iRxRead,1)
                    read(sFields(2),*) RxDC(iRxRead,2) 
 
                    if (len_trim( sFields(3) ) > 0 ) then
                        read(sFields(3),*) cRxNamesDC(iRxRead)    
                    else
                        cRxNamesDC(iRxRead)     = ' '
                    endif
                         
                ! end select                
                
                if (lprintData) then                   
                sFmt = '(2x,2(i12,2x),a12)'            
                    if ( iRxRead < 101 ) then
                        write(*,sFmt) RxDC(iRxRead,1:2), trim(adjustl(cRxNamesDC(iRxRead)))
                    elseif( iRxRead == 101 ) then
                        write(*,*) 'More than 101 DC receivers. Not writing the rest.'  
                    endif
                endif
                                                         
            enddo
            
            !
            ! Check and make sure all the receivers were read in:
            !
            if (iRxRead/=nRxDC) then
                write(*,*) ' Error reading DC receivers: iRxRead/=nRxDC:',iRxRead, nRxDC
                close(iof)
                call exitMARE2DEM            
            
            endif   
            
        case ('# dc transmitter electrodes')              
            
            read(sValue,*)   nTrodesTxDC
            write(cnum,'(i6)') nTrodesTxDC
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32(' # DC Transmitter Electrodes: ',cnum,6) 

            if (lprintData) write(*,'(a1,1x,3(a12,2x))') ' ','X','Y', 'Z' 
             
            allocate ( trodes_TxDC(nTrodesTxDC,3),  stat=nAllocErr )   
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many DC transmitter electrodes (', nTrodesTxDC, ')'
                call exitMARE2DEM 
            endif        
           
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            
            iTxRead = 0
            do while (iTxRead < nTrodesTxDC) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the parameters from sLine
                iTxRead  = iTxRead + 1    
 
                call parseFields( sLine, 3, sFields)
 
                read(sFields(1),*) trodes_TxDC(iTxRead,1)
                read(sFields(2),*) trodes_TxDC(iTxRead,2) 
                read(sFields(3),*) trodes_TxDC(iTxRead,3) 
  
                if (lprintData) then
                    sFmt = '(2x,3(f12.1,2x))'            
    
                    if ( iTxRead < 101 ) then
                        write(*,sFmt) trodes_TxDC(iTxRead,1:3) 
                    elseif( iTxRead == 101 ) then
                        write(*,*) 'More than 101 DC transmitter electrodes. Not writing the rest.'  
                    endif
                endif
                                                         
            enddo
            
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iTxRead/=nTrodesTxDC) then
                write(*,*) ' Error reading DC transmitters electrodes: iTxRead/=nTrodesTxDC:',iTxRead, nTrodesTxDC
                close(iof)
                call exitMARE2DEM            
            endif   
                        
        case ('# dc receiver electrodes')                     
            
            read(sValue,*)   nTrodesRxDC
            write(cnum,'(i6)')   nTrodesRxDC
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32(' # DC Receiver Electrodes: ',cnum,6) 
            
            write(*,'(a1,1x,3(a12,2x))') ' ','X','Y', 'Z' 
            
            allocate ( trodes_RxDC(nTrodesRxDC,3),  stat=nAllocErr )   
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many DC receiver electrode (', nTrodesRxDC, ')'
                call exitMARE2DEM 
            endif        
           
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            
            iRxRead = 0
            do while (iRxRead < nTrodesRxDC) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the parameters from sLine
                iRxRead  = iRxRead + 1    
 
                call parseFields( sLine, 3, sFields)
 
                read(sFields(1),*) trodes_RxDC(iRxRead,1)
                read(sFields(2),*) trodes_RxDC(iRxRead,2) 
                read(sFields(3),*) trodes_RxDC(iRxRead,3) 
                
                if (lprintData) then
                    sFmt = '(2x,3(f12.1,2x))'            

                    if ( iRxRead < 101 ) then
                        write(*,sFmt) trodes_RxDC(iRxRead,1:3) 
                    elseif( iRxRead == 101 ) then
                        write(*,*) 'More than 101 DC receiver electrodes. Not writing the rest.'  
                    endif
                endif                                        
            enddo
            
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iRxRead/=nTrodesRxDC) then
                write(*,*) ' Error reading DC receiver electrodes: iRxRead/=nTrodesRxDC:',iRxRead, nTrodesRxDC
                close(iof)
                call exitMARE2DEM            
            
            endif               
               
!----------------
! MT inputs:
!----------------                  
        case ('# mt frequencies')
            
            read(sValue,*) nFreqMT
            write(cnum,'(i3)') nFreqMT
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32('# MT Frequencies: ',cnum,6) 
            allocate (fTxMT(nFreqMT), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many MT frequencies (', fTxMT, ')'
                call exitMARE2DEM 
            endif   
            
            !
            ! Now read in block of frequencies, skipping any comment lines:
            !
            iFreqRead = 0
            do while (iFreqRead < nFreqMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iFreqRead  = iFreqRead + 1                  
                read(sLine,*)  fTxMT(iFreqRead)
                if (lprintData) write(*,'(g12.3)') fTxMT(iFreqRead)
               
        
            enddo
            !
            ! Check and make sure all the frequencies were read in:
            !
            if (iFreqRead/=nFreqMT) then
                write(*,*) ' Error reading MT frequencies: iFreqRead/=nFreqMT:',iFreqRead, nFreqMT
                close(iof)
                call exitMARE2DEM                        
            endif          
             
      
  
        case ('# mt receivers')
         
            read(sValue,*) nRxMT
            write(cnum,'(i6)') nRxMT
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32('# MT Receivers: ',cnum,6) 
            
 
            if (lprintData) write(*,'(a1,1x,9(a12,2x))') ' ','x','y', 'z','Theta','Alpha','Beta', 'Length', 'SolveStatic', 'Name'

                
            !
            ! Allocate the field arrays: 
            !
            allocate( xRxMT(nRxMT), yRxMT(nRxMT), zRxMT(nRxMT),ThetaRxMT(nRxMT), AlphaRxMT(nRxMT),BetaRxMT(nRxMT), &
                    &  cRxNamesMT(nRxMT), iSolveStatic(nRxMT),  lengthRxMT(nRxMT), stat=nAllocErr )
            iSolveStatic = 0 
            
            cRxNamesMT = ' ' 
            lengthRxMT = 0
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many MT receivers (', nRxMT, ')'
                call exitMARE2DEM 
            endif         
 
!         
            !
            ! Now read in block of receiver locations, skipping any comment lines:
            !
            iRxRead = 0
 
            do while (iRxRead < nRxMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1        
         
                 
                select case (trim(sDataFormat) )
            
                case ('emdata_2.0')  
                
                    call parseFields( sLine, 6, sFields)
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)         
                               
       
           
                case ('emdata_2.1')  
                
                    call parseFields( sLine, 8, sFields)
                     
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)
                     
                    read(sFields(7),*) iSolveStatic(iRxRead)
                     
                    if (len_trim( sFields(8) ) > 0 ) then
                        read(sFields(8),*) cRxNamesMT(iRxRead)
                    else
                        cRxNamesMT(iRxRead) = ' '
                    endif

                case ('emdata_2.2')  
                
                    call parseFields( sLine, 9, sFields)
                     
                    read(sFields(1),*) xRxMT(iRxRead)
                    read(sFields(2),*) yRxMT(iRxRead)
                    read(sFields(3),*) zRxMT(iRxRead)  
                    
                    read(sFields(4),*) ThetaRxMT(iRxRead)
                    read(sFields(5),*) AlphaRxMT(iRxRead)
                    read(sFields(6),*) BetaRxMT(iRxRead)
                    read(sFields(7),*) lengthRxMT(iRxRead)
                    read(sFields(8),*) iSolveStatic(iRxRead)
                     
                    if (len_trim( sFields(9) ) > 0 ) then
                        read(sFields(9),*) cRxNamesMT(iRxRead)
                    else
                        cRxNamesMT(iRxRead) = ' '
                    endif
               
                end select
                
                !
                ! Check rotations:
                !
                if ( (ThetaRxMT(iRxRead) /= 0) .or. (AlphaRxMT(iRxRead) /= 0) ) then 
                 	write(*,*) ' !! Non-zero MT receiver theta or alpha rotation angles detected for receiver ', iRxRead
                 	write(*,*) ' !! Sorry, but MARE2DEM currently only allows the beta (y-z plane) angle to be non-zero!'
                 	write(*,*) ' !! Please correct your data file. Stopping'
                	call exitMARE2DEM 
            	endif
            	if (lprintData) then
                if ( iRxRead < 101 ) then
                        write(*,'(2x,7(f12.1,2x),i12,2x,a12)')  xRxMT(iRxRead),yRxMT(iRxRead),zRxMT(iRxRead), & 
                             & ThetaRxMT(iRxRead),AlphaRxMT(iRxRead),BetaRxMT(iRxRead), lengthRxMT(iRxRead),iSolveStatic(iRxRead), &
                                trim(adjustl(cRxNamesMT(iRxRead)))  
                elseif( iRxRead == 101 ) then
                    write(*,*) 'More than 101 MT receivers. Not writing the rest.'  
                endif
                endif
                                          
              
            enddo
            !
            ! Check and make sure all the rx were read in:
            !
            if (iRxRead/=nRxMT) then
                write(*,*) ' Error reading MT receivers: iRxRead/=nRxMT:',iRxRead,nRxMT
                close(iof)
                call exitMARE2DEM                        
            endif                    
            
        case ('mt receiver names')             
               
            cRxNamesMT = ' '   
            iRxRead = 0
            if (lprintData) write(*,*) ' '
            if (lprintData) call print2col32('MT Receiver Names: ', ' ',6)
            do while (iRxRead < nRxMT) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                        
     
                cRxNamesMT(iRxRead) = trim(adjustl(sLine))
                if (lprintData) write(*,'(1x,a)') trim(cRxNamesMT(iRxRead)) 
 
         enddo  
   
    !
    ! Lastly we need to read in the DATA block.
    !
         case ('# data','#data','ndata')
         
            read(sValue,*) nd
            if (lprintData) write(*,*) ' ' 
            write(cnum,'(i6)') nd
            if (lprintData) call print2col32('Number of Data: ',cnum,6)
                      
            !
            ! Allocate the nd data arrays: 
            !
            allocate(dp(nd,cnDataParams), d(nd), d_wt(nd), sd(nd), dm(nd), stat=nAllocErr ) 
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many data points (', nd, ')'
                call exitMARE2DEM 
            endif        
            d_wt = 1. ! standard data weight, modified later on if joint MT + CSEM inversion
             
            iDataRead = 0
            if (lprintData) write(*,'(4(a12,1x),a12,2x,a12,2x)') 'Type','Time/Freq #','Tx #','Rx #', 'Data','Std Error'    
            
            do while (iDataRead < nd) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iDataRead  = iDataRead + 1     
                
                ! Read in 6 fields as character strings...
                call parseFields( sLine, 6, sFields)    
                
                ! Get the fields:
                read(sFields(1),*) dp(iDataRead,1)  
                read(sFields(2),*) dp(iDataRead,2)
                read(sFields(3),*) dp(iDataRead,3)
                read(sFields(4),*) dp(iDataRead,4)
                read(sFields(5),*) d(iDataRead)
                read(sFields(6),*) sd(iDataRead)

                
                ! DGM 10/2010 Don't spew everything to the screen - only need to
                ! check the first so many.
                if (lprintData) then
                    if( iDataRead < 101 ) then
                        write(*,fmt= '(4(i12,1x),ES12.3,2x,ES12.3,2x)') (dp(iDataRead,i),i=1,4),d(iDataRead),sd(iDataRead)
                    elseif( iDataRead == 101 ) then
                        write(*,*) 'More than 101 data. Not writing the rest.'  
                    endif
                endif
            enddo    
            if (lprintData) write(*,*) ' '    
            !
            ! Check and make sure all the data was read in:
            !
            if (iDataRead/=nd) then
                write(*,*) ' Error reading data: iDataRead/=nd:',iDataRead,nd
                close(iof)
                call exitMARE2DEM                        
            endif      
            !

            
            
        case default
                write(*,*) 'Error, unknown code in Data file: ',trim(sCode)
                write(*,*) 'skipping this, hopefully its not crucial!'
                !call exitMARE2DEM
            
        end select ! case (trim(sCode))
        
        
    enddo ! while (.true.) ! The main data reading loop
!
! Close the data file
!
    close(iof)


!
! Check that parameters have been defined:
!
    lerror = .false.
    if ((nTxCSEM < 1).and.(nRxMT<1).and.(nTxDC<1)) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of transmitters: nTxCSEM,nTxDC = ',nTxCSEM,nTxDC
    endif
    if ( (nRxCSEM < 1) .and. (nRxMT<1) .and.(nRxDC<1)) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of receivers: nRxCSEM, nRxMT,nRxDC = ',nRxCSEM, nRxMT,nRxDC
    endif
    if ((nFreqCSEM < 1) .and. (nFreqMT < 1).and.(nRxDC<1) ) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of frequencies: nFreqCSEM,nFreqMT = ',nFreqCSEM,nFreqMT
    endif
    if ( (nd < 1).and.(.not.lFwdFields) ) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of data: # Data = ',nd
    endif   
    
    if (lerror) then
        write(*,*)  ' Stopping ! '
        call exitMARE2DEM
    endif
 
!
! Count number of CSEM and MT data:
! 
    nMT   = 0
    nCSEM = 0 
    nDC   = 0 
    
    do i = 1,nd
        
        if (dp(i,3) <= 0) dp(i,3)  = dp(i,4)  ! if iTx not specified, make it the current receiver 
            
        select case (  dp(i,1) ) ! what data type is this?
                         
                         
            ! CSEM types:                                  
            case (indRealEx,indImagEx,indAmpEx,indPhsEx,indLog10AmpEx)   
                nCSEM = nCSEM + 1 
            case (indRealEy,indImagEy,indAmpEy,indPhsEy,indLog10AmpEy)  
                nCSEM = nCSEM + 1 
            case (indRealEz,indImagEz,indAmpEz,indPhsEz,indLog10AmpEz)  
                nCSEM = nCSEM + 1 
            case (indRealBx,indImagBx,indAmpBx,indPhsBx,indLog10AmpBx)  
                nCSEM = nCSEM + 1 
            case (indRealBy,indImagBy,indAmpBy,indPhsBy,indLog10AmpBy)  
                nCSEM = nCSEM + 1 
            case (indRealBz,indImagBz,indAmpBz,indPhsBz,indLog10AmpBz)                      
                nCSEM = nCSEM + 1    
            case (iPEmax,iPEmin)
                nCSEM = nCSEM + 1 
            case (iPBmax,iPBmin)
                nCSEM = nCSEM + 1 
                
             ! MT types:   
            case (indRhoZXY,indPhsZXY,indRealZXY,indImagZXY,indlog10RhoZXY, indRealMZY, indImagMZY,indAmpMZY,indPhsMZY )                 
                nMT = nMT + 1
            case (indRhoZYX,indPhsZYX,indRealZYX,indImagZYX,indlog10RhoZYX)    
                nMT = nMT + 1 
            case (indRealExMT,indImagExMT,indRealEyMT,indImagEyMT,indRealEzMT,indImagEzMT)    
                nMT = nMT + 1 
            case (indRealHxMT,indImagHxMT,indRealHyMT,indImagHyMT,indRealHzMT,indImagHzMT)    
                nMT = nMT + 1         
            case (indRhoZXX,indPhsZXX,indRealZXX,indImagZXX,indlog10RhoZXX,indRhoZYY,indPhsZYY,indRealZYY,indImagZYY,indlog10RhoZYY)
                write(*,*) ' Error in input data file, diagonal impedance elements not support in 2D. Stopping!'
                call exitMARE2DEM
                 
            ! DC resistivity:     
            case ( indAppRes_DC,indlog10AppRes_DC)
                nDC = nDC + 1 
                
            case default
                write(*,*) ' MARE2DEM: error in input data file, unknown data type!: ' ,  dp(i,1) 
                write(*,*) ' stopping! '  
                call exitMARE2DEM
                
        end select  ! case dp(i,1)
    enddo  ! loop over nd

!
! If lFwdFields then user wants all field components computed regardless of input data array:
!
    if (lFwdFields) then
        if ( nFreqCSEM > 0 ) then
            nCSEM       = nFreqCSEM*nRxCSEM*nTxCSEM*6
        endif
      
        if ( nFreqMT > 0 ) then
            nMT         = nFreqMT*nRxMT*6
        endif
        
        ! KWK debug February 27, 2019: lFwdFields not currently supported, but when resupported should reallocate data array here
        ! for all Tx,Rx,freq, component combinations...
    endif

    if (lprintData) write(*,*) 'Done reading the data file, here is what I found:'
    if (lprintData) write(*,*)  ' '
     
 
    write(cnum,'(i6)') nCSEM
    if (lprintData) call print2col32(' # CSEM Data: ',cnum,6)
    write(cnum,'(i6)') nMT
    if (lprintData) call print2col32(' # MT Data: ',cnum,6)
    write(cnum,'(i6)') nDC
    if (lprintData) call print2col32(' # DC Resistivity Data: ',cnum,6) 

!
! Allocate associated parameters that are carried through Occam for the best fitting model. Currently these are only used
! for the MT static shift parameter estimates:
!
    if (nRxMT > 0 ) then
        npm_assoc = nRxMT*2
        allocate( pm_assoc(npm_assoc) )
    endif 
    
!
! Finally, run a consistency check on the data parameter indices
!
! 'Type','Freq #','Tx #','Rx #', 'Data','Std Error' 
    do i = 1,nd
        select case (  dp(i,1) ) ! what data type is this?
                                          
        case (1:99)
            if (dp(i,2) > nFreqCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Frequency index exceeds # of CSEM frequencies'
                call exitMARE2DEM
            elseif (dp(i,3) > nTxCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Transmitter index exceeds # of CSEM transmitters'
                call exitMARE2DEM
            elseif (dp(i,4) > nRxCSEM) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Receiver index exceeds # of CSEM receivers'
                call exitMARE2DEM            
            endif
            
        case (101:199 )  
                 
            if (dp(i,2) > nFreqMT) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Frequency index exceeds # of MT frequencies'
                call exitMARE2DEM
           
            elseif (dp(i,4) > nRxMT) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Receiver index exceeds # of MT receivers'
                call exitMARE2DEM            
            endif
            
        case(201:299)
        
            if (dp(i,3) > nTxDC) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Transmitter index exceeds # of DC transmitters'
                call exitMARE2DEM
            elseif (dp(i,4) > nRxDC) then
                write(*,*) 'Error in data file on data matrix line: ',i
                write(*,*) 'Receiver index exceeds # of DC receivers'
                call exitMARE2DEM            
            endif       
            
        end select
    enddo
    write(*,*) ' '
    
    end subroutine readData       
!   
! !==================================================================================================================================! 
! !===================================================================================================================== writeJacobian
! !==================================================================================================================================! 
!     subroutine writeJacobian( nCurrentIter)
! 
!     use Occam 
!     use mare2dem_input_data_params
!     use mare2dem_global 
!  
!     
!     implicit none
!     
!     integer         :: nCurrentIter
!     character(50)   :: cNum
!     integer         :: lerr, i, j
!  
! 
!  
! !----------------------------
! !    
! ! Open Jacobian file
! !
!     write (cNum,*) nCurrentIter
!     open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.jacobian', iostat=lerr)
! 
! !
! ! Catch error opening file:
! !
!     if (lerr /= 0) then
!         write(*,*) ' Error opening Jacobian file, stopping!'
!         call exitMARE2DEM
!     end if
! !
! ! A-okay, write out the Jacobian
! ! 
! !
! ! Note that this routine was called from RunMARE2DEM and so WJ is just J (i.e. it has not been normalized by the data error bars):
! !
!     write(21,*)  nd, nParams
!  
!     do j = 1,nParams    
!         do i = 1,nd
!             write(21,'(g12.4)') wj(i,j)  
!         enddo
!     enddo
!      
!     close(21)
!       
!     end subroutine writeJacobian

      
          
!==================================================================================================================================! 
!======================================================================================================================= print2col32
!==================================================================================================================================! 
    subroutine print2col32(str1,str2,ioUnit) 
    
    implicit none
    
    character(len=*),intent(in)  :: str1, str2
    integer, intent(in)          :: ioUnit
    
 
    character(32) :: str1b
    
    str1b = str1 
    write(ioUnit,'(1x,a,a)') (adjustl(str1b)), trim(adjustl(str2))
     
    end subroutine print2col32
    
!-----------------------------------------------------------------------------------------------------------------------------------

    end module mare2dem_io
    

! These routines have to exist outside the module so that Occam can find it:    
!==================================================================================================================================! 
!===================================================================================================================== writeResponse
!==================================================================================================================================! 
    subroutine writeResponse( nCurrentIter)

    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global
 
    
    implicit none
    
    integer         :: nCurrentIter
    character(50)   :: cNum
    integer         :: lerr, i, j, ncolumns
    real(RealPrec)  :: resid
    

!    
! Open response file
!
    write (cNum,*) nCurrentIter
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.resp', iostat=lerr)

!
! Catch error opening response file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening response file, stopping!'
        call exitMARE2DEM
    end if
!
! No error, write out the responses:
!    
 
    write(21,'(a)') 'Format:    EMResp_2.2'
     
    ! UTM section:  
    if (len_trim(cUTMLine) > 0 ) then
        write(21,'(a,a)') 'UTM of x,y origin (UTM zone, N, E, 2D strike): ',trim(cUTMLine)
    endif   
    
    !
    ! CSEM section:
    !
    if (nTxCSEM>0) then
        write(21,'(a,a)') 'Phase convention: ',trim(adjustl(phaseConvention))
        write(21,'(a,a)') 'Reciprocity used: ',trim(adjustl(reciprocityUsed))
        write(cnum,'(i6)') nTxCSEM
        write(21,'(a,a)') '# Transmitters: ', adjustl(cnum)
        write(21,'(a1,8(a12,1x))') '!','X','Y','Z','Azimuth','Dip','Length','Type', 'Name'   
        do i = 1,nTxCSEM
           write(21,'(1x,6(f12.1,1x),a12,1x,a12)') xTxCSEM(i),yTxCSEM(i),zTxCSEM(i),azimuthTxCSEM(i),dipTxCSEM(i),lengthTxCSEM(i),&
                                                & cSourceType(i), trim(adjustl(cTxNamesCSEM(i)))
        enddo          
        write(cnum,'(i6)') nFreqCSEM
        write(21,'(a,a)') '# CSEM Frequencies: ', adjustl(cnum)    
        do i = 1,nFreqCSEM
            write(21,'(g13.5)') fTxCSEM(i)
        enddo        
        write(cnum,'(i6)') nRxCSEM
        write(21,'(a,a)') '# CSEM Receivers: ', adjustl(cnum)   
        write(21,'(a1,8(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta','Length', 'Name'  
        do i = 1,nRxCSEM
            write(21,'(1x,7(f12.1,1x),a12)') xRxCSEM(i),yRxCSEM(i),zRxCSEM(i), ThetaRxCSEM(i),AlphaRxCSEM(i),BetaRxCSEM(i), &
                                           & lengthRxCSEM(i), trim(adjustl(cRxNamesCSEM(i)))  
        enddo   
 
    endif
    
    !
    ! DC resistivity section:
    !
    if (nTxDC > 0) then
        write(cnum,'(i6)') nTxDC
        write(21,'(a,a)') '# DC Transmitters: ', adjustl(cnum)
        write(21,'(a1,3(a12,1x))') '!','Electrode A','Electrode B','Name'   
        do i = 1,nTxDC
           write(21,'(1x,2(i12,1x),1x,a12)') TxDC(i,1),TxDC(i,2), trim(adjustl(cTxNamesDC(i)))
        enddo          
    endif
    if (nRxDC > 0) then
        write(cnum,'(i6)') nRxDC
        write(21,'(a,a)') '# DC Receivers: ', adjustl(cnum)
        write(21,'(a1,3(a12,1x))') '!','Electrode M','Electrode N','Name'   
        do i = 1,nRxDC
           write(21,'(1x,2(i12,1x),1x,a12)') RxDC(i,1),RxDC(i,2), trim(adjustl(cRxNamesDC(i)))
        enddo          
    endif   
    if (nTrodesTxDC > 0) then
        write(cnum,'(i6)') nTrodesTxDC
        write(21,'(a,a)') '# DC Transmitter Electrodes: ', adjustl(cnum)
        write(21,'(a1,3(a12,1x))') '!','X','Y', 'Z' 
        do i = 1,nTrodesTxDC
           write(21,'(1x,3(f12.2,1x))') trodes_TxDC(i,1:3)
        enddo          
    endif    
    if (nTrodesRxDC > 0) then
        write(cnum,'(i6)') nTrodesRxDC
        write(21,'(a,a)') '# DC Receiver Electrodes: ', adjustl(cnum)
        write(21,'(a1,3(a12,1x))') '!','X','Y', 'Z' 
        do i = 1,nTrodesRxDC
           write(21,'(1x,3(f12.2,1x))') trodes_RxDC(i,1:3)
        enddo          
    endif     
    
    !
    ! MT section:
    !
    if (nRxMT>0) then    
        write(cnum,'(i6)') nFreqMT
        write(21,'(a,a)') '# MT Frequencies: ', adjustl(cnum)    
        do i = 1,nFreqMT
            write(21,'(g13.5)') fTxMT(i)
        enddo        
        write(cnum,'(i6)') nRxMT
        write(21,'(a,a)') '# MT Receivers: ', adjustl(cnum)   
        write(21,'(a1,9(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta','Length','SolveStatic','Name'
        do i = 1,nRxMT
            write(21,'(1x,7(f12.1,1x),i12,1x,a12)') xRxMT(i),yRxMT(i),zRxMT(i), ThetaRxMT(i),AlphaRxMT(i),BetaRxMT(i), &
                                       & lengthRxMT(i),iSolveStatic(i), trim(adjustl(cRxNamesMT(i)))
        enddo   
 
    endif    
    
    !
    ! Data and Model response section:
    !    
    write(cnum,'(i6)') nd      
    write(21,'(a,a)')'# Data: ', adjustl(cnum)   
    write(21,'(a1,4(a12,1x),4(a13,1x))') '!','Type','Freq#','Tx#','Rx#','Data','StdError','Response','Residual' 
    ncolumns = size(dp,2)
    do i = 1,nd
       if (sd(i) == 0) then  ! catch for dummy data file with 0 errors
            resid = 0
       else
            resid = (d(i)-dm(i))/sd(i)
       endif
       write(21,'(1x,4(i12,1x),3(g15.7,1x),g15.7)') (dp(i,j),j=1,ncolumns), d(i),sd(i), dm(i), resid
    enddo   
    close(21)
    
    end subroutine writeResponse
    
!==================================================================================================================================! 
!===================================================================================================================== writeJacobian
!==================================================================================================================================! 
    subroutine writeJacobian(nCurrentIter)

    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global 
       
    implicit none
    
    integer         :: nCurrentIter
    character(50)   :: cNum
    integer         :: lerr, i, j
 

 
!----------------------------
!    
! Open Jacobian file
!
    write (cNum,*) nCurrentIter-1  ! Note: since Jacobian is for starting model, we number if as current_iteration-1  

!
! Binary output:
!
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.jacobianBin', &
         &  form="unformatted", access="stream", status="replace", iostat=lerr)
!
! Catch error opening file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening Jacobian file, stopping!'
        call exitMARE2DEM
    end if
!
! A-okay, write out the Jacobian
! 
!
! Note that this routine was called from RunMARE2DEM and so WJ is just J (i.e. it has not been normalized by the data error bars):
!
    write(21)  nd, nParams
     
    do j = 1,nParams    
        do i = 1,nd
            write(21) wj(i,j)  
        enddo
    enddo
    close(21)


! Previous text file output version:
!
!     open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.jacobian', iostat=lerr)
! 
! !
! ! Catch error opening file:
! !
!     if (lerr /= 0) then
!         write(*,*) ' Error opening Jacobian file, stopping!'
!         call exitMARE2DEM
!     end if
! !
! ! A-okay, write out the Jacobian
! ! 
! !
! ! Note that this routine was called from RunMARE2DEM and so WJ is just J (i.e. it has not been normalized by the data error bars):
! !
!     write(21,*)  nd, nParams
!  
!     do j = 1,nParams    
!         do i = 1,nd
!             write(21,'(g12.4)') wj(i,j)  
!         enddo
!     enddo
!      
!     close(21)
      
    end subroutine writeJacobian
    
!==================================================================================================================================! 
!================================================================================================================= writeSensitivity
!==================================================================================================================================! 
    subroutine writeSensitivity(nCurrentIter)

    use Occam 
    use mare2dem_input_data_params
    use mare2dem_global 

    use fem2d_utilities  ! for generic FE operations
    use triangle_mesh
    use mare2dem_global
    use kx_io         ! for freeparams 
 
    implicit none
 
    character (24)  :: cend,tricommand        
    type(trimesh)   :: amesh ! see call_triangle.f90 for the definition of derived type

    
    integer         :: nCurrentIter
    character(50)   :: cNum
    integer         :: lerr, i, j
 
    real(8)                             :: sens, eleArea
    real(8), dimension(:), allocatable  :: area
    
     
        
    integer               :: e, n(3), iparamnum(4), iregion, nrho
    real(8), dimension(3) :: ye,ze,a,b,c
       
 
!----------------------------

  write (cNum,*) nCurrentIter-1  ! subtract one since Jacobian computed for starting model, which is current_iteration_number - 1


!
! Create a mesh:
!
    call copy_trimesh(inputmodel,amesh)
 
    write(cend,fmt='(F4.0)') minQangle
    tricommand = 'q'//trim(adjustl(cend))//'QpanjA'//CHAR(0)
   
    call call_triangle(tricommand,amesh)  

!
! Get number of resistivity components per model parameter:
!
    nrho = 0
    select case (trim(cAnisotropy))
    case ('isotropic')
        nrho = 1
    case ('isotropic_ip')
        nrho = 4    
    case ('isotropic_complex')
        nrho = 2               
    case ('triaxial')
        nrho = 3
    case ('tix','tiy','tiz')                               ! For transversely isotropic:
        nrho = 2          
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in subroutine writeSensitivity'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop
    end select  
 

!
! Get area of each parameter:
!    
    allocate(area(nParams))
    area = 0d0
    
    do e = 1,amesh%nele
   
        !
        ! Get the parameter number for each element:
        ! 
        iregion =  nint(abs(amesh%attr(e)))
       
        iparamnum = 0  
        do i = 1,nrho
            iparamnum(i) = iFreeParam(nrho*(iregion-1) + i)
        enddo  
        
        if ( all(iparamnum == 0) ) cycle  ! only include free parameters 
    
        !
        ! If free parameter, compute area:
        !

        n    = amesh%emap(1:3,e)
        ye   = amesh%y(n)
        ze   = amesh%z(n)   
        call get_abc_coeffs(ye,ze,a,b,c,eleArea)    
        
        do i = 1,nrho
        
            if (iparamnum(i)>0) area(iparamnum(i)) = area(iparamnum(i)) + eleArea
        
        enddo
    
    enddo
    
!    
! Open sensitivity file (binary version):
!
    open (unit=21, file=trim(outputFileRoot)//'.'//trim(adjustl(cNum))//'.sensitivity', &
         &  form="unformatted", access="stream", status="replace", iostat=lerr)

!
! Catch error opening file:
!
    if (lerr /= 0) then
        write(*,*) ' Error opening sensitivity file, stopping!'
        call exitMARE2DEM
    end if
 
!
! Compute sensitivity using area normalized sum of weighted Jacobian rows:
!
     
    do j = 1,nParams
    
        sens = 0d0    
    
        do i = 1,nd
            sens = sens + abs(wj(i,j))/sd(i) 
        enddo 
        
        write(21) sens/area(j)
        
    enddo
    close(21)

    deallocate(area)
    
    if ( allocated( inputmesh%attr ) )      call deallocate_trimesh(inputmesh,.false.)    
    
    end subroutine writeSensitivity   
