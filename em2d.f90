!-----------------------------------------------------------------------
!
!    Copyright 2008-2015
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
!======================================================================================================================== computeFwd
!==================================================================================================================================!  
    subroutine computeFwd( bDoPartials, currentMod )
!
! Routine to compute the 2D forward response and model Jacobian matrix
! This is called from Occam.
!
    use Occam           ! Passes out dm and wj here
    use kx_io           ! for rhoParams 
    use mare2dem_global
    
    implicit none

!
! Arguments
!
    logical, intent(in)        :: bDoPartials           ! logical flag for computing Jacobians
    real(RealPrec), intent(in) :: currentMod(nFree)   ! input model parameters

!
! Local variables:
!
    integer :: i 

!
! Use default amr setting for fwd call inside computeJacobian (bdopartials == .true.). No need to chance maxnadapt
! Use no amr (maxnadapt_noamr) otherwise and don't save the meshes
!
    ! non computeJacobian call
    if ( lReuseRefine .and. (bDoPartials == .false.)) then
       maxnadapt = maxnadapt_noamr
       lSaveMeshFiles = .false.
    endif
    !computeJacobian call
    if ( lReuseRefine .and. (bDoPartials == .true.)) then
       ! Check if this computeJacobian call will reuse refine
       if (mod(currentIteration,nbOccItToReuseRefine) > 0) then
           !disable amr and save meshes (enable intra occam iteration reuse refine)
           maxnadapt = maxnadapt_noamr
           lSaveMeshFiles = .false.
       else
           maxnadapt = maxnadapt_default
           lSaveMeshFiles = .true.
       endif
    endif
!
! Copy the model parameters (log10(resistivity)) in currentMod into the resistivity array used by the fwd code:
!   
    call insertFreeParams(currentMod)
      
! 
! Compute responses and sensitivities for 2D CSEM and MT data using the MARE2DEM parallel finite element kernel:
! 
    linversion = bDoPartials
    if (lPrintDebug) write(*,*) 'calling mpi_mare2dem...'
    
    call mpi_mare2dem

 !
 ! If inversion Jacobian call, there may be a few extra steps to perform:
 !
    if (bDoPartials) then

        
        ! Update data weights used for joint inversion (d_wt(:) = 1 if only CSEM or MT data present)
        call compute_DataWeights
        
        ! Check for nan's and infinity and replace them with 0's so the inversion can plow through:
        do i = 1,nd  
           where ( wj(i,:) /=  wj(i,:)   )  wj(i,:) = 0d0   ! replace nan's with 0 ! kwk debug: should put warning about this!
           where ( wj(i,:) >=  huge(1d0) )  wj(i,:) = 0d0    ! replace +infinity with 0
           where ( wj(i,:) <= -huge(1d0) )  wj(i,:) = 0d0    ! replace -infinity with 0
        enddo    
        
     endif

!
! Compute the MT static shifts, if requested:
!
    if (lPrintDebug) write(*,*) 'compute_MT_staticShifts...'
    
    call compute_MT_staticShifts    
    
    
    if (lPrintDebug) write(*,*) 'leaving computeFwd in EM2D.f90...'
    
    end subroutine computeFwd
 
!==================================================================================================================================! 
!============================================================================================================ displayJointInvMisfits
!==================================================================================================================================!     
    subroutine displayJointInvMisfits
!
! Prints joint inversion misfits (if joint inversion). Also shows estimated MT static shifts.
!    
    use Occam            
    use mare2dem_global 
    use mare2dem_input_data_params    
       
    implicit none
       
    integer         :: i     
    real(RealPrec)  :: residMT, residCSEM, residDC
    character(128)  :: cStr    

!----------------------------------------------
! Display MT and CSEM misfits:
!
    residMT  = 0
    residCSEM = 0
    do i = 1,nd
        if (sd(i) == 0 ) cycle ! skip if data uncertainty is 0
        
        if ( dp(i,1) < 100 ) then ! CSEM data:
               
             residCSEM = residCSEM + ( (d(i) - dm(i)) /sd(i))**2
             
        elseif ( (dp(i,1) > 100).and.(dp(i,1) < 200) ) then ! MT data: dp(i,1) > 100 
 
              residMT = residMT + ( (d(i) - dm(i)) /sd(i))**2

        elseif ( (dp(i,1) > 200).and.(dp(i,1) < 300) ) then ! DC resistivity
 
              residDC = residDC + ( (d(i) - dm(i)) /sd(i))**2
                            
        endif

    enddo
    
    ! Display misfit for each type separately:
    
    if ( (nCSEM > 0 ) .and.( (nMT > 0 ).or.(nDC > 0 ) ) )  then
        write(cStr,'(a32,g16.4)') 'CSEM Misfit: ', sqrt(residCSEM/nCSEM)
        call printOccamLog(cStr)   
    endif 
    if ( (nMT > 0 ) .and.( (nCSEM > 0 ).or.(nDC > 0 ) ) ) then
        write(cStr,'(a32,g16.4)') 'MT Misfit: ', sqrt(residMT/nMT)
        call printOccamLog(cStr)    
    endif
    if ( (nDC > 0 ) .and.( (nCSEM > 0 ).or.(nMT > 0 ) ) ) then
        write(cStr,'(a32,g16.4)') 'DC Resistivity Misfit: ', sqrt(residDC/nDC)
        call printOccamLog(cStr)    
    endif     
 
   
     ! Display any MT static shift solutions:
     
     if ( nRxMT > 0 ) then
     
         if (any(iSolveStatic > 0)) then
             
             write(cStr,'(a)') ' MT static shift estimates: '
             call printOccamLog(cStr)   
             write(cStr,'(4(a12,2x),a34)') ' Site #' ,'Name', 'TE Factor' , 'TM Factor', ' where linear ApRes = ApRes*Factor'
             call printOccamLog(cStr)   
             do i = 1,nRxMT
                 if (iSolveStatic(i) > 0) then
                    write(cStr,'(i12,2x,a12,2x,f12.3,2x,f12.3)') i, trim(cRxNamesMT(i)), 10**pm_assoc(2*i-1), 10**pm_assoc(2*i)
                    call printOccamLog(cStr) 
                 endif
             
             enddo
         endif
     
     endif
          
    end subroutine displayJointInvMisfits    
       
!==================================================================================================================================! 
!================================================================================================================== insertFreeParams
!==================================================================================================================================! 
    subroutine insertFreeParams(currentMod)
    
    use mare2dem_global
    use Occam           ! Passes out dm and wj here
    use kx_io           ! for rhoParams 
    
    real(RealPrec), intent(in) :: currentMod(nFree)   ! input model parameters
        
    integer :: i,j, ict
    
    ict = 0
    
    select case (trim(cAnisotropy ))
                    
    case ('isotropic','tix','tiy','tiz','triaxial','isotropic_complex')
        do i=1,nRegions
            do j = 1,nRhoPerRegion    
                if (iFreeParam((i-1)*nRhoPerRegion + j) > 0 ) then
                    ict = ict+1
                    rhoParams((i-1)*nRhoPerRegion + j) = 10d0**currentMod(ict) ! convert log10 to linear
                endif 
            enddo        
        enddo   
            
    case ('isotropic_ip') ! Cole-Cole Model
      
        do i=1,nRegions
            do j = 1,nRhoPerRegion    
                if (iFreeParam((i-1)*nRhoPerRegion + j) > 0 ) then
                    ict = ict+1

                  ! log:
                    if (j==1) then
                        rhoParams((i-1)*nRhoPerRegion + j) = 10d0**currentMod(ict) ! Rho
                    else
                        ! log:
                        !rhoParams((i-1)*nRhoPerRegion + j) = (10d0**currentMod(ict)) ! Eta, Tau, C
                        ! linear:
                         rhoParams((i-1)*nRhoPerRegion + j) = currentMod(ict) ! linear Eta, Tau, C 
                    endif
                     
                endif 
            enddo   
          
        enddo        
              
    end select  
    
    end subroutine insertFreeParams
!==================================================================================================================================! 
!=============================================================================================================== compute_DataWeights
!==================================================================================================================================!      
    subroutine compute_DataWeights 
    
! kwk debug, February 27 2019: this needs to be updated for DC resistivity data too...    
    
    use Occam
    use mare2dem_global
    use mare2dem_input_data_params
    use mare2dem_io
    
    implicit none
     
    integer         :: i
    character(32)   :: cNum
    real(RealPrec)  :: weightMT, weightCSEM   , chiMT, chiCS
 
    if (lPrintDebug) write(*,*) 'compute_DataWeights...'
    
 !
 ! Modify data weights if joint CSEM/MT inversion:
 !   
 
    d_wt = 1.0
    
    if ( (nFreqCSEM > 0 ) .and. (nFreqMT > 0 ) ) then  ! joint inversion, compute weights:       

!
! First get the Chi^2 misfit:
!           
        chiMT  = 0
        chiCS = 0
        do i = 1,nd
            if (dp(i,1) < 100)     then       
                chiCS = chiCS + ( (d(i) - dm(i)) /sd(i))**2
            else
                chiMT = chiMT + ( (d(i) - dm(i)) /sd(i))**2
            endif 
        enddo

! Unity weights:
!        weightMT   =  1d0 
!        weightCSEM =  1d0   
!        write(*,*) ' '
!        write(*,*) 'Joint inversion: not using any CSEM nor MT data weighting'
!        write(*,*) ' '
!  

! Normalized Chi^2:      
!        weightMT   =  1d0/sqrt(real(nMT))    
!        weightCSEM =  1d0 /sqrt(real(nCSEM))  
!        write(*,*) ' '
!        write(*,*) 'Using joint inversion data weights based on normalized Chi^2 '
!        write(*,*) ' '
! 
        
! Normalized Chi^2 with misfit balancing:         
        weightMT   =  1d0/sqrt(real(nMT))*sqrt(chiMT/nMT) 
        weightCSEM =  1d0/sqrt(real(nCSEM))*sqrt(chiCS/nCSEM)               
        write(*,*) ' '
        write(*,'(a)') 'Using joint inversion data weights based on normalized Chi^2 and misfit balancing'
        write(*,*) ' '
 
        !
        ! Insert appropriate weights into d_wt vector:
        !  
        write(cnum,'(g12.3)') weightCSEM
        call print2col32(' CSEM Relative Weights: ',cnum,6)
        write(cnum,'(g12.3)') weightMT
        call print2col32(' MT   Relative Weights: ',cnum,6)
        write(*,*) ' '
        do i = 1,nd
            select case (  dp(i,1) ) ! what data type is this?
        
            case (1:99) 
        
                d_wt(i) = weightCSEM
        
            case(100:199)
        
                d_wt(i) = weightMT              
        
            end select

        enddo    
    
    endif
    
    end subroutine compute_DataWeights 
        
!==================================================================================================================================! 
!=========================================================================================================== compute_MT_staticShifts
!==================================================================================================================================! 
    subroutine compute_MT_staticShifts
!
! A simplified static shift estimator based solely on the size of the mean residual. 
!
! ****** USE WITH CAUTION! *****
! 
    use Occam
    use mare2dem_input_data_params
    
    implicit none

    include 'em_parameters.inc'
    
    integer         :: iRx, i, nte, ntm
    
    real(RealPrec)  :: tesum, tmsum
    
!
! Solve for the static shift parameters at each site, if requested:
!
    do iRx = 1,nRxMT
        
        tesum     = 0
        tmsum     = 0
        pm_assoc(2*iRx-1:2*iRx) = 0
        
        if (iSolveStatic(iRx) /= 0) then
        
            nte = 0
            ntm = 0
            
            do i=1,nd

                if ( dp(i,4)  == iRx ) then
            
                    select case (  dp(i,1) )  ! what data type is this?   
                    
                    case (indRhoZXY)    ! TE Apparent resistivity
                        tesum = tesum + log10(d(i)) - log10(dm(i))         
                        nte   = nte + 1
                    case (indlog10RhoZXY)                                    
                        tesum = tesum + d(i) - dm(i)
                        nte   = nte + 1     
                    case (indRhoZYX)   ! TM Apparent resistivity
                        tmsum = tmsum + log10(d(i)) - log10(dm(i)) 
                        ntm   = ntm + 1
                    case (indlog10RhoZYX) 
                        tmsum = tmsum + d(i) - dm(i)
                        ntm   = ntm + 1                    
                    end select
                endif     
            enddo
            if (nte > 0 ) then
                tesum = tesum / nte
            endif
            if (ntm > 0 ) then
                tmsum = tmsum / ntm
            endif            
            
            ! If only TE or TM statics requested, zero the other static estimate:
            if (iSolveStatic(iRx) == 2) then ! TE only
                tmsum = 0
            elseif (iSolveStatic(iRx) == 3) then ! TM only
                tesum = 0
            endif 
            
            !write(*,*) 'tesum, tmsum: ',iRx, tesum, tmsum
            
            ! Store the static shifts:
            pm_assoc(2*iRx-1) = tesum
            pm_assoc(2*iRx  ) = tmsum    
            
            
            ! Finally do a second pass to apply the static shifts to the forward responses:
            do i=1,nd

                if ( dp(i,4)  == iRx ) then
            
                    select case (  dp(i,1) )  ! what data type is this?   
                    
                    case (indRhoZXY)    ! TE Apparent resistivity
                        dm(i) = dm(i)*10**tesum       
                    case (indlog10RhoZXY)                                    
                        dm(i) = dm(i) + tesum
                    case (indRhoZYX)   ! TM Apparent resistivity
                        dm(i) = dm(i)*10**tmsum 
                    case (indlog10RhoZYX) 
                        dm(i) = dm(i) + tmsum    
                                   
                    end select
                endif     
            enddo 
        endif
    enddo
 
    end subroutine compute_MT_staticShifts
    
 
