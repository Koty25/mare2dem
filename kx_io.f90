!-----------------------------------------------------------------------
!
!    Copyright 2018
!    Kerry Key
!    Lamont-Doherty Earth Observatory
!    Columbia University
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

module kx_io

    use triangle_mesh
    
    implicit none      
    
    
!--------  
! Inputs:
!--------
    integer,public :: myID  

!
! Resistivity parameters:  The code uses a lookup table for both fixed and free conductivity parameters.
!
! For cAnisotropy:  
!
!  'isotropic':      rhoParams(i)         is sig  for mesh%attr = i
!  'triaxial':       rhoParams(3*(i-1)+1) is sigx for mesh%attr = i
!                    rhoParams(3*(i-1)+2) is sigy for mesh%attr = i
!                    rhoParams(3*(i-1)+3) is sigz for mesh%attr = i   
! 'tix','tiy','tiz': 
!                    rhoParams(2*(i-1)+1) is sig in symmetry axis (sigx for tix, sigy for ti sigz for tiz) for mesh%attr = i, 
!                    rhoParams(2*(i-1)+2) is sig in transverse plane (sigyz for tix,sigyz for tiy, sigxy for tiz) for mesh%attr = i
! 'isotropic_ip' Cole-Cole model:
!                    rhoParams(4*(i-1)+1) is rho0 for mesh%attr = i
!                    rhoParams(4*(i-1)+2) is eta  for mesh%attr = i
!                    rhoParams(4*(i-1)+3) is tau  for mesh%attr = i  
!                    rhoParams(4*(i-1)+4) is c    for mesh%attr = i 
! 'isotropic_complex' IP:
!                    rhoParams(4*(i-1)+1) is rho real for mesh%attr = i
!                    rhoParams(4*(i-1)+2) is rho imag for mesh%attr = i
 
           
    character(24), public                       :: cAnisotropy  = 'isotropic'  !'isotropic','triaxial','tix','tiy','tiz', 'isotropic_ip','isotropic_complex' 
    integer, public                             :: nrhoParams = 0    ! # of resistivities in rhoParams (=nRegions*nRhoPerRegion)
    real(8), dimension(:), allocatable, public  :: rhoParams         ! array of conductivities for for fixed and free parameters    
    integer, dimension(:), allocatable, public  :: iFreeParam        ! Integer of free parameter # for each value in rhoParams

!
! Transmitter parameters:
!
    real(8), public :: kx   ! 2.5D Fourier transformation wavenumber, ignored for MT computations
    real(8), public :: w    ! Angular frequency 
    integer, public :: nTx  ! set to 1 for MT
    real(8), dimension(:), allocatable, public :: xTx, yTx, zTx, azimuthTx,dipTx, lengthTx ! x,y,z location of the transmitter 
    logical, dimension(:), allocatable, public :: lNeg_ikx ! Set to true for neg ikx transmitters.
    character(8), allocatable, public          :: TxType(:)  !'edipole' or 'bdipole',    

    character(2), public :: sType   ! 'cs','mt', or 'dc'
    logical, public      :: lMTscatteredField = .false.  
    logical, public      :: lUseReciprocity = .false.  ! if true, this properly scales J to B or M to E transmissions by ommu
    
!
! Receiver parameters:
!
    integer, public :: nRx                                     ! Number of receivers       
    real(8), dimension(:), allocatable, public :: xRx, yRx, zRx, azimuthRx,dipRx, lengthRx ! x,y,z location of the transmitter 
    character(8), allocatable, public          :: RxType(:)  !'edipole' or 'bdipole',   
    
    
!
! DC resistivity Rx and Tx parameters:
! 
    integer, public :: ntrodes_tx, ntrodes_rx
    real(8), dimension(:,:), allocatable, public ::  electrode_tx,electrode_rx ! x,y,z positions of electrodes for Tx and Rx.
                                           
!
! Quadrature weights used for finite dipole integrations:
!
    integer, public                     :: nquad_Tx = 0, nquad_Rx = 0

!
! Data masks. Set to true if EM fields to be computed at particular Rx-Tx pair.
!
    integer, dimension(:,:), allocatable, public :: iDataMask ! 0 = no data, 1 = no data but include for refinement, 2=data present
    logical :: l_has_TM_mode = .false. 
    logical :: l_has_TE_mode = .false.
    
!
! Mesh parameters are stored in the derived type 'trimesh' for ease of use
!
    type (trimesh), public :: mesh    
    character(256), public :: fileroot , cfileroot        ! base name for mesh files
    integer, public        :: meshnumber                  ! keeps track of mesh refinement number (1 = starting mesh)
    
!
! Local a priori refinement parameters:
!
    logical, public :: lLocalRefine  = .true.        ! Set to true to refine grid on input using skin depth rules
    real(8), public :: nTxEleArea = 1d0, nRxEleArea = 1d0,  nRxEleSkinDepths = 20., nSkinDepthsRegional = 1./8. ! Tolerances for a priori mesh refinement around Tx's and Rx's

    real(8), public :: nTxEleAreaDC = 1d-2, nRxEleAreaDC = 1d-2
                                                                                                                      
!
!  Adjustable parameters controlling the adaptive refinement and error estimation:
!
    real(8), public :: errortolerance    = 1.0       ! tolerance for relative error
    real(8), public :: minRangeProtector = 10        ! (m) protects agains over refinement close to the Tx
    real(8), public :: rxCloud           = 0         ! cloud of refinement around Rx and Tx. Just for testing, don't use this!
    real(8), public :: minArea           = 1d-4      ! minimum area to stop refinement from making very tiny elements  

    integer, public :: max_nsubrefine    = 1         ! max # of projection subrefinements for a given error estimator    
    real(8), public :: pct_refine        = 10.        ! percent of worst elements to refine per step (each of max_nsubrefine)    
    real(8), public :: ecutoff           = 1d-18     ! E(kx) error floor in V/m/Am  
    real(8), public :: hcutoff           = 1d-15     ! H(kx) error floor A/m/Am, note  B = 1.26x10^-6 * H
    integer, public :: maxRefinements    = 30        ! stop runaway refinements...    
    integer, public :: maxMeshNodes      = 100000    ! max# nodes in mesh, stop refinement if mesh grows this big. This will keep 
                                                     ! the code from running out of memory, but response accuracy may be compromised                        
!
! Model response derivatives:
!    
    logical, public :: lCompDerivs ! set =.true. to output the field derivatives with respect to conductivity
    
!
! Other adjustable parameters:
!
    real(8), public         :: minqangle          = 25.    ! quality angle for triangle.c (minimun inner angle for triangulation)
    character(32), public   :: linearSolver       = 'intelmkl'    ! linear solver for FE systems: 'superlu' or 'intelmkl'
    logical, public         :: lprintDebug_em2dkx = .false.
    logical, public         :: lprintDebug_dc2dkx = .false.
    logical, public         :: lDisplayRefinementStats = .true.   
    logical, public         :: lSaveLinearSystem  = .false. ! write out linear system Ax=b to files A = *.lhs, b=*.rhs, x=*.exhx
    logical, public         :: lSaveMeshFiles     = .false. ! set to .true. to write out the mesh to .poly, .ele and .node files   
    integer, public         :: idual_func         = 2       ! 0,1,2 for G0,G1,G2 error functionals from our paper
    logical, public         :: lSolveBumpGaussSeidel= .true.   
    logical, public         :: lprintTrace_em2dkx = .false. ! set to .true. to print tracing information of em2dkx
    real(8), public         :: t0_em2dkxTrace     = 0.       ! Initial time to be used during em2dkx tracing

!-----------
! Outputs:
! ---------
 
! Requested responses for various components: 

    integer, public                                 :: nComponents !kwk debug: these two are Inputs not outputs
    integer, dimension(:,:), allocatable, public    :: Components  ! [nComps x 2] for iRx,iTx, Rx and Tx types determine E or B.
    complex(8), dimension(:), allocatable, public   :: Fields      ! nComps array of field response for each requested component 
    complex(8), dimension(:,:), allocatable, public :: dFieldsdRho  ! sensitivity for each requested field component if inversion
    
    public :: deallocate_kx_io
    
    contains
    
!-----------------------------------------------------------------------------------------------------------------------------------
    subroutine deallocate_kx_io

!
! Deallocate public i/o variables in this module
!    
    if ( allocated( rhoParams ) )       deallocate( rhoParams )   
    if ( allocated( iFreeParam ) )      deallocate( iFreeParam )      
    if ( allocated(xTx) )               deallocate ( xTx,yTx,zTx,azimuthTx,dipTx,TxType, lengthTx ) 
    if ( allocated(lNeg_ikx) )          deallocate ( lNeg_ikx ) 

    if ( allocated(xRx) )               deallocate ( xRx,yRx,zRx,azimuthRx,dipRx,RxType, lengthRx ) 

    if ( allocated(electrode_rx) )     deallocate( electrode_tx,electrode_rx )
    
    if (allocated(Fields))             deallocate(Fields)
    if (allocated(dFieldsdRho))        deallocate(dFieldsdRho)
    if (allocated(Components))         deallocate(Components)
    if (allocated(iDataMask))          deallocate(iDataMask)
 
    
    end subroutine deallocate_kx_io   
    
end module kx_io
 
