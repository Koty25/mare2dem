# 1 "em2dkx.input.F90"
!-----------------------------------------------------------------------
!
!    Copyright 2008-2018
!    Kerry Key
!    Lamont-Doherty Earth Observatory
!    Columbia University
!
!    Formerly:
!    Scripps Institution of Oceanography
!    University of California San Diego
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
!
! Subroutine emd2kx compute the 2.5D CSEM responses at a given wavenumber (kx) using
! a goal-oriented adaptive finite element solution to the Ex-Hx CSEM formulation.
! Also can compute the 2D MT response using either a total or scattered field
! approach.
!
!==================================================================================================================================!
!======================================================================================================================== em2dkx_mod
!==================================================================================================================================!
# 1 "/home/users/fpmjunior/scorep/include/scorep/SCOREP_User.inc" 1 
!
!  This file is part of the Score-P software (http://www.score-p.org)
!
!  Copyright (c) 2009-2011,
!  RWTH Aachen University, Germany
!
!  Copyright (c) 2009-2011,
!  Gesellschaft fuer numerische Simulation mbH Braunschweig, Germany
!
!  Copyright (c) 2009-2011, 2014-2015, 2019,
!  Technische Universitaet Dresden, Germany
!
!  Copyright (c) 2009-2011,
!  University of Oregon, Eugene, USA
!
!  Copyright (c) 2009-2011, 2013, 2016,
!  Forschungszentrum Juelich GmbH, Germany
!
!  Copyright (c) 2009-2011,
!  German Research School for Simulation Sciences GmbH, Juelich/Aachen, Germany
!
!  Copyright (c) 2009-2011, 2015-2016,
!  Technische Universitaet Muenchen, Germany
!
!  Copyright (c) 2016,
!  Technische Universitaet Darmstadt, Germany
!
!  This software may be modified and distributed under the terms of
!  a BSD-style license.  See the COPYING file in the package base
!  directory for details.
!

!> @file
!! @ingroup    SCOREP_User_Interface
!!
!! @brief User interface for manual instrumentation of Fortran files.
!<



# 43




!************************************************************* Data types






!************************************************************* Constants



















!************************************************************* Regions

# 77


# 80


# 83


# 86


# 89


# 92


# 95


# 98












# 112




!************************************************************* Parameter

# 119


# 122


# 125


# 128


!************************************************************* User metric

# 133


# 136


# 139


# 142


# 145


!************************************************************* User topology macros

# 150


# 153


# 156


# 159


# 162


!************************************************************* User control







!************************************************************* Empty defines

# 213

# 38 "/home/users/fpmjunior/mare2dem-omeep-newest/em2dkx.f90" 2 
   
    module em2dkx_mod 

    use EM_constants
        
    use fem2d_utilities  ! for generic FE operations
    use MT1D_module      ! for the 1D boundary conditions for 2D MT computations
    use kdtree2_module   ! for fast searches
    use triangle_mesh    ! for trimesh structure
    
    use superlu_zsolver  ! for superlu sparse solver (optional)
    use intelmkl_solver
    
    use quad             ! for triquad integration of scattered field sources when jumps in sigma_0
    use sort             ! for quicksort routines
    
    use kx_io            ! holds i/o common to all 2D kernels (em2dkx,dc2dkx,etc)
    
    implicit none        ! I declare that you must declare everything
    
    private  ! The default is that all variables, subroutines and functions are private
    
!-----------------------------------------------------------------------
! Public interface:
!-----------------------------------------------------------------------
!
! Public subroutines:
!
    public :: em2dkx, displayRefinementStats_em2dkx


 
   
!-----------------------------------------------------------------------
! Private variables:
!-----------------------------------------------------------------------
!
! These are allocated and handled only inside em2dkx and the private attribute
! means you can't access these from outside the em2dkx module.
!
!
! A few adjustable parameters (don't change these unless you know what you are doing):
!

!
! Variables and arrays used across multiple subroutines:
!
    real(8), private                                 :: kx2
    complex(8), private                              :: ikx, ommu, iomeps 
    integer, dimension(:,:), allocatable, private    :: eRx 
    integer, dimension(:,:), allocatable, private    :: eTx
    logical, dimension(:), allocatable, private      :: lCloud
    integer, private                                 :: nedges 
    integer, dimension(:,:), allocatable, private    :: edges
    integer, dimension(:), allocatable, private      :: node2tri

!
! Sparse LHS and RHS arrays for linear FE's
!
    complex(8), dimension(:), allocatable, private   :: val
    integer,    dimension(:), allocatable, private   :: col, irw
    integer, private                                 :: nnz, info
    complex(8), dimension(:,:), allocatable, private :: rhs
    
    type(superlu_z)      :: superlu_p  ! Derived type that holds pointers to superlu factorization for primary/dual system
    type(mkl_dss_handle) :: intelmkl_p ! Derived type that holds pointers to intelmkl factorization for primary/dual system
    
!
! Arrays for uncoupled E and H MT solver:
!
    complex(8), dimension(:), allocatable, private   :: val_lE
    integer,    dimension(:), allocatable, private   :: col_lE, irw_lE
    integer, private                                 :: nnz_lE, info_lE 
    
    complex(8), dimension(:), allocatable, private   :: val_lH
    integer,    dimension(:), allocatable, private   :: col_lH, irw_lH 
    integer, private                                 :: nnz_lH, info_lH 
    
    complex(8), dimension(:,:), allocatable, private :: rhs_lE, rhs_lH
    
    type(superlu_z)      :: superlu_p_E, superlu_p_H   ! pointers to superlu factorization for primary/dual system
    type(mkl_dss_handle) :: intelmkl_p_E, intelmkl_p_H ! pointers to intelmkl factorization for primary/dual system
           
!
! Sparse LHS and RHS arrays for quadratic bump FE's
!
! Coupled anisotropic:
    complex(8), dimension(:), allocatable, private   :: val_q 
    integer,    dimension(:), allocatable, private   :: col_q, irw_q 
    integer, private                                 :: nnz_q, info_q

    complex(8), dimension(:,:), allocatable, private :: rhs_q 
    complex(8), dimension(:,:), allocatable, private :: rhs_wh 
    

    type(superlu_z)      :: superlu_q  ! pointers to superlu factorization for bump system
    type(mkl_dss_handle) :: intelmkl_q ! pointers to intelmkl factorization for bump system
    

    integer, dimension(:), allocatable, private      :: nodlab,bcnod   
    real(8), dimension(:), allocatable, private      :: errnrm          
    character(32) :: retricommand

    real(8), dimension(:), allocatable, private   :: jx,jy,jz  ! Electric transmitter  dipole moments for each nTx.
    real(8), dimension(:), allocatable, private   :: mx,my,mz  ! Magnetic transmitter dipole moments for each nTx.
    real(8), dimension(:,:), allocatable, private :: momentRx  ! 3 x nRx vector moments for each receiver component.

! quadrature weights and points:
    real(8), dimension(:), allocatable  :: quad_weights_Tx,     quad_xint_Tx  ! integration weights and local position
    real(8), dimension(:), allocatable  :: quad_weights_Rx,     quad_xint_Rx  ! integration weights and local position
    real(8), dimension(:,:), allocatable:: xTxQ, yTxQ, zTxQ, xRxQ, yRxQ, zRxQ                               
    
    
!
! Kd-tree pointer for fast element searches:
!
    type(kdtree2), pointer, private  :: tree    

!
! MT 1D fields for boundary conditions and possibly scattered field computations:
!
    type(MT1D)                  :: Left1DTE,Left1DTM,Right1DTE,Right1DTM,Scat1DTE,Scat1DTM
        
!
! Timers (only needed for debugging and performance optimizations)
!
    logical, private, parameter :: lprinttimers = .false.   
    real(8), private :: t0,t1,t2,t3,t4,t5,t6,t7,t7b,t8,t9
    real(8), private :: t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21    
    
!
! For adaptive refinement loop:
!
    integer :: newNumNodes, oldNumNodes, nRefinements  
    logical :: lConverged    
    complex(8), dimension(:), allocatable  ::  fieldsLast   ! used for convergence test
 
    real(8) :: maxerr           ! maximum error at any site
    logical :: lmadeNewMesh     ! true if a new mesh has been made
    type (trimesh) :: newmesh   ! temp storage for new adaptive meshes
    
    integer :: isubset,iRefinementGrp
    
                   
    contains


 
!==================================================================================================================================!
!============================================================================================================================ em2dkx
!==================================================================================================================================!
    subroutine em2dkx(iRefinementGrp_in,isubset_in,lrefine)   
!
! The main subroutine for computing the 2.5D EM fields in the kx domain
! for the input transmitter, receivers, frequency, wavenumber and mesh.
!
! This routine does both 2D MT and 2.5D CSEM computations.
!
! Updated 2018
!
! Kerry Key
! Lamont-Doherty Earth Observatory
! Columbia University
 
    implicit none
    
    integer, intent(in) :: isubset_in,iRefinementGrp_in
    logical             :: lrefine
    real(8)             :: maxRelDiff, time
        
    if (lprintDebug_em2dkx) write(*,*) myID,': entered em2dkx...'     

!
! Initialize a few things:
!
    allocate( fieldsLast(nComponents) ) 
    fieldsLast     = 0.
    lConverged     = .false.    
    isubset        = isubset_in
    iRefinementGrp = iRefinementGrp_in
    
!
! First do local refinement around Rx and Tx region:
!
       
    if ( lLocalRefine) call localRefinement_em2dkx()  
    
!
! Outer adaptive refinement loop:
!

    
    if (maxRefinements == 0) lrefine = .false.
    
    nRefinements = 0
    
    do while (.not.lConverged)
     
!
! Initializations:
!
        call cpu_time(t0)
    
        call em2dkx_initialize
 
        call cpu_time(t1)
    
!
! Set up dipoles
!
        call get_finite_dipoles
        
        call cpu_time(t2)
!
! Get element based lists (for node2tri, Rx and Tx locations):
!
        call get_element_lists
        
        call cpu_time(t3)
    
!
! Generate the node labels and boundary condition flags for nodes on mesh boundaries
!
    
        call gen_nodlab  
        
        call cpu_time(t4)

        call gen_bcnod   
        
        call cpu_time(t5)  

!
! Generate LHS matrix, rhs source vectors and solve the linear system(s):
!
        call gen_lhs 
    
        call cpu_time(t6)
    
        call gen_rhs
    
        call cpu_time(t7)
     
        call solve_primal
        
        call cpu_time(t8)
        
!
! Extract solution at receivers:
!
        call fill_solmtx  
        
        call cpu_time(t9)
        
!
! Check for solution convergence:
!
        call checkConvergence(maxRelDiff,lrefine)

        call cpu_time(t10)
                
!
! Compute error estimator and refine mesh if needed:
!
        if (.not.lConverged)  call estimateerror()

        call cpu_time(t20)         

!
! If requested, compute the adjoint partial derivatives:
!
        if ( lCompDerivs .and. lConverged )  call comp_adj_derivs      

        call cpu_time(t21)
!
! Display timing info:
!
        if (lprinttimers)  call em2dkx_printTimers  
        
!
! Display iteration stats:
!
        time = t21 - t0
        if (lDisplayRefinementStats) call displayRefinementStats_em2dkx(iRefinementGrp, isubset, meshnumber,lConverged, & 
                                         & oldNumNodes,newNumNodes,maxerr,maxRelDiff,time)
            
!
! Deallocate arrays:
!
        
        if (lmadeNewMesh) then
            call deallocate_trimesh(mesh,.false.)
            call copy_trimesh( newmesh, mesh )
        endif
        call deallocate_trimesh(newmesh,.false.)
        
        call em2dkx_deallocate  
        
            

    
    enddo ! do while (.not.lConverged)
       
    deallocate(fieldsLast) 
       
!
! All done, goodbye
!
    if (lprintDebug_em2dkx) write(*,*) myID,': leaving em2dkx...'
    
    end subroutine em2dkx

!==================================================================================================================================!
!===================================================================================================================== estimateerror
!==================================================================================================================================!
    subroutine estimateerror()
!
! The dual-weighted residual (DWR) error estimation subroutine.
!
 
!
! Local variables:
!
    integer         :: e, irefine,iTx 
    character(256)  :: filename, cadapt
    

    integer*8, save :: estimateerror = -1
	
    call SCOREP_F_Begin(estimateerror,"estimateerror",0,"",367) 

    
    if (lprintDebug_em2dkx) write(*,*) myID,': errorindicator...'
  
!
! First check all Rx elements to see if they are not too small. If they are, abort the error estimation:
!
!     ict = 0
!     do iRx = 1,nRx
!         do iTx = 1,nTx
!             call getWeight(iRx,iTx,weight)
!             if (weight > 0d0) ict = ict + 1
!         enddo
!     enddo
!     if (ict == 0) then
!         lConverged = .true.
!         return ! skip error estimation
!     endif
    
!
!  Allocate local arrays:
!
    allocate (rhs_q(2*nedges,nTx))  
    rhs_q = 0d0
!
! Generate RHS residual vector in bump space:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_F ...'
    call gen_rhs_bump_F 
    call cpu_time(t11)
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_B(0) ...'  
    call gen_rhs_bump_B(0)  
    
    call cpu_time(t12)    
    
!
! Generate LHS matrix in bump space:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_bump  ...'
 
    call gen_lhs_bump
    
    call cpu_time(t13)
!
! Solve the linear system using either SuperLU or Intel MKL solver (PARDISO)
!
    call solve_varepsilon(1)

    call cpu_time(t14)
 
!
! Compute goal oriented error estimator:
!
  
!
! 1. Solve dual problem in linear basis
!
    allocate(rhs_wh(2*mesh%nnod,nTx))

    call gen_rhs_dual  
    
    call cpu_time(t15)
    
    call solve_dual 


! Remember that for negative wavenumbers we solved [A B; B^T C][-u v] = [-c d] so take negative of u
    if (sType == 'cs') then
        do iTx = 1,size(rhs_wh,2)
            if (lNeg_ikx(iTx)) then
                rhs_wh(1:2*mesh%nnod:2,iTx) = -rhs_wh(1:2*mesh%nnod:2,iTx)  ! note Fortan increment is at end, not middle like matlab...
            endif        
        enddo        
    endif
        
    call cpu_time(t16)
 
!
! 2.  Approximate dual error in bump basis
!
    call gen_rhs_bump_G  

    call gen_rhs_bump_B(1)     

    call cpu_time(t17)   
       
!
! 3. Solve linear system for \delta_h:
!
    call solve_varepsilon(2) 

    call cpu_time(t18)  
   
!
! 4. Compute local error estimate
!
    call copy_trimesh( mesh, newmesh )
    
    do irefine = 1, max_nsubrefine 
        
        deallocate(errnrm)
        allocate(errnrm(newmesh%nele))
        
        call comp_drw_error_proj()  
        
        if ( (irefine == 1 ).and.lSaveMeshFiles) then
        
            write(cadapt,'(i6)') meshnumber  
            cadapt = adjustl(cadapt)    
            
            filename  = trim(fileroot)//'.'//trim(cadapt)//'.error'
            open(16,file=trim(filename),status='REPLACE') 
            
            do e = 1,mesh%nele
                write(16,*) errnrm(e)
            enddo
            close(16)
            
            filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_ex'
            open(16,file=trim(filename),status='REPLACE') 
            
            do e = 1,mesh%nnod
                write(16,*) abs(rhs(2*e-1,1))
            enddo
            close(16)
            
            filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_hx'
            open(16,file=trim(filename),status='REPLACE') 
            
            do e = 1,mesh%nnod
                write(16,*) abs(rhs(2*e  ,1))
            enddo
            close(16)
            
            filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_wh_ex'
            open(16,file=trim(filename),status='REPLACE') 
            
            do e = 1,mesh%nnod
                write(16,*) abs(rhs_wh(2*e-1,1))
            enddo
            close(16)
            
            filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs_wh_hx'
            open(16,file=trim(filename),status='REPLACE') 
            
            do e = 1,mesh%nnod
                write(16,*) abs(rhs_wh(2*e  ,1))
            enddo
            close(16)         
     
            filename  = trim(fileroot)//'.'//trim(cadapt)
            call write_triangle(newmesh,filename)
 
        endif

        maxerr = sqrt( maxval(errnrm) )
        
        if ( maxerr < errortolerance/1d2) then 
            lConverged = .true. 
        endif
        
        oldNumNodes = newmesh%nnod
        
        call refineMesh()
        
        newNumNodes = newmesh%nnod
        
        if (newNumNodes > oldNumNodes )  then   
            lmadeNewMesh= .true.
        else
            lConverged = .true. 
        endif 
        
        if ( newNumNodes > maxMeshNodes) then
            lConverged = .true. 
        endif        

    enddo
    call cpu_time(t19)

!
! Deallocate dual array
!
    deallocate (rhs_wh)

    
!
! Smooth the mesh:
!
    if (lmadeNewMesh)  then
        call smooth_TriMesh(newmesh,10,.false.)
        meshnumber = meshnumber + 1
    endif
 
!
! Free up Solver memory:
!
    if ( (.not.lSolveBumpGaussSeidel) ) then
        select case (trim(linearSolver))
        case ('superlu')
            call superlu_zfree (superlu_q)        
         case ('intelmkl')
            call intelmkl_zfree(intelmkl_q)   
        end select 
    endif   

!
! Deallocate arrays:
!
    if (allocated(rhs_q))    deallocate (rhs_q)
    if (allocated(val_q))    deallocate (val_q)
    if (allocated(col_q))    deallocate (col_q)
    if (allocated(irw_q))    deallocate (irw_q)
    

    call SCOREP_F_RegionEnd( estimateerror ) !end Instrumentation

    
    end subroutine estimateerror    
    
!==================================================================================================================================!
!=================================================================================================================== comp_adj_derivs
!==================================================================================================================================!
    subroutine comp_adj_derivs   
!
! Computes the model parameter sensitivity matrix using the adjoint method
!
    implicit none   
                               
    integer         ::  i,esite, e, n(3), iTx, iparamnum(4), iRx, etxr, iqTx,nqtx
    
    real(8)         :: ye(3),ze(3),a(3),b(3),c(3),area, l(3),dldy(3), dldz(3), yp,zp 
    complex(8)      :: fex(3), fhx(3), fey, fez, dexdz, dexdy, dhxdy, dhxdz, ee(3),hh(3)
    complex(8)      :: fex_adj(3), fhx_adj(3), fey_adj, fez_adj 
    
    complex(8)                  :: sxx, syy, szz, sss, zrho
    complex(8)                  :: sigx,sigy,sigz,dsds0,dsdeta,dsdtau,dsdc 
    complex(8)                  :: gammay2, gammaz2, expikdx     

    real(8), dimension(3,3)     :: M = [ 2d0, 1d0, 1d0,   1d0, 2d0, 1d0,   1d0, 1d0, 2d0] ! mass matrix, column major order
    
    integer                 :: iComp, isign, isignAdj(2), iq, nq, nrhs,irhs
    real(8)                 :: cosazi,sinazi,cosdip,sindip,sig0
    real(8)                 :: jax,jay,jaz,max,may,maz
 
    complex(8), dimension(:), allocatable :: ex1dmt,hx1dmt
    

    integer*8, save :: comp_adj_derivs = -1
	
    call SCOREP_F_Begin(comp_adj_derivs,"comp_adj_derivs",0,"",619) 


    if (lprintDebug_em2dkx) write(*,*) myID,': comp_adj_derivs...'     
!
! initialize to zero:
!
    dFieldsdRho  = 0d0          
    
    nrhs = 2
    if (sType == 'mt') then
        nrhs = 1
        isignAdj = 1
    else
! Adjoint field uses -ikx, unless primary field is for -ikx conjugate
        isignAdj(1) = -1  ! for +ikx primary
        isignAdj(2) =  1  ! for -ikx primary
    endif
 
    if ((sType == 'mt').and.(lMTscatteredField)) then  ! precompute 1D fields n times rather than 3 x nele (\approx 6n) if embedded in e loop
        allocate(ex1dmt(mesh%nnod),hx1dmt(mesh%nnod))
        do i = 1,mesh%nnod
            ex1dmt(i) = Scat1DTE%getEH1D(mesh%z(i),'e') 
            hx1dmt(i) = Scat1DTE%getEH1D(mesh%z(i),'h') 
        enddo
    endif
    
    allocate(rhs_wh(2*mesh%nnod,nrhs))   
            
!
! Loop over receivers and compute adjoint sources, then sensitivities:
!

    do iRx = 1,nRx

!
! Get dipole moments for source along component:
!
        jax = 0d0; jay = 0d0; jaz = 0d0
        max = 0d0; may = 0d0; maz = 0d0

        cosazi = cos(pi/180d0*azimuthRx(iRx) )
        sinazi = sin(pi/180d0*azimuthRx(iRx) )
        cosdip = cos(pi/180d0*dipRx(iRx))
        sindip = sin(pi/180d0*dipRx(iRx))   
               
        nq = 1       
        
        select case (trim(RxType(iRx)))
   
        case ('edipole')               
        
            jax = cosdip*cosazi
            jay = cosdip*sinazi
            jaz = sindip             
            
            if (lengthRx(iRx) > 0) nq = nquad_Rx    
             
        case ('bdipole')
            
            max = cosdip*cosazi
            may = cosdip*sinazi
            maz = sindip     
                                             
        end select 
        
        rhs_wh = 0d0
        
        
        do iq = 1,nq ! loop over receiver dipole's quadrature points
    
            yp = yRxQ(iRx,iq)
            zp = zRxQ(iRx,iq)                       
        
!
! Get anisotropic conductivity for element containing receiver
!
            esite = eRx(iRx,iq)
            call getSigs(esite,sigx,sigy,sigz) 
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)              

!
! Set up the basis functions and evaluate at the receiver location:
!
            n   = mesh%emap(1:3,esite)
            ye  = mesh%y(n)
            ze  = mesh%z(n)

            call get_linear_basis(ye,ze,yp,zp,l,dldy,dldz)               
            
!
! Compute adjoint fields for source at receiver:
!
            do irhs = 1,nrhs
                      
                ee =    -l*jax             - isignAdj(irhs)*ikx * ( dldy*jay/gammay2                  + dldz*jaz/gammaz2 )  & 
                    &                                    - ommu * ( dldy*maz*(sigy - iomeps)/gammay2  - dldz*may*(sigz - iomeps)/gammaz2 )  
                hh = -(                                  - ommu * ( dldy*jaz/gammaz2                  - dldz*jay/gammay2 )  &
                    &   -l*max*ommu - isignAdj(irhs)*ikx * ommu * ( dldy*may/gammaz2                  + dldz*maz/gammay2 )  )  


                ee = isignAdj(irhs)*ee ! for neg kx system

!
! Insert these into the correct rhs location
!
                if (nq > 1) then
                    expikdx = exp(-isignAdj(irhs)*ikx*(xRxQ(iRx,iq)-xRx(iRx)))
                    rhs_wh(2*n-1, irhs) = rhs_wh(2*n-1, irhs) + expikdx*ee*quad_weights_Rx(iq)/2.d0    ! Ex
                    rhs_wh(2*n  , irhs) = rhs_wh(2*n  , irhs) + expikdx*hh*quad_weights_Rx(iq)/2.d0    ! Hx
                else ! point dipole
                    rhs_wh(2*n-1, irhs) = ee ! Ex
                    rhs_wh(2*n  , irhs) = hh ! Hx
                endif                      
            
            enddo
            
        enddo  ! iq = 1,nq
        
!
! Solve the linear system:
!
        call solve_dual 
 
! Remember that for negative wavenumbers we solved [A B; B^T C][-u v] = [-c d] so take negative of u
        if (sType == 'cs') then
            do irhs = 1,nrhs
                rhs_wh(1:2*mesh%nnod:2,irhs) = isignAdj(irhs)*rhs_wh(1:2*mesh%nnod:2,irhs)  ! note Fortan increment is at end, not middle like matlab...
            enddo  
        endif
        
!
! Loop over all elements and integrate E^a(-ikx) dot E(ikx) :
!
        do e = 1, mesh%nele

            call getParamNums(e,iparamnum)

            if ( all(iparamnum == 0) ) cycle  ! only include free parameters

            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)

            call get_abc_coeffs(ye,ze,a,b,c,area)   
            dldy = b/2d0/area
            dldz = c/2d0/area           
                
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)                
            
!
! Loop over requested components:
!
            do iComp = 1,nComponents
            
                if ( Components(iComp,1) /= iRx ) cycle 
            
                iTx = Components(iComp,2)  

! Get fields from source at transmitter:
                fex  = rhs(2*n-1,iTx) 
                fhx  = rhs(2*n  ,iTx)   

                if ((sType == 'mt').and.(lMTscatteredField)) then
                    fex = fex + ex1dmt(n)
                    fhx = fhx + hx1dmt(n)
                endif
        
                dexdy =  sum(fex*dldy)
                dexdz =  sum(fex*dldz)
                dhxdy =  sum(fhx*dldy)
                dhxdz =  sum(fhx*dldz)
            
                isign =  1
                if (lNeg_ikx(iTx)) isign = -1
        
                fey  = (-ommu*dhxdz - isign*ikx*dexdy )/gammay2   ! constant across element since basis is linear
                fez  = ( ommu*dhxdy - isign*ikx*dexdz )/gammaz2  

! Get adjoint fields from source at receiver:
                if (isign == 1) then
                    irhs = 1 ! get -ikx adjoint fields or MT adjoint
                else                        
                    irhs = 2 ! get +ikx adjoint fields
                endif                               
            
                fex_adj  = rhs_wh(2*n-1, irhs) 
                fhx_adj  = rhs_wh(2*n  , irhs)       
            
                dexdy =  sum(fex_adj*dldy)
                dexdz =  sum(fex_adj*dldz)
                dhxdy =  sum(fhx_adj*dldy)
                dhxdz =  sum(fhx_adj*dldz)     

                fey_adj  = (-ommu*dhxdz - isignAdj(irhs)*ikx*dexdy )/gammay2   
                fez_adj  = ( ommu*dhxdy - isignAdj(irhs)*ikx*dexdz )/gammaz2  

! Sxx sensitivity can be computed using a mass matrix:
                sxx  = sum(fex_adj*matmul(M,fex))*area/12d0

! Syy and Szz sensitivities:
                syy  = fey*fey_adj*area
                szz  = fez*fez_adj*area

! if element contains the adjoint source, add on source delta terms (see eqs 14 and 15 in Key and Ovall, 2011)
                if (e==esite) then
                    syy = syy  + fey*( ommu*jay/gammay2  - isignAdj(irhs)*ikx*maz*ommu/gammay2  )   
                    szz = szz  + fez*( ommu*jaz/gammaz2  + isignAdj(irhs)*ikx*may*ommu/gammaz2  )
                endif 
                
! if element contains actual source, add on source delta terms (see eqs 14 and 15 in Key and Ovall, 2011)
                if (sType == 'cs') then
                    nqtx = 1
                    if (lengthTx(iTx) > 0 ) nqtx = nquad_Tx
                    do iqTx=1,nqtx
                        etxr = eTx(iTx,iqTx)
                        if (e==etxr) then
                            syy = syy  + fey_adj*( ommu*jy(iTx)/gammay2  - isign*ikx*mz(iTx)*ommu/gammay2  )   
                            szz = szz  + fez_adj*( ommu*jz(iTx)/gammaz2  + isign*ikx*my(iTx)*ommu/gammaz2  )                        
                        endif
                    
                    enddo
                endif

                select case (trim(RxType(iRx)))
   
                case ('edipole')               
! nothing to do
                 
                case ('bdipole')
            
! Scale magnetic terms by ommu since  M_s^a = ommu*H_s^a
                    sxx = sxx / ommu
                    syy = syy / ommu
                    szz = szz / ommu      
                                             
                end select 
                
! convert from d/dsig to d/dlog10rho:
                select case (trim(cAnisotropy))
                
                case ('isotropic_ip','isotropic_complex') ! do nothing here
                
                case default  ! scale everything else since its regular conductivity
                    sxx = -dlog(10.d0)*sigx*sxx
                    syy = -dlog(10.d0)*sigy*syy
                    szz = -dlog(10.d0)*sigz*szz
                end select
                
! Insert into output dFieldsdRho array:
                
                select case (trim(cAnisotropy))

                case ('isotropic')
                
                    dFieldsdRho(iComp,iparamnum(1)) =  dFieldsdRho(iComp,iparamnum(1)) + sxx + syy + szz
        
                case ('isotropic_ip')
                    
                    call getIPderivs(e,dsds0,dsdeta,dsdtau,dsdc,sig0)
                    
                    sss = sxx + syy + szz 
        
                    if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) -   dlog(10.d0)*sig0*dsds0*sss  ! convert from d/dsig to d/dlog10rho
                    if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) +   dsdeta*sss 
                    if (iparamnum(3)>0) dFieldsdRho(iComp,iparamnum(3)) = dFieldsdRho(iComp,iparamnum(3)) +   dsdtau*sss 
                    if (iparamnum(4)>0) dFieldsdRho(iComp,iparamnum(4)) = dFieldsdRho(iComp,iparamnum(4)) +   dsdc  *sss                    
                                        
                case ('isotropic_complex')
! rho = a+i*b
! dsig/dlog10(a) = dsig/drho*drho/da*da/dlog10(a)
!                = -1/rho^2*1*a/log10(e) = -1/rho^2* a log_e(10)
! dsig/dlog10(b) = dsig/drho*drho/db*db/dlog10(b)
!                = -1/rho^2*1i*b/log10(e) = -1i/rho^2* b log_e(10)
                    zrho = 1d0/sigx
! convert from d/dsig to d/dlog10rho
                    sss = -(dlog(10.d0)/zrho**2)*(sxx + syy + szz)
                    
                    if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) +     dble(zrho)*sss 
                    if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + ic*aimag(zrho)*sss  
                                        
                case ('triaxial')
                    if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) =  dFieldsdRho(iComp,iparamnum(1)) + sxx
                    if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) =  dFieldsdRho(iComp,iparamnum(2)) + syy 
                    if (iparamnum(3)>0) dFieldsdRho(iComp,iparamnum(3)) =  dFieldsdRho(iComp,iparamnum(3)) + szz 

                case ('tix') ! For transversely isotropic:
! First param is sig along symmetry axis; second param is along transverse plane
                    if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + sxx 
                    if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + syy  + szz

                case ('tiy')

                    if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + syy
                    if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + sxx + szz  

                case ('tiz')

                     if (iparamnum(1)>0) dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + szz 
                     if (iparamnum(2)>0) dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + syy + sxx

                end select              
            
            enddo ! iComp
            
            
        enddo ! e = 1 nele

    enddo  !  iRx
 
    deallocate(rhs_wh)
                
    if (allocated(ex1dmt)) deallocate(ex1dmt,hx1dmt)
    

    call SCOREP_F_RegionEnd( comp_adj_derivs ) !end Instrumentation


    end subroutine comp_adj_derivs
  
!==================================================================================================================================!
!================================================================================================================== checkConvergence
!==================================================================================================================================!
    subroutine checkConvergence(maxRelDiff,lrefine)
    
    real(8),intent(out) :: maxRelDiff
    logical,intent(in)  :: lrefine
    
    integer :: i,iRx,iTx
    real(8) :: wt,rTx 
    
    if (lprintDebug_em2dkx) write(*,*) myID,': checkConvergence...'  
    
    lConverged = .false.
    
    oldNumNodes = mesh%nnod
    newNumNodes = 0
    maxRelDiff  = 0d0    
    
    if (.not.lrefine)  then
        maxerr = -1 
        lConverged = .true.
        return
    endif
 
    if (nRefinements == 0) then
        fieldsLast = fields
        return
    endif
    
!
! Compute the numerical difference in the computed fields for each
! subsequent adaptive refinement iteration.
!
    if (sType == 'cs') then
 
        do i = 1,nComponents
        
            wt = 1.               
            iRx = Components(i,1)
            iTx = Components(i,2)
            rTx  = sqrt( ( yRx(iRx)- yTx(iTx))**2 + ( zRx(iRx)- zTx(iTx))**2 )
            if (rTx < minRangeProtector)  wt = 0.
  
            select case (Components(i,3))  ! KWK debug: need a total field option here to deal with singularities...
            case(1,2,3)
                fieldsLast(i) = wt*abs(fieldsLast(i) - fields(i))/(abs(fields(i))+ecutoff/kx)
            case(4,5,6)
                fieldsLast(i) = wt*abs(fieldsLast(i) - fields(i))/(abs(fields(i))+hcutoff/kx)
            end select
                  
        enddo
    
    else ! MT, separate TE and TM:
        
        
! Full MT responses:
        do i = 1,nComponents,6
            
! relative difference for impedance:
! TE:
            fieldsLast(i)  = abs( fields(i) / fields(i+4) -  fieldsLast(i)  / fieldsLast(i+4) ) / abs( fields(i) / fields(i+4) )

! TM:
            fieldsLast(i+3)  = abs( fields(i+1) / fields(i+3) -  fieldsLast(i+1)  / fieldsLast(i+3) ) / abs( fields(i+1) / fields(i+3) )

            fieldsLast(i+1) = 0d0
            fieldsLast(i+2) = 0d0
            fieldsLast(i+4) = 0d0
            fieldsLast(i+5) = 0d0
        
        enddo
           
    endif    
    maxRelDiff = maxval(abs(fieldsLast))

    lConverged = .true. 
    if (maxRelDiff > errortolerance/100.) lConverged = .false. 
          
    fieldsLast = fields

! Give up if we've already refined the mesh too many times:
    if (meshnumber > maxRefinements ) lConverged = .true.
    
! Stop refinement if the mesh has grown too large:
    if ( mesh%nnod > maxMeshNodes) lConverged = .true.
 
    end subroutine checkConvergence
        
!==================================================================================================================================!
!============================================================================================================ displayRefinementStats
!==================================================================================================================================!
    subroutine displayRefinementStats_em2dkx( iRefinementGrp,isubset, meshnumber,lconv, oldnodes,newnodes,sumerr,maxRelDiff,time)
    
    integer :: meshnumber, oldnodes, newnodes,isubset,iRefinementGrp
    real(8) :: sumerr, time,maxRelDiff
    logical :: lconv
    
    character(256) :: sFmt
    character(32)  :: snewnodes, status, stime, smaxRelDiff
    character(14)  :: smaxerr
    
    status = 'adaptive'
     
    if (sumerr >= 0) then         
        if ( sumerr > 0 ) then
             write(smaxerr,'(a2,2x,g9.2)') 'G:',sumerr
        else
            smaxerr = ' '
        endif
        if ( (maxRelDiff < 1d10 ) .and. (maxRelDiff > 0d0 ) ) then     
            write(smaxRelDiff,'(a9,2x,g9.2,1x,a1)') 'max(Rel):',maxRelDiff*1d2,'%'
        else
             smaxRelDiff = ' '
        endif
        if (lconv) then 
            status = 'converged'
            smaxerr = ' '
        endif
             
    elseif (sumerr == -2d0) then
        status = 'local' ! special flag for local refinement mode
        smaxerr = ' '
        smaxRelDiff = ' '
     elseif (sumerr == -1d0) then
        status = 'shared' ! special flag for local refinement mode
        smaxerr = ' '   
        smaxRelDiff = ' '
!lconv = .false.
    else
        smaxerr = ' '
        smaxRelDiff = ' '
        
    endif
    if (newnodes > 0) then
        write(snewnodes,'(a3,2x,i6)') '-->',newnodes
    else
        snewnodes   = ' '
    endif       
    
    write(stime,'(a8,f9.3,a2)') ' Timer: ',time, ' s'
    
    sFmt = '(a5,2x,i6,2x,a6,2x,i6,2x,a7,2x,i3,2x,a7,2x,i2,2x,a7,2x,a9,2x,a8,2x,i6,2x,a11,2x,a14,2x,a22,2x,a19)'     
    write(*,sFmt)  'Proc:',myID,'Group:',iRefinementGrp,'Subset:',isubset, 'Mesh #:',meshnumber,'Status:', trim(status), &
                & '# Nodes:',oldnodes,  trim(snewnodes), (smaxerr), trim(smaxRelDiff), trim(stime)
    
    end subroutine displayRefinementStats_em2dkx         

!==================================================================================================================================!
!================================================================================================================= em2dkx_initialize
!==================================================================================================================================!
    subroutine em2dkx_initialize   
!
! Initializes variables, checks input options and allocates some arrays
!
! Kerry Key
! Scripps Institution of Oceanography
!

    implicit none       

    character(32)  :: cend,cSide

    integer :: iTx, iRx    
    real(8) :: cosazi,sinazi,cosdip,sindip,left_conductance,right_conductance
!
! Check sType is valid:
!
    if ((sType /= 'mt').and.(sType/='cs')) then
        write(*,*) ' Error in em2dkx, unknown sType: ',trim(sType)
        write(*,*) ' Needs to be mt or cs ' 
        write(*,*) ' Stopping !'
        stop 
    endif
!
! Check kx:
!
    if (sType == 'mt')  kx = 0             
 
!
! Check input linear solver parameter:
!
    select case (trim(linearSolver))
        case ('superlu','intelmkl')
        case default
            write(*,*) ' Error in em2dkx, unknown linear solver: ',trim(linearSolver)
            write(*,*) ' Stopping !'
            stop
    end select
 
!
! Initialize a few things:
!
    kx2    = kx**2
    ikx    = ic*kx
    ommu   = ic*w*mu0
    iomeps = ic*w*EPSILON
    
    lmadeNewMesh = .false.
    
! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0) 
 
!
! Allocate private arrays used only in em2dkx:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': allocating in em2dkx...'
    
    allocate ( errnrm(mesh%nele) )
    allocate ( nodlab(mesh%nnod), bcnod(mesh%nnod)  )
    allocate ( rhs(2*mesh%nnod,nTx) )
    allocate ( eRx(nRx,nquad_Rx) )
    allocate ( lCloud(mesh%nele) )
    allocate ( edges(3,mesh%nele) )
    allocate ( eTx(nTx,nquad_Tx) )
    
    if (sType == 'cs') then
    
        allocate ( jx(nTx),jy(nTx),jz(nTx),mx(nTx),my(nTx),mz(nTx) )

!
! Assign dipole moments:
!
        jx = 0d0; jy = 0d0; jz = 0d0
        mx = 0d0; my = 0d0; mz = 0d0
    
        do iTx = 1,nTx    
 
             
            cosazi = cos(pi/180d0*azimuthTx(iTx) )
            sinazi = sin(pi/180d0*azimuthTx(iTx) )
            cosdip = cos(pi/180d0*dipTx(iTx))
            sindip = sin(pi/180d0*dipTx(iTx))   
                   
            select case (trim(TxType(iTx)))
       
            case ('edipole')               
            
                jx(iTx) = cosdip*cosazi
                jy(iTx) = cosdip*sinazi
                jz(iTx) = sindip             
            
            case ('bdipole')
                
                mx(iTx) = cosdip*cosazi
                my(iTx) = cosdip*sinazi
                mz(iTx) = sindip     
                                                 
            end select 
        
        enddo
        

!
! Get Tx quadrature points:
!
        allocate( quad_weights_Tx(nquad_Tx), quad_xint_Tx(nquad_Tx) )  
        call legendre_compute_dr( nquad_Tx,quad_xint_Tx,quad_weights_Tx )
       
! write(*,*) 'em2dkx_initialize, quadweights, nquad_Tx:',nquad_Tx
        
    endif

    allocate( quad_weights_Rx(nquad_Rx), quad_xint_Rx(nquad_Rx) )  
    call legendre_compute_dr( nquad_Rx,quad_xint_Rx,quad_weights_Rx )
 
    allocate ( momentRx(3,nRx) )

!
! Assign receiver dipole moments:
!
    momentRx = 0d0

    do iRx = 1,nRx    

        cosazi = cos(pi/180d0*azimuthRx(iRx) )
        sinazi = sin(pi/180d0*azimuthRx(iRx) )
        cosdip = cos(pi/180d0*dipRx(iRx))
        sindip = sin(pi/180d0*dipRx(iRx))   
                       
        momentRx(1,iRx) = cosdip*cosazi     ! x
        momentRx(2,iRx) = cosdip*sinazi     ! y
        momentRx(3,iRx) = sindip            ! z
 
    enddo
                
!
! If MT, get left and right side 1D layering models used for boundary conditions:
!
    if (sType == 'mt') then
         
! kwk debug note:  note that gen_1d_mt sets the frequency, i.e Left1DTE%omega
        call gen_1d_mt(Left1DTE , 'left ', 'x') 
        call gen_1d_mt(Right1DTE, 'right', 'x')           
        call gen_1d_mt(Left1DTM , 'left ', 'y')      
        call gen_1d_mt(Right1DTM, 'right', 'y')           
 

! if scattered field get 1D model for it and any subsequent modifications...
        
        if (lMTscatteredField) then
            
! For 1D background model, we will use the side that has the greatest integrated conductance.
! ignoring for now possible x vs y anisotropy effects...ideally for marine MT models, this picks
! the side that has the thick ocean in it.
            
            left_conductance  = Left1DTE%getTotalConductance()
            right_conductance = Right1DTE%getTotalConductance()
            
            if (left_conductance > right_conductance) then
                cSide = 'left'
            else
                cSide = 'right'
            endif            
            call gen_1d_mt(Scat1DTE ,cSide,'x')   
            call gen_1d_mt(Scat1DTM ,cSide,'y')  
 
        endif
        
    endif    
    
    if (lprintDebug_em2dkx) write(*,*) myID,': leaving em2dkx_initialize'
        
    end subroutine em2dkx_initialize
        
!==================================================================================================================================!
!================================================================================================================= em2dkx_deallocate
!==================================================================================================================================!
    subroutine em2dkx_deallocate   
!
! Deallocates all remaining internal arrays at end of em2dkx
!
! Kerry Key
! Scripps Institution of Oceanography
!
    implicit none   
    
!
! Deallocate all superlu memory and associated arrays:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': deallocating sparse matrices...'

    if (sType == 'cs') then

        select case (trim(linearSolver))
         
        case ('superlu')
            call superlu_zfree(superlu_p)   
        case ('intelmkl')
            call intelmkl_zfree(intelmkl_p)  
        end select  
        
        if ( allocated(val) )       deallocate(val)
        if ( allocated(col) )       deallocate(col)
        if ( allocated(irw) )       deallocate(irw)    
 
 
        if ( allocated(quad_weights_Tx) )       deallocate(quad_weights_Tx)   
        if ( allocated(quad_xint_Tx) )          deallocate(quad_xint_Tx)    
        if ( allocated(xTxQ) )                  deallocate(xTxQ, yTxQ, zTxQ)
        
    else  ! MT
        
        select case (trim(linearSolver))
        
        case ('superlu')
            call superlu_zfree(superlu_p_E)   
            call superlu_zfree(superlu_p_H)  
        case ('intelmkl')
            if (l_has_TE_mode) call intelmkl_zfree(intelmkl_p_E)   
            if (l_has_TM_mode) call intelmkl_zfree(intelmkl_p_H)  
        end select 
        
        if ( allocated(rhs_lE) )    deallocate(rhs_lE) 
        if ( allocated(val_lE) )    deallocate(val_lE)
        if ( allocated(col_lE) )    deallocate(col_lE)
        if ( allocated(irw_lE) )    deallocate(irw_lE) 
        
        if ( allocated(rhs_lH) )    deallocate(rhs_lH)                
        if ( allocated(val_lH) )    deallocate(val_lH)
        if ( allocated(col_lH) )    deallocate(col_lH)
        if ( allocated(irw_lH) )    deallocate(irw_lH)                                           
 
    end if 



    if (lprintDebug_em2dkx) write(*,*) myID,': deallocating other arrays...'
    deallocate ( errnrm ) 
    deallocate ( nodlab, bcnod ) 
    deallocate ( rhs )  
    deallocate ( eRx )
    deallocate( lCloud)
    deallocate ( edges )
    deallocate ( eTx )
    call kdtree2_destroy(tree)  
    deallocate( node2tri )
 
    if ( allocated(jx) )                deallocate ( jx,jy,jz )  
    if ( allocated(mx) )                deallocate ( mx,my,mz )
    if ( allocated(momentRx) )          deallocate ( momentRx )
    
    if ( allocated(quad_weights_Rx) )   deallocate(quad_weights_Rx)    
    if ( allocated(quad_xint_Rx) )      deallocate(quad_xint_Rx)     
    if ( allocated(xRxQ) )              deallocate(xRxQ, yRxQ, zRxQ)    
    
    if (sType == 'mt') then    
     
! Deallocate 1D arrays:
        call delete(Left1DTE)   
        call delete(Right1DTE)
        call delete(Left1DTM)   
        call delete(Right1DTM)
        
        if (lMTscatteredField) then
            call delete(Scat1DTE)
            call delete(Scat1DTM)               
        endif
    
    endif    
          

    if ( allocated( newmesh%attr ) )    call deallocate_trimesh(newmesh,.false.)
    if (allocated(newmesh%area))        deallocate(newmesh%area)   
          
    end subroutine em2dkx_deallocate  
    
!==================================================================================================================================!
!================================================================================================================ em2dkx_printTimers
!==================================================================================================================================!
    subroutine em2dkx_printTimers 
    
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx em2dkx_initialize:    ',t1-t0
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx get_finite_dipoles:   ',t2-t1
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx get_element_lists:    ',t3-t2
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_nodlab:           ',t4-t3
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_bcnod             ',t5-t4        
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_lhs:              ',t6-t5
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_rhs:              ',t7-t6
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx solve_primal:         ',t8-t7
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx fill_solmtx:          ',t9-t8
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx checkConvergence:     ',t10-t9
    
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_rhs_bump_F:       ',t11-t10
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_rhs_bump_B:       ',t12-t11
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_lhs_bump:         ',t13-t12
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx solve_varepsilon(1):  ',t14-t13
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_rhs_dual:         ',t15-t14
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx solve_dual:           ',t16-t15        
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx gen_rhs_bump_G:       ',t17-t16    
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx solve_varepsilon(2):  ',t18-t17 
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx estimateerror dealloc:',t20-t19   
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx comp_adj_derivs:      ',t21-t20   
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx total:                ',t21-t0

    end subroutine em2dkx_printTimers

!==================================================================================================================================!
!============================================================================================================ localRefinement_em2dkx
!==================================================================================================================================!
    subroutine localRefinement_em2dkx

    implicit none
    
    integer                 :: iTx, iRx, i,e, n(3), inele(1), newnodes, oldnodes
    real(8)                 :: area, skinDepth,r, rmin, scaling, rminRx,rminTx 
    real(8), dimension(3)   :: ye,ze 
    character(32)           :: cend  
    logical                 :: lMakeNewMesh 
    real(kdkind)            :: qv(2) 
    type(kdtree2_result)    :: results(1)   
    type(kdtree2), pointer  :: treeTx, treeRx, treeNodes      
    real(8), dimension(:), allocatable  :: sigTemp
    complex(8)                          :: sigx,sigy,sigz
    character(256)      :: cadapt, filename          
    real(8), dimension(:,:), allocatable :: yz
    real(8) :: tstart,tend
            

    integer*8, save :: localrefinement = -1
	
    call SCOREP_F_Begin(localrefinement,"localrefinement",0,"",1424) 

            
    if (lprintDebug_em2dkx) write(*,*) myID,': entered localRefinement_em2dkx...'            
    
    call cpu_time(tstart)
 

! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0) 
    
    oldnodes = mesh%nnod  
    
!
! Create some kdtrees of the Tx's and Rx's for fast closest point searches:
!
    if (nRx > 1) then
        allocate(yz(2,nRx))  
        yz(1,:) = yRx
        yz(2,:) = zRx 
        treeRx => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
        deallocate(yz)
        nullify(treeRx%the_data)  
    endif    
    if (nTx > 1) then
        allocate(yz(2,nTx))   
        yz(1,:) = yTx
        yz(2,:) = zTx   
        treeTx => kdtree2_create( yz, sort=.true.,rearrange=.true.)  
        deallocate(yz)
        nullify(treeTx%the_data)  
    endif    
 
!
! Loop on a priori refinement until all elements are smaller than requested tolerances:
!
    do  
    
        lMakeNewMesh = .false.

        allocate ( mesh%area(mesh%nele) )
        
        mesh%area  = -1d0
        
        if (lprintDebug_em2dkx) write(*,*) 'local refinement loop, mesh%nnod: ',mesh%nnod
        
      
!
! Refine elements containing receivers:
!
        allocate (sigTemp(mesh%nele))
                    
        do e = 1,mesh%nele
            call getSigs(e,sigx,sigy,sigz)
            sigTemp(e) = dble(sigx+sigy+sigz)/3d0
        enddo

        allocate (node2tri(mesh%nnod))
        call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)                
 
        allocate(yz(2,size(mesh%y)))   
        yz(1,:) = mesh%y
        yz(2,:) = mesh%z    
        treeNodes => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
        deallocate(yz)        
        nullify(treeNodes%the_data) 
                  
        do iRx = 1,nRx
  
            call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
            & sigTemp, 1, yRx(iRx), zRx(iRx), inele)
 
            call getArea(inele(1),area,0) 
 
            call getSkinDepth(inele(1),skinDepth,0,.false.)
            
            if ( area >  nRxEleArea ) then  
!if ( sqrt(area) > skinDepth/nRxEleSkinDepths ) then
                mesh%area(inele) = area/4d0  
                lMakeNewMesh= .true.  
                
            endif
            
        enddo      
    
        deallocate (sigTemp,node2tri)
        call kdtree2_destroy(treeNodes)
                   
          
        do e = 1,mesh%nele
  
            call getArea(e,area,0) 
            call getSkinDepth(e,skinDepth,0,.false.)
 
       
            rmin = 1d100
            
            rminRx = 1d100
            rminTx = 1d100
            
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)    

            do i = 1,3

                qv(1) = ye(i)
                qv(2) = ze(i)
                            
!
! Find closest Rx:
!
                if (nRx > 1) then
                    call kdtree2_n_nearest(tp=treeRx,qv=qv,nn=1,results=results) 
                else
                    results(1)%idx = 1
                endif
                
                r = sqrt( (yRx(results(1)%idx) - ye(i))**2 + (zRx(results(1)%idx) - ze(i))**2)
            
                if (r < rminRx) rminRx = r  
            
!
! Find closest Tx:
!
                if (sType == 'cs') then
                    if (nTx > 1) then
                        call kdtree2_n_nearest(tp=treeTx,qv=qv,nn=1,results=results) 
                    else
                        results(1)%idx = 1
                    endif                
                    r = sqrt( (yTx(results(1)%idx) - ye(i))**2 + (zTx(results(1)%idx) - ze(i))**2)
            
                    if (r < rminTx) rminTx = r  
                
                 endif                  
                 
            enddo
            
!
! Refinement in a cloud surrounding the Rx's and Tx's:
!
           rmin = min(rminRx,rminTx)
           
           if  ( (rmin/(1.0*skinDepth)) > 50.)  then  ! catch to avoid overflow
               scaling = exp(50.)
           elseif (rmin < (1.0*skinDepth)) then
                scaling = 1d0
           else
               scaling = exp(rmin/(1.0*skinDepth)-1d0)
           endif

           if ( sqrt(area) > (  skinDepth*nSkinDepthsRegional*scaling) ) then
               mesh%area(e) = area/4d0
               lMakeNewMesh= .true.  
           endif  
            
              
        enddo ! Loop over elements

        if (sType == 'cs') then
        
            allocate(yz(2,size(mesh%y)))   
            yz(1,:) = mesh%y
            yz(2,:) = mesh%z    
            treeNodes => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
            deallocate(yz)       
            nullify(treeNodes%the_data) 
            
            allocate (sigTemp(mesh%nele))
             
            do e = 1,mesh%nele
                call getSigs(e,sigx,sigy,sigz)
                sigTemp(e) = dble(sigx+sigy+sigz)/3d0
            enddo
    
            allocate (node2tri(mesh%nnod))
    
            call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)
    
            do iTx = 1,nTx
 
                call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
                & sigTemp, 1, yTx(iTx), zTx(iTx), inele)
                
                call getArea(inele(1),area,0) 
                if ( area > nTxEleArea) then
                
                    mesh%area(inele) = area/4.0
                    lMakeNewMesh= .true.  
                    
                endif
                
            enddo
            
            deallocate(node2tri,sigTemp)           
            call kdtree2_destroy(treeNodes)
 
        
        endif
             
!
! Refine the mesh:
!
        if (lMakeNewMesh ) then        
            call call_triangle(retricommand, mesh )
            if (allocated(mesh%area))   deallocate(mesh%area)     
            
! Smooth the new mesh:
            call smooth_TriMesh(mesh,3,.false.)

        endif
         
        if (.not.lMakeNewMesh) then
            if (lprintDebug_em2dkx) write(*,*) 'a priori refinement: ',mesh%nnod, mesh%nele
            exit       
        endif    
    enddo
    newnodes    = mesh%nnod  
    if (newnodes > oldnodes) meshnumber  = meshnumber + 1       
               
   
    if (lSaveMeshFiles)   then  
        write(cadapt,'(i6)') meshnumber  
        cadapt = adjustl(cadapt)    
        filename  = trim(fileroot)//'.'//trim(cadapt)
        call write_triangle(mesh,filename)
    endif   

    
    if (allocated(mesh%area))   deallocate(mesh%area) 
 
    if (nRx > 1) call kdtree2_destroy(treeRx)  
    if (nTx > 1) call kdtree2_destroy(treeTx)
    
    call cpu_time(tend)
    if (lprinttimers)  write(*,'(i4,1x,a32,1x,f8.4)') myID,' em2dkx localRefinement_em2dkx:    ',tend-tstart
 
    if (lDisplayRefinementStats) call displayRefinementStats_em2dkx(iRefinementGrp,isubset,meshnumber,.false., &
                                             & oldnodes,newnodes,-2d0,-2d0, tend-tstart)
                                             

    call SCOREP_F_RegionEnd( compute ) !end Instrumentation

                                        
    end subroutine localRefinement_em2dkx
  
!==================================================================================================================================!
!====================================================================================================================== getSkinDepth
!==================================================================================================================================!
    subroutine getSkinDepth(e,skinDepth,iopt,lTx)
    
    integer, intent(in)   :: e
    integer, intent(in)   :: iopt
    logical, intent(in), optional :: lTx
    real(8), intent(out)  :: skinDepth 
    
    complex(8)            :: sig,sigx,sigy,sigz
         
!
! Get conductivity of element e:
!
    call getSigs(e,sigx,sigy,sigz,iopt)
    sig = maxval([dble(sigx),dble(sigy),dble(sigz)])
       
    skinDepth = sqrt(2d0 / (sig*MU0*w) )  
    
    end subroutine getSkinDepth

!==================================================================================================================================!
!=========================================================================================================================== getArea
!==================================================================================================================================!
    subroutine getArea(e,area,iopt)
    
    integer, intent(in)   :: e
    real(8), intent(out)  :: area
    integer, intent(in)    :: iopt
        
    integer               :: n(3)  
    real(8), dimension(3) :: ye,ze,a,b,c
    
!
! Area of element e:
!
    if (iopt == 1) then
        n    = newmesh%emap(1:3,e)
        ye   = newmesh%y(n)
        ze   = newmesh%z(n)              
    else
        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)   
    endif  
    call get_abc_coeffs(ye,ze,a,b,c,area)    
    
    
    end subroutine getArea
           
!==================================================================================================================================!
!================================================================================================================ get_finite_dipoles
!==================================================================================================================================!
    subroutine get_finite_dipoles 
!
! Creates array of point dipole located at the Gauss Legendre integration points to approximate a finite dipole
!
    integer                     :: iTx,iRx,j
    real(8)                     :: dip,azm,xx
 
        
    if (lprintDebug_em2dkx) write(*,*) myID,': get_finite_dipoles....', nTx,nquad_Tx
    
    if (sType == 'cs') then
        
        allocate(xTxQ(nTx,nquad_Tx), yTxQ(nTx,nquad_Tx), zTxQ(nTx,nquad_Tx))
    
        do iTx = 1,nTx 
        
            if (lengthTx(iTx) == 0) then    ! only a point dipole
                xTxQ(iTx,1) = xTx(iTx)
                yTxQ(iTx,1) = yTx(iTx)
                zTxQ(iTx,1) = zTx(iTx)
            
            else ! a finite dipole:
 
                dip = dipTx(iTx) * deg2rad
                azm = azimuthTx(iTx) * deg2rad

                do j = 1,nquad_Tx

                    xx = quad_xint_Tx(j)*lengthTx(iTx)/2d0 ! xint goes from -1 to 1, so scale to position point dipole array from -l/2 to +l/2

                    xTxQ(iTx,j)       = xTx(iTx) + (xx*cos(dip) * cos(azm))   
                    yTxQ(iTx,j)       = yTx(iTx) + (xx*cos(dip) * sin(azm))
                    zTxQ(iTx,j)       = zTx(iTx) + (xx*sin(dip))
                
                enddo
            
            endif
        enddo           
    endif
    
    allocate(xRxQ(nRx,nquad_Rx), yRxQ(nRx,nquad_Rx), zRxQ(nRx,nquad_Rx))
    
    do iRx = 1,nRx 
        
     if (lengthRx(iRx) == 0) then    ! only a point dipole
            xRxQ(iRx,1) = xRx(iRx)
            yRxQ(iRx,1) = yRx(iRx)
            zRxQ(iRx,1) = zRx(iRx)
            
        else ! a finite dipole:
 
            dip = dipRx(iRx) * deg2rad
            azm = azimuthRx(iRx) * deg2rad

            do j = 1,nquad_Rx

                xx = quad_xint_Rx(j)*lengthRx(iRx)/2d0 ! xint goes from -1 to 1, so scale to position point dipole array from -l/2 to +l/2

                xRxQ(iRx,j)       = xRx(iRx) + (xx*cos(dip) * cos(azm))   
                yRxQ(iRx,j)       = yRx(iRx) + (xx*cos(dip) * sin(azm))
                zRxQ(iRx,j)       = zRx(iRx) + (xx*sin(dip))
                
            enddo
            
        endif
    
    enddo     
    
    end subroutine get_finite_dipoles
!==================================================================================================================================!
!================================================================================================================= get_element_lists
!==================================================================================================================================!
    subroutine get_element_lists
!
! Computes the edge list, site element list and the transmitter element lists
!
! Kerry Key
! Scripps Institution of Oceanography
!
    implicit none
    
!
! Local variables:
!
    integer                             :: e, iTx, iRx, n
    real(8), dimension(:), allocatable  :: sigTemp
    complex(8)                          :: sigx,sigy,sigz
 
     
    real(8), dimension(:,:),allocatable :: yz
     
!
! Create element edge list:
!
    if (lprintDebug_em2dkx) write(*,*) myID,': computeEdges...'
    
    
    call computeEdges( mesh%nele, mesh%emap, mesh%neighborlist, edges, nedges )      
 
        
!
! Get site element lists
!
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...'
    
!
! First make the kdtree:
!
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...kdtree2',size(mesh%y),size(mesh%z)
    
    allocate(yz(2,size(mesh%y)))   
    yz(1,:) = mesh%y
    yz(2,:) = mesh%z    
    tree => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
    deallocate(yz)
    nullify(tree%the_data)
 
    
! Get node2tri array:
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...getNode2Tri'
    
    allocate (node2tri(mesh%nnod))
    
    call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)
 
!
! Now call the new jump-and-walk findElement routine:
!
    
    if (lprintDebug_em2dkx) write(*,*) myID,': getting site element lists...findElement'
    
    allocate (sigTemp(mesh%nele))
    
! We will use the average conductivity to define the element for the searches
! since if the search finds the point is incident on a node or edge, it chooses the most conductive element.
    do e = 1,mesh%nele
        call getSigs(e,sigx,sigy,sigz)
        sigTemp(e) = dble(sigx+sigy+sigz)/3d0
    enddo
    
    do iRx = 1,nRx
        n = 1
        if (lengthRx(iRx) > 0) n = nquad_Rx     
            
        call findElement(tree, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
        & sigTemp, n, yRxQ(iRx,1:n), zRxQ(iRx,1:n), eRx(iRx,1:n))
        
! write(*,*) 'iRx,eRx(iTx,1): ',iRx,eRx(iTx,1)
        
    enddo
    
    deallocate(sigTemp)    
    
!
! If CSEM computation, get transmitter element list eTx:
!
    if (sType == 'cs') then 
 
        allocate (sigTemp(mesh%nele)) 
    
        do e = 1,mesh%nele
          call getSigs(e,sigx,sigy,sigz)
          sigTemp(e) = dble(sigx+sigy+sigz)/3d0
        enddo
        
        do iTx = 1,nTx
            n = 1
            if (lengthTx(iTx) > 0) n = nquad_Tx     
            
!write(*,'(2(i3,1x),50(e8.1,1x))') iTx,n,yTxQ(iTx,1:n)
                      
            call findElement(tree,mesh%nnod,mesh%y,mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
            & sigTemp, n,yTxQ(iTx,1:n),zTxQ(iTx,1:n),eTx(iTx,1:n)) 
               
        enddo
    
        deallocate(sigTemp)
    
    endif
    
    end subroutine get_element_lists
           
!==================================================================================================================================!
!======================================================================================================================== gen_nodlab
!==================================================================================================================================!
    subroutine gen_nodlab
!
! Subroutine to label the mesh boundaries.  Assumes that the mesh
! has an exact rectangular outer boundary aligned along the coordinate
! axes!
!
!  nodlab array flags:
!   0       Mesh interior
!   1       Top edge, including left and right most corners
!   2       Bottom edge, including left and right most corners
!   3       Left edge
!   4       Right edge
!
! y is positive right, z is positive down
!
!
    implicit none
    
    real(8)  :: ymin,ymax,zmin,zmax 

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_nodlab...' 
      
    nodlab = 0 ! default is interior node
    
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)
    zmin = minval(mesh%z)
    zmax = maxval(mesh%z)
    
! First set left and right side flags:
    where (mesh%y == ymax) nodlab = 4  ! right
    where (mesh%y == ymin) nodlab = 3  ! left
    where (mesh%z == zmax) nodlab = 2  ! bottom
    where (mesh%z == zmin) nodlab = 1  ! top
   
    end subroutine gen_nodlab      


!==================================================================================================================================!
!========================================================================================================================= gen_bcnod
!==================================================================================================================================!
    subroutine gen_bcnod
!
!  Routine to remap the vector of node labels in 'nodlab' to a
!  new vector of labels bcnod that is appropriate for particular boundary
!  conditions.
!
! Mapping:
!
!  nodlab   bcnod
!   0           0   Mesh interior
!   1           1   Top edge, including left and right most corners
!   2           1   Bottom edge, including left and right most corners
!   3           1   Left edge
!   4           1   Right edge
!
!
!  Values within bcnod indicate the fate of the node they represent:
!
!    0  ...  interior node where FE solution is obtained.
!    1  ...  boundary node where Dirichlet condition is applied.
!
!
    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_bcnod...'  
      
    bcnod = 0
!
!  Remap the node labels
!
    where ( nodlab > 0) bcnod = 1
    
    end subroutine gen_bcnod
    
!==================================================================================================================================!
!=========================================================================================================================== getSigs
!==================================================================================================================================!

    subroutine getSigs(e,sigx,sigy,sigz,iopt) 
!
! Gets anisotropic conductivity parameters for element e
!

    implicit none
 
    integer, intent(in)             :: e
    integer, intent(in), optional   :: iopt
    complex(8), intent(out)         :: sigx,sigy,sigz
    
! Local variables:
    integer     :: iparamnum
    real(8)     :: sig0, eta, tau, c
    complex(8)  :: p, zrho
      
!
! Get the anisotropy type and its associated conductivities
!
    if (present(iopt)) then
        if (iopt==1) then 
            iparamnum = nint(abs(newmesh%attr(e)))
        else
            iparamnum = nint(abs(mesh%attr(e)))
        endif
    else
        iparamnum =  nint(abs(mesh%attr(e)))
    endif
    
    select case (trim(cAnisotropy))
    
    case ('isotropic')
        sigx = 1d0/rhoParams(iparamnum)
        sigy = 1d0/rhoParams(iparamnum) 
        sigz = 1d0/rhoParams(iparamnum)
        
    case ('isotropic_ip')  ! Pelton Cole-Cole:
    
        sig0 = 1d0/rhoParams(4*(iparamnum-1)+1)
        eta  = rhoParams(4*(iparamnum-1)+2)
        tau  = rhoParams(4*(iparamnum-1)+3)  
        c    = rhoParams(4*(iparamnum-1)+4)    
        
!write(*,'(4(g12.2,1x))') sig0,eta,tau,c
        p = 0d0
        if (c > 0d0) p = (-ic*w*tau)**c ! conditional to avoid case when  p is undefined due to tau=0 and c=0, i.e.  0**0
                 
        sigx = sig0*(1d0+p)/(1d0 + (1d0-eta)*p) 
        sigy = sigx
        sigz = sigx
        
!write(*,*) 'aimag(sigx): ', e, aimag(sigx)
        
    case ('isotropic_complex')  ! Complex resistivity
            
        zrho = cmplx(rhoParams(2*(iparamnum-1)+1), rhoParams(2*(iparamnum-1)+2))
        sigx = 1d0/zrho
        sigy = sigx
        sigz = sigx 
        
!write(*,*) e,dble(zrho),aimag(zrho)
                     
    case ('triaxial')
        sigx = 1d0/rhoParams(3*(iparamnum-1)+1)
        sigy = 1d0/rhoParams(3*(iparamnum-1)+2)
        sigz = 1d0/rhoParams(3*(iparamnum-1)+3)
    case ('tix')                               ! For transversely isotropic:
        sigx = 1d0/rhoParams(2*(iparamnum-1)+1)     ! First param is sig along symmetry axis
        sigy = 1d0/rhoParams(2*(iparamnum-1)+2)     ! Second param is along transverse plane
        sigz = 1d0/rhoParams(2*(iparamnum-1)+2)         
    case ('tiy')
        sigx = 1d0/rhoParams(2*(iparamnum-1)+2)
        sigy = 1d0/rhoParams(2*(iparamnum-1)+1)
        sigz = 1d0/rhoParams(2*(iparamnum-1)+2)         
    case ('tiz')
        sigx = 1d0/rhoParams(2*(iparamnum-1)+2)
        sigy = 1d0/rhoParams(2*(iparamnum-1)+2)
        sigz = 1d0/rhoParams(2*(iparamnum-1)+1)         
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in em2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
    
!write(*,*) iparamnum, 1./real(sigx),1./real(sigy),1./real(sigz)
     
    end subroutine getSigs
    
    
!==================================================================================================================================!
!======================================================================================================================= getIPderivs
!==================================================================================================================================!
    subroutine getIPderivs(e,dsds0,dsdeta,dsdtau,dsdc,sig0)
 
    implicit none
 
    integer, intent(in)           :: e
    complex(8), intent(out)       :: dsds0,dsdeta,dsdtau,dsdc
    real(8), intent(out)           :: sig0
    
! Local variables:
    integer     :: iregion
    real(8)     :: eta, tau, c
    complex(8)  :: p
      
    iregion =  nint(abs(mesh%attr(e)))
 
! Pelton Cole-Cole:
    
    sig0 = 1d0/rhoParams(4*(iregion-1)+1)
    eta  = rhoParams(4*(iregion-1)+2)
    tau  = rhoParams(4*(iregion-1)+3)  
    c    = rhoParams(4*(iregion-1)+4)   

    dsds0   = 0d0 
    dsdeta  = 0d0
    dsdtau  = 0d0 
    dsdc    = 0d0
    
    p = 1d0
    if (c > 0d0) p = (-ic*w*tau)**c   
!    sig = sig0*(1d0+p)/(1d0 + (1d0-eta)*p)
         
! sig0:
    dsds0 = (1d0+p)/(1d0 + (1d0-eta)*p) 
    
    if (tau > 0d0 ) then
        
! eta:
        dsdeta = sig0*(1d0+p)*p/(1d0 + (1d0-eta)*p)**2 
        
! tau:
        dsdtau = sig0*eta*c*(-ic*w)*(-ic*w*tau)**(c-1d0)/(1d0 + (1d0-eta)*p)**2    
        
! c:
        dsdc   = sig0*p*eta*log(-ic*w*tau)/(1d0 + (1d0-eta)*p)**2
    
    endif
  
    end subroutine getIPderivs
        
!==================================================================================================================================!
!====================================================================================================================== getParamNums
!==================================================================================================================================!
    subroutine getParamNums(e,iparamnum) 
!
! Gets parameter indices for element e
!

    implicit none
 
    integer, intent(in)  :: e
    integer, intent(out) :: iparamnum(4)
    
    integer :: iregion
    
    iparamnum = 0
    
!
! Get parameter numbers, depending on anisotropy type:
!
    iregion =  nint(abs(mesh%attr(e)))
       
    select case (trim(cAnisotropy))
    
    case ('isotropic')
        iparamnum(1) = iFreeParam(iregion)
    case ('triaxial')
        iparamnum(1) = iFreeParam(3*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(3*(iregion-1) + 2)
        iparamnum(3) = iFreeParam(3*(iregion-1) + 3)
    case ('tix','tiy','tiz','isotropic_complex')  
        iparamnum(1) = iFreeParam(2*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(2*(iregion-1) + 2)    
    case ('isotropic_ip')   ! Cole-Cole IP
        iparamnum(1) = iFreeParam(4*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(4*(iregion-1) + 2)
        iparamnum(3) = iFreeParam(4*(iregion-1) + 3)    
        iparamnum(4) = iFreeParam(4*(iregion-1) + 4)                     
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in em2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
     
     
    end subroutine getParamNums
        
!==================================================================================================================================!
!========================================================================================================================= gen_1d_mt
!==================================================================================================================================!
      subroutine gen_1d_mt(model1D,cSide,cSigDirection,maxLayers)
!
! Routine to generate 1D layered models based on the conductivity
! profile of the right side of the mesh (y-max) and another along the
! left side (y-min).  Now also computes the 1D MT layer coefficients
! at the current angular frequency "w".
!


    implicit none
  
    character(5), intent(in)  :: cSide         ! 'left' or 'right'
    character(1), intent(in)  :: cSigDirection ! 'x', 'y' for TE or TM anisotropic 1D conductivity
    type(MT1D), intent(inout) :: model1D     
    integer, optional         :: maxLayers     ! limit number of 1D layers

    real(8)                  :: yposition  ! position where 1D model to be generated (needs to be leftmost or rightmost side)
    integer                  :: nSideNodes,n(3)
    real(8), dimension(:), allocatable :: sigSide, zSide 
    
    integer     :: i,e, eSide
    real(8)     :: d1,d2
    complex(8)  :: sigx,sigy,sigz  

    
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_1d_mt ...'  
!
! First count the number of nodes on the requested side
!
    select case (trim(cSide))
    case ('left')
        yposition =  minval(mesh%y)
    case('right')
        yposition =  maxval(mesh%y)  
    case default
        write(*,*) 'Error in gen_1d_mt, undefined cSide: ', trim(cSide) 
        write(*,*) 'stopping!'
        stop
    end select
    
    nSideNodes  = count(mesh%y == yposition)
    
    if (nSideNodes <= 0) then
        write(*,*) 'Error in gen_1d_mt, number of side nodes: ', nSideNodes, ' y position for profile: ', yposition
        write(*,*) 'stopping!'
        stop
    endif
 
!
! Allocate buffer arrays:
!
    allocate(sigSide(nSideNodes), zSide(nSideNodes) )
    

!
!  Extract an unsorted list of element conductivities along the requested side of the mesh:
!
    eSide = 0
    do e=1,mesh%nele  
        n = mesh%emap(1:3,e)
        
        do i = 1,3
        
            d1 = abs(mesh%y(n(eps(i)))   - yposition) 
            d2 = abs(mesh%y(n(eps(i+1))) - yposition) 
            
            if ( (d1 < 1d-6) .and. (d2 < 1d-6) ) then
                eSide = eSide + 1
                
                call getSigs(e,sigx,sigy,sigz)
                
                select case (trim(cSigDirection))
                case('x')
                    sigSide(eSide) = dble(sigx)
                case('y')
                    sigSide(eSide) = dble(sigy)
                case default
                    write(*,*) 'Error in gen_1d_mt, undefined cSide: ', trim(cSigDirection) 
                    write(*,*) 'stopping!'
                    stop                   
                end select 
                
                zSide(eSide)   = min(mesh%z(n(eps(i))),mesh%z(n(eps(i+1))))            ! Stores top depth
            
            endif
       
        enddo
    enddo
!
! Sort the lists on z:
!
    call quick_sort_d(zSide(1:eSide),sigSide(1:eSide))
  
!
! Find layers of constant conductivity:
!
     nSideNodes = 1
     do i=2,eSide
        if ( sigSide(i) /= sigSide(nSideNodes)  ) then
            nSideNodes = nSideNodes + 1 
! Record new layer top:
            zSide(nSideNodes) = zSide(i)
            sigSide(nSideNodes) = sigSide(i)
        endif
     enddo
  
!
! Now insert these into the PW1D derived type:
!
    if (present(maxLayers))   nSideNodes = min(maxLayers,nSideNodes)
    
    model1D = new_MT1D(zSide(1:nSideNodes),sigSide(1:nSideNodes),w)
 
!call print(model1D)
           
!
! Deallocate local arrays:
!
    deallocate( zSide, sigSide )
      
    end subroutine gen_1d_mt
!==================================================================================================================================!
!====================================================================================================================== getStiff_lin
!==================================================================================================================================!
    subroutine getStiff_lin(e,mtx)
!
! Generates the 6 x 6 nodal linear FE stiffness matrix for the 2.5D EM problem.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 21, 2011: Modified to handle triaxial anisotropy, woo-hoo!
!
      
    implicit none
    
    integer, intent(in) :: e
    complex(8),dimension(21), intent(out) :: mtx   

! Local variables:
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze
    complex(8)              :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2
    real(8), dimension(3)   :: a,b,c
    real(8)                 :: ar
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms
  
!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / (gammay2*4d0*ar)
    aez = (sigz - iomeps) / (gammaz2*4d0*ar)
    bex = (sigx - iomeps) * ar / 6d0
    cey =  ikx / (gammay2*4d0*ar)
    cez =  ikx / (gammaz2*4d0*ar)

    ahz = -ommu / (gammay2*4d0*ar)    ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / (gammaz2*4d0*ar)    ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 6d0
    chy = -cez 
    chz = -cey 

!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
! E1
    mtx(1)  = aey*b(1)*b(1) + aez*c(1)*c(1) + bex      ! e1 e1
    mtx(2)  = cez*b(1)*c(1) - cey*c(1)*b(1)            ! e1 h1 not zero when anisotropic!
    mtx(3)  = aey*b(2)*b(1) + aez*c(2)*c(1) + bex/2d0  ! e1 e2
    mtx(4)  = cez*b(2)*c(1) - cey*c(2)*b(1)            ! e1 h2
    mtx(5)  = aey*b(3)*b(1) + aez*c(3)*c(1) + bex/2d0  ! e1 e3
    mtx(6)  = cez*b(3)*c(1) - cey*c(3)*b(1)            ! e1 h3
! H1
    mtx(7)  = ahy*b(1)*b(1) + ahz*c(1)*c(1) + bhx      ! h1 h1
    mtx(8)  = chz*b(2)*c(1) - chy*c(2)*b(1)            ! h1 e2
    mtx(9)  = ahy*b(2)*b(1) + ahz*c(2)*c(1) + bhx/2d0  ! h1 h2
    mtx(10) = chz*b(3)*c(1) - chy*c(3)*b(1)            ! h1 e3
    mtx(11) = ahy*b(3)*b(1) + ahz*c(3)*c(1) + bhx/2d0  ! h1 h3
! E2
    mtx(12) = aey*b(2)*b(2) + aez*c(2)*c(2) + bex      ! e2 e2
    mtx(13) = cez*b(2)*c(2) - cey*c(2)*b(2)            ! e2 h2 not zero when anisotropic!
    mtx(14) = aey*b(3)*b(2) + aez*c(3)*c(2) + bex/2d0  ! e2 e3
    mtx(15) = cez*b(3)*c(2) - cey*c(3)*b(2)            ! e2 h3
!H2
    mtx(16) = ahy*b(2)*b(2) + ahz*c(2)*c(2) + bhx      ! h2 h2
    mtx(17) = chz*b(3)*c(2) - chy*c(3)*b(2)            ! h2 e3
    mtx(18) = ahy*b(3)*b(2) + ahz*c(3)*c(2) + bhx/2d0  ! h2 h3
!E3
    mtx(19) = aey*b(3)*b(3) + aez*c(3)*c(3) + bex      ! e3 e3
    mtx(20) = cez*b(3)*c(3) - cey*c(3)*b(3)            ! e3 h3 not zero when anisotropic!
!H3
    mtx(21) = ahy*b(3)*b(3) + ahz*c(3)*c(3) + bhx      ! h3 h3

    end subroutine getStiff_lin

!==================================================================================================================================!
!===================================================================================================================== getStiff_bump
!==================================================================================================================================!
    subroutine getStiff_bump(e,mtx)
!
! Routine to compute the element stiffness matrix for the coupled 2.5D
! CSEM system in the space of quadratic edge bump functions.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! June 2008:        Created
! April 21, 2011:   Modified to handle triaxial anisotropy. Now has cross
!                   coupling terms between E-H that vanish if isotropic.
       
    implicit none
    
    integer, intent(in)     :: e
    complex(8),intent(out)  :: mtx(21)   

! Local variables:
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze
    complex(8)              :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2
    real(8), dimension(3)   :: a,b,c
    real(8)                 :: ar, ar1o3,ar2o3
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms

!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / gammay2
    aez = (sigz - iomeps) / gammaz2
    bex = (sigx - iomeps) * ar / 45d0
    cey =  ikx / gammay2  
    cez =  ikx / gammaz2

    ahz = -ommu / gammay2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / gammaz2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 45d0
    chy = -cez 
    chz = -cey 
   
    ar2o3 = 2d0/(3d0*ar)
    ar1o3 = 1d0/(3d0*ar)
    
!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
! E1
    mtx(1)  = ar2o3*(aey*(b(1)*b(1) - b(2)*b(3)) + aez*(c(1)*c(1) - c(2)*c(3)) ) + bex*8d0 ! e1 e1
    mtx(2)  = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e1 h1, not zero when anisotropic!
    mtx(3)  = ar2o3*(aey*b(2)*b(1) + aez*c(2)*c(1) )                             + bex*4d0 ! e1 e2
    mtx(4)  = ar1o3*(cez-cey)*(b(1)*c(2) + b(2)*c(1))                                      ! e1 h2
    mtx(5)  = ar2o3*(aey*b(1)*b(3) + aez*c(1)*c(3) )                             + bex*4d0 ! e1 e3
    mtx(6)  = ar1o3*(cez-cey)*(b(1)*c(3) + b(3)*c(1))                                      ! e1 h3
! H1
    mtx(7)  = ar2o3*(ahy*(b(1)*b(1) - b(2)*b(3)) + ahz*(c(1)*c(1) - c(2)*c(3)) ) + bhx*8d0 ! h1 h1
    mtx(8)  = ar1o3*(chz-chy)*(b(2)*c(1) + b(1)*c(2))                                      ! h1 e2
    mtx(9)  = ar2o3*(ahy*b(2)*b(1) + ahz*c(2)*c(1) )                             + bhx*4d0 ! h1 h2
    mtx(10) = ar1o3*(chz-chy)*(b(3)*c(1) + c(3)*b(1))                                      ! h1 e3
    mtx(11) = ar2o3*(ahy*b(3)*b(1) + ahz*c(3)*c(1) )                             + bhx*4d0 ! h1 h3
! E2
    mtx(12) = ar2o3*(aey*(b(2)*b(2) - b(1)*b(3)) + aez*(c(2)*c(2) - c(1)*c(3)) ) + bex*8d0 ! e2 e2
    mtx(13) = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e2 h2, not zero when anisotropic!
    mtx(14) = ar2o3*(aey*b(3)*b(2) + aez*c(3)*c(2))                              + bex*4d0 ! e2 e3
    mtx(15) = ar1o3*(cez-cey)*(b(3)*c(2) + c(3)*b(2))                                      ! e2 h3
!H2
    mtx(16) = ar2o3*(ahy*(b(2)*b(2) - b(1)*b(3)) + ahz*(c(2)*c(2) - c(1)*c(3)) ) + bhx*8d0 ! h2 h2
    mtx(17) = ar1o3*(chz-chy)*(b(3)*c(2) + c(3)*b(2) )                                     ! h2 e3
    mtx(18) = ar2o3*(ahy*b(3)*b(2) + ahz*c(3)*c(2) )                             + bhx*4d0 ! h2 h3
!E3
    mtx(19) = ar2o3*(aey*(b(3)*b(3) - b(1)*b(2)) + aez*(c(3)*c(3) - c(1)*c(2)) ) + bex*8d0 ! e3 e3
    mtx(20) = ar1o3*(cez-cey)*(b(1)*c(1) + b(2)*c(2) + b(3)*c(3) )                         ! e3 h3, not zero when anisotropic!
!H3
    mtx(21) = ar2o3*(ahy*(b(3)*b(3) - b(1)*b(2)) + ahz*(c(3)*c(3) - c(1)*c(2)) ) + bhx*8d0 ! h3 h3
        
 
    end subroutine getStiff_bump

!==================================================================================================================================!
!================================================================================================================  getStiff_lin_bump
!==================================================================================================================================!
    subroutine getStiff_lin_bump(e,mtx)
!
! Routine to compute the element stiffness matrix for the coupled 2.5D
! CSEM system for B(u,v) where u is linear and v is bump space.
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011:   Created
       
    implicit none
     
    integer, intent(in)    :: e
    complex(8),intent(out) :: mtx(21)   

! Local variables:
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze,a,b,c
    real(8)                 :: ar, aro3 
    complex(8)              :: gammay2,gammaz2, sigx,sigy,sigz
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms

!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
    gammay2 = kx2-ommu*(sigy - iomeps)
    gammaz2 = kx2-ommu*(sigz - iomeps)
        
!
! Get Linear basis Coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D EM  coefficients for triaxial conductivity:
!
    aey = (sigy - iomeps) / gammay2
    aez = (sigz - iomeps) / gammaz2
    bex = (sigx - iomeps) * ar / 15d0
    cey =  ikx / gammay2  
    cez =  ikx / gammaz2

    ahz = -ommu / gammay2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    ahy = -ommu / gammaz2     ! Note the negative sign for H gives E-H coupling terms matrix symmetry
    bhx = -ommu * ar / 15d0
    chy = -cez 
    chz = -cey 
   
    aro3 = 1d0/(3d0*ar)
    
!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
! E1
    mtx(1)  = -aro3*(aey*b(1)*b(1) + aez*c(1)*c(1))  + bex      ! e1 e1
    mtx(2)  =  aro3*(cey*b(1)*c(1) - cez*b(1)*c(1))             ! e1 h1, not zero when anisotropic!
    mtx(3)  = -aro3*(aey*b(1)*b(2) + aez*c(2)*c(1))  + bex*2d0  ! e1 e2
    mtx(4)  =  aro3*(cey*b(1)*c(2) - cez*b(2)*c(1))             ! e1 h2
    mtx(5)  = -aro3*(aey*b(1)*b(3) + aez*c(3)*c(1))  + bex*2d0  ! e1 e3
    mtx(6)  =  aro3*(cey*b(1)*c(3) - cez*b(3)*c(1))             ! e1 h3
! H1
    mtx(7)  = -aro3*(ahy*b(1)*b(1) + ahz*c(1)*c(1))  + bhx      ! h1 h1
    mtx(8)  =  aro3*(chy*b(1)*c(2) - chz*b(2)*c(1))             ! h1 e2
    mtx(9)  = -aro3*(ahy*b(1)*b(2) + ahz*c(2)*c(1))  + bhx*2d0  ! h1 h2
    mtx(10) =  aro3*(chy*b(1)*c(3) - chz*b(3)*c(1))             ! h1 e3
    mtx(11) = -aro3*(ahy*b(1)*b(3) + ahz*c(3)*c(1))  + bhx*2d0  ! h1 h3
! E2
    mtx(12) = -aro3*(aey*b(2)*b(2) + aez*c(2)*c(2))  + bex      ! e2 e2
    mtx(13) =  aro3*(cey*b(2)*c(2) - cez*b(2)*c(2))             ! e2 h2, not zero when anisotropic!
    mtx(14) = -aro3*(aey*b(2)*b(3) + aez*c(3)*c(2))  + bex*2d0  ! e2 e3
    mtx(15) =  aro3*(cey*b(2)*c(3) - cez*b(3)*c(2))             ! e2 h3
!H2
    mtx(16) = -aro3*(ahy*b(2)*b(2) + ahz*c(2)*c(2))  + bhx      ! h2 h2
    mtx(17) =  aro3*(chy*b(2)*c(3) - chz*b(3)*c(2))             ! h2 e3
    mtx(18) = -aro3*(ahy*b(2)*b(3) + ahz*c(3)*c(2))  + bhx*2d0  ! h2 h3
!E3
    mtx(19) = -aro3*(aey*b(3)*b(3) + aez*c(3)*c(3))  + bex      ! e3 e3
    mtx(20) =  aro3*(cey*b(3)*c(3) - cez*b(3)*c(3))             ! e3 h3, not zero when anisotropic!
!H3
    mtx(21) = -aro3*(ahy*b(3)*b(3) + ahz*c(3)*c(3))  + bhx      ! h3 h3
        
    end subroutine getStiff_lin_bump
     
 
!==================================================================================================================================!
!=========================================================================================================================== gen_rhs
!==================================================================================================================================!
    subroutine gen_rhs 
!
! Generates the rhs vectors for CSEM and MT problems using optional
! methods

    if (sType == 'mt') then ! MT using decoupled linear systems:
        call gen_rhs_mt  
    else
        call gen_rhs_dipole
    endif    
 

         
    end subroutine gen_rhs          
!==================================================================================================================================!
!==================================================================================================================== gen_rhs_dipole
!==================================================================================================================================!
    subroutine  gen_rhs_dipole
!
! Routine to compute the RHS of the FEM system of equations using
! a total field formulation for an electric dipole.
!
! Kerry Key
! Scripps Institution of Oceanography
!
    implicit none 
    
    integer                 :: i,n(3),e, iTx, isign, iq, nq
    complex(8)              :: sigx,sigy,sigz
    complex(8)              :: gammay2,gammaz2,expikdx        
    real(8), dimension(3)   :: l,dldy,dldz, ye,ze   
    complex(8),dimension(3) :: ee,hh  

    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_dipole ...'  

!
!  Initialize the RHS
!
    rhs = (0d0,0d0)     
   
!
! Loop over transmitters:
!
    do iTx = 1,nTx
    
        nq = 1
        if (lengthTx(iTx) > 0 ) nq = nquad_Tx
        
        do iq = 1,nq
    
!
! Get a single element that has the transmitter:
!
            e = eTx(iTx,iq)
        
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)              
    
!
! Set up the basis functions:
!
            n   = mesh%emap(1:3,e)
            ye  = mesh%y(n)
            ze  = mesh%z(n)

!
! Evaluate the basis functions and their derivatives at the source location:
!
            call get_linear_basis(ye,ze,yTxQ(iTx,iq),zTxQ(iTx,iq),l,dldy,dldz)    
        
            isign =  1
            if (lNeg_ikx(iTx)) isign = -1
!
! Apply coefficients:
!
            ee =    -l*jx(iTx) - isign*ikx  * ( dldy*jy(iTx)/gammay2      + dldz*jz(iTx)/gammaz2 )  & 
               &               - ommu * ( dldy*mz(iTx)*(sigy - iomeps)/gammay2 - dldz*my(iTx)*(sigz - iomeps)/gammaz2 )  !ommu for H src
            hh = -(                              - ommu * ( dldy*jz(iTx)/gammaz2      - dldz*jy(iTx)/gammay2 )  &
               &    -l*mx(iTx)*ommu -  isign*ikx * ommu * ( dldy*my(iTx)/gammaz2      + dldz*mz(iTx)/gammay2 )  )  
! Note minus sign applied to hh in order to make the matrix symmetric

            ee = isign*ee ! for neg kx system
!
! Insert these into the correct rhs location
!
            if (nq > 1) then
                expikdx = exp(-isign*ikx*(xTxQ(iTx,iq)-xTx(iTx))) 
                rhs(2*n-1, iTx) = rhs(2*n-1, iTx) + expikdx*ee*quad_weights_Tx(iq)/2.d0    ! Ex
                rhs(2*n  , iTx) = rhs(2*n  , iTx) + expikdx*hh*quad_weights_Tx(iq)/2.d0    ! Hx
            else ! point dipole
                rhs(2*n-1, iTx) = ee ! Ex
                rhs(2*n  , iTx) = hh ! Hx
            endif           
        enddo
    
    enddo

! Lastly, insert the dirichlet boundary condition into the entries
! of the RHS whose nodes lie on the border.
!
! In the total field formulation this isn't necessary unless somebody
! mistakenly places a Tx on an exterior node boundary, but then the solution
! is going to be incorrect anyway...
!
    do i=1,mesh%nnod
        if(bcnod(i).ne.0) then
          rhs(2*i-1,1:nTx) = (0d0, 0d0)
          rhs(2*i  ,1:nTx) = (0d0, 0d0)
        endif
    enddo    
    
    end subroutine gen_rhs_dipole
!==================================================================================================================================!
!======================================================================================================================== gen_rhs_mt
!==================================================================================================================================!
    subroutine  gen_rhs_mt
 
    implicit none 
    
    integer                     :: i,e, n(3)
    integer,dimension(3,2)      :: imtxE, imtxH ! maps for BCs
    real(8)                     :: ymin, ymax, ediff, eta, y,z     
    complex(8)                  :: sigx,sigy,sigz, dsigTE,dsigTM
    real(8)                     :: ye(3),ze(3),a(3),b(3),c(3),area,l(3),dldy(3),dldz(3)
    complex(8)                  :: elTE,erTE,hlTM,hrTM,comboTE,comboTM
    complex(8)                  :: exTE, hxTM 
    complex(8),dimension(21)    :: mtx   
  
    integer                  :: nqp, k,nq_order
    real(8), dimension(3,50) :: qp  
    real(8), dimension(50)   :: qwt 
    complex(8), dimension(3) :: qsum
    real(8)                  :: zq,bary(3)
 
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_mt ...'
    
!
! Set up mapping arrays for pulling out parts of 6 x 6 symmetric stiffness matrix stored as 21 upper triangular coeffs:
!
    imtxE(1,1) = 3   ! node 1 on boundary, select inner product with nodes 2 and 3
    imtxE(1,2) = 5
    imtxE(2,1) = 14  ! node 2 on boundary, select inner product with nodes 3 and 1
    imtxE(2,2) = 3 
    imtxE(3,1) = 5   ! node 3 on boundary, select inner product with nodes 1 and 2
    imtxE(3,2) = 14    
    imtxH(1,1) = 9
    imtxH(1,2) = 11
    imtxH(2,1) = 18
    imtxH(2,2) = 9 
    imtxH(3,1) = 11
    imtxH(3,2) = 18       
  
!
! Get quadrature points for scattered field integration:
!
    if (lMTscatteredField) then   
        nq_order = 10
        call triquad(nq_order,qp,qwt,nqp)  
    endif   
!
!  Initialize the RHS
!
    rhs = (0d0,0d0)     
!
!  Compute the min and max x-values.
!
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)       
      
!
!  Loop over elements.
!
    do e=1,mesh%nele
    
        n = mesh%emap(1:3,e)
        
        if (lMTscatteredField) then
  
! Get basis coefficients and area of element e:
            ye   = mesh%y(n)
            ze   = mesh%z(n)
            call get_abc_coeffs(ye,ze,a,b,c,area)   
                      
! Locate corresponding 1D layer for TE and TM models:
            z = sum(ze)/3d0
            y = sum(ye)/3d0
            
! Get scattering conductivities:
            call getSigs(e,sigx,sigy,sigz) ! get anisotropic conductivity terms for element e
       
            
!
! If scattered field formulation, insert current sources:
!
            
! here i'm using quadrature to deal with non-constant dsig if
! 1D layer cuts through element e so dsig is a 2D step function.
! order 10 quadrature seems to work really well and doesn't take much time.
        
! TE
            qsum = 0d0  
            do k = 1,nqp
                bary = qp(1:3,k)
                zq   = sum(bary*ze) 
                dsigTE = sigx - Scat1DTE%getSig(zq)
                if (dsigTE /= 0d0 ) then
                    qsum = qsum + dsigTE*Scat1DTE%getEH1D(zq,'e')*bary*qwt(k)
                endif
            enddo 
            rhs(2*n-1,1) = rhs(2*n-1,1) - qsum*area ! TE is easy since n(1:3) get same source contributions
            
! TM:
            qsum = 0d0  
            do k = 1,nqp
                bary = qp(1:3,k)
                zq   = sum(bary*ze) 
                dsigTM = sigy - Scat1DTM%getSig(zq) 
                if (dsigTM /= 0d0 ) then
                    qsum = qsum + dsigTM/sigy*Scat1DTM%getEH1D(zq,'e')*qwt(k)   ! note using -ey here, but sign cancels with coeffs
                endif
            enddo 
                
            call get_linear_basis(ye,ze,y,z,l,dldy,dldz)            
            rhs(2*n  ,1) = rhs(2*n  ,1) - qsum*dldz*area ! TM part, note dldz(1:3) for n(1:3) ! negative sign to match -sign in system matrix for symmetry in CSEM coupled system
         
    
        endif ! lMTscatteredField
        
!
! Boundary conditions:  left and right side solution fixed to 1D solutions:
!
!    Total field formulation:
!
!       u(y_left,z)  = e1dleft(z)
!       u(y_right,z) = e1dright(z)
!       u(y,z_top) and  u(y,z_bot) get a cosine taper between the left and right side values at the top and bottom.
!
!    Scattered field formulation where u_total = u1d + u_scat :
!
!       u_total(y_left,z)  =  e1dleft(z)  = u_1d(z) + u_scat(y_left,z)
!       u_total(y_right,z) =  e1dright(z) = u_1d(z) + u_scat(y_right,z)
!
!       so therefore the boundary conditions for the scattered field are then:
!
!       u_scat(y_left,z)  = e1dleft(z)  - u_1d(z)
!       u_scat(y_right,z) = e1dright(z) - u_1d(z)
!
!       and along the top and bottom we use a cosine taper between the left and right side boundary conditions:
!
!       eta    = (y-ymin)/(ymax-ymin)
!       ediff = (1d0-cos(PI*eta))/2d0     ! ediff = 0 at y_left and ediff = 1 at y_right
!
!       u_scat(y,z_top)  = (1-ediff)*u_scat(y_left,z_top) + ediff*u_scat(y_right,z_top)
!                        = (1-ediff)*(e1dleft(z)  - u_1d(z))) + ediff*(e1dright(z) - u_1d(z))
!                        = (1-ediff)*e1dleft(z) + ediff*e1dright(z)  - u_1d(z)
!                        = ediff*(e1dright(z)-e1dleft(z)) + e1dleft(z)  - u_1d(z)
!
!       and similar for z_bottom. Note that for the total field formulation, u_1d = 0 in the above equation.
!
!      That was easy!
!
        
            
!
! Loop over nodes, if BC(n(i)), then subtract off contribution to n(eps(i+1)) and n(eps(i+2)):
!
        do i=1,3
        
            if (bcnod(n(i)) == 1) then 
            
! Get stiffness matrix:
                call getStiff_lin(e, mtx)
                
                z =  mesh%z(n(i))
                
! Get the 1D MT fields
                elTE = Left1DTE%getEH1D( z,'e')
                erTE = Right1DTE%getEH1D(z,'e')
                hlTM = Left1DTM%getEH1D( z,'h')
                hrTM = Right1DTM%getEH1D(z,'h')
                                
                exTE = 0d0 
                hxTM = 0d0
                
                if (lMTscatteredField) then                 
                    exTE = Scat1DTE%getEH1D(z,'e')            
                    hxTM = Scat1DTM%getEH1D(z,'h')
                endif

                eta = (mesh%y(n(i))-ymin)/(ymax-ymin)            
                ediff = (1d0-dcos(PI*eta))/2d0    
                comboTE = (erTE - elTE)*ediff + elTE - exTE
                comboTM = (hrTM - hlTM)*ediff + hlTM - hxTM

                
! Now the hard part is subtracting this off the RHS appropriately:
! remembering that there is no Ex-Hx coupling for 2D MT...
! TE Ex part:
                rhs(2*n(eps(i+1))-1,1) = rhs(2*n(eps(i+1))-1,1) - mtx(imtxE(i,1))*comboTE
                rhs(2*n(eps(i+2))-1,1) = rhs(2*n(eps(i+2))-1,1) - mtx(imtxE(i,2))*comboTE

! TM Hx part:
                rhs(2*n(eps(i+1))  ,1) = rhs(2*n(eps(i+1)) ,1 ) - mtx(imtxH(i,1))*comboTM
                rhs(2*n(eps(i+2))  ,1) = rhs(2*n(eps(i+2)) ,1 ) - mtx(imtxH(i,2))*comboTM           
    
            endif
        enddo  ! i = 1,3
    enddo  ! e = 1,nele

!
! Lastly, we need to set the boundary node values:
!
    do i = 1,mesh%nnod
        if (bcnod(i) == 1) then
            
            z = mesh%z(i)
            
! Get the 1D MT fields:
            elTE = Left1DTE%getEH1D( z,'e')
            erTE = Right1DTE%getEH1D(z,'e')
            hlTM = Left1DTM%getEH1D( z,'h')
            hrTM = Right1DTM%getEH1D(z,'h')            

            exTE = 0d0 
            hxTM = 0d0
            
            if (lMTscatteredField) then
                exTE = Scat1DTE%getEH1D(z,'e')
                hxTM = Scat1DTM%getEH1D(z,'h')
            endif
         
! Cosine taper between left and right sides:
            eta = (mesh%y(i)-ymin)/(ymax-ymin)            
            ediff = (1d0-dcos(PI*eta))/2d0        
            comboTE = (erTE - elTE)*ediff + elTE - exTE  
            comboTM = (hrTM - hlTM)*ediff + hlTM - hxTM
            
            rhs(2*i - 1 ,1) = comboTE
            rhs(2*i     ,1) = comboTM
                
        endif
    enddo
    
    end subroutine gen_rhs_mt

!==================================================================================================================================!
!==================================================================================================================== gen_rhs_bump_f
!==================================================================================================================================!
    subroutine  gen_rhs_bump_f
!
! Routine to compute F(v)  for all v in W_h (bump space)
!
! Kerry Key
! Scripps Institution of Oceanography
!
! April 2011    Added support for triaxial anisotropy
! Sept 2009     Added support for 2DMT total field rhs, super difficult to code, see below...
! June 2008
!

    implicit none

!
! Local variables:
!
    integer               :: n(3),e,iTx,j, iedge, isign, iq, nq
    real(8), dimension(3) :: q,dqdy,dqdz,ye,ze  
    complex(8)            :: sigx,sigy,sigz
    complex(8)            :: gammay2,gammaz2,expikdx,ee(3),hh(3)       
 
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_f ...'  
 
    do iTx = 1,nTx
    
        if (sType == 'cs') then ! CSEM E dipole, x,y,z or yz

            nq = 1
            if (lengthTx(iTx) > 0 ) nq = nquad_Tx
        
            do iq = 1,nq
                
            isign =  1
            if (lNeg_ikx(iTx)) isign = -1
               
!
! Compute F(v) terms using a delta function for the source:
!
         
!
! Get a single element that has the transmitter:
!
            e = eTx(iTx,iq)
            
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)
 
!
! Set up the basis functions:
!
            n   = mesh%emap(1:3,e)
            ye  = mesh%y(n)
            ze  = mesh%z(n)
            
!
! Evaluate bump basis at transmitter location:
!
            call get_bump_basis(ye,ze,yTxQ(iTx,iq),zTxQ(iTx,iq),q,dqdy,dqdz)  
            
!
! Apply coefficients:
!

            ee = -q*jx(iTx)            - isign*ikx  * ( dqdy*jy(iTx)/gammay2                 + dqdz*jz(iTx)/gammaz2 )   & 
               &                             - ommu * ( dqdy*mz(iTx)*(sigy - iomeps)/gammay2 - dqdz*my(iTx)* (sigz - iomeps)/gammaz2 ) 
            hh = -(                          - ommu * ( dqdy*jz(iTx)/gammaz2                 - dqdz*jy(iTx)/gammay2 ) &
               & -q*mx(iTx)*ommu - isign*ikx * ommu * ( dqdy*my(iTx)/gammaz2                 + dqdz*mz(iTx)/gammay2 ) )
            
            ee = isign*ee ! for neg kx system

!
! Insert these into the correct rhs location
!
            if (nq > 1) then
            
                expikdx = exp(-isign*ikx*(xTxQ(iTx,iq)-xTx(iTx))) ! Fourier phase shift
                
                rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + expikdx*ee*quad_weights_Tx(iq)/2.d0    ! Ex
                rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + expikdx*hh*quad_weights_Tx(iq)/2.d0    ! Hx
            
            else ! point dipole
                rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + ee ! Ex
                rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + hh ! Hx
            
            endif                                            
               
            enddo ! iq = 1,nq
        else! MT  ! this is essentially the same as gen_rhs_mt, except that now we have a bump basis function
! Since for total field MT F = 0, all we have to do here is account for the boundary conditions.
! Tried using homogeneous Dirichlet and it works okay, but now we account for small curvature error along edges...
            call gen_rhs_bump_f_MT

        endif
        
    enddo

!
! Set 0 Dirichlet condition for CSEM:
!
    if (sType == 'cs') then
        do e = 1,mesh%nele
             do j = 1,3
!
! ex and hx on boundaries are fixed to 0 so boundaries are fixed to 0 for error terms as well:
                 if ( (bcnod( mesh%emap(eps(j+1),e) ) /= 0) .and. (bcnod( mesh%emap(eps(j+2),e) ) /= 0) ) then
                    iedge = edges(j,e) ! Global edge number:
                    rhs_q(2*iedge-1,1:nTx) = 0d0
                    rhs_q(2*iedge  ,1:nTx) = 0d0
                endif
            enddo
        enddo
    endif
    
    
    end subroutine gen_rhs_bump_f
!==================================================================================================================================!
!================================================================================================================= gen_rhs_bump_f_MT
!==================================================================================================================================!
    subroutine  gen_rhs_bump_f_MT 
    
    implicit none 
    
    integer                     :: i,e, n(3)
    integer,dimension(3,2)      :: imtxE, imtxH ! maps for BCs
    real(8)                     :: ymin, ymax, z1, z2, z0, y2, y1  
    complex(8)                  :: ex0,ex1,ex2,hx0,hx1,hx2, errTE, errTM
    complex(8),dimension(21)    :: mtx   
 
    complex(8)                  :: sigx,sigy,sigz, dsigTE,dsigTM
    real(8)                     :: ye(3),ze(3),a(3),b(3),c(3),area, z

    integer                  :: nqp, k,nq_order
    real(8), dimension(3,50) :: qp
    real(8), dimension(50)   :: qwt 
    complex(8), dimension(3) :: qsum
    real(8) :: yq,zq,bary(3),q(3),dqdy(3),dqdz(3)

       
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_f_MT ...'
    
!
! Set up mapping arrays for pulling out parts of 6 x 6 symmetric stiffness matrix stored as 21 upper triangular coeffs:
!
    imtxE(1,1) = 3   ! edge 1 on boundary, select inner product with edges 2 and 3
    imtxE(1,2) = 5
    imtxE(2,1) = 14  ! edge 2 on boundary, select inner product with edges 3 and 1
    imtxE(2,2) = 3 
    imtxE(3,1) = 5   ! edge 3 on boundary, select inner product with edges 1 and 2
    imtxE(3,2) = 14    
    imtxH(1,1) = 9
    imtxH(1,2) = 11
    imtxH(2,1) = 18
    imtxH(2,2) = 9 
    imtxH(3,1) = 11
    imtxH(3,2) = 18       
  

!
! If scattered field, insert scattering sources in f:
!
   if (lMTscatteredField) then   

!
! Get quadrature points for scattered field integration:
!
        nq_order = 10
        call triquad(nq_order,qp,qwt,nqp)  
 
        do e=1,mesh%nele
    
            n = mesh%emap(1:3,e)
        
! Locate corresponding 1D layer for TE and TM models:
            z = sum(mesh%z(n))/3d0

! Get scattering conductivities:
            call getSigs(e,sigx,sigy,sigz) ! get anisotropic conductivity terms for element e
            
            
! Get basis coefficients and area of element e:
            ye   = mesh%y(n)
            ze   = mesh%z(n)
            call get_abc_coeffs(ye,ze,a,b,c,area)               
                     
!
! If scattered field formulation, insert current sources:
!
 
! TE:
            
            qsum = 0d0  
            do k = 1,nqp
                bary = qp(1:3,k)
                yq   = sum(bary*ye)
                zq   = sum(bary*ze) 
                call get_bump_basis(ye,ze,yq,zq,q,dqdy,dqdz)    
                dsigTE = sigx - Scat1DTE%getSig(zq)
                if (dsigTE /= 0d0 ) then    
                    qsum = qsum + dsigTE*Scat1DTE%getEH1D(zq,'e')*q*qwt(k)
                endif
            enddo 
            rhs_q(2*edges(1:3,e)-1,1) =  rhs_q(2*edges(1:3,e)-1,1) - qsum*area  ! TE is easy since n(1:3) get same source contributions
                        
 
! TM:
            
            qsum = 0d0  
            do k = 1,nqp
                bary = qp(1:3,k)
                yq   = sum(bary*ye)
                zq   = sum(bary*ze)                 
                call get_bump_basis(ye,ze,yq,zq,q,dqdy,dqdz)    
                dsigTM = sigy - Scat1DTM%getSig(zq) 
                if (dsigTM /= 0d0 ) then
                    qsum = qsum - dsigTM/sigy*Scat1DTM%getEH1D(zq,'e')*dqdz*qwt(k)
                endif
            enddo      
            rhs_q(2*edges(1:3,e)  ,1) = rhs_q(2*edges(1:3,e)  ,1)  + qsum*area ! TM part
            
                
        enddo      
    endif
   
!
! Apply boundary condition:
!
!
!  Compute the min and max x-values.
!
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)       

            
    do e = 1,mesh%nele

        n    = mesh%emap(1:3,e)
        
        do i = 1,3
    
! Apply boundary conditions for i-j pairs (i-i terms are set to 1 afterwards):
        if ( ( bcnod(n(eps(i+1))) /= 0 ) .and. ( bcnod(n(eps(i+2))) /= 0 )  ) then ! we got a boundary:
 
! Ge bounding node depths and the midpoint depth:
            z1 = mesh%z(n(eps(i+1)))
            z2 = mesh%z(n(eps(i+2)))

            z0 = (z1 + z2)/2d0  ! depth midpoint
            
            y1 = mesh%y(n(eps(i+1)))
            y2 = mesh%y(n(eps(i+2)))
               
            if (y1 /= y2) cycle !write(*,*) 'error,y1/=y1!!!, element: ',e,i, y1, y2
            if (z1 == z2) cycle ! we only care about vertical boundaries
                     
! Get the 1D MT fields:
            if (y1 == ymin) then ! left side:
                ex0 = Left1DTE%getEH1D(z0,'e')
                ex1 = Left1DTE%getEH1D(z1,'e')
                ex2 = Left1DTE%getEH1D(z2,'e')
                hx0 = Left1DTM%getEH1D(z0,'h')
                hx1 = Left1DTM%getEH1D(z1,'h')
                hx2 = Left1DTM%getEH1D(z2,'h')
            else ! right side
                ex0 = Right1DTE%getEH1D(z0,'e')
                ex1 = Right1DTE%getEH1D(z1,'e')
                ex2 = Right1DTE%getEH1D(z2,'e')
                hx0 = Right1DTM%getEH1D(z0,'h')
                hx1 = Right1DTM%getEH1D(z1,'h')
                hx2 = Right1DTM%getEH1D(z2,'h')                          
            endif
            
! The error in the bump basis is then:
            errTE = ex0 - (ex1 + ex2 ) / 2d0               
            errTM = hx0 - (hx1 + hx2 ) / 2d0      
    
! Set RHS equal to this:
            rhs_q(2*edges(i,e)-1,1) = errTE
            rhs_q(2*edges(i,e)  ,1) = errTM
            
! Subtract off RHS of other edges:
            call getStiff_bump(e, mtx)
            
! Note that E-H coupling is 0 with anisotropic MT since kx=0, so only use decoupled terms for TE and TM modes:
            rhs_q(2*edges(eps(i+1),e)-1,1) = rhs_q(2*edges(eps(i+1),e)-1,1) - mtx(imtxE(i,1))*errTE
            rhs_q(2*edges(eps(i+2),e)-1,1) = rhs_q(2*edges(eps(i+2),e)-1,1) - mtx(imtxE(i,2))*errTE
            
            rhs_q(2*edges(eps(i+1),e)  ,1) = rhs_q(2*edges(eps(i+1),e)  ,1) - mtx(imtxH(i,1))*errTM   
            rhs_q(2*edges(eps(i+2),e)  ,1) = rhs_q(2*edges(eps(i+2),e)  ,1) - mtx(imtxH(i,2))*errTM
            
        endif
       
       enddo !i
                               
    enddo ! e

        
    end subroutine gen_rhs_bump_f_MT
!==================================================================================================================================!
!==================================================================================================================== gen_rhs_bump_b
!==================================================================================================================================!
    subroutine  gen_rhs_bump_b(imode)
!
! Routine to compute F(v) - B(u_h,v) for all v in W_h (bump space)
! F(v) contribution is handled in another subroutine.
! Also computes G(v) - B(v,wh) for dual error.
!
    integer, intent(in) :: imode

!
! Local variables:
!
    integer         :: i, j, k, e, iedge, iTx, n(3), ict, ij(6,6), isign
    complex(8)      :: ul(6), mtx(21), smtx(21)
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_b ...'   
!
! Create 6 x 6 lookup table for indices to symmetric upper triangular stiffness matrix:
!
    ict = 0
    do i = 1,6    
        do k = i,6
            ict = ict + 1
            ij(i , k ) = ict
            ij(k , i ) = ict
       enddo
    enddo
  
!
! Create int( -B(u_h,v) ) where u_h is in V_h and v is in W_h
!
    
!
! Loop over elements:
!
    do e = 1,mesh%nele
    
!
! Get B(u,v) lin-bump stiffness matrix for this element:
!
        call getStiff_lin_bump(e,mtx)
 
!
! Compute B(u,v) for this element, subtract from RHS vector (rhs = F - B(u,v)):
!
! B(u,v) =  MTX*u
!
        n  = mesh%emap(1:3,e)        
        do iTx = 1,nTx
            if (sType == 'mt') then            
                isign = 1
            else
                isign =  1
                if (lNeg_ikx(iTx)) isign = -1
            endif    
            smtx = mtx
            smtx([2,4,6,8,10,13,15,17,20]) = isign*smtx([2,4,6,8,10,13,15,17,20])  ! 3 x 3 off diagonal block
                    
            do i = 1,3 ! nodes
                if (imode == 0 ) then ! use rhs from primal problem:
                    ul(2*i-1) = (rhs( 2*n(i)-1, iTx ))
                    ul(2*i)   = (rhs( 2*n(i)  , iTx ))
                else
                    ul(2*i-1) = (rhs_wh( 2*n(i)-1, iTx ))
                    ul(2*i)   = (rhs_wh( 2*n(i)  , iTx ))
                endif  
            enddo
            
            do j = 1,3 ! edges
                iedge = edges(j,e) ! Global edge number
                rhs_q(2*iedge-1, iTx) = rhs_q(2*iedge-1, iTx) - sum(smtx(ij(2*j-1,1:6))*ul)
                rhs_q(2*iedge  , iTx) = rhs_q(2*iedge  , iTx) - sum(smtx(ij(2*j  ,1:6))*ul)
            enddo              
    
        enddo ! iTx

    enddo ! e = 1,nele
         
! Remember that for negative wavenumbers we solved [A B; B^T C][-u v] = [-c d] so take negative of u
    if (sType == 'cs') then ! CSEM dipole
        do iTx = 1,nTx
            if (lNeg_ikx(iTx)) then
                 rhs_q(1:2*nedges:2,iTx)  = -rhs_q(1:2*nedges:2,iTx)  ! for neg kx system  ! note Fortan increment is at end
            endif        
        enddo                             
    endif
          
    do e = 1,mesh%nele
        do j = 1,3
!
! ex and hx on boundaries are fixed to 0, so same for error terms on boundaries:
             if ( (bcnod( mesh%emap(eps(j+1),e) ) /= 0) .and. (bcnod( mesh%emap(eps(j+2),e) ) /= 0) ) then
                iedge = edges(j,e) ! Global edge number:
                rhs_q(2*iedge-1,1:nTx) = 0d0
                rhs_q(2*iedge  ,1:nTx) = 0d0
            endif
        enddo
    enddo

    
    end subroutine gen_rhs_bump_b         
    
!==================================================================================================================================!
!==================================================================================================================== gen_rhs_bump_g
!==================================================================================================================================!
    subroutine  gen_rhs_bump_g
!
! Computes G(v)  for dual error.
! Since E and H are decoupled on LHS (in bump basis), we generate
! separate RHS vectors so that the linear systems are smaller.
 
!
! Local variables:
!
    integer                  :: n(3),e, iTx, iRx,i, iq,nq
    real(8)                  :: smtx(3,3), mmtx(3,3),denomEx,denomHx,denomEg,denomHg
    real(8)                  :: wtEx,wtHx,wtEg,wtHg, ye(3),ze(3),a(3),b(3),c(3),area,yp,zp,qwt
    complex(8), dimension(3) :: fel,fhl,feq,fhq
    complex(8), dimension(3) :: exrn,etrn
    complex(8), dimension(3) :: hxrn,htrn
    complex(8)               :: exrd, etrd, hxrd, htrd
 
    complex(8) , dimension(:,:), allocatable:: rhs_q0
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_bump_g ...'  
        
 
    allocate(rhs_q0(2*nedges,1:nTx)   )
    rhs_q0 = rhs_q
    
! Reset to zero:
    rhs_q = 0
     
!
! Loop over transmitters and receivers:

    do iTx = 1,nTx
        
        do iRx = 1,nRx

            nq = 1
            if (lengthRx(iRx) > 0) nq = nquad_Rx     
                    
            do iq= 1,nq
        
            call markCloud(iRx,iTx,iq)

            qwt = 1d0
            if (nq > 1) qwt = quad_weights_Rx(iq)/2.d0
         
            denomEx = 0.
            denomHx = 0.
            denomEg = 0.
            denomHg = 0. 
        
!
! Integrate denominator:
!
            do e = 1,mesh%nele

                if (.not.lCloud(e)) cycle
            
                n    = mesh%emap(1:3,e)
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
 
                fel(1:3) = rhs(2*n-1,iTx)
                fhl(1:3) = rhs(2*n  ,iTx)         
            
! l_i l_j
                mmtx      =  area/12d0
                mmtx(1,1) =  area/6d0
                mmtx(2,2) =  area/6d0
                mmtx(3,3) =  area/6d0                  
              
                exrd = (sum(matmul(mmtx,fel)*conjg(fel)))
                hxrd = (sum(matmul(mmtx,fhl)*conjg(fhl)))
                       
            
! grad l_i dot grad l_j
        
                do i = 1,3
                    smtx(i,:) = (b*b(i)+c*c(i))/(4d0*area)
                enddo
 
                etrd = dble(sum(matmul(smtx,fel)*conjg(fel)))
                htrd = dble(sum(matmul(smtx,fhl)*conjg(fhl)))

                if (sType == 'mt') then
                    denomEx = denomEx + exrd + abs(ecutoff)**2*area
                    denomHx = denomHx + hxrd + abs(hcutoff)**2*area

                    denomEg = denomEg + etrd + (ecutoff**2)*area  
                    denomHg = denomHg + htrd + (hcutoff**2)*area  
                
                else
                    denomEx = denomEx + exrd + abs(ecutoff/kx)**2*area
                    denomHx = denomHx + hxrd + abs(hcutoff/kx)**2*area

                    denomEg = denomEg + etrd + abs(ecutoff/kx)**2*area
                    denomHg = denomHg + htrd + abs(hcutoff/kx)**2*area
                           
                endif
                

            enddo

!
! Integrate numerator:
!
                            
            do e = 1,mesh%nele
        
                if (.not.lCloud(e)) cycle
   
                n    = mesh%emap(1:3,e)   
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
!                yp = sum(ye)/3.
!                zp = sum(ze)/3.
                yp = yRxQ(iRx,iq)
                zp = zRxQ(iRx,iq)

                fel(1:3) = rhs(2*n-1,iTx)
                fhl(1:3) = rhs(2*n  ,iTx)    
                feq(1:3) = rhs_q0( 2*edges(1:3,e)-1,iTx )
                fhq(1:3) = rhs_q0( 2*edges(1:3,e)  ,iTx )
 
 
! Compute weights:
                call getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)
       
! Numerator:
! q_i * q_j:
                mmtx      = area*4d0/45d0
                mmtx(1,1) = area*8d0/45d0
                mmtx(2,2) = area*8d0/45d0
                mmtx(3,3) = area*8d0/45d0      
                exrn = matmul(mmtx,(feq))
                hxrn = matmul(mmtx,(fhq))
        
! grad_q_i dot grad_q_j:
                do i = 1,3
                    smtx(i,:) = (b*b(i)+c(i)*c)*2d0/(3d0*area)
                enddo           
                smtx(1,1) = (b(1)*b(1)-b(2)*b(3) + c(1)*c(1) - c(2)*c(3) )*2d0/(3d0*area)
                smtx(2,2) = (b(2)*b(2)-b(3)*b(1) + c(2)*c(2) - c(3)*c(1) )*2d0/(3d0*area)
                smtx(3,3) = (b(3)*b(3)-b(1)*b(2) + c(3)*c(3) - c(1)*c(2) )*2d0/(3d0*area)
                etrn = matmul(smtx,feq)
                htrn = matmul(smtx,fhq)     
         
                exrn = abs(exrn)
                hxrn = abs(hxrn)
                etrn = abs(etrn)
                htrn = abs(htrn)
 
!
! User selected functional:
!
                select case (idual_func)
            
! G0:  residual error only:
                case (0)
                         
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEx*( exrn/denomEx )*qwt
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHx*( hxrn/denomHx )*qwt
 
! G1:  gradient of residual error only:
                case(1)
                
                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEg*( etrn/denomEg )*qwt
                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHg*( htrn/denomHg )*qwt
                
! G2: residual error plus gradient of residual error:
                case(2)
   
!                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEx*( exrn/denomEx )*qwt
!                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHx*( hxrn/denomHx )*qwt
!
!                    rhs_q(2*edges(1:3,e)-1,iTx) = rhs_q(2*edges(1:3,e)-1,iTx) + wtEg*( etrn/denomEg )*qwt
!                    rhs_q(2*edges(1:3,e)  ,iTx) = rhs_q(2*edges(1:3,e)  ,iTx) + wtHg*( htrn/denomHg )*qwt
                    do i =1,3
                    rhs_q(2*edges((i),e)-1,iTx) = max(abs(rhs_q(2*edges((i),e)-1,iTx)) , abs(wtEx*( exrn(i)/denomEx )*qwt))
                    rhs_q(2*edges((i),e)  ,iTx) = max(abs(rhs_q(2*edges((i),e)  ,iTx)) , abs(wtHx*( hxrn(i)/denomHx )*qwt))
 
                    rhs_q(2*edges((i),e)-1,iTx) = max(abs(rhs_q(2*edges((i),e)-1,iTx)) , abs(wtEg*( etrn(i)/denomEg )*qwt))
                    rhs_q(2*edges((i),e)  ,iTx) = max(abs(rhs_q(2*edges((i),e)  ,iTx)) , abs(wtHg*( htrn(i)/denomHg )*qwt))
                    enddo
                                                           
                end select
  
            enddo ! e
            enddo ! iq
        enddo  ! iRx
    enddo  ! iTx
 
    deallocate (  rhs_q0) 
        
    end subroutine gen_rhs_bump_g
 
!==================================================================================================================================!
!====================================================================================================================== gen_rhs_dual
!==================================================================================================================================!
    subroutine gen_rhs_dual
!
!  Computes RHS vector for dual system used for DRW
!
    implicit none

    integer                  :: i,iRx,e,n(3),iTx,iq,nq
    real(8)                  :: smtx(3,3), mmtx(3,3), denomEx, denomHx, denomEg, denomHg, qwt
    real(8)                  :: wtEx,wtHx,wtEg,wtHg,ye(3),ze(3), a(3),b(3),c(3),area,yp,zp 
    complex(8), dimension(3) :: fel,fhl,feq,fhq 
    complex(8), dimension(3) :: exrn,etrn
    complex(8), dimension(3) :: hxrn,htrn
    complex(8)               :: exrd, etrd, hxrd, htrd
    
 
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_rhs_dual ...'     
    
!
!  Initialize the RHS
!
    rhs_wh = 0d0
     
!
! Loop over transmitters and receivers:
!
    do iTx = 1,nTx
    
        do iRx = 1,nRx

            nq = 1
            if (lengthRx(iRx) > 0) nq = nquad_Rx     
                    
            do iq= 1,nq
        
            call markCloud(iRx,iTx,iq)

            qwt = 1d0
            if (nq > 1) qwt = quad_weights_Rx(iq)/2.d0

        
               
            denomEx = 0.
            denomHx = 0.
            denomEg = 0.
            denomHg = 0. 
        
!
! Integrate denominator:
!
            do e = 1,mesh%nele

                if (.not.lCloud(e)) cycle
            
                n    = mesh%emap(1:3,e)
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
 
                fel(1:3) = rhs(2*n-1,iTx)  
                fhl(1:3) = rhs(2*n  ,iTx)    
            
! l_i l_j
                mmtx      =  area/12d0
                mmtx(1,1) =  area/6d0
                mmtx(2,2) =  area/6d0
                mmtx(3,3) =  area/6d0                  
              
                exrd = (sum(matmul(mmtx,fel)*conjg(fel)))
                hxrd = (sum(matmul(mmtx,fhl)*conjg(fhl)))                     
            
! grad l_i dot grad l_j
                do i = 1,3
                    smtx(i,:) = (b*b(i)+c*c(i))/(4d0*area)
                enddo
 
                etrd = dble(sum(matmul(smtx,fel)*conjg(fel)))
                htrd = dble(sum(matmul(smtx,fhl)*conjg(fhl)))


            if (sType == 'mt') then
                denomEx = denomEx + exrd + abs(ecutoff)**2*area
                denomHx = denomHx + hxrd + abs(hcutoff)**2*area

                denomEg = denomEg + etrd + (ecutoff**2)*area  
                denomHg = denomHg + htrd + (hcutoff**2)*area   
            else
                denomEx = denomEx + exrd + abs(ecutoff/kx)**2*area
                denomHx = denomHx + hxrd + abs(hcutoff/kx)**2*area

                denomEg = denomEg + etrd + abs(ecutoff/kx)**2*area
                denomHg = denomHg + htrd + abs(hcutoff/kx)**2*area
            endif

            
            enddo
        
!
! Integrate numerator:
!
           do e = 1,mesh%nele
       
                if (.not.lCloud(e)) cycle
            
!
! Set up the basis functions:
!
                n    = mesh%emap(1:3,e)
                ye   = mesh%y(n)
                ze   = mesh%z(n)
                call get_abc_coeffs(ye,ze,a,b,c,area)      
 
!                yp = sum(ye)/3.
!                zp = sum(ze)/3.
                yp = yRxQ(iRx,iq)
                zp = zRxQ(iRx,iq)
 
                fel(1:3) = rhs(2*n-1,iTx)
                fhl(1:3) = rhs(2*n  ,iTx)
                feq(1:3) = rhs_q(2*edges(1:3,e)-1,iTx)
                fhq(1:3) = rhs_q(2*edges(1:3,e)  ,iTx)   
            
! Compute weights:
                call getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)

! Numerator:
       
! q_i * l_j:
                mmtx      =  area*2d0/15d0
                mmtx(1,1) =  area/15d0   
                mmtx(2,2) =  area/15d0
                mmtx(3,3) =  area/15d0      
            
                exrn = matmul(mmtx,(feq))
                hxrn = matmul(mmtx,(fhq))
            
! grad_q_i dot grad_l_i:
                do i = 1,3
                    smtx(i,:) = -(b*b(i)+c(i)*c)/(3d0*area)
                enddo            
                etrn = matmul(smtx,feq)
                htrn = matmul(smtx,fhq)
 
                exrn = abs(exrn)
                hxrn = abs(hxrn)
                etrn = abs(etrn)
                htrn = abs(htrn)
 
!
! User selected functional:
!
                select case (idual_func)
            
! G0:  residual error only:
                case (0)
                
                    rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEx*( exrn/denomEx )*qwt
                    rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHx*( hxrn/denomHx )*qwt
                       
! G1:  gradient of residual error only:
                case(1)

                    rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEg*( etrn/denomEg )*qwt
                    rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHg*( htrn/denomHg )*qwt              
                
! G2: residual error plus gradient of residual error:
                case(2)
 
!                    rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEx*( exrn/denomEx )*qwt
!                    rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHx*( hxrn/denomHx )*qwt
!
!                    rhs_wh(2*n-1,iTx) = rhs_wh(2*n-1,iTx) + wtEg*( etrn/denomEg )*qwt
!                    rhs_wh(2*n  ,iTx) = rhs_wh(2*n  ,iTx) + wtHg*( htrn/denomHg )*qwt
                    
!             write(*,'(i5,1x,4(g12.4,1x))') iRx,  wtEx*( exrn/denomEx ),wtHx*( hxrn/denomHx ),wtEg*( etrn/denomEg ) ,  wtHg*( htrn/denomHg )
          
          
                    do i =1,3
                    rhs_wh(2*n(i)-1,iTx) = max(abs(rhs_wh(2*n(i)-1,iTx)) , abs(wtEx*( exrn(i)/denomEx )*qwt ))
                    rhs_wh(2*n(i)  ,iTx) = max(abs(rhs_wh(2*n(i)  ,iTx)) , abs(wtHx*( hxrn(i)/denomHx )*qwt ))
 
                    rhs_wh(2*n(i)-1,iTx) = max(abs(rhs_wh(2*n(i)-1,iTx)) , abs(wtEg*( etrn(i)/denomEg )*qwt ))
                    rhs_wh(2*n(i)  ,iTx) = max(abs(rhs_wh(2*n(i)  ,iTx)) , abs(wtHg*( htrn(i)/denomHg )*qwt ))  
                    enddo      
                end select
 
            enddo    ! loop over elements
            enddo ! iq
        enddo ! loop over receivers
    enddo ! loop over transmitters
    
! Set sign for negative kx:
    if (sType == 'cs') then ! CSEM dipole
        do iTx = 1,nTx
            if (lNeg_ikx(iTx)) then
                rhs_wh(1:2*mesh%nnod:2,iTx) = -rhs_wh(1:2*mesh%nnod:2,iTx)  ! note Fortan increment is at end
            endif        
        enddo         
    endif 
         
!
!  Set sourcing function to 0 on outer boundaries
!
    do i=1,mesh%nnod
        if (bcnod(i).ge.1) then
            rhs_wh(2*i-1,1:nTx) = 0d0
            rhs_wh(2*i  ,1:nTx) = 0d0
        endif
    enddo
               

    end subroutine gen_rhs_dual      
       
!==================================================================================================================================!
!========================================================================================================================= markCloud
!==================================================================================================================================!
 
 subroutine markCloud(iRx,iTx,iq)

    integer :: e,i,n(3), iTx, iRx,iq
    real(8) :: r,ye(3),ze(3)
 
      
    lCloud = .false.
!
! Mark all elements that are close to the receiver:
!
    do e = 1,mesh%nele

        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)    

        do i = 1,3

            r = sqrt( (yRxQ(iRx,iq) - ye(i))**2 + (zRxQ(iRx,iq) - ze(i))**2)

            if (r < rxCloud)  lCloud(e) = .true.

        enddo
    enddo            

!
! Also mark all the element containing the receiver:
!

     if ( iDataMask(iRx,iTx) > 0 ) then    
! write(*,'(3(i3,1x),f12.5,1x,f12.5,1x,i6)') iRx,iTx,iq,yRxQ(iRx,iq),zRxQ(iRx,iq),eRx(iRx,iq)
        lCloud(eRx(iRx,iq)) = .true.
     endif   

 
 end subroutine markCloud
                
                
!==================================================================================================================================!
!======================================================================================================================== getWeights
!==================================================================================================================================!
    subroutine getWeights(iRx,iTx,yp,zp,area,ye,ze,fel,feq,fhl,fhq,wtEx,wtHx,wtEg,wtHg)

!
! This subroutine sets some damping weights for the dual source functions used to guide adaptive refinement.
!
    
    integer, intent(in)                     :: iTx, iRx
    real(8), intent(in)                     :: area,yp,zp 
    real(8), dimension(3), intent(in)       :: ye,ze
    complex(8), dimension(3), intent(in)    :: fel,feq,fhl,fhq 
    real(8), intent(out)                    :: wtEx,wtHx,wtEg,wtHg 
  
    real(8)                     :: rTx,wrTx
    complex(8), dimension(3)    :: fiex,fihx,fieq,fihq 
       
!
! First get a range dependent weight, to avoid over-refinement where receivers are located close to the source:
!
    wrTx = 1d0
        
    if (sType == 'cs') then
        if (lengthTx(iTx) > 0) then
            rTx  = minval(sqrt( ( yp - yTxQ(iTx,1:nquad_Tx))**2 + ( zp - zTxQ(iTx,1:nquad_Tx)))**2 )
        else ! point dipole
            rTx  = sqrt( ( yp - yTx(iTx))**2 + ( zp - zTx(iTx))**2 )
        endif 
        if (rTx < minRangeProtector)   wrTx = 0  
   
    endif           

! Data mask weights. If isite,iTx combo not needed, set weight to zero
    if ( iDataMask(iRx,iTx) == 0 ) wrTx = 0d0  
    
    if ( area < minArea ) wrTx = 0d0
 
! map to component weights:
    wtEx = wrTx
    wtHx = wrTx
    wtEg = wrTx
    wtHg = wrTx

!
! Interpolate field to yp,zp location:
!
    call fget_interp_lin(yp,zp,fel,ye,ze,fiex)
    call fget_interp_lin(yp,zp,fhl,ye,ze,fihx)
    call fget_interp_bump(yp,zp,feq,ye,ze,fieq)
    call fget_interp_bump(yp,zp,fhq,ye,ze,fihq)        
            
! Set weights to zero if fields or errors are below cutoff thresholds:
    

    if (sType == 'mt') then
  
        if ( ( abs(fiex(1)+fieq(1)) < ecutoff )) then  
            wtEx = 0d0
            wtEg = wtEx   
        endif
        
        if ( ( abs(fihx(1)+fihq(1)) < hcutoff )) then  
            wtHx = 0d0
            wtHg = wtHx
        endif            
            
    else
    
        if ( ( abs(fiex(1)+fieq(1)) < ecutoff/kx )) then  
            wtEx = 0d0   
            wtEg = wtEx   
        endif
  
        if ( ( abs(fihx(1)+fihq(1)) < hcutoff/kx )) then  
            wtHx = 0d0 
            wtHg = wtHx
        endif
     
    endif 
          
    end subroutine getWeights
     
!==================================================================================================================================!
!=========================================================================================================================== gen_lhs
!==================================================================================================================================!
    subroutine gen_lhs 
!
! Generates the LHS matrix for CSEM and MT problems using optional
! methods
!
! Kerry Key
! Scripps Institution of Oceanography
!


    integer*8, save :: gen_lhs = -1
	
    call SCOREP_F_Begin(gen_lhs,"gen_lhs",0,"",3894) 


    if (sType == 'mt') then ! MT using decoupled linear systems:
    
        kx = 0d0 ! Make sure kx=0
        call gen_lhs_mt
    else
        call gen_lhs_csem
    endif
    

    call SCOREP_F_RegionEnd( gen_lhs ) !end Instrumentation


    end subroutine gen_lhs        

!==================================================================================================================================!
!====================================================================================================================== gen_lhs_csem
!==================================================================================================================================!
      subroutine gen_lhs_csem 
!
! Generates the lhs matrix for the coupled system in compressed sparse row format.
!
! Kerry Key
! Scripps Institution of Oceanography
 
    
!
! Local parameters:
!
    integer                             :: i,j,e, n(3),v1,v2, indi, indj, indm, ind, i0,i1,ii,jj
    complex(8), dimension(21)           :: mtx
    integer,dimension(:), allocatable   :: nAdj
        
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs ...'
 
!
! Allocate the sparse matrix arrays:
!
    nnz = 4*mesh%nnod + 8*nedges        
    allocate ( val(nnz), col(nnz), irw(2*mesh%nnod+1) )
    allocate (   nAdj(mesh%nnod) )  ! temporary for CSR matrix
    
!
! Initialize a few things
!
    val  = (0d0,0d0)     ! Initialize the lhs matrix in compressed sparse row format:
    col = 0
    irw = 0
 
!
! Loop through all elements and count node adjacencies:
!
    nAdj = 1  ! for the self adjacency
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)   
        do j = 1,3
            v1 = n(eps(j))
            v2 = n(eps(j+1))       
            if (v1 > v2) then ! use direction of edge to ensure we count each connection only once
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1              
             elseif (mesh%neighborlist(eps(j+2),e) == -1) then ! this is a boundary edge with reverse order, add it:
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1     
            endif
        enddo
    enddo   

!
! Make irw row pointer array:
!
    irw(1) =  1
   
    do e =  1,mesh%nnod
        irw(2*e      )  = irw(2*e-1) + 2*nAdj(e)  ! Hx that goes with last Ex
        irw(2*(e+1)-1)  = irw(2*e  ) + 2*nAdj(e)  ! Ex for next row
    enddo
   
!
! Loop over all elements and insert entries into the sparse matrix:
!
    do e = 1,mesh%nele
    
! Get stiffness matrix for this element:
        call getStiff_lin(e, mtx)
            
        n = mesh%emap(1:3,e)    

        do i = 1,3
        
            indi = n(i)
                
            do j = 1,3
                
                indj = n(j)                
                
                indm = getUpperTriIndex(6,2*i-1,2*j-1)  ! Ex-Ex
                
                do ind = irw(2*indi-1),irw(2*indi)-1                 
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj-1) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj-1
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i-1,2*j  )  ! Ex-Hx
                
                do ind = irw(2*indi-1),irw(2*indi)-1
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj
                        exit
                    endif 
                enddo
 
                indm = getUpperTriIndex(6,2*i  ,2*j  )  ! Hx-Hx
                
                do ind = irw(2*indi),irw(2*indi+1)-1                    
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj  ) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj  
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i  ,2*j-1 )  ! Hx-Ex
                
                do ind = irw(2*indi),irw(2*indi+1)-1                  
                    if ( (col(ind) == 0) .or. (col(ind) == 2*indj-1  ) ) then  
                        val(ind) = val(ind) + mtx(indm)    
                        col(ind) = 2*indj-1  
                        exit
                    endif 
                enddo
            
 
            enddo ! j
            
        enddo  ! i
        
   enddo  ! e

!
! Apply boundary conditions:
!
    do i = 1,2*mesh%nnod

        do i1 = irw(i),irw(i+1)-1
        
             j  = col(i1) 
             ii = floor(dble(i+1.)/2.)
             jj = floor(dble(j+1.)/2.) 
             
             if (   (  (bcnod(ii)>0).or.(bcnod(jj)>0) ) .and. (i /= j ) ) then
                val(i1) = 0d0    
             elseif (  (bcnod(ii)>0 ) .and. (i == j) ) then
                val(i1) = 1d0   
!write(*,*)  i,j,ii,jj,i1,val(i1)
             endif 
 
       enddo
    enddo
 
!
! Finish up by sorting each row by column number:
!
    do i = 1,2*mesh%nnod
        i0 = irw(i)
        i1 = irw(i+1)-1
        call quick_sort(col(i0:i1),zlist1=val(i0:i1)) 
    enddo
   
    deallocate (nAdj) 
 
         
    end subroutine gen_lhs_csem 
!==================================================================================================================================!
!======================================================================================================================== gen_lhs_mt
!==================================================================================================================================!
    subroutine gen_lhs_mt
!
! Generates the lhs matrix in compressed sparse row format for the
! uncoupled MT E and H systems.
!
! Kerry Key
! Scripps Institution of Oceanography
!
    
    implicit none
    
    integer                             :: nnz, e, indm, n(3), i, j, ind, i0, i1, indi, indj, v1, v2
    complex(8), dimension(21)           :: mtx
    integer,dimension(:), allocatable   :: nAdj 
       
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_mt...'

!
! Initialize a few things
!
    nnz =  mesh%nnod + 2*nedges  
    
    allocate ( rhs_lE(mesh%nnod,nTx), rhs_lH(mesh%nnod,nTx)  ) 
    allocate ( val_lE(nnz), col_lE(nnz), irw_lE(mesh%nnod+1) )
    allocate ( val_lH(nnz), col_lH(nnz), irw_lH(mesh%nnod+1) )
    
    allocate (   nAdj(mesh%nnod) )  ! temporary for CSR matrix
    
         
    val_lE = (0d0,0d0)     ! Initialize the lhs matrices in compressed sparse row format:
    col_lE = 0
    irw_lE = 0
    nnz_lE = nnz
    
    val_lH = (0d0,0d0)     ! Initialize the lhs matrices in compressed sparse row format:
    col_lH = 0
    irw_lH = 0
    nnz_lH = nnz  
 
!
! Loop through all elements and count node adjacencies:
!
 
    nAdj = 1  ! for the self adjacency
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)   
        do j = 1,3
            v1 = n(eps(j))
            v2 = n(eps(j+1))       
            if (v1 > v2) then ! use direction of edge to ensure we count each connection only once
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1              
             elseif (mesh%neighborlist(eps(j+2),e) == -1) then ! this is a boundary edge with reverse order, add it:
                nAdj(v1) = nAdj(v1) + 1
                nAdj(v2) = nAdj(v2) + 1     
            endif
        enddo
    enddo   
 

!
! Make irw row pointer arrays:
!
    irw_lE(1) =  1
    irw_lH(1) =  1 
    do e =  1,mesh%nnod
         irw_lE(e+1)  = irw_lE(e) + nAdj(e)
         irw_lH(e+1)  = irw_lH(e) + nAdj(e)
    enddo
   
!
! Loop over all elements and insert entries into the sparse matrix:
!
    do e = 1,mesh%nele
    
! Get stiffness matrix for this element:
        call getStiff_lin(e, mtx)

        n = mesh%emap(1:3,e)    

        do i = 1,3
        
            indi = n(i)
                
            do j = 1,3
                
                indj = n(j)                
                
! Ex:
                indm = getUpperTriIndex(6,2*i-1,2*j-1) 
                do ind = irw_lE(indi),irw_lE(indi+1)-1
                    if ( (col_lE(ind) < 1) .or. (col_lE(ind) == indj) ) then  
                        val_lE(ind) = val_lE(ind) + mtx(indm)    
                        col_lE(ind) = indj
                        exit
                    endif 
                enddo
 
                
! Hx:
                indm = getUpperTriIndex(6,2*i,2*j) 
                do ind = irw_lH(indi),irw_lH(indi+1)-1
                    if ( (col_lH(ind) < 1) .or. (col_lH(ind) == indj) ) then  
                        val_lH(ind) = val_lH(ind) + mtx(indm)    
                        col_lH(ind) = indj
                        exit
                    endif 
                enddo
 
            enddo ! j
            
        enddo  ! i
        
   enddo  ! e

  
!
! Apply boundary conditions:
!
    do i = 1,mesh%nnod

        do i1 = irw_lE(i),irw_lE(i+1)-1
        
             j = col_lE(i1)

             if (   (  (bcnod(i)>0).or.(bcnod(j)>0) ) .and. (i /= j ) ) then
                val_lE(i1) = 0d0  
                val_lH(i1) = 0d0      
             elseif (  (bcnod(i)>0 ) .and. (i == j) ) then
                val_lE(i1) = 1d0  
                val_lH(i1) = 1d0     
             endif   
 
       enddo
    enddo    

!
! Finish up by sorting each row by column number:
!
    do i = 1,mesh%nnod
        i0 = irw_lE(i)
        i1 = irw_lE(i+1)-1
        call quick_sort(col_lE(i0:i1),zlist1=val_lE(i0:i1))
        call quick_sort(col_lH(i0:i1),zlist1=val_lH(i0:i1)) 
    enddo
    
    deallocate (  nAdj)  
        
    end subroutine gen_lhs_mt

!==================================================================================================================================!
!================================================================================================================== getUpperTriIndex
!==================================================================================================================================!
    integer function getUpperTriIndex(n,i,j) result(ind)
    
    integer :: n,i,j
    
    if (i<=j) then
        ind  = ((i-1)*(2*n-i))/2 + j
    else
        ind  = ((j-1)*(2*n-j))/2 + i
    endif
    
    end function getUpperTriIndex

!==================================================================================================================================!
!====================================================================================================================== gen_lhs_bump
!==================================================================================================================================!
      subroutine gen_lhs_bump
!
! Generates the lhs matrix for the bump space in compressed sparse row
! format.  Coupled anisotropic matrix is 2nedges x 2nedges.
!
! Kerry Key
! Scripps Institution of Oceanography
 
!
! Local parameters for element stiffness comps:
!
    integer                             :: i,j,e, nnz, indm, indi,indj, n(3),  i1, i0, ind 
    integer, dimension(:),allocatable   :: bcedge, eAdj
    complex(8), dimension(21)           :: mtx
    
    if (lprintDebug_em2dkx) write(*,*) myID,': gen_lhs_bump ...'   
 
!
! Allocate the sparse matrix arrays:
!
    nnz = 4*nedges + 24*mesh%nele          
    allocate (  bcedge(2*nedges)   )  ! temporary for CSR matrix
    allocate ( col_q(nnz), val_q(nnz),irw_q(2*nedges+1)  )

!
! Initialize a few things
!
    val_q  = (0d0,0d0)     ! Initialize the lhs matrix in compressed sparse row format:
    col_q = 0
    irw_q = 0
    nnz_q = nnz
 
    bcedge = 0

! create flag for edges on boundary:
    do e = 1,mesh%nele
        n = mesh%emap(1:3,e)       
        do i = 1,3
            if ( (bcnod(n(eps(i+1))) > 0) .and. (bcnod(n(eps(i+2))) > 0) ) then
                bcedge(2*edges(i,e)-1) = 1
                bcedge(2*edges(i,e)  ) = 1
            endif
        enddo    
    enddo  
    
    
!
! Loop through all elements and count edge adjacencies:
!
    allocate(eAdj(2*nedges))
    eAdj = 2  ! for the self connections
    do e = 1,mesh%nele
        do i = 1,3
            indi = edges(i,e)            
            eAdj(2*indi-1) = eAdj(2*indi-1) + 4  ! Ex edges
            eAdj(2*indi  ) = eAdj(2*indi  ) + 4  ! Hx edges
        enddo
     
    enddo

!
! Make irw_q row pointer array:
!
    irw_q(1) =  1 
    do e = 1,2*nedges
         irw_q(e+1)  = irw_q(e) + eAdj(e)
    enddo
 
!
! Loop over all elements and insert stiffness matrix components:
!
        
    do e = 1,mesh%nele
    
! Get stiffness matrix for this element:
        call getStiff_bump(e,mtx)     
       

! Loop over rows and columns and insert 2 x 2 entries:
            
        do i = 1,3  ! row
            
            indi = edges(i,e)
            
            do j = 1,3  ! column
                
                indj = edges(j,e)
                             
                indm = getUpperTriIndex(6,2*i-1,2*j-1)  ! Ex-Ex
                
                do ind = irw_q(2*indi-1),irw_q(2*indi)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj-1) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj-1
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i-1,2*j  )  ! Ex-Hx
                
                do ind = irw_q(2*indi-1),irw_q(2*indi)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj
                        exit
                    endif 
                enddo
 
                indm = getUpperTriIndex(6,2*i  ,2*j  )  ! Hx-Hx
                
                do ind = irw_q(2*indi),irw_q(2*indi+1)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj  ) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj  
                        exit
                    endif 
                enddo
                
                indm = getUpperTriIndex(6,2*i  ,2*j-1 )  ! Hx-Ex
                
                do ind = irw_q(2*indi),irw_q(2*indi+1)-1
                    
                    if ( (col_q(ind) == 0) .or. (col_q(ind) == 2*indj-1  ) ) then  
                        val_q(ind) = val_q(ind) + mtx(indm)    
                        col_q(ind) = 2*indj-1  
                        exit
                    endif 
                enddo
                                   
            enddo
            
        enddo
    enddo
 
!
! Apply boundary conditions:
!
    do i = 1,2*nedges

        do i1 = irw_q(i),irw_q(i+1)-1
        
             j = col_q(i1)
             
             if ( ( (bcedge(i)>0) .or. (bcedge(j)>0) ) .and. (i /= j ) ) then
                val_q(i1) = 0.    
             elseif ( (bcedge(i)>0 ) .and. (i == j) ) then
                val_q(i1) = 1.   
             endif   
 
       enddo
    enddo
 
 
!
! Finish up by sorting each row by column number:
!
 
    do i = 1,2*nedges
        i0 = irw_q(i)
        i1 = irw_q(i+1)-1
        call quick_sort(col_q(i0:i1),zlist1=val_q(i0:i1)) 
    enddo
 
    deallocate (bcedge,eAdj)
 
         
    end subroutine gen_lhs_bump

!==================================================================================================================================!
!======================================================================================================================= fill_solmtx
!==================================================================================================================================!
    subroutine fill_solmtx
 
!
! Computes the kx domain fields at the receivers for the current wavenumber
! using a total field formulation.
!
    implicit none
    
! local variables
    real(8)                     :: yp,zp, qwt 
    integer                     :: iRx,e, iTx, isign, iComp, iq, nq
    complex(8)                  :: exx,hxx,dexdy,dexdz,dhxdy,dhxdz, efield(3), hfield(3), expikdx   
    
    complex(8)                  :: sigx,sigy,sigz,dsig
    complex(8)                  :: gammay2,gammaz2  
    
        
! FE Nodal quantities:
    integer,   dimension(3) :: n 
    complex(8),dimension(3) :: fex,fhx,fiex,fihx 
    real(8),   dimension(3) :: ye,ze  
 
    
    if (lprintDebug_em2dkx) write(*,*) myID,': fill_solmtx ...'  
        
!
! Loop over the sites and fill in the kx domain solutions:
!
    do iComp = 1,nComponents
        
        Fields(iComp) =  0d0  ! initialize
 
        iRx = Components(iComp,1) 
        iTx = Components(iComp,2)

        nq = 1
        if (lengthRx(iRx) > 0) nq = nquad_Rx     
        
        if (Components(iComp,3)>3) nq = 1 ! don't integrate mags
        
        do iq = 1,nq ! loop over quadrature points
        
            yp = yRxQ(iRx,iq)
            zp = zRxQ(iRx,iq)
    
!
! Get the site element index and material params:
!
            e = eRx(iRx,iq)
        
            call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
            gammay2 = kx2-ommu*(sigy - iomeps)
            gammaz2 = kx2-ommu*(sigz - iomeps)         

!
! Extract solution and gradient vectors:
!
! Nodes positions:
            n    = mesh%emap(1:3,e)
            ye   = mesh%y(n)
            ze   = mesh%z(n)     
      
! Get linear fields:
            fex = rhs(2*n-1,iTx) 
            fhx = rhs(2*n  ,iTx)     
        
! Interpolate field to site location:
            call fget_interp_lin(yp,zp,fex,ye,ze,fiex)
            call fget_interp_lin(yp,zp,fhx,ye,ze,fihx)

!
! Remap:
!
            exx   = fiex(1) 
            dexdy = fiex(2) 
            dexdz = fiex(3) 
            hxx   = fihx(1) 
            dhxdy = fihx(2) 
            dhxdz = fihx(3) 


            isign =  1
            if ((sType == 'cs') .and. (lNeg_ikx(iTx)) ) isign = -1

!
! Compute the total fields at current wavenumber, assumes Rx is not located at location of point sources, so those terms
! ignored here
!
! Ex:
            efield(1) = exx
      
! Ey
            efield(2) = (-ommu*dhxdz - isign*ikx*dexdy )/gammay2               
    
! Ez
            efield(3) = ( ommu*dhxdy - isign*ikx*dexdz )/gammaz2
                    
! Hx
            hfield(1) = hxx             

! Hy
            hfield(2) = (-isign*ikx*dhxdy - (sigz - iomeps)*dexdz )/gammaz2                
                                                    
! Hz
            hfield(3) = (-isign*ikx*dhxdz + (sigy - iomeps)*dexdy )/gammay2  
      
      
! Add on primary fields if using scattered field formulation:
            if ( (sType == 'mt').and.(lMTscatteredField)) then
        
! add on additional term to get scattered field for ey component:
                dsig      = sigy - Scat1DTM%getSig(zp)   
                efield(2) = efield(2) + (dsig/sigy)*Scat1DTM%getEH1D(zp,'e')  

                efield(1) = efield(1) + Scat1DTE%getEH1D(zp,'e')
                hfield(2) = hfield(2) + Scat1DTE%getEH1D(zp,'h')

                hfield(1) = hfield(1) + Scat1DTM%getEH1D(zp,'h')
                efield(2) = efield(2) - Scat1DTM%getEH1D(zp,'e') ! note minus sign for hx/ey

            endif
            
! For CSEM data, apply reciprocity scaling for mag Tx with elec Rx, or elec Tx with mag Rx:
            if ((lUseReciprocity).and.(sType == 'cs')) then
            
                select case (trim(TxType(iTx)))       
                case ('edipole')               
                    select case (trim(RxType(iRx)))
                    case('bdipole')
                         hfield =  ic*w*hfield    
                    case default ! nothing to do
                    end select                 
            
                case ('bdipole')
                    select case (trim(RxType(iRx)))
                    case('edipole')
                        efield =  -ic/w*efield    
                    case default ! nothing to do
                    end select               
                                                 
                end select             
   
            endif
            
!
! Get desired component at requested angle using the momentRx vector:
!
            if (nq > 1) then
                expikdx = exp(+isign*ikx*(xRxQ(iRx,iq)-xRx(iRx))) 
                qwt     = quad_weights_Rx(iq)/2.d0
            else
                expikdx = 1d0
                qwt     = 1d0
            endif  
            
            select case (RxType(iRx))
            case('edipole')
                Fields(iComp) = Fields(iComp) + sum(momentRx(1:3,iRx)*efield)*qwt*expikdx
            case('bdipole')
                Fields(iComp) = Fields(iComp) + sum(momentRx(1:3,iRx)*hfield)*qwt*expikdx
                
            end select        
           
        enddo ! iq = 1,nq loop over quadrature integration points
                        
    enddo  ! loop over nComponents
    
    end subroutine fill_solmtx
    
!==================================================================================================================================!
!=============================================================================================================== comp_drw_error_proj
!==================================================================================================================================!
    subroutine  comp_drw_error_proj( )
!
! Routine to compute  | F(delta_h) - B(u_h,delta_h) |_{\tau}
! where u_h is linear FE solution and delta_h is dual error in quadratic
! bump space.
!
  
    integer                 :: nele, e, n(3), iTx,k
    real(8)                 :: isign
    complex(8)              :: gammay2,gammaz2,sigx,sigy,sigz        
    complex(8)              :: aey,aez,bex,cey,cez,ahy,ahz,bhx,chy,chz ! pde coefficients for e and h terms
    real(8)                 :: areanew, erre,errh

    real(8), dimension(3)                :: ye,ze,yeq,zeq,a,b,c 
    complex(8), dimension(3)             :: ex,hx,deltaE,deltaH,feli,fhli,feqi,fhqi 
    
    real(8), dimension(:), allocatable    :: yp,zp
    integer, dimension(:), allocatable    :: inme
    complex(8), dimension(:), allocatable :: ee,hh 

    real(8), dimension(:), allocatable :: sigTemp
    
    complex(8), dimension(4) :: quade,quadh 
    real(8), dimension(3)    :: bary
  
    integer :: nqp 
    real(8), dimension(4,3) :: qp
    real(8), dimension(4)   :: qwt 
        
    real(8) :: yq,zq
 
    nqp = 4
    qp(1,:) = [ 1d0/3d0, 1d0/3d0, 1d0/3d0 ]
    qp(2,:) = [ 3d0/5d0, 1d0/5d0, 1d0/5d0 ]
    qp(3,:) = [ 1d0/5d0, 3d0/5d0, 1d0/5d0 ]
    qp(4,:) = [ 1d0/5d0, 1d0/5d0, 3d0/5d0 ]  
    qwt = [-0.5625d0, 0.520833333333333d0, 0.520833333333333d0, 0.520833333333333d0]

!    nqp = 1
!    qp(1,:) = [ 1d0/3d0, 1d0/3d0, 1d0/3d0 ]
!    qwt = [1d0,0d0,0d0,0d0]
            
    if (lprintDebug_em2dkx) write(*,*) myID,': comp_drw_error_proj ...'  
       
!
! Get list of elements in old mesh that contain centroids of new elements:
!
    nele = newmesh%nele
    
    allocate( yp(nele), zp(nele), inme(nele))
 
    yp = ( newmesh%y(newmesh%emap(1,:)) +  newmesh%y(newmesh%emap(2,:)) +  newmesh%y(newmesh%emap(3,:)) ) /3d0
    zp = ( newmesh%z(newmesh%emap(1,:)) +  newmesh%z(newmesh%emap(2,:)) +  newmesh%z(newmesh%emap(3,:)) ) /3d0
    
    allocate (sigTemp(mesh%nele))
    do e = 1,mesh%nele
        call getSigs(e,sigx,sigy,sigz)
        sigTemp(e) = dble(sigx+sigy+sigz)/3d0
    enddo
    
    call findElement(tree, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri,  &
    & sigTemp, newmesh%nele,yp, zp,inme)
    
    deallocate(sigTemp)
   
    allocate( ee(nTx),hh(nTx)  )
 
!
! Loop over number of elements in newmesh:
!
    do e = 1,nele 
            
!
! Generate the element coefficients from the old mesh:
!
        call getSigs(inme(e),sigx,sigy,sigz) ! gets anisotropic conductivity components for element e
        gammay2 = kx2-ommu*(sigy - iomeps)
        gammaz2 = kx2-ommu*(sigz - iomeps) 
        
        aey = (sigy - iomeps) / gammay2
        aez = (sigz - iomeps) / gammaz2
        bex = (sigx - iomeps)  
        cey =  ikx / gammay2  
        cez =  ikx / gammaz2
        
        ahy =  ommu / gammay2  
        ahz =  ommu / gammaz2  
        bhx =  ommu  
        chy =  ikx / gammay2
        chz =  ikx / gammaz2         
                             
        n  = newmesh%emap(1:3,e)
        ye = newmesh%y(n)        
        ze = newmesh%z(n)                   
        call get_abc_coeffs(ye,ze,a,b,c,areanew)  
        
        n  = mesh%emap(1:3,inme(e))
        yeq = mesh%y(n)
        zeq = mesh%z(n)                       

        do iTx = 1,nTx
 
            ex     = (rhs( 2*n-1 ,iTx)) ! linear
            hx     = (rhs( 2*n   ,iTx)) 
                                 
            deltaE = (rhs_q(2*edges(1:3,inme(e))-1,iTx))   ! bump
            deltaH = (rhs_q(2*edges(1:3,inme(e))  ,iTx))     
 
            isign =  1d0
            if (lNeg_ikx(iTx))   isign = -1d0
        
            do k = 1,nqp
   
! Get quadrature point (in element of the new mesh):
                bary = qp(k,1:3)
                yq   = sum(bary*ye)
                zq   = sum(bary*ze)
           
! Get linear fields:
                call fget_interp_lin(yq,zq,ex,yeq,zeq,feli)   
                call fget_interp_lin(yq,zq,hx,yeq,zeq,fhli)
                    
! Get bump fields:
                call fget_interp_bump(yq,zq,deltaE,yeq,zeq,feqi)
                call fget_interp_bump(yq,zq,deltaH,yeq,zeq,fhqi)                                     
      
! Evaluate form:
                quade(k) = aey*feli(2)*feqi(2) + aez*feli(3)*feqi(3) + bex*feli(1)*feqi(1) + isign*(cez*fhli(2)*feqi(3)-cey*fhli(3)*feqi(2))  
                quadh(k) = ahz*fhli(2)*fhqi(2) + ahy*fhli(3)*fhqi(3) + bhx*fhli(1)*fhqi(1) + isign*(chy*feli(2)*fhqi(3)-chz*feli(3)*fhqi(2))   
                                     
                quade(k) =  (quade(k))*qwt(k)
                quadh(k) =  (quadh(k))*qwt(k)
         
            enddo
           
            ee(iTx) = sum((quade))*areanew  
            hh(iTx) = sum((quadh))*areanew
         
        enddo ! iTx = 1,nTx
 
        erre =  maxval((abs(ee)))       ! ignoring F(delta) term since point sources only contribute as delta functions
        errh =  maxval((abs(hh))) 
        errnrm(e) =  max(erre,errh)  
              
    enddo ! do e = 1,nele   in newmesh
 
    deallocate(ee,hh,yp,zp,inme)
    
    end subroutine comp_drw_error_proj
   
!==================================================================================================================================!
!======================================================================================================================== refineMesh
!==================================================================================================================================!
    subroutine refineMesh()  
!
! Refines mesh associated with iTx, updates edges and site element lists.
!
    implicit none
    
    integer               :: i,e
    real(8)               :: myarea

    integer, dimension(:),allocatable  ::   indx    

    if (lprintDebug_em2dkx) write(*,*) myID,': refineMesh...'   
    
!
! Compute area of each element and set to 1/2 if element needs refinement
!
 
    allocate(indx(newmesh%nele))
    call indexx(newmesh%nele, errnrm,indx)    
    
    allocate (newmesh%area(newmesh%nele))
    newmesh%area = -1d0 
 
    do i = newmesh%nele,nint(newmesh%nele*(1d0 - pct_refine/100d0)),-1
 
        e = indx(i)    
        
        if (sqrt(errnrm(e)) > errortolerance/1d2/2d0) then  

            call getArea(e,myarea,1)
 
            if ( ( myarea > minArea ) ) then ! in case refinement goes bonkers, limit smallest elements possible
                  
                newmesh%area(e) = myarea/min(max(4d0, (sqrt(errnrm(e))/(errortolerance/1d2) )),6d0) 
                
            endif

        endif
    enddo

!
! Refine the mesh:
!
    call call_triangle(retricommand, newmesh )
    
    nRefinements = nRefinements + 1

    deallocate(indx)

    end subroutine refineMesh     
    
!==================================================================================================================================!
!=========================================================================================================== factor_lhs_bump_superlu
!==================================================================================================================================!
    subroutine factor_lhs_bump_superlu()
!
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
    
    implicit none    

    if (lprintDebug_em2dkx) write(*,*) myID,': factor_lhs_bump_superlu ...' 
    
    call superlu_zfactor(superlu_q,val_q,col_q,irw_q)
      
    if (superlu_q%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr factor_lhs_bump_superlu' 
        write(*,*) ' Error code = ',superlu_q%error
        write(*,*) ' Stopping !!!'
        stop
    endif
    
    end subroutine factor_lhs_bump_superlu
!==================================================================================================================================!
!========================================================================================================== factor_lhs_bump_intelmkl
!==================================================================================================================================!
    subroutine factor_lhs_bump_intelmkl()

    
    implicit none    

    if (lprintDebug_em2dkx) write(*,*) myID,': factor_lhs_bump_intelmkl ...' 
    
    call intelmkl_zfactor(intelmkl_q,val_q,col_q,irw_q)
    
    end subroutine factor_lhs_bump_intelmkl
        
!==================================================================================================================================!
!========================================================================================================================== write_Ab
!==================================================================================================================================!
    subroutine write_Ab()
!
! Writes A and b from Ax=b out to file.
! Matrix A is saved as:  irow,jcol, real, imag
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
        
        
    if (sType == 'cs') then
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) 2*mesh%nnod+1,nnz
        do i = 1,2*mesh%nnod+1
            write(16,*) irw(i)
        enddo                 
        do j = 1,nnz
            write(16,*) col(j), dble(val(j)), aimag(val(j)) 
        enddo
        close(16)       
                  
        
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) 2*mesh%nnod,  nTx
        do j = 1,nTx
            do i = 1,2*mesh%nnod
                write(16,*)   dble(rhs(i,j)), aimag(rhs(i,j)) 
            enddo
        enddo
        close(16)                      
        
    else  ! MT using decoupled linear systems:
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx_lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*)nnz_lH
        do i = 1,mesh%nnod
            do j = irw_lH(i),irw_lH(i+1)-1
                write(16,*) i, col_lH(j), dble(val_lH(j)), aimag(val_lH(j)) 
            enddo
        enddo
        close(16)           
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx_rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i,1)), aimag(rhs(2*i,1)) 
        enddo
        close(16)    
            
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex_lhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) nnz_lE
        do i = 1,mesh%nnod
            do j = irw_lE(i),irw_lE(i+1)-1
                write(16,*) i, col_lE(j), dble(val_lE(j)), aimag(val_lE(j)) 
            enddo
        enddo
        close(16)           
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex_rhs'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i-1,1)), aimag(rhs(2*i-1,1)) 
        enddo
        close(16)    
                               
    end if
               
    end subroutine write_Ab
!==================================================================================================================================!
!=========================================================================================================================== write_x
!==================================================================================================================================!
    subroutine write_x()
!
! Writes x from Ax=b out to file.
! Matrix A is saved as:  irow,jcol, real, imag
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
        
        
    if (sType == 'cs') then  ! CSEM dipole
    
        
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.x'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) 2*mesh%nnod,nTx
        do j = 1,nTx
            do i = 1,2*mesh%nnod
                write(16,*)   dble(rhs(i,j)), aimag(rhs(i,j)) 
            enddo
        enddo
        close(16)                      
        
    else ! MT using decoupled linear systems:
    
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.hx'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i,1)), aimag(rhs(2*i,1)) 
        enddo
        close(16)    
            
        
    
        filename  = trim(fileroot)//'.'//trim(cadapt)//'.ex'
        open(16,file=trim(filename),status='REPLACE') 
        write(16,*) mesh%nnod
        do i = 1,mesh%nnod
            write(16,*) dble(rhs(2*i-1,1)), aimag(rhs(2*i-1,1)) 
        enddo
        close(16)    
                               
    end if         
               
    end subroutine write_x
                    
!==================================================================================================================================!
!====================================================================================================================== solve_primal
!==================================================================================================================================!
    subroutine solve_primal 
!
! Solves the linear system AX=B for CSEM and MT problems using optional
! methods
!
! Kerry Key
! Scripps Institution of Oceanography
!

    integer :: iTx
    

    integer*8, save :: solve_primal = -1
	
    call SCOREP_F_Begin(solve_primal,"solve_primal",0,"",4988) 

    
    if (lSaveLinearSystem) call write_Ab
             
    if (sType == 'cs') then ! CSEM dipole
    
        select case (trim(linearSolver))
        case ('superlu')
            call solve_primal_superlu
        case('intelmkl')
            call solve_primal_intelmkl  
        end select
         
! Matrix has been factored, so we can deallocate the CSR arrays:
        deallocate(val,col,irw)       
        
        if (lSaveLinearSystem) call write_x
        
! Remember that for negative wavenumbers we solved [A B; B^T C][-u v] = [-c d] so take negative of u
        do iTx = 1,nTx
            if (lNeg_ikx(iTx)) then
                rhs(1:2*mesh%nnod:2,iTx) = -rhs(1:2*mesh%nnod:2,iTx)  ! note Fortan increment is at end, not middle like matlab...
            endif        
        enddo 
        
    else ! MT using decoupled linear systems:
    
        select case (trim(linearSolver)) 
        case ('superlu')
            call solve_primal_superlu_mt
        case('intelmkl')
            call solve_primal_intelmkl_mt  
        end select  
        
! Matrices have been factored, so we can deallocate the CSR arrays:
        deallocate(val_lE,col_lE,irw_lE)  
        deallocate(val_lH,col_lH,irw_lH)  
        
        if (lSaveLinearSystem) call write_x
                      
    end if


    call SCOREP_F_RegionEnd( solve_primal ) !end Instrumentation

                 
    end subroutine solve_primal   
                    

!==================================================================================================================================!
!============================================================================================================== solve_primal_superlu
!==================================================================================================================================!
    subroutine solve_primal_superlu()
!
! Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].
! The lhs A is passed in in compressed sparse row format.
! SuperLU takes the lhs and factors it, then  solves for all vectors
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the SuperLU factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new superlu_solver module
!

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_superlu ...'
          
!
! Factor the matrix first:
!
    call superlu_zfactor(superlu_p,val,col,irw)
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! Now solve  AX=B for all column vectors in B:
!
    call superlu_zsolve (superlu_p,rhs )
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif
 
    end subroutine solve_primal_superlu
    
!==================================================================================================================================!
!============================================================================================================= solve_primal_intelmkl
!==================================================================================================================================!
    subroutine solve_primal_intelmkl()

    implicit none

    integer :: j
    complex(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_intelmkl ...'
          
!
! Factor the matrix first:
!
    call intelmkl_zfactor(intelmkl_p,val,col,irw)
    
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs,1)) )
    
    do j = 1,size(rhs,2)  
 
        x = rhs(:,j)     
        call intelmkl_zsolve(intelmkl_p,x )
        rhs(:,j) = x
 
    enddo
 
    deallocate(x)
 
    end subroutine solve_primal_intelmkl
        
!==================================================================================================================================!
!=========================================================================================================== solve_primal_superlu_mt
!==================================================================================================================================!
    subroutine solve_primal_superlu_mt()
!
! Solves AX=B for all vectors in matrix B=rhs=[rhs1, rhs2,...].
! The lhs A is passed in in compressed sparse row format.
! SuperLU takes the lhs and factors it, then  solves for all vectors
! in X corresponding to vectors in the rhs. X is returned in rhs to
! conserve memory.  I leave the SuperLU factors in memory since the error
! estimator will use these for solving the dual problem later on.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! May 2011: Updated to use my new superlu_solver module
!

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_superlu_mt ...'
        
!
! Decoupled Matrix for MT problem:
!
! Factor  E matrix
!
    call superlu_zfactor(superlu_p_E,val_lE,col_lE,irw_lE)
    
    if (superlu_p_E%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! Factor  H matrix
!
    call superlu_zfactor(superlu_p_H,val_lH,col_lH,irw_lH)
    
    if (superlu_p_H%error /= 0) then
        write(*,*) ' Error in superlu_zfactor from calling in subr  solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif    
!
! Now solve  AX=B for all column vectors in B:
!
! E system:
!
    rhs_lE = rhs(1:2*mesh%nnod:2,1:1) 
    
    call superlu_zsolve (superlu_p_E,rhs_lE )
    
    if (superlu_p_E%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_E%error
        write(*,*) ' Stopping !!!'
        stop
    endif
!
! H system:
!
    rhs_lH = rhs(2:2*mesh%nnod:2,1:1) 
    
    call superlu_zsolve (superlu_p_H,rhs_lH )  
    if (superlu_p_H%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_primal_superlu()' 
        write(*,*) ' Error code = ',superlu_p_H%error
        write(*,*) ' Stopping !!!'
        stop
    endif 
   
! Push these back into the full RHS vector:
    rhs(1:2*mesh%nnod:2,1) = rhs_lE(1:mesh%nnod,1)
    rhs(2:2*mesh%nnod:2,1) = rhs_lH(1:mesh%nnod,1) 
    
    
    end subroutine solve_primal_superlu_mt

!==================================================================================================================================!
!========================================================================================================== solve_primal_intelmkl_mt
!==================================================================================================================================!
    subroutine solve_primal_intelmkl_mt()

    implicit none

    complex(8) , dimension(:), allocatable :: x
    
    integer :: n
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_primal_intelmkl_mt ...'
 
 
!
! Factor and solve  AX=B for all column vectors in B:
!
    n = size(rhs,1)/2
    
    allocate ( x(n) )
    
! TE Ex:
    if (l_has_TE_mode) then
        x = rhs(1:2*mesh%nnod:2,1)     
        call intelmkl_zfactor(intelmkl_p_E,val_lE,col_lE,irw_lE)  
        call intelmkl_zsolve(intelmkl_p_E,x )
        rhs(1:2*mesh%nnod:2,1)  = x
    endif
    
! TM Hx:
    if (l_has_TM_mode) then
         x = rhs(2:2*mesh%nnod:2,1) 
        call intelmkl_zfactor(intelmkl_p_H,val_lH,col_lH,irw_lH)   
        call intelmkl_zsolve(intelmkl_p_H,x )
        rhs(2:2*mesh%nnod:2,1)  = x        
    endif 
    
    deallocate(x)
 
    end subroutine solve_primal_intelmkl_mt
    
    
!==================================================================================================================================!
!============================================================================================================= solve_varepsilon
!==================================================================================================================================!
    subroutine solve_varepsilon(imode)
 
    integer,intent(in) :: imode ! 1 factor and solve, otherwise just solve
    
    integer :: iTx
    
    if (lSolveBumpGaussSeidel) then
        
        call solve_varepsilon_GaussSeidel
        
    else
    
        select case (trim(linearSolver))
    
        case ('superlu')
           if (imode == 1) call factor_lhs_bump_superlu
           call solve_varepsilon_superlu
      
        case ('intelmkl')
           if (imode == 1) call factor_lhs_bump_intelmkl
           call solve_varepsilon_intelmkl     

        case default 
            write(*,*) 'Error, no linear system solver parameter specified, stopping!'
            stop
              
        end select    
        
    endif
 
! Remember that for negative wavenumbers we solved [A B; B^T C][-u v] = [-c d] so take negative of u
    if (sType == 'cs') then ! CSEM dipole
        do iTx = 1,nTx
            if (lNeg_ikx(iTx)) then
                rhs_q(1:2*nedges:2,iTx) = -rhs_q(1:2*nedges:2,iTx)  ! note Fortan increment is at end, not middle like matlab...
            endif        
        enddo         
    endif 
                
    end subroutine solve_varepsilon     
    
!==================================================================================================================================!
!====================================================================================================== solve_varepsilon_GaussSeidel
!==================================================================================================================================!
    subroutine solve_varepsilon_GaussSeidel
 
    implicit none

    integer :: i,j,iTx,iter,ival
    
    complex(8) , dimension(:), allocatable :: x,b,res
    
    complex(8)  :: diag,xk
    real(8)     :: maxRel !, omega
    
    integer, parameter :: maxNumIter = 5
    real(8) :: tolGS 
    
    tolGS = errortolerance/1d2
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_GaussSeidel nhrs = ...',size(rhs_q,2)
 
    allocate ( x(size(rhs_q,1)) , b(size(rhs_q,1)),res(size(rhs_q,1)))
    
    do iTx = 1,size(rhs_q,2)  
      
        b       = rhs_q(:,iTx)     
        x       = 0d0
            
        do iter = 1, maxNumIter
            
            maxRel = -1
            
            do i = 1, 2*nedges
 
                xk = b(i)

                do ival = irw_q(i), irw_q(i+1) - 1
        
                    j = col_q(ival) 
        
                    if (j /= i) then
                        xk = xk - val_q(ival) * x(col_q(ival))
                    else
                        diag = val_q(ival) 
                    endif
        
                enddo
                
                xk = xk / diag 
                
                maxRel = max( (abs(x(i) - xk)/(max(abs(x(i)),1d-80))),maxRel)
              
! Gauss-s=Seidel:
                 x(i) = xk
                     
! SOR
!                omega = 1.1d0
!                x(i) = omega*xk +(1d0-omega)*x(i)
                
            enddo
!
!             do i = 1, 2*nedges
!                xk = 0
!                do ival = irw_q(i), irw_q(i+1) - 1
!                    j = col_q(ival)
!                    xk  =  xk  + val_q(ival) * x(col_q(ival))
!                enddo
!                res(i) = b(i) - xk
!            enddo
!            write(*,*) iter, maxRel, sqrt(dble(dot_product(res,res)))
!!
            if (maxRel < tolGS) exit
            
        enddo ! iter = 1, maxNumIter

        rhs_q(:,iTx) = x
 
    enddo
    
    deallocate(x,b,res)

    end subroutine solve_varepsilon_GaussSeidel

!==================================================================================================================================!
!========================================================================================================== solve_varepsilon_superlu
!==================================================================================================================================!
    subroutine solve_varepsilon_superlu
!
!  Written by Kerry Key
!  Scripps Institution of Oceanography
!
! May 2011:         Updated to use my new superlu_solver module
! Summer 2008:      Implemented
 

    implicit none

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_superlu ...'
    
!
! Solves AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_q derived type.
!

    call superlu_zsolve (superlu_q,rhs_q)
    
    if (superlu_q%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_varepsilon_superlu()' 
        write(*,*) ' Error code = ',superlu_q%error
        write(*,*) ' Stopping !!!'
        stop
    endif

    end subroutine solve_varepsilon_superlu

!==================================================================================================================================!
!========================================================================================================= solve_varepsilon_intelmkl
!==================================================================================================================================!
    subroutine solve_varepsilon_intelmkl
 
    implicit none

    integer :: j
    complex(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_varepsilon_intelmkl ...'
 
    allocate ( x(size(rhs_q,1)) )
    
    do j = 1,size(rhs_q,2)  
 
        x = rhs_q(:,j)     
        
        call intelmkl_zsolve(intelmkl_q,x )
        
        rhs_q(:,j) = x
 
    enddo
    
!    write(*,*) rhs_q(1:3,1)
    deallocate(x)
    
    end subroutine solve_varepsilon_intelmkl
    
    
    
!==================================================================================================================================!
!======================================================================================================================== solve_dual
!==================================================================================================================================!
    subroutine solve_dual 
!
!  Solve dual problem
!
    select case (trim(linearSolver))
    
    case ('superlu')
        
        if (sType == 'mt') then
            call solve_dual_superlu_mt 
        else
            call solve_dual_superlu
        endif

    case('intelmkl')
        
        if (sType == 'mt') then
            call solve_dual_intelmkl_mt  
        else
            call solve_dual_intelmkl
        endif    
                       
    end select    
  

    end subroutine solve_dual
   
!==================================================================================================================================!
!================================================================================================================ solve_dual_superlu
!==================================================================================================================================!
    subroutine solve_dual_superlu
!
!  Solve dual problem for error weighting with superlu.
!
! Kerry Key
! Scripps Institution of Oceanography
!
    
    implicit none
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_superlu ...'
     
!
! Solve AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_p derived type.
!
 
    call superlu_zsolve (superlu_p,rhs_wh)
    
    if (superlu_p%error /= 0) then
        write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
        write(*,*) ' Error code = ',superlu_p%error
        write(*,*) ' Stopping !!!'
        stop
    endif

    end subroutine solve_dual_superlu
    
!==================================================================================================================================!
!=============================================================================================================== solve_dual_intelmkl
!==================================================================================================================================!
    subroutine solve_dual_intelmkl()

    implicit none

    integer :: j
    complex(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_intelmkl ...'
 
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs_wh,1)) )
    
    do j = 1,size(rhs_wh,2)  
 
        x = rhs_wh(:,j)     
        
        call intelmkl_zsolve(intelmkl_p,x )
        
        rhs_wh(:,j) = x
 
    enddo
 
    deallocate(x)
 
    end subroutine solve_dual_intelmkl
            
!==================================================================================================================================!
!============================================================================================================= solve_dual_superlu_mt
!==================================================================================================================================!
    subroutine solve_dual_superlu_mt
!
! Solve dual problem for error weighting with superlu.
!
! Kerry Key
! Scripps Institution of Oceanography
!
    
    implicit none
    integer :: j  

    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_superlu_mt ...'
  
!
! Solve AX=B for all column vectors in B:
! Assumes factors have already be computed and stored in superlu_p derived type.
!
    do j = 1,size(rhs_wh,2)  
    
        rhs_lE(:,1) = rhs_wh(1:2*mesh%nnod:2,j) 
        call superlu_zsolve (superlu_p_E,rhs_lE)
        
        if (superlu_p_E%error /= 0) then
            write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
            write(*,*) ' Error code = ',superlu_p_E%error
            write(*,*) ' Stopping !!!'
            stop
        endif
    
        rhs_lH(:,1) = rhs_wh(2:2*mesh%nnod:2,j) 
        call superlu_zsolve (superlu_p_H,rhs_lH)
        
        if (superlu_p_H%error /= 0) then
            write(*,*) ' Error in superlu_zsolve from calling in subr solve_dual_superlu()' 
            write(*,*) ' Error code = ',superlu_p_H%error
            write(*,*) ' Stopping !!!'
            stop
        endif
    
!
! Fill back into rhs_wh array
!
        rhs_wh(1:2*mesh%nnod:2,j) = rhs_lE(:,1)
        rhs_wh(2:2*mesh%nnod:2,j) = rhs_lH(:,1)  
    
    enddo 
    
    end subroutine solve_dual_superlu_mt

!==================================================================================================================================!
!============================================================================================================ solve_dual_intelmkl_mt
!==================================================================================================================================!
    subroutine solve_dual_intelmkl_mt()

    implicit none

    integer :: j, n
    complex(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_em2dkx) write(*,*) myID,': solve_dual_intelmkl_mt ...'
 
!
! Now solve  AX=B for all column vectors in B:
!
    n = size(rhs_wh,1)/2
    
    allocate ( x(n) )
    
    do j = 1,size(rhs_wh,2)  
 
! TE Ex:
        if (l_has_TE_mode) then           
            x = rhs_wh(1:2*mesh%nnod:2,j)     
            call intelmkl_zsolve(intelmkl_p_E,x )
            rhs_wh(1:2*mesh%nnod:2,j) = x
        endif 
        
! TM Hx:
        if (l_has_TM_mode) then  
            x = rhs_wh(2:2*mesh%nnod:2,j)     
            call intelmkl_zsolve(intelmkl_p_H,x )
            rhs_wh(2:2*mesh%nnod:2,j) = x
        endif
    enddo
 
    deallocate(x)
 
    if (lprintDebug_em2dkx) write(*,*) myID,': leaving solve_dual_intelmkl_mt ...'
     
    end subroutine solve_dual_intelmkl_mt

!==================================================================================================================================!
!=============================================================================================================== legendre_compute_dr
!==================================================================================================================================!
    subroutine legendre_compute_dr ( order, xtab, weight )

!*****************************************************************************80
!
!! LEGENDRE_COMPUTE_DR computes a Gauss-Legendre rule, Davis-Rabinowitz method.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ORDER, the order of the rule.
!    ORDER must be greater than 0.
!
!    Output, real ( kind = 8 ) XTAB(ORDER), the abscissas.
!
!    Output, real ( kind = 8 ) WEIGHT(ORDER), the weights.
!    The weights are positive, symmetric, and should sum to 2.
!
    implicit none
    
    integer ( kind = 4 ) order
    
    real    ( kind = 8 ) d1
    real    ( kind = 8 ) d2pn
    real    ( kind = 8 ) d3pn
    real    ( kind = 8 ) d4pn
    real    ( kind = 8 ) dp
    real    ( kind = 8 ) dpn
    real    ( kind = 8 ) e1
    real    ( kind = 8 ) fx
    real    ( kind = 8 ) h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) iback
    integer ( kind = 4 ) k
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mp1mi
    integer ( kind = 4 ) ncopy
    integer ( kind = 4 ) nmove
    real    ( kind = 8 ) p
    real    ( kind = 8 ) :: pi = 3.141592653589793D+00
    real    ( kind = 8 ) pk
    real    ( kind = 8 ) pkm1
    real    ( kind = 8 ) pkp1
    real    ( kind = 8 ) t
    real    ( kind = 8 ) u
    real    ( kind = 8 ) v
    real    ( kind = 8 ) x0
    real    ( kind = 8 ) xtab(order)
    real    ( kind = 8 ) xtemp
    real    ( kind = 8 ) weight(order)
    
    if ( order < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_COMPUTE - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
    stop
    end if
    
    e1 = real ( order * ( order + 1 ), kind = 8 )
    
    m = ( order + 1 ) / 2
    
    do i = 1, m
    
    mp1mi = m + 1 - i
    
    t = real ( 4 * i - 1, kind = 8 ) * pi &
      / real ( 4 * order + 2, kind = 8 )
    
    x0 = cos ( t ) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 &
      / real ( order, kind = 8 ) ) &
      / real ( 8 * order * order, kind = 8 ) )
    
    pkm1 = 1.0D+00
    pk = x0
    
    do k = 2, order
      pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) &
        / real ( k, kind = 8 )
      pkm1 = pk
      pk = pkp1
    end do
    
    d1 = real ( order, kind = 8 ) * ( pkm1 - x0 * pk )
    
    dpn = d1 / ( 1.0D+00 - x0 * x0 )
    
    d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0 * x0 )
    
    d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) &
      / ( 1.0D+00 - x0 * x0 )
    
    d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) / &
      ( 1.0D+00 - x0 * x0 )
    
    u = pk / dpn
    v = d2pn / dpn
!
!  Initial approximation H:
!
    h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v * v - d3pn / &
      ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
    p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00 &
      * ( d3pn + 0.25D+00 * h * d4pn ) ) )
    
    dp = dpn + h * ( d2pn + 0.5D+00 * h * ( d3pn + h * d4pn / 3.0D+00 ) )
    
    h = h - p / dp
    
    xtemp = x0 + h
    
    xtab(mp1mi) = xtemp
    
    fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00 &
      * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )
    
    weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp * xtemp ) / ( fx * fx )
    
    end do
    
    if ( mod ( order, 2 ) == 1 ) then
    xtab(1) = 0.0D+00
    end if
!
!  Shift the data up.
!
    nmove = ( order + 1 ) / 2
    ncopy = order - nmove
    
    do i = 1, nmove
    iback = order + 1 - i
    xtab(iback) = xtab(iback-ncopy)
    weight(iback) = weight(iback-ncopy)
    end do
!
!  Reflect values for the negative abscissas.
!
    do i = 1, order - nmove
    xtab(i) = - xtab(order+1-i)
    weight(i) = weight(order+1-i)
    end do
    
    end subroutine legendre_compute_dr


!
!real(8) function dist3D_Segment_to_Segment(s1p0,s1p1,s2p0,s2p1)
!
!!/ dist3D_Segment_to_Segment(): get the 3D minimum distance between 2 segments
!!//    Input:  two 3D line segments S1 and S2
!!//    Return: the shortest distance between S1 and S2
!!float
!!dist3D_Segment_to_Segment( Segment S1, Segment S2)
!!{
!    real(8), dimension(2), intent(in) :: s1p0,s1p1,s2p0,s2p1
!
!    real(8), dimension(2) :: u,v,w, dp
!    real(8)               :: a,b,c,d,e,dd,sc, sN, sD,tc, tN, tD
!
!    real(8), parameter :: small_num = 1d-3
!
!    u = s1p1 - s1p0
!    v = s2p1 - s2p0
!    w = s1p0 - s2p0
!    a = dot_product(u,u)
!    b = dot_product(u,v)
!    c = dot_product(v,v)
!    d = dot_product(u,w)
!    e = dot_product(v,w)
!    D = a*c - b*b        !// always >= 0
!
!    sD = D      ! // sc = sN / sD, default sD = D >= 0
!    tD = D      ! // tc = tN / tD, default tD = D >= 0
!
!   ! // compute the line parameters of the two closest points
!    if (D < small_num) then !{ // the lines are almost parallel
!        sN = 0.0           !// force using point P0 on segment S1
!        sD = 1.0           !// to prevent possible division by 0.0 later
!        tN = e
!        tD = c
!
!    else                  !// get the closest points on the infinite lines
!        sN = (b*e - c*d)
!        tN = (a*e - b*d)
!
!        if (sN < 0.0) then  ! // sc < 0 => the s=0 edge is visible
!            sN = 0.0
!            tN = e
!            tD = c
!
!        elseif (sN > sD) then !{  // sc > 1  => the s=1 edge is visible
!            sN = sD
!            tN = e + b
!            tD = c
!        endif
!    endif
!
!    if (tN < 0.0) then ! {            // tc < 0 => the t=0 edge is visible
!        tN = 0.0
!        !// recompute sc for this edge
!        if (-d < 0.0) then
!            sN = 0.0
!        elseif (-d > a) then
!            sN = sD
!        else
!            sN = -d
!            sD = a
!        endif
!    elseif (tN > tD) then ! {      // tc > 1  => the t=1 edge is visible
!        tN = tD
!        !// recompute sc for this edge
!        if ((-d + b) < 0.0) then
!            sN = 0
!        elseif ((-d + b) > a) then
!            sN = sD
!        else
!            sN = (-d +  b)
!            sD = a
!        endif
!    endif
!    !// finally do the division to get sc and tc
!!    sc = (abs(sN) < SMALL_NUM ? 0.0 : sN / sD);
!!    tc = (abs(tN) < SMALL_NUM ? 0.0 : tN / tD);
!    sc =   sN / sD
!    tc =   tN / tD
!
!    !// get the difference of the two closest points
!    dP = w + (sc * u) - (tc * v) !  // =  S1(sc) - S2(tc)
!
!    dist3D_Segment_to_Segment =  norm2(dP) !   // return the closest distance
!
!    end function dist3D_Segment_to_Segment
    
    end module em2dkx_mod
    

