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
! Subroutine dx2dkx computes the 2.5D DC resistivity potentials at a given wavenumber (kx) using 
! a goal-oriented adaptive finite element solution to the scalar potential equation for triaxial anisotropy.
! Tx electrodes are modeled as separate linear systems and thus they can
! have arbitrary locations in x,y,z since we don't assume x=0 here. 
!
! !!! KWK note: a posteriori error estimation and adaptive refinement not yet implemented! 
!               code currently does brute force a priori local mesh refinement around 
!               electrode spread.
!
!==================================================================================================================================! 
!======================================================================================================================== dc2dkx_mod
!==================================================================================================================================!     
    module dc2dkx_mod 

    use EM_constants
        
    use fem2d_utilities  ! for generic FE operations
    use kdtree2_module   ! for fast searches
    use triangle_mesh    ! for trimesh structure 
    
    use intelmkl_solver
    
    use quad 
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
    public :: dc2dkx, displayRefinementStats_dc2dkx


   
!-----------------------------------------------------------------------
! Private variables:
!-----------------------------------------------------------------------
!
! These are allocated and handled only inside dc2dkx and the private attribute 
! means you can't access these from outside the dc2dkx module.
!
!
! A few adjustable parameters (don't change these unless you know what you are doing): 
!

!
! Variables and arrays used across multiple subroutines:
!
    real(8), private                                 :: kx2
    integer, dimension(:), allocatable, private      :: eRx 
    integer, dimension(:), allocatable, private      :: eTx
    logical, dimension(:), allocatable, private      :: lCloud
    integer, private                                 :: nedges 
    integer, dimension(:,:), allocatable, private    :: edges
    integer, dimension(:), allocatable, private      :: node2tri

!
! Sparse LHS and RHS arrays for linear FE's
!
    real(8), dimension(:), allocatable, private      :: val
    integer, dimension(:), allocatable, private      :: col, irw
    integer, private                                 :: nnz, info
    real(8), dimension(:,:), allocatable, private    :: rhs, rhs_wh 
 
    type(mkl_dss_handle) :: intelmkl_p ! Derived type that holds pointers to intelmkl factorization for primary/dual system  
    
 
    logical, dimension(:),   allocatable, private      :: lBoundary 
    real(8), dimension(:),   allocatable, private      :: errnrm  
    real(8), dimension(:,:), allocatable, private      :: localError, dualError        
    character(32) :: retricommand
 
    
!
! Kd-tree pointer for fast element searches:
!   
    type(kdtree2), pointer, private  :: tree    

    
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
    real(8), dimension(:), allocatable  ::  fieldsLast   ! used for convergence test
 
    real(8) :: maxerr           ! maximum error at any site       
    logical :: lmadeNewMesh     ! true if a new mesh has been made       
    type (trimesh) :: newmesh   ! temp storage for new adaptive meshes 
    
    integer :: isubset,iRefinementGrp
    
                   
    contains
 
 
!==================================================================================================================================! 
!============================================================================================================================ dc2dkx
!==================================================================================================================================!
    subroutine dc2dkx(iRefinementGrp_in,isubset_in,lrefine)   
!
! The main subroutine for computing the 2.5D DC fields in the kx domain
! for the input transmitter, receivers, wavenumber of mesh

! Created April-May 2018
!
! Kerry Key
! Lamont-Doherty Earth Observatory
! Columbia University 
! 
    
    integer, intent(in) :: isubset_in,iRefinementGrp_in
    logical             :: lrefine
    real(8)             :: maxRelDiff, time
        
    if (lprintDebug_dc2dkx) write(*,*) myID,': entered dc2dkx...maxRefinements:'  ,maxRefinements   

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
    if (lLocalRefine) call localRefinement_dc2dkx()  

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
    
        call dc2dkx_initialize
 
        call cpu_time(t1)
    
    !
    ! Get element based lists (for node2tri, Rx and Tx electrode locations):
    !        
        call get_element_lists
        
        call cpu_time(t2)
    !
    ! Generate the node labels and boundary condition flags for nodes on mesh boundaries
    !
        call gen_lBoundary
                
        call cpu_time(t3)
        
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
    ! Extract solution at receiver electrodes:
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
    !    if (.not.lConverged)  call estimateerror()

        call cpu_time(t20)         

    !
    ! If requested, compute the adjoint partial derivatives:
    !    
        if ( lCompDerivs .and. lConverged )  call comp_adj_derivs      

        call cpu_time(t21)
    !
    ! Display timing info:
    !
        if (lprinttimers)  call dc2dkx_printTimers  
        
    !
    ! Display iteration stats:
    !  
        time = t21 - t0
        if (lDisplayRefinementStats) call displayRefinementStats_dc2dkx(iRefinementGrp, isubset, meshnumber,lConverged, & 
                                         & oldNumNodes,newNumNodes,maxerr,maxRelDiff,time)
            
    !
    ! Deallocate arrays:
    !   
        
        if (lmadeNewMesh) then
            call deallocate_trimesh(mesh,.false.)
            call copy_trimesh( newmesh, mesh )
        endif
        call deallocate_trimesh(newmesh,.false.)
        
        call dc2dkx_deallocate  
        
        
    enddo ! do while (.not.lConverged)  
       
    deallocate(fieldsLast) 
       
!
! All done, goodbye
!    
    if (lprintDebug_dc2dkx) write(*,*) myID,': leaving dc2dkx...'
    
    end subroutine dc2dkx
    
!==================================================================================================================================! 
!=================================================================================================================== comp_adj_derivs
!==================================================================================================================================! 
    subroutine comp_adj_derivs   
!
! Computes the model parameter sensitivity matrix using the adjoint method
!
       
    integer         ::  e, n(3), iTx, iparamnum(3), iRx, iComp
    
    real(8)         :: ye(3),ze(3),a(3),b(3),c(3),area,l(3),dldy(3), dldz(3) 
    real(8)         :: sxx, syy, szz, vp(3),dvpdy,dvpdz,va(3),dvady,dvadz, sigx,sigy,sigz

    real(8), dimension(3,3)     :: M != [ 2d0, 1d0, 1d0,   1d0, 2d0, 1d0,   1d0, 1d0, 2d0] ! mass matrix, column major order
    
     
    if (lprintDebug_dc2dkx) write(*,*) myID,': comp_adj_derivs...'
!
! fill matrix M
!

M(1,1) = 2d0
M(1,2) = 1d0
M(1,3) = 1d0
M(2,1) = 1d0
M(2,2) = 2d0
M(2,3) = 1d0
M(3,1) = 1d0
M(3,2) = 1d0
M(3,3) = 2d0
         
!
! initialize to zero:
!    
    dFieldsdRho  = 0d0          
 
    allocate(rhs_wh(mesh%nnod,1))   
            
    !
    ! Loop over receivers and compute adjoint sources, then sensitivities:
    !

    do iRx = 1,ntrodes_rx
            
        rhs_wh = 0d0
              
        !
        ! Set up the basis functions and evaluate at the receiver location:
        !
        e   = eRx(iRx)
        n   = mesh%emap(1:3,e)
        ye  = mesh%y(n)
        ze  = mesh%z(n)                           
        
        call get_linear_basis(ye,ze,electrode_rx(iRx,2),electrode_rx(iRx,3),l,dldy,dldz)               

        !
        ! Insert into RHS:
        ! 
        
        rhs_wh(n, 1) = l  
            
        !
        ! Solve the linear system:
        !       
        call solve_dual 

 
        !
        ! Loop over all elements and integrate V^a(-ikx) dot V(ikx) :
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
     
            !
            ! Loop over requested components:
            !
            do iComp = 1,nComponents
            
                if ( Components(iComp,1) /= iRx ) cycle 
            
                iTx = Components(iComp,2)  

                ! Get potential from source at transmitter:
                vp  = rhs(n,iTx) 
                
                ! Get potential from source at receiver:
                va = -rhs_wh(n,1) 
                      
                ! Get gradients of potential:      
  
                dvpdy =  sum(vp*dldy)
                dvpdz =  sum(vp*dldz)
                dvady =  sum(va*dldy)
                dvadz =  sum(va*dldz)

                ! Sxx sensitivity can be computed using a mass matrix:
                sxx  = kx2*sum(va*matmul(M,vp))*area/12d0

                ! Syy and Szz sensitivities:
                syy  = dvpdy*dvady*area
                szz  = dvpdz*dvadz*area
                
                ! Convert from dF/dsig to dF/dlog10(rho):
                call getSigs(e,sigx,sigy,sigz) ! gets anisotropic conductivity elements for element e
                sxx = -dlog(10.d0)*sigx*sxx
                syy = -dlog(10.d0)*sigy*syy
                szz = -dlog(10.d0)*sigz*szz

                ! insert into sensitivity array: 
                select case (trim(cAnisotropy))

                case ('isotropic')
                    dFieldsdRho(iComp,iparamnum(1)) =  dFieldsdRho(iComp,iparamnum(1)) + ( sxx + syy + szz )
        
                case ('triaxial')
                    dFieldsdRho(iComp,iparamnum(1)) =  dFieldsdRho(iComp,iparamnum(1))  + ( sxx )
                    dFieldsdRho(iComp,iparamnum(2)) =  dFieldsdRho(iComp,iparamnum(2))  + ( syy ) 
                    dFieldsdRho(iComp,iparamnum(3)) =  dFieldsdRho(iComp,iparamnum(3))  + ( szz ) 

                case ('tix') ! For transversely isotropic:
                     ! First param is sig along symmetry axis; second param is along transverse plane
                    dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + ( sxx ) 
                    dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + ( syy  + szz ) 

                case ('tiy')

                    dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + ( syy )
                    dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + ( sxx + szz )  

                case ('tiz')

                     dFieldsdRho(iComp,iparamnum(1)) = dFieldsdRho(iComp,iparamnum(1)) + ( szz ) 
                     dFieldsdRho(iComp,iparamnum(2)) = dFieldsdRho(iComp,iparamnum(2)) + ( syy + sxx ) 

                end select              
            
            enddo ! iComp
            
            
        enddo ! e = 1 nele      

    enddo  !  iRx

! print for debugging:
!     do iComp = 1,nComponents
!         do e=1,size(dFieldsdRho,2)
!             write(*,*) i,e,dFieldsdRho(iComp,e) 
!         enddo
!     enddo
    
    deallocate(rhs_wh)

    end subroutine comp_adj_derivs
  
!==================================================================================================================================! 
!================================================================================================================== checkConvergence
!==================================================================================================================================!
    subroutine checkConvergence(maxRelDiff,lrefine)
    
    real(8),intent(out) :: maxRelDiff
    logical,intent(in)  :: lrefine
    
    integer :: i 
    real(8) :: wt 
    
    lConverged = .false.
    
    oldNumNodes = mesh%nnod
    newNumNodes = 0
    maxRelDiff  = 0d0    
    
    !write(*,*) 'checkConvergence: lrefine',lrefine
    
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
    ! Compute the numerical difference in the computed potentials for each
    ! subsequent adaptive refinement iteration.
    ! 
    do i = 1,nComponents
    
        wt = 1.               
!         iRx = Components(i,1)
!         iTx = Components(i,2)
!         rTx  = sqrt( ( yRx(iRx)- yTx(iTx))**2 + ( zRx(iRx)- zTx(iTx))**2 )
!         if (rTx < minRangeProtector)  wt = 0.
        
        fieldsLast(i) = wt*abs(fieldsLast(i) - fields(i))/(abs(fields(i)))
              
    enddo
  
             
    maxRelDiff = maxval(abs(fieldsLast))
    
    !write(*,*) 'maxRelDiff:',maxRelDiff

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
    subroutine displayRefinementStats_dc2dkx( iRefinementGrp,isubset, meshnumber,lconv, oldnodes,newnodes,sumerr,maxRelDiff,time)
    
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
    
    end subroutine displayRefinementStats_dc2dkx         

!==================================================================================================================================! 
!================================================================================================================= dc2dkx_initialize
!==================================================================================================================================! 
    subroutine dc2dkx_initialize   
!
! Initializes variables, checks input options and allocates some arrays
! 
    character(32)  :: cend   

! 
! Check input linear solver parameter:
!
    select case (trim(linearSolver))
        case ('intelmkl')
         ! a-okay
        case default
            write(*,*) ' Error in dc2dkx, unknown linear solver: ',trim(linearSolver)
            write(*,*) ' Stopping !'
            stop
    end select
 
! 
! Initialize a few things:
!        
    kx2          = kx**2    
    lmadeNewMesh = .false.
    
    ! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0) 
 
!
! Allocate private arrays used only in dc2dkx:
!
    if (lprintDebug_dc2dkx) write(*,*) myID,': allocating in dc2dkx...'
    
    allocate ( rhs(mesh%nnod,ntrodes_tx) )
    allocate ( lBoundary(mesh%nnod)  )
    
    allocate ( eRx(ntrodes_rx),eTx(ntrodes_tx) )
    allocate ( lCloud(mesh%nele) )
    allocate ( edges(3,mesh%nele) )
 
    if (lprintDebug_dc2dkx) write(*,*) myID,': leaving dc2dkx_initialize'
        
    end subroutine dc2dkx_initialize
        
!==================================================================================================================================! 
!================================================================================================================= dc2dkx_deallocate
!==================================================================================================================================! 
    subroutine dc2dkx_deallocate   
!
! Deallocates all remaining internal arrays at end of dc2dkx
!
 
    if (lprintDebug_dc2dkx) write(*,*) myID,': deallocating arrays...'

!     select case (trim(linearSolver))
!      
!     case ('superlu')
!         call superlu_zfree(superlu_p)   
!    case ('intelmkl')
        call intelmkl_zfree(intelmkl_p)  
!    end select  
    
    if ( allocated(val) )       deallocate(val)
    if ( allocated(col) )       deallocate(col)
    if ( allocated(irw) )       deallocate(irw)    
 
    if (allocated(errnrm)) deallocate ( errnrm ) 
    deallocate ( lBoundary ) 
    deallocate ( rhs )  
    deallocate ( eRx, eTx )
    deallocate( lCloud)
    deallocate ( edges )  !kwk is this used at all?
    if (allocated(node2tri)) deallocate( node2tri ) !kwk is this used at all?
    !kwk debug what about localError, dualError  ?
    
    call kdtree2_destroy(tree)  
    
       
    end subroutine dc2dkx_deallocate  
    
!==================================================================================================================================! 
!================================================================================================================ dc2dkx_printTimers
!==================================================================================================================================! 
    subroutine dc2dkx_printTimers 
    
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx dc2dkx_initialize:    ',t1-t0
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx get_finite_dipoles:   ',t2-t1
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx get_element_lists:    ',t3-t2
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_nodlab:           ',t4-t3
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_bcnod             ',t5-t4        
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_lhs:              ',t6-t5
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_rhs:              ',t7-t6
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx solve_primal:         ',t8-t7
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx fill_solmtx:          ',t9-t8
    write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx checkConvergence:     ',t10-t9
    
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_rhs_bump_F:       ',t11-t10
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_rhs_bump_B:       ',t12-t11
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_lhs_bump:         ',t13-t12
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx solve_varepsilon(1):  ',t14-t13
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_rhs_dual:         ',t15-t14
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx solve_dual:           ',t16-t15        
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx gen_rhs_bump_G:       ',t17-t16    
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx solve_varepsilon(2):  ',t18-t17 
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx estimateerror dealloc:',t20-t19   
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx comp_adj_derivs:      ',t21-t20   
!     write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx total:                ',t21-t0

    end subroutine dc2dkx_printTimers

!==================================================================================================================================! 
!============================================================================================================ localRefinement_dc2dkx
!==================================================================================================================================! 
    subroutine localRefinement_dc2dkx

    implicit none
    
    integer                 :: iTx, iRx, i,e, n(3), inele(1), newnodes, oldnodes
    real(8)                 :: area,r, rmin, rminRx,rminTx, areaScaled, a1,a2,a0,r1,r2,frac,ae1,ae2
    real(8)                 :: d1,d2,c    
    real(8), dimension(3)   :: ye,ze 
    
    logical                 :: lMakeNewMesh 
    real(kdkind)            :: qv(2) 
    type(kdtree2_result)    :: results(1)   
    type(kdtree2), pointer  :: treeTx, treeRx, treeNodes     
     
    real(8), dimension(:), allocatable  :: sigTemp
    real(8)                             :: sigx,sigy,sigz
    
    character(256)                       :: cadapt, filename,cend  
    real(8), dimension(:,:), allocatable :: yz
 
    real(8) :: tstart,tend
           
    if (lprintDebug_dc2dkx) write(*,*) myID,': entered localRefinement_dc2dkx...'            
    
    call cpu_time(tstart)
 
    oldnodes = mesh%nnod  
    
    ! Triangle command:
    write(cend,fmt='(F4.0)') minQangle
    retricommand = 'q'//trim(adjustl(cend))//'rpnajQ'//CHAR(0)     
    
!
! Create some kdtrees of the Tx's and Rx's for fast closest point searches:
! 
    if (ntrodes_rx > 1) then
        allocate(yz(2,ntrodes_rx))  
        yz(1,:) = electrode_rx(:,2)
        yz(2,:) = electrode_rx(:,3) 
        treeRx => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
        deallocate(yz)
        nullify(treeRx%the_data)  
    endif    
    if (ntrodes_tx > 1) then
        allocate(yz(2,ntrodes_tx))   
        yz(1,:) = electrode_tx(:,2) 
        yz(2,:) = electrode_tx(:,3)   
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
        
        if (lprintDebug_dc2dkx) write(*,*) 'local refinement loop, mesh%nnod: ',mesh%nnod
      
        !
        ! Refine elements containing receivers:
        !
        allocate (sigTemp(mesh%nele))
                    
        do e = 1,mesh%nele
            call getSigs(e,sigx,sigy,sigz)
            sigTemp(e) = abs(sigx+sigy+sigz)/3d0
        enddo
        
        allocate (node2tri(mesh%nnod))
        call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)                
 
        treeNodes => kdtree2_create( transpose(reshape( [mesh%y, mesh%z]  , [size(mesh%y),2]) ), sort=.true.,rearrange=.true.) 
        nullify(treeNodes%the_data) 
                   
        do iRx = 1,ntrodes_rx
  
            call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
            & sigTemp, 1, electrode_rx(iRx,2), electrode_rx(iRx,3), inele)
 
            call getArea(inele(1),area,0) 
            
            !write(*,'(a,i,1x,g,1x,g,1x,i)') ' iRx,area,nRxEleArea,inele:',iRx,area,nRxEleArea,inele
 
            if ( area > nRxEleAreaDC ) then  
            
                mesh%area(inele) = area/2.0
                lMakeNewMesh= .true.  
                
            endif
            
        enddo      
        
        deallocate (sigTemp,node2tri)
        call kdtree2_destroy(treeNodes)
        
        !
        ! Refine elements containing transmitters:
        !
   
        treeNodes => kdtree2_create( transpose(reshape( [mesh%y, mesh%z]  , [size(mesh%y),2]) ), sort=.true.,rearrange=.true.) 
        nullify(treeNodes%the_data) 

        allocate (sigTemp(mesh%nele))

        do e = 1,mesh%nele
            call getSigs(e,sigx,sigy,sigz)
            sigTemp(e) = abs(sigx+sigy+sigz)/3d0
        enddo

        allocate (node2tri(mesh%nnod))

        call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)

        do iTx = 1,ntrodes_tx

            call findElement(treeNodes, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
            & sigTemp, 1, electrode_tx(iTx,2), electrode_tx(iTx,3), inele)

            call getArea(inele(1),area,0) 
            
            !write(*,'(a,i,1x,g,1x,g,1x,i)') ' iTx,area,nTxEleArea,inele:',iTx,area,nTxEleArea,inele
            
            if ( area > nTxEleAreaDC) then
                if ( ( mesh%area(inele(1)) < 0 )  .or. (mesh%area(inele(1)) >  area/2d0) ) then
                    mesh%area(inele) = area/2d0
                    lMakeNewMesh= .true.  
                endif
            endif

        enddo

        deallocate(node2tri,sigTemp)           
        call kdtree2_destroy(treeNodes)
 
                 
        !
        ! Radial and elliptical refinement around Rx-Tx pairs
        !   

        if (maxRefinements == 0) then ! only do this when adaptive refinement is not being done                        
          
        do e = 1,mesh%nele
  
            call getArea(e,area,0)  
       
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
                
                r = sqrt( (electrode_rx(results(1)%idx,2) - ye(i))**2 + (electrode_rx(results(1)%idx,3) - ze(i))**2)
            
                if (r < rminRx) rminRx = r  
            
                !
                ! Find closest Tx:
                !
                 
                if (nTx > 1) then
                    call kdtree2_n_nearest(tp=treeTx,qv=qv,nn=1,results=results) 
                else
                    results(1)%idx = 1
                endif                
                r = sqrt( (electrode_tx(results(1)%idx,2) - ye(i))**2 + (electrode_tx(results(1)%idx,3) - ze(i))**2)
        
                if (r < rminTx) rminTx = r  
                
                             
                 
            enddo

            rmin = min(rminRx,rminTx)
  
            !
            ! Radial refinement around closest Rx or Tx to current element:
            !
                       
!             r1 = 40.  
!             r2 = 300.
            r1 = 5.  
            r2 = 500.            
            !a0 = .1**2.
             a0 = 1.**2.
!             a1 = 20.**2.
!             a2 = 200.**2.
!             a1 = 5.**2.
!             a2 = 20.**2.
            a1 = 20.**2.
            a2 = 80.**2.        
            areaScaled = 1d100
            
            if (rmin <= r1)  then
                frac = (rmin/r1)
                areaScaled = a0*(1.-frac) + a1*frac                        
            elseif (rmin <= r2)  then 
                frac = (rmin-r1)/(r2-r1)
                areaScaled = a1*(1.-frac) + a2*frac       
            endif                  
            if ( areaScaled < area) then
                if ( ( mesh%area(e) < 0 )  .or. (mesh%area(e) >  area/2d0) ) then
                    mesh%area(e) = area/2d0
                    lMakeNewMesh= .true.
                endif
            endif           

            !
            ! a1 size area refinement in ellipses connection all Rx-Tx pairs:
            !               
            do iRx = 1,ntrodes_rx
                do iTx = 1,ntrodes_tx

                    if (iDataMask(iRx,iTx)>0.) then
        
                        do i = 1,3
 
                            d1 = sqrt( (ye(i) - electrode_rx(iRx,2))**2 + (ze(i) - electrode_rx(iRx,3))**2 )
                            d2 = sqrt( (ye(i) - electrode_tx(iTx,2))**2 + (ze(i) - electrode_tx(iTx,3))**2 )
    
                            c =  sqrt( (electrode_rx(iRx,2) - electrode_tx(iTx,2))**2  + &
                                &      (electrode_rx(iRx,3) - electrode_tx(iTx,3))**2     )/2d0
                                
                            ae1 = c + r1        
                            ae2 = c + r2
                        
                            if ( d1 + d2  <= 2.*ae1 ) then
         
                                areaScaled = a1
  
                            elseif ( d1 + d2  <= 2.*ae2 ) then 
                            
                                areaScaled = a2
                            
                            else ! shrink outside ellipse
                                            
                                areaScaled = a2*( (d1 + d2 )/(2.*ae2) )**5
             
                            endif  

                            if ( areaScaled < area) then
                                if ( ( mesh%area(e) < 0 )  .or. (mesh%area(e) >  area/2d0) ) then
                                    mesh%area(e) = area/2d0
                                    lMakeNewMesh= .true.
                                endif
                            endif   

                        enddo
    
                    endif
                enddo
            enddo         


        enddo ! Loop over elements

        endif
             
        ! 
        ! Refine the mesh:
        !
        if (lMakeNewMesh ) then        
    
            call call_triangle(retricommand, mesh )           
            if (allocated(mesh%area))   deallocate(mesh%area)     

            ! Smooth the new mesh:
            call smooth_TriMesh(mesh,10,.false.)
          
        endif
 
        if (.not.lMakeNewMesh) then
            if (lprintDebug_dc2dkx) write(*,*) 'a priori refinement: ',mesh%nnod, mesh%nele
            exit    
        endif    
    enddo
    
    newnodes    = mesh%nnod  
    newnodes    = mesh%nnod  
    if (newnodes > oldnodes) meshnumber  = meshnumber + 1  

    if (lSaveMeshFiles)   then  
        write(cadapt,'(i6)') meshnumber  
        cadapt = adjustl(cadapt)    
        filename  = trim(fileroot)//'.'//trim(cadapt)
        call write_triangle(mesh,filename)
    endif   
    
    if (allocated(mesh%area))   deallocate(mesh%area) 
 
    if (ntrodes_rx > 1) call kdtree2_destroy(treeRx)  
    if (ntrodes_tx > 1) call kdtree2_destroy(treeTx)
    
    call cpu_time(tend)
    if (lprinttimers)  write(*,'(i4,1x,a32,1x,f8.4)') myID,' dc2dkx localRefinement_dc2dkx:    ',tend-tstart
 
    if (lDisplayRefinementStats) call displayRefinementStats_dc2dkx(iRefinementGrp,isubset,meshnumber,.false., &
                                             & oldnodes,newnodes,-2d0,-2d0, tend-tstart)
                                        
    end subroutine localRefinement_dc2dkx
  
!==================================================================================================================================! 
!=========================================================================================================================== getArea
!==================================================================================================================================! 
    subroutine getArea(e,area,iopt)
    
    integer, intent(in)   :: e
    real(8), intent(out)  :: area
    integer, intent(in)   :: iopt
        
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
!================================================================================================================= get_element_lists
!==================================================================================================================================! 
    subroutine get_element_lists
!
! Compute receiver and transmitter element lists
!    
    
!
! Local variables:
!
    integer                             :: e 
    real(8), dimension(:), allocatable  :: sigTemp
    real(8)                             :: sigx,sigy,sigz     
    real(8), dimension(:,:),allocatable :: yz
    
    !
    ! Create element edge list:
    !
    if (lprintDebug_dc2dkx) write(*,*) myID,': computeEdges...'
    
    call computeEdges( mesh%nele, mesh%emap, mesh%neighborlist, edges, nedges )      
            
            
    if (lprintDebug_dc2dkx) write(*,*) myID,': getting site element lists...'
    
    !
    ! First make the kdtree:
    !  
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': getting site element lists...kdtree2',size(mesh%y),size(mesh%z)
    
    allocate(yz(2,size(mesh%y)))   
    yz(1,:) = mesh%y
    yz(2,:) = mesh%z    
    tree => kdtree2_create( yz, sort=.true.,rearrange=.true.) 
    deallocate(yz)
    nullify(tree%the_data)
    
    ! Get node2tri array:
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': getting site element lists...getNode2Tri'
    
    allocate (node2tri(mesh%nnod))
    
    call getNode2Tri( mesh%nele, mesh%emap,mesh%nnod, node2tri)
 
    !
    ! Jump-and-walk findElement routine:
    !
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': getting Rx and Tx element lists...findElement'
    
    allocate (sigTemp(mesh%nele))
    
    ! We will use the average conductivity to define the element for the searches
    ! since if the search finds the point is incident on a node or edge, it chooses the most conductive element.
    do e = 1,mesh%nele
        call getSigs(e,sigx,sigy,sigz)
        sigTemp(e) = (sigx+sigy+sigz)/3d0
    enddo
    
    call findElement(tree, mesh%nnod, mesh%y, mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
    & sigTemp, ntrodes_rx, electrode_rx(:,2), electrode_rx(:,3), eRx)
    
    call findElement(tree,mesh%nnod,mesh%y,mesh%z, mesh%nele, mesh%emap, mesh%neighborlist, node2tri, &
    & sigTemp, ntrodes_Tx, electrode_tx(:,2), electrode_tx(:,3), eTx)
    
    deallocate(sigTemp)    
  
    end subroutine get_element_lists
    
!==================================================================================================================================! 
!===================================================================================================================== gen_lBoundary
!==================================================================================================================================! 
    subroutine gen_lBoundary
!
! Subroutine to label the mesh boundaries.  Assumes that the mesh
! has an exact rectangular outer boundary aligned along the coordinate 
! axes
! 
    real(8)  :: ymin,ymax,zmin,zmax 

    if (lprintDebug_dc2dkx) write(*,*) myID,': gen_lBoundary...' 
 
    ymin = minval(mesh%y)
    ymax = maxval(mesh%y)
    zmin = minval(mesh%z)
    zmax = maxval(mesh%z)
    
    lBoundary = .false.
    where (mesh%y == ymax) lBoundary = .true.   ! right
    where (mesh%y == ymin) lBoundary = .true.   ! left
    where (mesh%z == zmax) lBoundary = .true.   ! bottom 
    where (mesh%z == zmin) lBoundary = .true.   ! top

    end subroutine gen_lBoundary      
    
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
    real(8), intent(out)            :: sigx,sigy,sigz
    
    ! Local variables:
    integer     :: iparamnum

      
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
        
    case ('triaxial')
        sigx = 1d0/rhoParams(3*(iparamnum-1)+1)
        sigy = 1d0/rhoParams(3*(iparamnum-1)+2)
        sigz = 1d0/rhoParams(3*(iparamnum-1)+3)
    case ('tix')                                    ! For transversely isotropic:
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
        write(*,*) 'Error decoding anisotropy in getSigs in dc2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
  
    end subroutine getSigs
    
!==================================================================================================================================! 
!====================================================================================================================== getParamNums
!==================================================================================================================================!  
    subroutine getParamNums(e,iparamnum) 
!
! Gets parameter indices for element e
! 

    implicit none
 
    integer, intent(in)  :: e
    integer, intent(out) :: iparamnum(3)
    
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
    case ('tix','tiy','tiz')                               ! For transversely isotropic:
        iparamnum(1) = iFreeParam(2*(iregion-1) + 1)
        iparamnum(2) = iFreeParam(2*(iregion-1) + 2)    
        
    case default
        write(*,*) 'Error decoding anisotropy in getSigs in dc2dkx.f90'
        write(*,*) 'Unknown or unsupported anisotropy code: ',trim(cAnisotropy)
        stop

    end select
     
     
    end subroutine getParamNums
        
!==================================================================================================================================! 
!====================================================================================================================== getStiff_lin 
!==================================================================================================================================! 
    subroutine getStiff_lin(e,mtx)
!
! Generates the 3 x 3 nodal linear FE stiffness matrix for the 2.5D triaxially anisotropic DC resistivity problem.
! 
    integer, intent(in)                 :: e
    real(8),dimension(3,3), intent(out) :: mtx   

! Local variables:  
    integer, dimension(3)   :: n
    real(8), dimension(3)   :: ye,ze
    real(8)                 :: sigx,sigy,sigz
    real(8), dimension(3)   :: a,b,c
    real(8)                 :: ar
    real(8)                 :: aey,aez,bex 
  
!
! Get triaxial conductivity for element e
!
    call getSigs(e,sigx,sigy,sigz) 
  
!
! Get linear basis coefficients:
!
    n    = mesh%emap(1:3,e)
    ye   = mesh%y(n)
    ze   = mesh%z(n)   
    call get_abc_coeffs(ye,ze,a,b,c,ar)
 
!
! 2.5D DC coefficients for triaxial conductivity:
!
    aey =  sigy  / (4d0*ar)
    aez =  sigz /  (4d0*ar)
    bex =  kx2*sigx * ar / 6d0
    
!
! Elements of the stiffness matrix (upper triangular part only since its symmetric):
!
    ! upper triangle with diag:
    mtx(1,1) = aey*b(1)*b(1) + aez*c(1)*c(1) + bex      ! e1 e1
    mtx(1,2) = aey*b(2)*b(1) + aez*c(2)*c(1) + bex/2d0  ! e1 e2
    mtx(1,3) = aey*b(3)*b(1) + aez*c(3)*c(1) + bex/2d0  ! e1 e3    
    mtx(2,2) = aey*b(2)*b(2) + aez*c(2)*c(2) + bex      ! e2 e2
    mtx(2,3) = aey*b(3)*b(2) + aez*c(3)*c(2) + bex/2d0  ! e2 e3 
    mtx(3,3) = aey*b(3)*b(3) + aez*c(3)*c(3) + bex      ! e3 e3
    
    !fill in lower triangle:
    mtx(2,1) = mtx(1,2)
    mtx(3,1) = mtx(1,3)
    mtx(3,2) = mtx(2,3)
    
    end subroutine getStiff_lin
 
!==================================================================================================================================! 
!=========================================================================================================================== gen_rhs 
!==================================================================================================================================! 
    subroutine gen_rhs 
! 
! Generates the 2.5D FE rhs array for DC resistivity point sources
!
    implicit none 
    
    integer                 :: i,n(3),e, iTx 
    real(8), dimension(3)   :: l,dldy,dldz, ye,ze   
 

    if (lprintDebug_dc2dkx) write(*,*) myID,': gen_rhs_dipole ...'  

!
!  Initialize the RHS 
!
    rhs = 0d0
   
!
! Loop over transmitters:
!
    do iTx = 1,ntrodes_tx
    
        !
        ! Set up the basis functions:
        !
            e   = eTx(iTx)
            n   = mesh%emap(1:3,e)
            ye  = mesh%y(n)
            ze  = mesh%z(n)

        !
        ! Evaluate the basis functions at the source location:
        !
            call get_linear_basis(ye,ze,electrode_tx(iTx,2),electrode_tx(iTx,3),l,dldy,dldz)    
 
        !
        ! Insert into RHS:
        ! 
            rhs(n, iTx) = l  
  
    enddo

!
! Insert the dirichlet boundary condition into the entries 
! of the RHS whose nodes lie on the border:  
!
! In the total field formulation this isn't necessary unless somebody 
! mistakenly places a Tx on an exterior node boundary, but then the solution
! is going to be incorrect anyway...
!
    do i=1,mesh%nnod
        if( lBoundary(i) ) rhs(i,1:ntrodes_tx) = 0d0
    enddo    
 
     
        
    end subroutine gen_rhs

!==================================================================================================================================! 
!=========================================================================================================================== gen_lhs 
!==================================================================================================================================! 
     subroutine gen_lhs 
! 
! Generates the lhs matrix for the DC resistivity with triaxial anisotropy 
!
 
    integer                             :: e,i,j,v1,v2,ind,ierr,i0,i1  
    real(8), dimension(3,3)             :: mtx
    integer, dimension(:), allocatable  :: row  
     
! 
! Number of non-zeros in the LHS matrix:
! 
    nnz = mesh%nnod + 2*nedges  !  self node components + node-node components   
 
!  
! Allocate CSR arrays val,col,irw:
!    
    allocate ( val(nnz), col(nnz), row(nnz), irw(mesh%nnod + 1),stat=ierr) 
    if (ierr .ne. 0) then
        write(*,*) ' dc2kdkx.f subroutine gen_lhs: Error allocating, too many non-zeros elements to fit into memory: ', nnz 
        stop
    end if  
    val = 0d0 


!
! Loop over each element, compute element stiffness matrix and insert into the global CSR matrix:
!   
    do e = 1,mesh%nele
    

        !
        ! Get element stiffness matrix:
        !
        call getStiff_lin(e, mtx) ! code below assumes full matrix returned...
        
        !
        ! Insert entries into the (unsorted) CSR matrix:
        !
        
        !
        ! self node:                                    [ nnod ]
        !
        do i = 1,3
            ind = mesh%emap(i,e) 
            val(ind) = val(ind)  + mtx(i,i)  
            row(ind) = ind
            col(ind) = ind
        enddo
        
        !
        ! node-node:                                    [ 2*nedges ]
        !
        do i = 1,3
            v1 = mesh%emap(eps(i+1),e) 
            v2 = mesh%emap(eps(i+2),e)
            
            ind = mesh%nnod + 2*edges(i,e)-1               ! note that v1,v2 order doesn't matter since matrix is symmetric    
            val(ind) = val(ind)  + mtx(eps(i+1),eps(i+2))  
            row(ind) = v1
            col(ind) = v2
            
            ind = mesh%nnod + 2*edges(i,e)
            val(ind) = val(ind)  + mtx(eps(i+2),eps(i+1))  
            row(ind) = v2
            col(ind) = v1            
            
        enddo        
        
 
       
    enddo  ! e = 1,nele
     

 
    !
    ! Apply boundary conditions (easy homogeneous Dirichlet)
    !

    !
    ! First pass implements the diagonal entries:
    !
    do e = 1,nnz
        i = row(e)
        j = col(e)
        if ( ( i == j ) .and. (lBoundary(i)) ) then
            val(e)   = 1d0
!            bval     = 0d0 ! kwk debug: insert call to a subroutine boundary(x)  here
!            rhs(i,1) = bval
        endif
        
    enddo
 
    
    !
    ! Second pass implements the off-diagonal entries:
    !
    do e = 1,nnz
        i = row(e)
        j = col(e)
        
        if ( ((lBoundary(i)) .or. (lBoundary(j)) ) .and. ( i /= j ) )  then 
            
!            if (lBoundary(j)) then 
!                rhs(i,1) = rhs(i,1) - val(e)*rhs(j,1) 
!            endif
            val(e) = 0d0
        endif
        
    enddo   

!
! Now sort on the row and then column indices:
!

   ! 
   ! First sort based on the row:
   !
    call quick_sort(row,col,dlist1=val) 

   !
   ! Then sort each row's columns, and set the CSR row index:
   !     
    i0 = 1
    irw(1) = 1 
    do i = 1,mesh%nnod 
    
        ! Find the end of this row's entires:
        do j = i0,nnz
       
            if (row(j) /= i) then
                exit
            else
                i1 = j
            endif
              
        enddo
        
        if (i1>i0) call quick_sort(col(i0:i1),dlist1=val(i0:i1)) 
        i0 = i1+1
        irw(i+1) = i0
    enddo

!     write(*,*) 'post sort:'    
!     do i = 1,nnz
!         write(*,'(2(i,1x),g)') row(i),col(i),val(i)
!     enddo
 

!
! Deallocate local arrays:
!
    deallocate(row)
    
         
    end subroutine gen_lhs

!==================================================================================================================================! 
!======================================================================================================================= fill_solmtx
!==================================================================================================================================! 
    subroutine fill_solmtx
!
! Extracts potentials at receiver electrodes
!
    implicit none
    
  ! local variables  
    real(8)               :: yp,zp
    integer               :: iRx,e,iTx,iComp
    integer, dimension(3) :: n 
    real(8), dimension(3) :: ye,ze,fv,l,dldy,dldz 
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': fill_solmtx ... nComponents:'  ,nComponents
        
!
! Loop over the sites and fill in the kx domain solutions:
!
    do iComp = 1,nComponents
        
        ! Get rx and tx electrode indices for this component:
        iRx = Components(iComp,1) 
        iTx = Components(iComp,2)
                
        ! get element arrays:
        e    = eRx(iRx)
        n    = mesh%emap(1:3,e)
        ye   = mesh%y(n)
        ze   = mesh%z(n)     
      
        ! Get linear potential coefficients for this element and transmitter:
        fv = rhs(n,iTx) 

        ! Interpolate potential to receiver electrode location: 
        yp = electrode_rx(iRx,2)
        zp = electrode_rx(iRx,3)   
        
        call get_linear_basis(ye,ze,yp,zp,l,dldy,dldz)
        
        !write(*,'(i,1x,i,1x,i,1x,g,1x,5(g,1x))') iRx,iTx,iComp,kx,electrode_rx(iRx,2),electrode_tx(iTx,2),sum(l*fv) 
        
        ! Save into Fields array:
        Fields(iComp) = sum(l*fv)   
                   
    enddo  ! loop over nComponents    
    
    end subroutine fill_solmtx

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

    if (lprintDebug_dc2dkx) write(*,*) myID,': refineMesh...'   
    
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
!========================================================================================================================== write_Ab 
!==================================================================================================================================!    
    subroutine write_Ab()
!
! Writes A and b from Ax=b out to file. 
! Matrix A is saved as:  irow,jcol, real  
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
       
    
    filename  = trim(fileroot)//'.'//trim(cadapt)//'.lhs'
    open(16,file=trim(filename),status='REPLACE') 
    write(16,*) mesh%nnod+1,nnz
    do i = 1,mesh%nnod+1
        write(16,*) irw(i)
    enddo                 
    do j = 1,nnz
        write(16,*) col(j), dble(val(j)) !, aimag(val(j)) 
    enddo
    close(16)       
              
    
    filename  = trim(fileroot)//'.'//trim(cadapt)//'.rhs'
    open(16,file=trim(filename),status='REPLACE') 
    write(16,*) mesh%nnod,  ntrodes_tx
    do j = 1,ntrodes_tx
        do i = 1,mesh%nnod
            write(16,*)   dble(rhs(i,j)) !, aimag(rhs(i,j)) 
        enddo
    enddo
    close(16)                      

               
    end subroutine write_Ab
!==================================================================================================================================! 
!=========================================================================================================================== write_x
!==================================================================================================================================!    
    subroutine write_x()
!
! Writes x from Ax=b out to file. 
! Matrix A is saved as:  irow,jcol, real 
!
    implicit none
    
    integer         :: i,j
    
    character(256)  :: filename, cadapt 
        
    write(cadapt,'(i6)') meshnumber  
    cadapt = adjustl(cadapt)       
        
    filename  = trim(fileroot)//'.'//trim(cadapt)//'.x'
    open(16,file=trim(filename),status='REPLACE') 
    write(16,*) mesh%nnod,ntrodes_tx
    do j = 1,ntrodes_tx
        do i = 1,mesh%nnod
            write(16,*)   dble(rhs(i,j)) !, aimag(rhs(i,j)) 
        enddo
    enddo
    close(16)                      
        
 
    end subroutine write_x
                    
!==================================================================================================================================! 
!====================================================================================================================== solve_primal 
!==================================================================================================================================! 
    subroutine solve_primal 
!
! Solves the linear system AX=B using one of various solvers:
!

    if (lSaveLinearSystem) call write_Ab
             
!     select case (trim(linearSolver))
!     case ('superlu')
!         call solve_primal_superlu
!     case('intelmkl')
        call solve_primal_intelmkl  
!    end select
     
    ! Matrix has been factored, so we can deallocate the CSR arrays:
    deallocate(val,col,irw)       
        
    if (lSaveLinearSystem) call write_x
            
    end subroutine solve_primal   
                    
!==================================================================================================================================! 
!============================================================================================================= solve_primal_intelmkl 
!==================================================================================================================================!    
    subroutine solve_primal_intelmkl()

    implicit none

    integer :: j
    real(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': solve_primal_intelmkl ...'
          
!
! Factor the matrix first:
!
    call intelmkl_dfactor(intelmkl_p,val,col,irw)
    
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs,1)) )
    
    do j = 1,size(rhs,2)  
 
        x = rhs(:,j)     
        call intelmkl_dsolve(intelmkl_p,x )
        rhs(:,j) = x
 
    enddo
 
    deallocate(x)
 
    end subroutine solve_primal_intelmkl
 
!==================================================================================================================================! 
!======================================================================================================================== solve_dual
!==================================================================================================================================!     
    subroutine solve_dual 
!
!  Solve dual problem
!                                                                      
    select case (trim(linearSolver))
    
    case('intelmkl')
  
        call solve_dual_intelmkl
          
    end select   
    
    end subroutine solve_dual
 
!==================================================================================================================================! 
!=============================================================================================================== solve_dual_intelmkl 
!==================================================================================================================================!    
    subroutine solve_dual_intelmkl()

    implicit none

    integer :: j
    real(8) , dimension(:), allocatable :: x
    
    if (lprintDebug_dc2dkx) write(*,*) myID,': solve_dual_intelmkl ...'
 
!
! Now solve  AX=B for all column vectors in B:
!
    allocate ( x(size(rhs_wh,1)) )
    
    do j = 1,size(rhs_wh,2)  
 
        x = rhs_wh(:,j)     
        
        call intelmkl_dsolve(intelmkl_p,x )
        
        rhs_wh(:,j) = x
 
    enddo
 
    deallocate(x)
 
    end subroutine solve_dual_intelmkl


 
    end module dc2dkx_mod
    

