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
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------  module MT1D
!----------------------------------------------------------------------------------------------------------------------------------!

! Fortran module for 1D EM plane-wave computations
!
 

    module MT1D_module
    
    use EM_constants, only : ic, mu0, pi
    
    implicit none
    
    private
    
    type, public :: MT1D 
        
        ! Public fields: 
        integer, private                             :: nlayer
        real(8), dimension(:), allocatable, private  :: zlay, sig  ! top depth and conductivity   
        real(8), private                             :: omega      ! angular frequency
        
        ! Private fields:
        complex(8), dimension(:), allocatable, private :: a,b  ! layer coefficients
       
      contains ! type bound procedures:
      
        procedure, public :: getEH1D
        procedure, public :: getSig  ! getSig(z)
        procedure, public :: setSig  ! setSig(ilayer,sig)
        procedure, public :: getLayer
        procedure, public :: getTotalConductance  

       
    end type MT1D

    ! public procedures:
    public :: new_MT1D, delete, print 

    ! generic interface for delete_M1D  
    interface delete                ! e.g., if Layers1D is of type MT1D, this allows you to use
        module procedure delete_MT1D       !   delete(Layers1D)  to deallocate and close the object
    end interface
    interface print
        module procedure print_MT1D
    end interface    

    contains 

!----------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------ print_MT1D
!----------------------------------------------------------------------------------------------------------------------------------!
     subroutine print_MT1D(this)  
           
    
    class (MT1D), intent(in) :: this
     
    integer :: i    
    
    write(*,'(3(a12,1x))') 'Layer','Depth','S/m'
    
    do i = 1,this%nlayer   
        write(*,'(i12,1x,f12.1,1x,g12.3,1x)') i, this%zlay(i), this%sig(i)
    enddo
    write(*,'(a,1x,g12.3)') ' Frequency(rad/s): ', this%omega

    end subroutine print_MT1D 
        
!----------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------- getTotalConductance
!----------------------------------------------------------------------------------------------------------------------------------!
    real(8) function getTotalConductance(this) result(conductance)  
           
    class (MT1D)  ,intent(inout)        :: this   
                    
    integer :: i    
    
    !write(*,'(5(a12,1x))') 'Layer','Depth','S/m','thick (m)','conductance'
    
    conductance = 0d0
    do i = 1,this%nlayer-1     
        conductance = conductance + this%sig(i)*(this%zlay(i+1)-this%zlay(i))
        !write(*,'(i12,1x,f12.1,1x,3(g12.3,1x))') i, this%zlay(i), this%sig(i),(this%zlay(i+1)-this%zlay(i)),this%sig(i)*(this%zlay(i+1)-this%zlay(i))
    enddo
    !write(*,'(a,1x,1(g12.3,1x))') 'total conductance: ',conductance
    
    end function getTotalConductance 
!----------------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------- new_MT1D
!----------------------------------------------------------------------------------------------------------------------------------!
    type(MT1D) function new_MT1D(zlay,sig,omega)  result(this)
           
    real(8), dimension(:), intent(in)  :: zlay,sig
    real(8), intent(in)                :: omega 
 
    integer :: nlayer
    
    !
    ! Set the model parameters:
    !
    nlayer  = size(zlay,1)
    
    this%nlayer = nlayer
    allocate (this%zlay(nlayer),this%sig(nlayer))
    this%zlay  = zlay(1:nlayer)
    this%sig   = sig(1:nlayer)

    !
    ! Set the angular frequency:
    !
    this%omega   = omega
 
    !
    ! Precompute the layer coefficients for this frequency:
    !
    call getMT1Dcoeffs(this)    
        
    end function new_MT1D

!----------------------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------------- getSig
!----------------------------------------------------------------------------------------------------------------------------------!
    real(8) function getSig(this,z) result(sig)
!
! Returns 1D conductivity at depth z:
! Note that layer depths are for the top of the layer and note the equality convention used:
!  z_top(ilay) <= z < z_top(ilay+1) 
!    
    
    class (MT1D), intent(in) :: this
    real(8), intent(in)      :: z  
    
    integer :: ilay    
    
    ilay = this%getLayer(z)
    sig  = this%sig(ilay)
    
    end function getSig

!----------------------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------------- setSig
!----------------------------------------------------------------------------------------------------------------------------------!
    subroutine setSig(this,ilay,sig)  
!
! Sets 1D conductivity of layer ilay
!
    class (MT1D), intent(inout) :: this
    integer, intent(in)         :: ilay    
    real(8), intent(in)         :: sig 
    

    ! Update conductivity at depth z
    this%sig(ilay) = sig
 
    ! Re-compute the layer coefficients for the updated model:
    call getMT1Dcoeffs(this)   
 
    end subroutine setSig    
        
!----------------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------- getLayer
!----------------------------------------------------------------------------------------------------------------------------------!
    integer function getLayer(this,z) result(ilay)
!
! Returns index of the layer containing z.
! Note that layer depths are for the top of the layer and note the equality convention used:
!  z_top(ilay) <= z < z_top(ilay+1) 
!    
    
    class (MT1D), intent(in) :: this
    real(8), intent(in)     :: z  
    
    integer :: i    
    
    ilay = 1
    do i = 1,this%nlayer
        if (z >= this%zlay(i))  ilay = i
    enddo

    end function getLayer
    
!----------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------- getMT1Dcoeffs
!----------------------------------------------------------------------------------------------------------------------------------!
    subroutine getMT1Dcoeffs(this)
!
! Computes the layer propagation coefficients for plane-wave magnetotelluric electric and magnetic fields as
! a function of depth in a 1D model.  Uses a stable recursion similar to the ones I used in Dipole1D, but modified for 
! the simpler case of vertically incident plane wave. 
!
! The routine finds E field coefficients a,b for each layer i, so that fields 
! at depth z have the form:
!
!   E_i(z) = a_i * exp( +k_i * (z - z_i+1) ) + b_i * exp( -k_i * (z - z_i) )
!
!   where k_i = sqrt( sqrt(-1)*omega*mu0*sig_i ) and an exp(+iwt) time dependence has been used.
!
! Field is normalized so that E at the top of the model has unit magnitude.
!

    
    ! I/O:
    type (MT1D), intent(inout) :: this

    ! Local variables:
    integer    :: i
    complex(8) :: gmogp, rjexp, Atop 
    complex(8), dimension(:), allocatable :: kk, Rp, expmgh
    
!
! Allocate storage:
!
    allocate( kk(0:this%nlayer), Rp(0:this%nlayer), expmgh(0:this%nlayer) )

!
! Initialize a few variables:
!        
    kk(1:this%nlayer)  = sqrt(ic*this%omega*mu0*this%sig(1:this%nlayer))
    kk(0) = sqrt(ic*this%omega*mu0*this%sig(1))    
    expmgh = 0d0    
    do i = 1,this%nlayer-1
        expmgh(i) = exp( -kk(i)*( this%zlay(i+1) - this%zlay(i) ) )
    enddo
   
    Rp = 0d0

!
! Compute recursion:
!    
    do i = this%nlayer-1,1,-1
       gmogp   = (ic*this%omega*mu0*(this%sig(i)-this%sig(i+1) )) / (kk(i)+kk(i+1))**2
       rjexp   = Rp(i+1) * expmgh(i+1)
       Rp(i)   =  (gmogp + rjexp)   / ( 1.d0 + gmogp*rjexp)   * expmgh(i)
    enddo
    
!
! Back propagate after setting top layer to boundary condition (a=1,b=0):
!   
    if ( (.not.allocated(this%a) ) .or. ( size(this%a) /= this%nlayer+1 ) ) then
        allocate( this%a(0:this%nlayer), this%b(0:this%nlayer))
    endif
    this%a    = 0
    this%b    = 0 
    this%a(0) = 1 
 
    ! Propagate the solution across each interface via field continuity conditions:
    do i = 1,this%nlayer
    
         ! E coefficients
        Atop = (this%a(i-1) + this%b(i-1)*expmgh(i-1) )   ! Potential at layer top
        this%b(i)  = Atop / (1.d0 + Rp(i)*expmgh(i) )
        this%a(i)  = this%b(i)*Rp(i)
    
    enddo
!
! Deallocate temporary storage:
!
    deallocate ( kk, Rp, expmgh )
    
    end subroutine getMT1Dcoeffs
    
!----------------------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------- getMT1Dfields
!----------------------------------------------------------------------------------------------------------------------------------!
    complex(8) function getEH1D(this,z,component) result(field)  
!
!
! Outputs the E and H fields at absolute depth z using:
!
!   E_i(z) = a_i * exp( +k_i * (z - z_i+1) ) + b_i * exp( -k_i * (z - z_i) )
!   H_i(z) = d(E_i(z))/dz * 1/(-sqrt(-1)*omega*mu) = k_i*(a_i*exp(+k_i*(z - z_i+1)) - b_i * exp(-k_i*(z - z_i)) )/(-ic*omega*mu)
!
!   where k_i = sqrt( sqrt(-1)*omega*mu0*sig_i ) and an exp(+iwt) time dependence has been used.
!
!   
! Fields are normalized so that E at the top of the this (i.e., z = z_1)
! has unit magnitude.
!

    class (MT1D)  ,intent(inout)        :: this   
    real(8)       ,intent(in)           :: z          ! depth of field evaluation
    character(*)  ,intent(in)           :: component  ! 'e' or 'h'                       

    integer    :: ilay
    complex(8) :: kk, expp, expm
    
!
! Find the layer at depth z:
!
    ilay = this%getLayer(z)

!
! Compute the fields at depth z:
!

    kk = sqrt(ic*this%omega*mu0*this%sig(ilay))
    
    expp = 0
    expm = 0
    
    expm = exp( -kk * (z - this%zlay(ilay) ) )
    if (ilay /= this%nlayer) expp = exp( kk * (z - this%zlay(ilay+1) ))
    
    select case (trim(component))
    
    case ('e','E')
        field = this%a(ilay) * expp  + this%b(ilay) * expm  ! E   
    case('h','H')
        field = kk/(-ic*this%omega*mu0) * ( this%a(ilay) * expp  - this%b(ilay) * expm )  ! H
    case default 
        write(*,*) 'Error in mt1d.f90 function getEH1D, requested component unrecognized: ',trim(component)
        stop
    end select
   
    ! tweak for exp(-iwt) dependence in mare2dem
    field = conjg(field)
    
    ! reality check: note that if there is only one layer in the model then b = 0 and then
    ! Ex/Hz = Z = omega*mu/k, as required by the halfspace MT theory.

    ! Zero out low values to avoid numerical issues in other codes using this module:
    if (abs(field) < 1d-80) field = 0d0
 
    
    end function getEH1D
    
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------- delete_MT1D
!----------------------------------------------------------------------------------------------------------------------------------!
    subroutine delete_MT1D(this)
 !
 ! Deallocates fields of MT1D derived type 
 !
    type (MT1D), intent(inout)  :: this   
    
    deallocate(this%a,this%b,this%sig,this%zlay)
    
    end subroutine delete_MT1D
    
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!

    end module MT1D_module