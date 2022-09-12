! intelmkl_solver (double precision complex and real options)
! A Fortran module for calling the sparse LU solver included in the Intel MKL (pardiso)
! Separate subroutines are available for factoring and solving a system
! so that the factors can be reused for multiple right-hand-sides.
!
!  Kerry Key 
!  Scripps Institution of Oceanography
!
 
include '/home/intel/oneapi/mkl/latest/include/mkl_dss.f90'
    
module intelmkl_solver

    use mkl_dss
 
    implicit none
!
! Derived type to store intelmkl variables for a given matrix factorization
!
 
 
    public :: intelmkl_zfactor  ! Step 1. Factor a complex sparse matrix: A = LU.
    public :: intelmkl_zsolve   ! Step 2. Solve the linear system: A x = LU x = b.
    public :: intelmkl_zfree    ! Step 3. Deallocate the memory for the factors LU.
 
    public :: intelmkl_dfactor  ! Step 1. Factor a complex sparse matrix: A = LU.
    public :: intelmkl_dsolve   ! Step 2. Solve the linear system: A x = LU x = b.
    public :: intelmkl_dfree    ! Step 3. Deallocate the memory for the factors LU.
        
    contains

!------------------------------------------------------------------------- 
!----------------------------------------------------------intelmkl_zfactor 
!------------------------------------------------------------------------- 
    subroutine intelmkl_zfactor(this,Az,Ai,Ap)
!
! Factors a complex matrix A, with pointers to solution being stored in 
! intelmkl derived type variable "this".
!

!
! Arguments:
!
    TYPE(MKL_DSS_HANDLE), intent(inout)   :: this  
     
    complex(8), dimension(:), intent(in)  :: Az    
    integer, dimension(:), intent(inout)  :: Ai, Ap     ! modified for c indexing here, corrected before returing
    
    integer ::   opt
    
    
!
! Local variables:
!
    integer    :: n, nz, info, perm(1)
    
!
! Get matrix size
!
    n = ubound(Ap,1) - 1
    nz = Ap(n+1) - 1

! 
! Initialize the solver:
!
    opt  = MKL_DSS_DEFAULTS  
    info = dss_create( this, opt )
    
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_CREATE, info = ',info
        stop
    endif
 
!
! Define the non-zero structure of the matrix:
!

    info = dss_define_structure( this, MKL_DSS_NON_SYMMETRIC_COMPLEX, Ap, n, n, Ai, nz )
    !info = dss_define_structure( this, MKL_DSS_SYMMETRIC_STRUCTURE_COMPLEX, Ap, n, n, Ai, nz )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_DEFINE_STRUCTURE, info = ',info
        stop
    endif
 
!
! Reorder the matrix:
!
    perm(1) = 0
    opt = MKL_DSS_AUTO_ORDER
    info = dss_reorder( this, opt, perm )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_REORDER, info = ',info
        stop
    endif 

      
! 
! Factor the matrix:
!
    opt = MKL_DSS_POSITIVE_DEFINITE   !MKL_DSS_INDEFINITE   
    info = dss_factor_complex( this, opt, Az )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_FACTOR_REAL, info = ',info
        stop
    endif 
 
    end subroutine intelmkl_zfactor

!------------------------------------------------------------------------- 
!-----------------------------------------------------------intelmkl_zsolve
!------------------------------------------------------------------------- 
    subroutine intelmkl_zsolve (this,b)
    
    implicit none
    
!
! Arguments:
!
    TYPE(MKL_DSS_HANDLE), intent(inout)      :: this  
    complex(8), dimension(:), intent(inout)  :: b   
    
!
! Local variables:
!
    integer :: n, nrhs, info, opt

        
    complex(8), dimension(:), allocatable :: sol
 
!
! Allocate temporary solution array 
!
    n = ubound(b,1)    
    nrhs =  1 ! !ubound(b,2)
    
    allocate(sol(n)) 
 
!
! Solve the linear system Ax=b:
!
    opt = MKL_DSS_DEFAULTS + MKL_DSS_REFINEMENT_OFF 
    info = dss_solve_complex(this, opt, b, nrhs, sol )
    
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_SOLVE_COMPLEX, info = ',info
        stop
    endif 
 
    b = sol
    
    deallocate(sol)
  
    end subroutine intelmkl_zsolve    
    
!------------------------------------------------------------------------- 
!------------------------------------------------------------intelmkl_zfree
!------------------------------------------------------------------------- 
    subroutine intelmkl_zfree(this)
    
    implicit none
        
    TYPE(MKL_DSS_HANDLE), intent(inout)    :: this  
    
    integer :: info
    
!
! Free numeric factorization memory:
!
    info = dss_delete( this, MKL_DSS_DEFAULTS )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_DELETE, info = ',info
        stop
    endif 
        
    end subroutine intelmkl_zfree
 
!------------------------------------------------------------------------- 
!----------------------------------------------------------intelmkl_dfactor 
!------------------------------------------------------------------------- 
    subroutine intelmkl_dfactor(this,Az,Ai,Ap)
!
! Factors a complex matrix A, with pointers to solution being stored in 
! intelmkl derived type variable "this".
!

!
! Arguments:
!
    TYPE(MKL_DSS_HANDLE), intent(inout)   :: this  
     
    real(8), dimension(:), intent(in)     :: Az    
    integer, dimension(:), intent(inout)  :: Ai, Ap     ! modified for c indexing here, corrected before returing
    
    integer ::   opt
    
    
!
! Local variables:
!
    integer    :: n, nz, info, perm(1)
    
!
! Get matrix size
!
    n = ubound(Ap,1) - 1
    nz = Ap(n+1) - 1

! 
! Initialize the solver:
!
    opt  = MKL_DSS_DEFAULTS  
    info = dss_create( this, opt )
    
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_CREATE, info = ',info
        stop
    endif
 
!
! Define the non-zero structure of the matrix:
!

    info = dss_define_structure( this, MKL_DSS_NON_SYMMETRIC, Ap, n, n, Ai, nz )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_DEFINE_STRUCTURE, info = ',info
        stop
    endif
 
!
! Reorder the matrix:
!
    perm(1) = 0
    opt = MKL_DSS_AUTO_ORDER
    info = dss_reorder( this, opt, perm )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_REORDER, info = ',info
        stop
    endif 

      
! 
! Factor the matrix:
!
    opt = MKL_DSS_POSITIVE_DEFINITE   !MKL_DSS_INDEFINITE   
    info = dss_factor_real( this, opt, Az )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_FACTOR_REAL, info = ',info
        stop
    endif 
 
    end subroutine intelmkl_dfactor

!------------------------------------------------------------------------- 
!-----------------------------------------------------------intelmkl_dsolve
!------------------------------------------------------------------------- 
    subroutine intelmkl_dsolve (this,b)
    
    implicit none
    
!
! Arguments:
!
    TYPE(MKL_DSS_HANDLE), intent(inout)      :: this  
    real(8), dimension(:), intent(inout)  :: b   
    
!
! Local variables:
!
    integer :: n, nrhs, info, opt

        
    real(8), dimension(:), allocatable :: sol
 
!
! Allocate temporary solution array 
!
    n = ubound(b,1)    
    nrhs =  1 ! !ubound(b,2)
    
    allocate(sol(n)) 
 
!
! Solve the linear system Ax=b:
!
    opt = MKL_DSS_DEFAULTS + MKL_DSS_REFINEMENT_OFF 
    info = dss_solve_real(this, opt, b, nrhs, sol )
    
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_SOLVE_COMPLEX, info = ',info
        stop
    endif 
 
    b = sol
    
    deallocate(sol)
  
    end subroutine intelmkl_dsolve    
    
!------------------------------------------------------------------------- 
!------------------------------------------------------------intelmkl_dfree
!------------------------------------------------------------------------- 
    subroutine intelmkl_dfree(this)
    
    implicit none
        
    TYPE(MKL_DSS_HANDLE), intent(inout)    :: this  
    
    integer :: info
    
!
! Free numeric factorization memory:
!
    info = dss_delete( this, MKL_DSS_DEFAULTS )
    if (info /= MKL_DSS_SUCCESS) then
        write(*,*) 'Error calling DSS_DELETE, info = ',info
        stop
    endif 
        
    end subroutine intelmkl_dfree

 

 

end module intelmkl_solver