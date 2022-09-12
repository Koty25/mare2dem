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
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
    module mare2dem_scalapack
!
! Module for scalapack interface used for Gauss-Newton inversion steps including
! dense matrix multiplication (A^T*A) and dense matrix inversion via Cholesky factorization.
!    
! Kerry Key
! Scripps Institution of Oceanography
!
    
    implicit none
    
    ! Variables shared by common blacs context:
    integer             :: myid, ierr, nprocs, ictxt, nprow, npcol, myrow, mycol
    integer,parameter   :: descriptor_len=9


    contains

    
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
    subroutine setupBlacsGrid
      
    call blacs_pinfo(myid, nprocs )                 ! Initialize the blacs.  Note: processors are counted starting at 0.
    call gridsetup(nprocs,nprow,npcol)              ! get dimensions of 2D processor grid 
    call blacs_get(-1,0,ictxt)                      !  Initialize a single blacs context.
    call blacs_gridinit(ictxt,'r',nprow,npcol)      !  Initialize a single blacs context.
    call blacs_gridinfo(ictxt, nprow,npcol,myrow,mycol)  ! determine row and column of this processor on the processor grid


!    write(*,*) 'myid, nprocs: ' ,myid, nprocs 
!    write(*,*) 'myid,nprow,npcol: ',myid,nprow,npcol
!    write(*,*) 'myid, ictxt: ',myid, ictxt
!    write(*,*) 'myid,myrow,mycol: ',myid,myrow,mycol
    
    end subroutine setupBlacsGrid
    
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
    subroutine allocateDistributedMatrix(m,n,desca,a_loc)
     
    integer, intent(in)                              :: m,n
    integer, dimension(descriptor_len), intent(in)   :: desca 
    real(8), dimension(:,:), allocatable,intent(out) :: a_loc
    
    integer             :: nb,l_m,l_n,nb_c
    
    integer :: istat, info
    
    integer :: numroc ! external function in scalapack library
 
    !
    ! Broadcast size of matrix from root process:
    !    
    if ( (myrow == 0) .AND. (mycol == 0) ) then
        call igebs2d(ictxt, 'All', 'i-ring', 1, 1, m, 1 )
    else
        call igebr2d(ictxt, 'All', 'i-ring', 1, 1, m, 1, 0, 0 )
    endif
    if ( (myrow == 0) .AND. (mycol == 0) ) then
        call igebs2d(ictxt, 'All', 'i-ring', 1, 1, n, 1 )
    else
        call igebr2d(ictxt, 'All', 'i-ring', 1, 1, n, 1, 0, 0 )
    endif    
    
    !
    ! Get block size for distributed matrix:
    !   
    if (n == 1) then
        call blockset( nb, 64, m, m, nprow, npcol)  
    else
        call blockset( nb, 64, m, n, nprow, npcol)  
    endif
    
    !
    ! Get local number of rows and columns for the distributed matrix:
    !
    l_m = numroc(m,nb,myrow,0,nprow)
    l_n = numroc(n,nb,mycol,0,npcol)
  
    if (n == 1) then
        l_n = 1
        nb_c = 1
    else
        nb_c = nb
    endif

    
    !
    ! Create the matrix description:
    !
    call descinit( desca, m, n, nb, nb_c, 0, 0, ictxt, l_m, info )  
    
    !
    ! Allocate storage for the distributed matrix:
    !
    allocate (a_loc(l_m, l_n ), stat=istat)
      
    a_loc = 0d0
    
    !write(*,'(a,6(i6,1x))') 'myid,m,n,nb,nb_c,l_m,l_n:',myid,m,n,nb,nb_c, l_m,l_n

    ! kwk debug, need to add istat error check here
    
    end subroutine allocateDistributedMatrix
    
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
    subroutine allocateDistributedVector(m,desca,a_loc)
     
    integer, intent(in)                            :: m
    integer, dimension(descriptor_len), intent(in) :: desca 
    real(8), dimension(:), allocatable,intent(out) :: a_loc
    
    integer :: nb,l_m 
    
    integer :: istat, info
    
    integer :: numroc ! external function in scalapack library
 
    !
    ! Broadcast size of matrix from root process:
    !    
    if ( (myrow == 0) .AND. (mycol == 0) ) then
        call igebs2d(ictxt, 'All', 'i-ring', 1, 1, m, 1 )
    else
        call igebr2d(ictxt, 'All', 'i-ring', 1, 1, m, 1, 0, 0 )
    endif
 
    !
    ! Get block size for distributed matrix:
    !   
    call blockset( nb, 64, m, m, nprow, npcol)  
    
    !
    ! Get local number of rows and columns for the distributed matrix:
    !
    l_m = numroc(m,nb,myrow,0,nprow)
    
    !
    ! Create the matrix description:
    !
    call descinit( desca, m, 1, nb, 1, 0, 0, ictxt, l_m, info )  
    
    !
    ! Allocate storage for the distributed matrix:
    !
    allocate (a_loc(l_m), stat=istat)
      
    a_loc = 0d0
 
    end subroutine allocateDistributedVector
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
   subroutine scalapack_scatterMatrix(a, b, descb)
!
! Scatters matrix A from root node to distributed matrix B  
! 
! Routine that mimic's ScaLAPACK's pdgemr2d but that works for arrays with more than 2^31 elements.
!
! Based on the tutorial C code found at:
! https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
!    


    integer, dimension(:), intent(in)   :: descb(descriptor_len)
    real(8), intent(in)                 :: a(:,:),b(:,:)
    
    integer :: sendr,sendc,recvr,recvc,r,c,nr,nc,ldb,ncol_b,m,n,nbr,nbc
          
    !
    ! Scatter matrix A onto B:
    !
    sendr = 0
    sendc = 0
    recvr = 0
    recvc = 0
    
    ldb = descb(9)
    
    m   = descb(3)
    n   = descb(4)
    nbr = descb(5)
    nbc = descb(6)
    
    ncol_b =  size(b,dim=2)
    
    !write(*,'(a,8(i6,1x))') 'myid,m,n,nbr,nbc,nldb,ncol_b:',myid,m,n,nbr,nbc,ldb,ncol_b

    ! loop over rows to send:
    do r = 1,m,nbr
        sendc = 0
        
        ! number of rows to send, check for last row:
        nr = nbr
        if (m-r+1 < nbr) nr = m-r+1
        
        ! loop over columns to send:
        do c = 1,n,nbc
            
            ! number of cols to send, check for last column:
            nc = nbc
            if (n-c+1 < nbc) nc = n-c+1
            
            ! send it if root process:
            if (myid == 0) then
                call dgesd2d( descb(2),nr,nc, a(r,c), m, sendr, sendc)
            endif
                
            ! receive if it other
            if (myrow ==  sendr .and. mycol == sendc ) then
            
                call  dgerv2d( descb(2), nr, nc, b(recvr+1,recvc+1), ldb, 0, 0)
                recvc = mod(recvc+nc,ncol_b)
 
            endif
            sendc = mod(sendc+1,npcol)
        enddo
        if (myrow == sendr) recvr = mod(recvr+nr,ldb)
        sendr = mod(sendr+1,nprow)
    enddo
    

    end subroutine scalapack_scatterMatrix
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
   subroutine scalapack_gatherMatrix(a, desca, b )
!
! Gathers distributed matrix A to matrix B on root node
!
! Routine that mimic's ScaLAPACK's pdgemr2d but that works for arrays with more than 2^31 elements.
!
! Based on the tutorial C code found at:
! https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
!

    integer, dimension(:), intent(in)   :: desca(descriptor_len)
    real(8), intent(in)                 :: a(:,:),b(:,:)
      
    integer :: sendr,sendc,recvr,recvc,r,c,nr,nc,lda,ncol_a,m,n,nbr,nbc
          
    !
    ! Gather matrix A onto B:
    !       
    sendr = 0
    sendc = 0
    recvr = 0
    recvc = 0
    
    lda = desca(9) ! rows in local matrix storage
    
    m   = desca(3) ! rows in global matrix
    n   = desca(4) ! cols in global matrix
    nbr = desca(5) ! row block factor used to distribute matrix
    nbc = desca(6) ! column block factor used to distribute matrix
    
    ncol_a =  size(a,dim=2)  
    
 
    ! loop over rows to send:
    do r = 1,m,nbr
        sendc = 0
        
        ! number of rows to send, check for last row:
        nr = nbr
        if (m-r+1 < nbr) nr = m-r+1
        
        ! loop over columns to send:
        do c = 1,n,nbc
            
            ! number of cols to send, check for last column:
            nc = nbc
            if (n-c+1 < nbc) nc = n-c+1

            ! send it if it other
            if (myrow ==  sendr .and. mycol == sendc ) then
        
                call dgesd2d( desca(2),nr,nc, a(recvr+1,recvc+1), lda, 0, 0)
                recvc = mod(recvc+nc,ncol_a)
                
            endif
                
            ! receive if it root
             if (myid == 0) then         
                call  dgerv2d( desca(2), nr, nc, b(r,c), m, sendr, sendc)
            endif

            sendc = mod(sendc+1,npcol)
            
            call blacs_barrier(desca(2),'all')  ! this helps keep the root process from being 
            ! overwhelmed by all the dgesd2d data send operations. Sync the processes here so they 
            ! send data one at a time...otherwise the gather will bottleneck and take FOREVER 
            ! compared to the dense matrix operations (as tested on both laptop and TSCC).
            
        enddo
        if (myrow == sendr) recvr = mod(recvr+nr,lda)
        sendr = mod(sendr+1,nprow)
    enddo
    
    end subroutine scalapack_gatherMatrix
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
   subroutine scalapack_scatterVector(a, b, descb)
!
! Scatters vector A from root node to distributed vector B  
! 
! Routine that mimic's ScaLAPACK's pdgemr2d but that works for arrays with more than 2^31 elements.
!
! Based on the tutorial C code found at:
! https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
!    


    integer, dimension(:), intent(in)   :: descb(descriptor_len)
    real(8), intent(in)                 :: a(:),b(:)
    
    integer :: sendr,sendc,recvr,recvc,r,c,nr,nc,ldb,ncol_b,m,n,nbr,nbc
          
    !
    ! Scatter matrix A onto B:
    !
    sendr = 0
    sendc = 0
    recvr = 0
    recvc = 0
    
    ldb = descb(9)
    
    m   = descb(3)
    n   = descb(4)
    nbr = descb(5)
    nbc = descb(6)
    
    ncol_b = 1
    
    !write(*,'(a,8(i6,1x))') 'myid,m,n,nbr,nbc,nldb,ncol_b:',myid,m,n,nbr,nbc,ldb,ncol_b

    ! loop over rows to send:
    do r = 1,m,nbr
        sendc = 0
        
        ! number of rows to send, check for last row:
        nr = nbr
        if (m-r+1 < nbr) nr = m-r+1
        
        ! loop over columns to send:
        do c = 1,n,nbc
            
            ! number of cols to send, check for last column:
            nc = nbc
            if (n-c+1 < nbc) nc = n-c+1
            
            ! send it if root process:
            if (myid == 0) then
                call dgesd2d( descb(2),nr,nc, a(r), m, sendr, sendc)
            endif
                
            ! receive if it other
            if (myrow ==  sendr .and. mycol == sendc ) then
            
                call  dgerv2d( descb(2), nr, nc, b(recvr+1), ldb, 0, 0)
                recvc = mod(recvc+nc,ncol_b)
 
            endif
            sendc = mod(sendc+1,npcol)
        enddo
        if (myrow == sendr) recvr = mod(recvr+nr,ldb)
        sendr = mod(sendr+1,nprow)
    enddo
    

    end subroutine scalapack_scatterVector
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
   subroutine scalapack_gatherVector(a, desca, b )
!
! Gathers distributed vector A to vector B on root node
!
! Routine that mimic's ScaLAPACK's pdgemr2d but that works for arrays with more than 2^31 elements.
!
! Based on the tutorial C code found at:
! https://andyspiros.wordpress.com/2011/07/08/an-example-of-blacs-with-c/
!

    integer, dimension(:), intent(in)   :: desca(descriptor_len)
    real(8), intent(in)                 :: a(:),b(:)
      
    integer :: sendr,sendc,recvr,recvc,r,c,nr,nc,lda,ncol_a,m,n,nbr,nbc
          
    !
    ! Gather matrix A onto B:
    !       
    sendr = 0
    sendc = 0
    recvr = 0
    recvc = 0
    
    lda = desca(9) ! rows in local matrix storage
    
    m   = desca(3) ! rows in global matrix
    n   = desca(4) ! cols in global matrix
    nbr = desca(5) ! row block factor used to distribute matrix
    nbc = desca(6) ! column block factor used to distribute matrix
    
    ncol_a =  1 
    
 
    ! loop over rows to send:
    do r = 1,m,nbr
        sendc = 0
        
        ! number of rows to send, check for last row:
        nr = nbr
        if (m-r+1 < nbr) nr = m-r+1
        
        ! loop over columns to send:
        do c = 1,n,nbc
            
            ! number of cols to send, check for last column:
            nc = nbc
            if (n-c+1 < nbc) nc = n-c+1

            ! send it if it other
            if (myrow ==  sendr .and. mycol == sendc ) then
        
                call dgesd2d( desca(2),nr,nc, a(recvr+1), lda, 0, 0)
                recvc = mod(recvc+nc,ncol_a)
                
            endif
                
            ! receive if it root
             if (myid == 0) then         
                call  dgerv2d( desca(2), nr, nc, b(r), m, sendr, sendc)
            endif

            sendc = mod(sendc+1,npcol)
        enddo
        if (myrow == sendr) recvr = mod(recvr+nr,lda)
        sendr = mod(sendr+1,nprow)
    enddo
    
    end subroutine scalapack_gatherVector
        
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
    subroutine gridsetup(nproc,nprow,npcol)
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol,  that are the sizes of the 2d processors mesh.
!
! Written by Carlo Cavazzoni
!
    integer, intent(in)  :: nproc
    integer, intent(out) :: nprow,npcol
    
    integer sqrtnp,i

    sqrtnp = int( sqrt( dble(nproc) ) + 1 )
    do i=1,sqrtnp
    if(mod(nproc,i).eq.0) nprow = i
    end do
    npcol = nproc/nprow

    end subroutine gridsetup

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
    subroutine blockset( nb, nbuser, m,n, nprow, npcol)
!
!     This subroutine try to choose an optimal block size
!     for the distributd matrix.
!
!     Written by Carlo Cavazzoni, CINECA
!
    integer, intent(in)  :: m, n
    integer, intent(in)  :: nprow, npcol,nbuser
    integer, intent(out) :: nb

    nb = min ( m/nprow, n/npcol )
    if(nbuser.gt.0) then
    nb = min ( nb, nbuser )
    endif
    nb = max(nb,1)

    end subroutine blockset   

!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------  
      
    end module mare2dem_scalapack