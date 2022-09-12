  module sort
  
    implicit none
    
    private
    
!
! Public subroutines:
!       
    public ::   quick_sort     ! call quick_sort(ilist1,ilist2,zlist1) 
    public ::   quick_sort_d   ! call quick_sort(dlist1,dlist2,ilist2) 

    contains


!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------            
   
    recursive subroutine quick_sort(ilist1,ilist2,zlist1,dlist1) 
    !
    ! Does quicksort based on order of integers in ilist1, applies same ordering to optionals ilist2, zlist1 and dlist1.  
    !
    
! KWK May 2018: might be faster to just allocate another index array outside and use that to post-sort real/complex arrays 
! since the code below is getting quite heinous...
    
    !
    ! Modified from Brainerd's F2003 book
    !
       integer,    dimension(:), intent(in out)             :: ilist1
       integer,    dimension(:), intent(in out), optional   :: ilist2
       complex(8), dimension(:), intent(in out), optional   :: zlist1
       real(8),    dimension(:), intent(in out), optional   :: dlist1
       
       integer              :: i, j, n 
       integer              :: chosen, temp
       complex(8)           :: ztemp
       real(8)              :: dtemp
       integer, parameter   :: max_simple_sort_size = 6
       
       n = size(ilist1) 
       
       if (n <= max_simple_sort_size) then 
          ! Use interchange sort for small lists 
          if ( (present(ilist2)) .and. (present(zlist1)) .and. (present(dlist1)) ) then
                call interchange_sort(ilist1,ilist2=ilist2,zlist1=zlist1,dlist1=dlist1) 
          elseif ( (present(zlist1)) .and. (present(dlist1)) ) then
                call interchange_sort(ilist1,zlist1=zlist1,dlist1=dlist1) 
          elseif ( (present(ilist2)) .and. (present(zlist1)) ) then
                call interchange_sort(ilist1,ilist2,zlist1) 
          elseif ( (present(ilist2)) .and. (present(dlist1)) ) then
                call interchange_sort(ilist1,ilist2,dlist1=dlist1)                 
          elseif ( (present(zlist1)) ) then
                call interchange_sort(ilist1,zlist1=zlist1)
          elseif ( (present(dlist1)) ) then
                call interchange_sort(ilist1,dlist1=dlist1)        
          elseif ( (present(ilist2)) ) then
                call interchange_sort(ilist1,ilist2=ilist2)            
          else
                call interchange_sort(ilist1) 
          endif
       else 
          ! Use partition (“quick”) sort 
          chosen = ilist1(n/2) 
          i = 0 
          j = n + 1 
          do 
             ! Scan list from left end 
             ! until element >= chosen is found 
             do 
                i = i + 1 
                if (ilist1(i) >= chosen) exit 
             end do 
             ! Scan list from right end 
             ! until element <= chosen is found 
             do 
                j = j - 1 
                if (ilist1(j) <= chosen) exit 
             end do    
             if (i < j) then 
             
                ! Swap two out of place elements 
                temp      = ilist1(i) 
                ilist1(i) = ilist1(j) 
                ilist1(j) = temp                 
                 
                if (present(ilist2)) then
                    temp      = ilist2(i) 
                    ilist2(i) = ilist2(j) 
                    ilist2(j) = temp                       
                endif
                 
                if  (present(zlist1)) then  
                    ztemp     = zlist1(i) 
                    zlist1(i) = zlist1(j) 
                    zlist1(j) = ztemp 
                endif
                
                if  (present(dlist1)) then  
                    dtemp     = dlist1(i) 
                    dlist1(i) = dlist1(j) 
                    dlist1(j) = dtemp 
                endif
                
             else if (i == j) then 
                i = i + 1 
                exit 
             else 
                exit 
             end if 
          end do 
          if ( (present(ilist2)).and.(present(zlist1)) .and. (present(dlist1)) ) then 
                if (1 < j) call quick_sort(ilist1(:j),ilist2=ilist2(:j),zlist1=zlist1(:j),dlist1=dlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),ilist2=ilist2(i:),zlist1=zlist1(i:),dlist1=dlist1(i:)) 
          elseif ( (present(zlist1)) .and. (present(dlist1)) ) then 
                if (1 < j) call quick_sort(ilist1(:j),zlist1=zlist1(:j),dlist1=dlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),zlist1=zlist1(i:),dlist1=dlist1(i:)) 
          elseif ( (present(ilist2)) .and. (present(zlist1)) ) then
                if (1 < j) call quick_sort(ilist1(:j),ilist2(:j),zlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),ilist2(i:),zlist1(i:)) 
          elseif ( (present(ilist2)) .and. (present(dlist1)) ) then
                if (1 < j) call quick_sort(ilist1(:j),ilist2(:j),dlist1=dlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),ilist2(i:),dlist1=dlist1(i:))                 
          elseif ( (present(zlist1)) ) then
                if (1 < j) call quick_sort(ilist1(:j),zlist1=zlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),zlist1=zlist1(i:))  
          elseif ( (present(dlist1)) ) then
                if (1 < j) call quick_sort(ilist1(:j),dlist1=dlist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:),dlist1=dlist1(i:))        
          elseif ( (present(ilist2)) ) then
                if (1 < j) call quick_sort(ilist1(:j),ilist2=ilist2(:j)) 
                if (i < n) call quick_sort(ilist1(i:),ilist2=ilist2(i:))           
          else
                if (1 < j) call quick_sort(ilist1(:j)) 
                if (i < n) call quick_sort(ilist1(i:))      
          endif          
          
 
       end if  ! test for small array 
    end subroutine quick_sort 

!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------            
 
    subroutine interchange_sort(ilist1,ilist2,zlist1,dlist1)
     
       integer,    dimension(:), intent(in out)             :: ilist1
       integer,    dimension(:), intent(in out), optional   :: ilist2
       complex(8), dimension(:), intent(in out), optional   :: zlist1
       real(8),    dimension(:), intent(in out), optional   :: dlist1
       
       integer      :: i, j 
       integer      :: temp 
       complex(8)   :: ztemp
       real(8)      :: dtemp
       
       do i = 1, size(ilist1) - 1 
       
          do j = i + 1, size(ilist1) 
          
             if (ilist1(i) >  ilist1(j)) then 
             
                temp      = ilist1(i) 
                ilist1(i) = ilist1(j) 
                ilist1(j) = temp 
                             
                if (present(ilist2)) then
                    temp      = ilist2(i) 
                    ilist2(i) = ilist2(j) 
                    ilist2(j) = temp                       
                endif
                 
                if  (present(zlist1)) then  
                    ztemp     = zlist1(i) 
                    zlist1(i) = zlist1(j) 
                    zlist1(j) = ztemp 
                endif
        
                if  (present(dlist1)) then  
                    dtemp     = dlist1(i) 
                    dlist1(i) = dlist1(j) 
                    dlist1(j) = dtemp 
                endif
                                
             end if 
          end do 
       end do 
    end subroutine interchange_sort     
    
!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------            
   
    recursive subroutine quick_sort_d(dlist1,dlist2,ilist2) 
    !
    ! Does quicksort based on order of integers in ilist1, applies same ordering to optionals ilist2 and zlist1.  
    !
    
    !
    ! Modified from Brainerd's F2003 book
    !
       real(8),    dimension(:), intent(in out)             :: dlist1
       real(8),    dimension(:), intent(in out), optional   :: dlist2
       integer,    dimension(:), intent(in out), optional   :: ilist2
   
       
       integer              :: i, j, n 
       integer              :: temp
       real(8)              :: dtemp, chosen
       integer, parameter   :: max_simple_sort_size = 6
       
       n = size(dlist1) 
       
       if (n <= max_simple_sort_size) then 
          ! Use interchange sort for small lists 
           if ( (present(ilist2)).and.(present(dlist2)) ) then
                call interchange_sort_d(dlist1,dlist2=dlist2,ilist2=ilist2)            
          elseif ( (present(ilist2)) ) then
                call interchange_sort_d(dlist1,ilist2=ilist2)         
          elseif ( (present(dlist2)) ) then
                call interchange_sort_d(dlist1,dlist2=dlist2)          
          else
                call interchange_sort_d(dlist1) 
          endif
       else 
          ! Use partition (“quick”) sort 
          chosen = dlist1(n/2) 
          i = 0 
          j = n + 1 
          do 
             ! Scan list from left end 
             ! until element >= chosen is found 
             do 
                i = i + 1 
                if (dlist1(i) >= chosen) exit 
             end do 
             ! Scan list from right end 
             ! until element <= chosen is found 
             do 
                j = j - 1 
                if (dlist1(j) <= chosen) exit 
             end do    
             if (i < j) then 
             
                ! Swap two out of place elements 
                dtemp     = dlist1(i) 
                dlist1(i) = dlist1(j) 
                dlist1(j) = dtemp                 
                 
                if (present(ilist2)) then
                    temp      = ilist2(i) 
                    ilist2(i) = ilist2(j) 
                    ilist2(j) = temp                       
                endif
                
                if (present(dlist2)) then
                    dtemp     = dlist2(i) 
                    dlist2(i) = dlist2(j) 
                    dlist2(j) = dtemp                       
                endif
                
             else if (i == j) then 
                i = i + 1 
                exit 
             else 
                exit 
             end if 
          end do 
          
          if ( (present(ilist2)).and.(present(dlist2)) ) then
                if (1 < j) call quick_sort_d(dlist1(:j),dlist2=dlist2(:j),ilist2=ilist2(:j)) 
                if (i < n) call quick_sort_d(dlist1(i:),dlist2=dlist2(i:),ilist2=ilist2(i:))      
          elseif ( (present(ilist2)) ) then
                if (1 < j) call quick_sort_d(dlist1(:j),ilist2=ilist2(:j)) 
                if (i < n) call quick_sort_d(dlist1(i:),ilist2=ilist2(i:))   
          elseif ( (present(dlist2)) ) then
                if (1 < j) call quick_sort_d(dlist1(:j),dlist2=dlist2(:j)) 
                if (i < n) call quick_sort_d(dlist1(i:),dlist2=dlist2(i:))                           
          else
                if (1 < j) call quick_sort_d(dlist1(:j)) 
                if (i < n) call quick_sort_d(dlist1(i:))      
          endif          
          
 
       end if  ! test for small array 
    end subroutine quick_sort_d 

!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------        
!-----------------------------------------------------------------------------------------------------------------------------------            
 
    subroutine interchange_sort_d(dlist1,dlist2,ilist2)
     
       real(8),    dimension(:), intent(in out)             :: dlist1
       real(8),    dimension(:), intent(in out), optional   :: dlist2
       integer,    dimension(:), intent(in out), optional   :: ilist2
        
       
       integer      :: i, j 
       integer      :: temp 
       real(8)      :: dtemp
       
       do i = 1, size(dlist1) - 1 
       
          do j = i + 1, size(dlist1) 
          
             if (dlist1(i) >  dlist1(j)) then 
             
                dtemp     = dlist1(i) 
                dlist1(i) = dlist1(j) 
                dlist1(j) = dtemp 
                
                if (present(dlist2)) then
                    dtemp     = dlist2(i) 
                    dlist2(i) = dlist2(j) 
                    dlist2(j) = dtemp                       
                endif          
                 
                if (present(ilist2)) then
                    temp      = ilist2(i) 
                    ilist2(i) = ilist2(j) 
                    ilist2(j) = temp                       
                endif
                
             end if 
          end do 
       end do 
    end subroutine interchange_sort_d         
    
    
    end   module sort