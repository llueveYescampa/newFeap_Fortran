subroutine prtstr(dt,st,numnp,n1,n2,n3)
implicit none
  integer          :: numnp,n1,n2,n3
  double precision :: dt(numnp),st(numnp,7)

!  Purpose: Output projected nodal stress values

!  Inputs:
!     dt(*)       - Value of 8-th nodal value
!     st(numnp,*) - Value of 1-7  nodal values
!     numnp       - Number of nodes in mesh
!     n1          - First node to output
!     n2          - Last  node to output
!     n3          - Increment to n1 for next node to output, etc.

   integer   i, kount, n

   include 'iofile.h'
   kount = 0
   do n = n1,n2,n3
     kount = kount - 1
     if(kount.le.0) then
       call prthed(iow)
       write(iow,2000)
       if(ior.lt.0) then 
         write(*,2000)
       end if  
       kount = 17
     end if
     write(iow,'(i5,4e13.5/5x,4e13.5/1x)') n,(st(n,i),i=1,7),dt(n)
     if(ior.lt.0) then 
       write(*,'(i5,4e13.5/5x,4e13.5/1x)') n,(st(n,i),i=1,7),dt(n)
     end if  
   end do  

2000 format('   N o d al   S t r e s s e s'//' node',4x,'11-stress',4x, &
     '12-stress',4x,'22-stress',4x,'33-stress'/10x,'1-stress',5x,      &
     '2-stress',4x,'max-shear',8x,'angle')

end
