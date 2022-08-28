subroutine prtrea(r,ndf,numnp,n1,n2,n3)
implicit  none
integer   ndf,numnp, n1,n2,n3
double precision    r(ndf,1)

!  Purpose: Print nodal reactions

!  Inputs:
!     r(ndf,*)  - Nodal reactions
!     ndf       - Number dof/node
!     numnp     - Number of nodes in mesh
!     n1        - First node to output
!     n2        - Last  node to output
!     n3        - Increment to n1 for next node to output, etc.

!  Outputs:
!     none


   integer   i, k, n, kount

   double precision    rr(6),rsum(6),asum(6),psum(6)

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   call pconsd(rsum,ndf,0.0d0)
   call pconsd(psum,ndf,0.0d0)
   call pconsd(asum,ndf,0.0d0)

   do i = 1,numnp
     do k = 1,ndf
       rsum(k) = rsum(k) - r(k,i)
       asum(k) = asum(k) + abs(r(k,i))
     end do
   end do  

   kount = 0
   do n = n1,n2,n3
     kount = kount - 1
     if(kount.le.0) then
       call prthed(ioWrite)
       write(ioWrite,2000) (k,k=1,ndf)
       if(ioRead.lt.0) write(*,2000) (k,k=1,ndf)
       kount = 50
     end if

     do k = 1,ndf
       rr(k) = -r(k,n)
       psum(k) = psum(k) + rr(k)
     end do  

     write(ioWrite,'(i10,6e13.4)') n,(rr(k),k=1,ndf)
     if(ioRead.lt.0) then 
       write(*,'(i10,6e13.4)') n,(rr(k),k=1,ndf)
     end if  
   end do  

!  Print statics check

   write(ioWrite,'(/7x,a,6e13.4)') 'sum', (rsum(k),k=1,ndf)
   write(ioWrite,'( 3x,a,6e13.4)') 'prt sum',(psum(k),k=1,ndf)
   write(ioWrite,'( 3x,a,6e13.4)') 'abs sum',(asum(k),k=1,ndf)
   if(ioRead.lt.0) then
     write(*,'(/7x,a,6e13.4)') 'sum',(rsum(k),k=1,ndf)
     write(*,'( 3x,a,6e13.4)') 'prt sum', (psum(k),k=1,ndf)
     write(*,'( 3x,a,6e13.4)') 'abs sum', (asum(k),k=1,ndf)
   end if

2000  format('  N o d a l    R e a c t i o n s'//6x,'node',6(i9,' dof'))

end
