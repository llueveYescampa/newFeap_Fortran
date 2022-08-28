subroutine prtdis(x,b,ttim,prop,ndm,ndf,n1,n2,n3)
implicit  none
integer   ndm,ndf,n1,n2,n3
double precision    x(ndm,*),b(ndf,*),ttim,prop

!  Purpose: Output nodal displacement values

!  Inputs:
!     x(ndm,*)  - Nodal coordinates for mesh
!     b(ndf,*)  - solution vector for nodes
!     ttim      - Current time in solution
!     prop      - Current proportional load value
!     ndm       - Spatial dimension of mesh
!     ndf       - Number dof/node
!     n1        - First node to output
!     n2        - Last  node to output
!     n3        - Increment to n1 for next node to output, etc.

!  Outputs:
!     none

   character cd*6,di*6
   integer   i, n, kount

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   data cd/' coord'/,di/' displ'/

   kount = 0
   do n = n1,n2,n3
     kount = kount - 1
     if(kount.le.0) then
       call prthed(ioWrite)
       write(ioWrite,2000) ttim,prop,(i,cd,i=1,ndm),(i,di,i=1,ndf)
       if(ioRead.lt.0) then
	 write(*,2000) ttim,prop,(i,cd,i=1,ndm),(i,di,i=1,ndf)
       end if
       kount = 48
     endif
     if(x(1,n).ne. -999.) then
       write(ioWrite,'(i6,1p9e13.6)') n,(x(i,n),i=1,ndm),(b(i,n),i=1,ndf)
       if(ioRead.lt.0) then
         write(*,'(i6,1p9e13.6)') n,(x(i,n),i=1,ndm),(b(i,n),i=1,ndf)
       end if
     else
       write(ioWrite,'(i6,a)') n,' not input.'
       if(ioRead.lt.0) then
         write(*,'(i6,a)') n,' not input.'
       end if  
     end if
   end do  

2000  format('  N o d a l   D i s p l a c e m e n t s',5x,       &
             'time',e18.5/31x,'prop. ld. (eigenvalue)',e13.5//   &
             '  node',9(i7,a6))

end
