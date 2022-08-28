subroutine matl05(d)
implicit none
double precision d(13)

   logical         pcomp
   character       typ*5
   integer iocheck
   
   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   character*80    yyy
   common /ydata/  yyy

   do while (.true.)
     if(ioRead.lt.0) then
       write(*,3000)
       read(*,'(a5)') typ
     else
       read(ioRead,'(a5)') typ
     end if
     call pintio(yyy,10)
     read(yyy,'(7f10.0)',IOSTAT=iocheck) d
     if (iocheck .eq. 0) then 
       if(pcomp(typ,'beam')) then
         d(2)  = 0.0d0
         d(13) = 0.0d0
       else
         d(13) = 1.0d0
       end if
       write(ioWrite,2000) typ,d
       if(ioRead.lt.0) then 
         write(*,2000) typ,d
       end if  
   
!      Set beam/shell in-plane and bending stiffness values
   
       d(7)  = d(1)*d(3)/(1.0d0 - d(2)*d(2))
       d(8)  = d(2)*d(7)
       d(9)  = d(7)*d(3)*d(3)/12.d0
       d(10) = d(2)*d(9)
       d(11) = d(6)*d(3)
       d(12) = d(3)*d(3)/12.d0
       exit
     else
       call pperror('PCELM5',yyy)
     end if
   end do
   return

2000  format(5x,'Rectangular Beam/Axisymmetric Shell Model'/  &
       5x,'Type: ',a5/5x,'Elastic Modulus',g15.4/             &
       5x,'Poisson Ratio',g17.4/5x,'Thickness  ',g19.4/       &
       5x,'Normal Load  ',g17.4/                              &
       5x,'Linear(0=l 1=nl)',f9.1/5x,'Mass Density',g18.4/1x)

3000  format(' Input: 1. Type (beam or shell)'/               &
             '        2. E nu h press lin rho'/3x,'>',$)
end
