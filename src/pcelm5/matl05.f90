subroutine matl05(d)
implicit none
  double precision :: d(*)

   logical         pcomp
   character       typ*5
   integer iocheck
   character :: yyy*80
   
   include 'iofile.h'

   do while (.true.)
     if(ior.lt.0) then
       write(*,3000)
       read(*,'(a5)') typ
     else
       read(ior,'(a5)') typ
     end if
     call pintio(yyy,10)
     read(yyy,'(7f10.0)',IOSTAT=iocheck) d(1:6)
     if (iocheck .eq. 0) then 
       if(pcomp(typ,'beam')) then
         d(2)  = 0.0d0
         d(13) = 0.0d0
       else
         d(13) = 1.0d0
       end if
       write(iow,2000) typ,d(1),d(2),d(3),d(4),d(5),d(6)
       if(ior.lt.0) then 
         write(*,2000) typ,d(1),d(2),d(3),d(4),d(5),d(6)
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
       call myPerror('matl05',yyy)
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
