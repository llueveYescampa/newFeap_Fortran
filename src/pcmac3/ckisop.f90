subroutine ckisop(ix,xl,shp,ndm)
implicit  none
integer   ix(*),ndm
double precision    xl(ndm,*),shp(3,*)
      
!  Purpose: Check isoparametric elements for errors

!  Inputs:
!     ix(*)     - Node numbers on element
!     xl(ndm,*) - Nodal coordinates for element
!     shp(3,*)  - Shape function storage
!     ndm       - Spatial dimension of mesh

!  Outputs:
!     none

   integer   i,ineg, l
   integer   xn(9),yn(9),ic(18)
   double precision  ss,tt,xsj

   double precision dm
   integer             n,ma,mct,iel,nel
   common /eldata/  dm,n,ma,mct,iel,nel

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   data xn/-1,1,1,-1,0,1,0,1,0/
   data yn/-1,-1,1,1,-1,0,1,0,0/

!  check the element for input errors

   ineg = 0
   do l = 1,nel
     if(xl(1,l).eq. -999.0 .and. ix(l).ne.0) then
       ic(ineg+1) = l
       ic(ineg+2) = abs(ix(l))
       ineg = ineg + 2
     end if
   end do  
   if(ineg.gt.0) then
     write(ioWrite,2000) n,(ic(i),i=1,ineg)
     if(ioRead.lt.0) then 
       write(*,2000) n,(ic(i),i=1,ineg)
     end if  
   else
     do l = 1,nel
       ss = xn(l)
       tt = yn(l)
       call  shapeFunc(ss,tt,xl,shp,xsj,ndm,nel,ix,.false.)
       if(xsj.le.0.0d0) then
         ic(ineg+1) = l
         ic(ineg+2) = abs(ix(l))
         ineg = ineg + 2
       end if
     end do  
     if(ineg.gt.0) then
       write(ioWrite,2001) n,(ic(i),i=1,ineg)
       if(ioRead.lt.0) then 
         write(*,2001) n,(ic(i),i=1,ineg)
       end if  
     end if
   end if
      
2000  format(' >Element',i4,' coordinates not input for nodes:'/ &
            ('                Local =',i3,' Global =',i4))
2001  format(' >Element',i4,' has negative jacobian at nodes:'/  &
            ('                Local =',i3,' Global =',i4))
end
