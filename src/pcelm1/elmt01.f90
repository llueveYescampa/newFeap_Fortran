subroutine elmt01(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
implicit  none
integer   ix(*),                    ndf,ndm,nst,isw
double precision    d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,*),p(*)

!  Purpose: Plane linear elastic displacement model element routine

!  Inputs:
!     d(*)      - Material set parameters
!     ul(ndf,*) - Solution variables for element
!     xl(ndf,*) - Nodal coordinates for element
!     ix(*)     - Nodal connections for element
!     tl(*)     - Nodal temperatures for element
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nst       - Size of element matrix
!     isw       - Switch parameter to control action taken by routine

!  Outputs:
!     s(nst,*)  - Element matrix
!     p(*)      - Element vector

   integer   ityp, j,j1, k,k1, l,lint
   
   double precision a11,a21,a31,a41,a12,a32,a42, w11,w12,w22
   double precision e,xnu,alp,t0, xx,yy,zz, dv, xsj, sigr4

   double precision eps(4),sigr(6),shp(3,9),sg(16),tg(16),wg(16),ang
   character wd(3)*12 !,yyy*80

   integer iocheck1, iocheck2


   include 'maxa.h'      

   include 'adata.h'

   include 'cdata.h'
   include 'eldata.h'
   include 'iofile.h'
   include 'ydata.h'

   data wd/'Plane Stress','Plane Strain','Axisymmetric'/

!  Go to correct array processor

   l    = int(d(5))
   k    = int(d(6))
   ityp = int(d(15))
   lint = 0
   select case (isw)
   case(1)
!    Input material properties
     do while (.true.)
       if(ioRead.lt.0) then
         write(*,5000)
       end if
       call pintio(yyy,10)
       read(yyy,'(3f10.0,3i10)',IOSTAT=iocheck1) e,xnu,d(4),l,k,ityp
       if (iocheck1 .eq. 0) then
         do while (.true.)
           if(ioRead.lt.0) then 
             write(*,5001)
           end if  
           call pintio(yyy,10)
           read(yyy,'(8f10.0)',IOSTAT=iocheck2) d(14),d(11),d(12),alp,t0
           if (iocheck2 .eq. 0) then
!            Set material parameter type and flags
             ityp = max(1,min(ityp,3))
             j = min(ityp,2)
             d(1) = e*(1.+(1-j)*xnu)/(1.+xnu)/(1.-j*xnu)
             d(2) = xnu*d(1)/(1.+(1-j)*xnu)
             d(3) = e/2./(1.+xnu)
             d(13)= d(2)*(j-1)
             if(d(14).le.0.0d0 .or. ityp.ge.2) then 
               d(14) = 1.0
             end if  
             d(15) = ityp
             d(16) = e
             d(17) = xnu
             d(18) = -xnu/e
             l = min(4,max(1,l))
             k = min(4,max(1,k))
             d(5) = l
             d(6) = k
             d(9) = t0
             d(10)= e*alp/(1.-j*xnu)
             lint = 0
             write(ioWrite,2000) wd(ityp),d(16),d(17),d(4),l,k, &
                           d(14),d(11),d(12),alp,t0
             if(ioRead.lt.0) then
               write(*,2000) wd(ityp),d(16),d(17),d(4),l,k, &
                           d(14),d(11),d(12),alp,t0
             end if
             d(4) = d(4)*d(14)
             return
           else
!            Read error messages
             call pperror('PCELM1',yyy)
           end if
         end do
       else
!        Read error messages
         call pperror('PCELM1',yyy)
       end if
     end do
   case(2)
!    Check mesh
     call ckisop(ix,xl,shp,ndm)
     return
   case(3,4,5,6)
!    Stiffness/residual computation
     
     if(isw.eq.4) l = k
     if(l*l.ne.lint) call pgauss(l,lint,sg,tg,wg)
     
!    Compute integrals of shape functions
     
     do l = 1,lint
       call shapeFunc(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
       call stre01(d,xl,ul,tl,shp,eps,sigr,xx,yy,ndm,ndf,nel,ityp)
       xsj = xsj*wg(l)*d(14)
     
!      Compute jacobian correction
     
       if(ityp.le.2) then
         dv  = xsj
         xsj = 0.0
         zz  = 0.0
         sigr4 = -d(11)*dv
       else
         dv  = xsj*xx
         zz  = 1./xx
         sigr4 = sigr(4)*xsj - d(11)*dv
       endif
       j1 = 1
     
!      Compute the mass term
     
       if(isw.eq.5) then
         dv = dv*d(4)
         do j = 1,nel
           p(j1  ) = p(j1) + shp(3,j)*dv
           p(j1+1) = p(j1)
           j1 = j1 + ndf
         end do  
       else if(isw.eq.4) then
         call pstres(sigr,sigr(5),sigr(6),ang)
     
!        Output stresses and strains
     
         mct = mct - 2
         if(mct.le.0) then
           call prthed(ioWrite)
           write(ioWrite,2001)
           if(ioRead.lt.0) write(*,2001)
           mct = 50
         end if
         write(ioWrite,2002) n,xx,sigr,ma,yy,eps,ang
         if(ioRead.lt.0) write(*,2002) n,xx,sigr,ma,yy,eps,ang
       else
     
!        Loop over rows
     
         do j = 1,nel
           w11 = shp(1,j)*dv
           w12 = shp(2,j)*dv
           w22 = shp(3,j)*xsj
     
!          Compute the internal forces
     
           p(j1  ) = p(j1  ) - (shp(1,j)*sigr(1)+shp(2,j)*sigr(2))*dv &
                               -  shp(3,j)*sigr4
           p(j1+1) = p(j1+1) - (shp(1,j)*sigr(2)+shp(2,j)*sigr(3))*dv &
                               + d(12)*shp(3,j)*dv
       
!          Loop over columns (symmetry noted)
     
           if(isw.eq.3) then
             k1 = j1
             a11 = d(1)*w11 + d(2)*w22
             a21 = d(2)*w11 + d(1)*w22
             a31 = d(2)*(w11+w22)
             a41 = d(3)*w12
             a12 = d(2)*w12
             a32 = d(1)*w12
             a42 = d(3)*w11
             do k = j,nel
               w11 = shp(1,k)
               w12 = shp(2,k)
               w22 = shp(3,k)*zz
               s(j1  ,k1  ) = s(j1  ,k1  ) + w11*a11+w22*a21+w12*a41
               s(j1+1,k1  ) = s(j1+1,k1  ) + (w11 + w22)*a12+w12*a42
               s(j1  ,k1+1) = s(j1  ,k1+1) + w12*a31 + w11*a41
               s(j1+1,k1+1) = s(j1+1,k1+1) + w12*a32 + w11*a42
               k1 = k1 + ndf
           end do
           end if
           j1 = j1 + ndf
         end do  
       end if
     end do  
       
!    Make stiffness symmetric and compute a residual
     
     if(isw.eq.3) then
       do j = 1,nst
         do k = j,nst
           s(k,j) = s(j,k)
         end do
       end do  
     end if
     return
   case(7)
!    Compute the stress errors
     if(l*l.ne.lint) then 
       call pgauss(l,lint,sg,tg,wg)
     end if  
     call ster01(ix,d,xl,ul,tl,shp,aa(numnp+1),ndf,ndm, &
                numnp,numel,sg,tg,wg,sigr,eps,lint,ityp)
     return
   case(8)
!      Compute the nodal stress values

     if(l*l.ne.lint) then
       call pgauss(l,lint,sg,tg,wg)
     end if  
     call stcn01(ix,d,xl,ul,tl,shp,aa,aa(numnp+1),ndf,ndm,nel, &
                numnp,sg,tg,sigr,eps,lint,ityp)
     return
   end select 
      


!     Formats for input-output

2000  format(/5x,a12,' Linear Elastic Element'//                       &
     10x,'Modulus',e18.5/10x,'Poisson ratio',f8.5/10x,'Density',e18.5/ &
     10x,'Gauss pts/dir',i3/10x,'Stress pts',i6/10x,'Thickness',e16.5/ &
     10x,'1-gravity',e16.5/10x,'2-gravity',e16.5/10x,'Alpha',e20.5/    &
     10x,'Base temp',e16.5/)
2001  format(5x,'Element Stresses'//' elmt  1-coord',2x,'11-stress',2x,&
     '12-stress',2x,'22-stress',2x,'33-stress',3x,'1-stress',3x,       &
     '2-stress'/' matl  2-coord',2x,'11-strain',2x,'12-strain',2x,     &
     '22-strain',2x,'33-strain',6x,'angle'/39(' -'))
2002  format(i5,0p1f9.3,1p6e11.3/i5,0p1f9.3,1p4e11.3,0p1f11.2/)
5000  format(' Input: E, nu, rho, pts/stiff, pts/stre',                &
    ', type(1=stress,2=strain,3=axism)',/3x,'>',$)
5001  format(' Input: Thickness, 1-body force, 2-body force, alpha,'   &
           ,' Temp-base'/3x,'>',$)

end
