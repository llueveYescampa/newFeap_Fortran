subroutine elmt03(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
implicit none
  integer          :: ix(*), ndf,ndm,nst,isw
  double precision :: d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,*),p(*)

!  Purpose: Plane stress/strain Pian Sumihara Element

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


   integer   i,is,j
   double precision ssa,tta,ssg,ttg,xx,yy
   double precision sig(4),p1,p2,p3
   character*4  wd(2)
   character yyy*80
   
   integer iocheck

   include 'adata.h'
   include 'cdata.h'
   include 'eldata.h'
   include 'elcom3.h'
   include 'iofile.h'

   data wd/'ress','rain'/

!  Go to correct array processor


!  Input/output material properties
   select case (isw)
   case (1)
     do while (.true.)
       if(ior.lt.0) then
         write(*,3000)
       end if  
       call pintio(yyy,10)
       read(yyy,'(4f10.0,i10)',IOSTAT=iocheck) (d(i),i=8,11),is
       if (iocheck .eq. 0) then
         is = max(1,min(is,2))
         d(12)= is
         d(2) = d(9)*d(8)/(1.+d(9))/(1.-is * d(9))
         d(4) = d(8)/(1.+d(9))
         d(3) = d(4)/2.
         d(1) = d(2) + d(4)
         if(is.eq.1) then
!          Set parameters for plane stress (is = 1)
     
           if(d(11).le.0.0d0) d(11) = 1.0
           d(5) =  1./d(8)
           d(6) = -d(9)/d(8)
           d(13)=  0.0
         else
!          Set parameters for plane strain (is = 2)
     
           d(11)=  1.0d0
           d(5) =  (1.-d(9))/d(4)
           d(6) = -d(9)/d(4)
           d(13)=  d(9)
         end if
         write(iow,2000) wd(is),(d(i),i=8,11)
         if(ior.lt.0) then
           write(*,2000) wd(is),(d(i),i=8,11)
         end if  
         d(10) = d(10)*d(11)
         return
       else
         call myPerror('elmt03',yyy)
       end if
     end do
   case (2)
!    Check mesh
     call ckisop(ix,xl,aa,ndm)
     return
   case (3,4,5,6,8)
!    Compute Pian-Sumihara arrays for elastic: compute jacobian
     xs = (-xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4))/4.
     ys = (-xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4))/4.
     xt = (-xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4))/4.
     yt = (-xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4))/4.
     xh = ( xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4))/4.
     yh = ( xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4))/4.
     xj0 = xs*yt - xt*ys
     xj1 = xs*yh - xh*ys
     xj2 = xh*yt - xt*yh
     if(isw.eq.5) then 
!      Compute a lumped mass matrix
       
       p(      1) = (xj0-(xj1+xj2)/3.)*d(10)
       p(  ndf+1) = (xj0+(xj1-xj2)/3.)*d(10)
       p(2*ndf+1) = (xj0+(xj1+xj2)/3.)*d(10)
       p(3*ndf+1) = (xj0-(xj1-xj2)/3.)*d(10)
       do i = 2,nst,ndf
         p(i) = p(i-1)
       end do ! i
       return
     end if  
     ssa = xj1/xj0/3.d0
     tta = xj2/xj0/3.d0
!    Form stiffness for elastic part and compute the beta parameters
     
     call pian03(d,ul,s,p,nst,ndf,isw)
     if(isw.eq.4) then
!      Compute the stresses

       is = d(12)
!      Compute the stresses at the center and the specified points
       
       ssg    = -ssa*beta(5)
       ttg    = -tta*beta(4)
       sig(1) =  beta(1) + a1(1)*ttg + a2(1)*ssg
       sig(2) =  beta(2) + a1(2)*ttg + a2(2)*ssg
       sig(3) =  beta(3) + a1(3)*ttg + a2(3)*ssg
       sig(4) =  d(13)*(sig(1)+sig(2))
       call pstres(sig,p1,p2,p3)
       xx = (xl(1,1)+xl(1,2)+xl(1,3)+xl(1,4))/4.0
       yy = (xl(2,1)+xl(2,2)+xl(2,3)+xl(2,4))/4.0
       mct = mct - 1
       if(mct.le.0) then 
         call prthed(iow)
         write(iow,2001) wd(is)
         if(ior.lt.0) then
           write(*,2001) wd(is)
         end if  
         mct = 25
       end if  
       write(iow,'(2i9,2f12.4,2e12.4,f8.2/30x,4e12.4)') &
              n,ma,xx,yy,p1,p2,p3,(sig(i),i=1,4)
       if(ior.lt.0) then 
         write(*,'(2i9,2f12.4,2e12.4,f8.2/30x,4e12.4)') &
              n,ma,xx,yy,p1,p2,p3,(sig(i),i=1,4)
       end if  
       return
     else if(isw.eq.8) then 
!      Compute the nodal stress values
     
       call stcn03(ix,d,ssa,tta,aa,aa(numnp+1),numnp)
       return
     end if  
     
!    Compute symetric part of s
     
       do i = 1,nst
         do j = i,nst
           s(j,i) = s(i,j)
         end do ! j
       end do ! i
     return
   
   case (7)
!    Error estimator goes here!
     return
   end select

!  Format statements


2000  format(5x,'Plane St',a4,' Element'//10x,'modulus      =',e12.5/ & 
        10x,'poisson ratio=', f8.5/10x,'mass density =',e12.5/        &
        10x,'thickness    =', e12.5)

2001  format(5x,'Plane St',a4,' Stresses'//'  element material',   &
        5x,'1-coord',5x,'2-coord',8x,'sig1',8x,'sig2',3x,'angle'/  &
        38x,'s-11',8x,'s-12',8x,'s-22',8x,'s-33'/1x)

3000  format(' Input: e, nu, rho, th, is (1=pl.stress,2=pl.strain)'/ &
          3x,'mate>',$)

end
