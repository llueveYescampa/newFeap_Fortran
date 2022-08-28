subroutine elmt06(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
implicit  none
integer  ix(*), ndf,ndm,nst,isw
double precision  d(*),ul(ndf,*),xl(ndm,*),tl(*),s(nst,nst),p(nst)

!  Purpose: Two dimensional laplace equation with a reaction term
  
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
  
  
   character wlab(2)*12
   integer   i,i1,j,j1,kat,l,lint,nn
   double precision    q1,q2,qm,qq,dq,uu, a1,a2,a3, rr,zz, xsj
   double precision    coord ,shp(3,9)
   double precision    sg(9),tg(9),wg(9)
   
   integer iocheck
   logical test
  
   integer maxa
include 'maxa.h'     
   double precision  aa
   common /adata/    aa(maxa)

   character*4     head
   common /bdata/  head(20)

   integer         numnp,numel,nummat,nen,neq
   common /cdata/  numnp,numel,nummat,nen,neq

   double precision dm
   integer             n,ma,mct,iel,nel
   common /eldata/  dm,n,ma,mct,iel,nel

   integer         ioRead,ioWrite
   common /iofile/ ioRead,ioWrite

   character       yyy*80
   common /ydata / yyy

   data wlab/'  p l a n e ','axisymmetric'/

!  Compute the quadrature points

   if(isw .gt. 2) then
     l = max(1,min(nel/2,3))
     call pgauss(l,lint,sg,tg,wg)
   end if


   if(isw.eq.1) then
!    Input material properties
     test = .true.
     do while (test)
       if(ioRead.lt.0) then 
         write(*,3000)
       end if  
       call pintio(yyy,10)
       read(yyy,'(6f10.0,2i10)',IOSTAT=iocheck) (d(i),i=1,6),nn,kat
       if (iocheck .eq. 0) then
         if(kat .ne. 2) then 
           kat=1
         end if  
         nn = max(-1,min(1,nn))
         write(ioWrite,2000) (d(i),i=1,6),nn,wlab(kat)
         if(ioRead.lt.0) then 
           write(*,2000) (d(i),i=1,6),nn,wlab(kat)
         end if  
         d(7) = nn
         d(8) = kat
         d(9) = d(2)*d(3)
         test = .false.
         return
       else  
         call pperror('PCELM6',yyy)
       end if
     end do
   else if(mod(isw,3).eq.0) then
!    Compute conductivity matrix and residual
     nn  = d(7)
     kat = d(8)
     do l=1,lint
       call shapeFunc(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
       xsj = xsj*wg(l)
       if(kat.eq.2) xsj = xsj*coord(xl,shp,ndm,nel)
       call flo06(1.0d0,shp,ul,q1,q2,qm,uu,ndf,nel)
       if(nn.eq.0) then
         qq = dm*d(4)
         dq = 0.0
       else if(nn.eq.1) then
         qq = dm*d(4)*exp(d(5)*uu)
         dq = d(5)
       else if(uu.ne.0.0d0) then
         qq = dm*d(4)*exp(d(5)-d(5)*d(6)/uu)
         dq = d(5)*d(6)/uu**2
       else
         write(*,*) ' ** ELMT06 ERROR ** T = 0.0: stop'
         stop
       end if
       j1 = 1
       do j=1,nel
         a1    = d(1)*shp(1,j)*xsj
         a2    = d(1)*shp(2,j)*xsj
         a3    = qq*shp(3,j)*xsj
!        Compute residual
         p(j1) = p(j1) + a1*q1 + a2*q2 - a3
         a3    = a3*dq
!        Compute tangent
         if(isw .eq. 3) then
           i1  = 1
           do i=1,nel
             s(i1,j1) = s(i1,j1)+a1*shp(1,i)+a2*shp(2,i)+a3*shp(3,i)
             i1 = i1 + ndf
           end do ! i
         end if
         j1 = j1 + ndf
       end do ! j
     end do ! l
   else if(isw .eq. 4) then
!    Compute the flows in each element
     call shapeFunc(0.0d0,0.0d0,xl,shp,xsj,ndm,nel,ix,.false.)
     rr = d(1)
     call flo06(rr,shp,ul,q1,q2,qm,uu,ndf,nel)
     rr  = coord(xl(1,1),shp,ndm,nel)
     zz  = coord(xl(2,1),shp,ndm,nel)
     mct = mct - 1
     if(mct.lt.0) then
       mct = 50
       write(ioWrite,2001) head
       if(ioRead.lt.0) write(*,2001) head
     end if
     write(ioWrite,'(2i5,0p2f11.3,1p4e11.3)') n,ma,rr,zz,q1,q2,qm,uu
     if(ioRead.lt.0) then 
       write(*,'(2i5,0p2f11.3,1p4e11.3)') n,ma,rr,zz,q1,q2,qm,uu
     end if  
   else if(isw .eq. 5) then
!    Compute heat capacity (mass) matrix
     do l=1,lint
       call shapeFunc(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.true.)
       xsj = xsj*wg(l)
       if(kat.eq.2) xsj = xsj*coord(xl,shp,ndm,nel)
       j1 = 1
       do j=1,nel
         p(j1) = p(j1) + d(9)*shp(3,j)*xsj
         j1 = j1 + ndf
       end do ! j
     end do ! l
   else if(isw .eq. 8) then
!    Compute the nodal flow values
     call stcn06(ix,d,xl,ul,shp,aa,aa(numnp+1),ndf,ndm,nel,numnp,sg,&
                 tg,wg,lint)
   end if

!  Formats


2000  format(3x,'Laplace Element with Reaction Loading'//       &
        4x,'Conductivity  ',e12.4/4x,'Specific Heat ',e12.4/    &
        4x,'Density       ',e12.4/4x,'Load Amplitude',e12.4/    &
        4x,'Reaction Exp. ',e12.4/4x,'Ambient Temp. ',e12.4/    &
        4x,'n - (Temp**n) ',i5/   4x,a12,' analysis')

2001  format(1x,20a4//'  L a p l a c e   E q u a t i o n'//     &
       ' elem  mat    1-coord    2-coord    1-flow     2-flow', &
       '   max flow     U-value'/)


3000  format(' Input:K, c, rho, Q, r, Ta, nn, geom(1=plane,2=axisym)' &
             /3x,'>',$)

end

