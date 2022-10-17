subroutine stcn01(ix,d,xl,ul,tl,shp,dt,st,ndf,ndm,nel,numnp,sg,tg, &
                  sig,eps,lint,ityp)
implicit none
integer  ix(*),                 lint,ityp,ndf,ndm,nel,numnp
double precision d(*),xl(ndm,*),ul(ndf,*),tl(*),shp(3,4),dt(numnp), &
                 st(numnp,*),sg(9),tg(9),sig(4),eps(4)
      
!  Purpose: Project stresses onto nodes

!  Inputs
!     ix(*)     - Nodes connected to element
!     d(*)      - Material set parameters
!     xl(ndf,*) - Nodal coordinates for element
!     ul(ndf,*) - Solution variables for element
!     ix(*)     - Nodal connections for element
!     tl(*)     - Nodal temperatures for element
!     shp(3,*)  - Shape function and derivatives
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nel       - Number of nodes on element
!     numnp     - Number of nodes in mesh
!     sg(*)     - Quadrature points
!     tg(*)     - Quadrature points
!     sig(4)    - Stresses at quadrature points
!     eps(4)    - Stresses at quadrature points
!     lint      - number of quadrature points
!     ityp      - Problem type

!  Outputs:
!     dt(*)     - Diagonal projection matrix
!     st(*,*)   - Integral of FEM stresses over element
   
   integer  i, l,ll
   double precision   gr,xsj,xsji, xx,yy

   include 'errind.h'
   
   gr = (1.+d(17))/d(16)
   do l = 1,lint
     call shapeFunc(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
     call stre01(d,xl,ul,tl,shp,eps,sig,xx,yy,ndm,ndf,nel,ityp)
     enerr = enerr+((sig(1)+sig(3)+sig(4))**2)*d(18)*xsj +       &
                  gr*xsj*sig(2)**2
     do i = 1,4
       enerr  = enerr  + gr*sig(i)**2*xsj
     end do  
     do i = 1,nel
       ll = abs(ix(i))
       if(ll.gt.0) then
         xsji = xsj*shp(3,i)
         dt(ll) = dt(ll) + xsji
         st(ll,1) = st(ll,1) + sig(1)*xsji
         st(ll,2) = st(ll,2) + sig(2)*xsji
         st(ll,3) = st(ll,3) + sig(3)*xsji
         st(ll,4) = st(ll,4) + sig(4)*xsji
       end if
     end do  
   end do  
end
