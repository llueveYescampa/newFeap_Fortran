subroutine stcn02(ix,d,xl,ul,shp,dt,st,ndf,ndm,nnp,sg,tg,sig,eps,lint,type,ib)
implicit  none
integer   ix(*),ndf,ndm,nnp,lint,ib
double precision d(*), xl(ndm,*), ul(ndf,*),shp(3,4),dt(nnp), st(nnp,*),sg(*),&
                 tg(*),sig(*),eps(*),type
      
!  Purpose: Project stresses onto nodes

!  Inputs
!     ix(*)     - Nodes connected to element
!     d(*)      - Material set parameters
!     xl(ndf,*) - Nodal coordinates for element
!     ul(ndf,*) - Solution variables for element
!     shp(3,*)  - Shape function and derivatives
!     ndf       - Number dof/node
!     ndm       - Spatial dimension of FEM mesh
!     nnp       - Number of nodes in mesh
!     sg(*)     - Quadrature points
!     tg(*)     - Quadrature points
!     sig(4)    - Stresses at quadrature points
!     eps(4)    - Stresses at quadrature points
!     lint      - number of quadrature points
!     type      - Problem type
!     ib        - Formulation: =0 B-bar; else displacement model

!  Outputs:
!     dt(*)     - Diagonal projection matrix
!     st(*,*)   - Integral of FEM stresses over element


   integer   i,ii,jj,ll
   double precision    shpj,rr,zz,xsj(1)


!  Compute stresses at nodes from history terms

   do jj = 1,lint
     call shapeFunc(sg(jj),tg(jj),xl,shp,xsj,ndm,4,ix,.false.)
     call strn02(shp,xl,ul,type,rr,zz,ndm,ndf,eps)
     call modl02(d,ul,eps,sig,1.0d0,ndf,ib)
     do ii = 1,4
       ll = abs(ix(ii))
       if(ll.gt.0) then
         shpj = shp(3,i)*xsj(1)
         dt(ll) = dt(ll) + shpj
         do i = 1,4
           st(ll,i) = st(ll,i) + sig(i)*shpj
         end do ! i
       end if
     end do ! ii
   end do ! jj

end
