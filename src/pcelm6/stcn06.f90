subroutine stcn06(ix,d,xl,ul,shp,dt,st,ndf,ndm,nel,numnp,sg,tg,wg,lint)
implicit none
integer          :: ix(*),ndf,ndm,nel,numnp,lint
double precision :: d(*),xl(ndm,*),ul(ndf,*),shp(3,*),dt(numnp), &
                    st(numnp,*),sg(*),tg(*),wg(*)
      
!  Project values to nodes
   
   integer ii,l,ll
   double precision xsj,xsji, q1,q2,qm, uu

   do l = 1,lint
      call shapeFunc(sg(l),tg(l),xl,shp,xsj,ndm,nel,ix,.false.)
      xsji = d(1)
      call flo06(xsji,shp,ul,q1,q2,qm,uu,ndf,nel)
      do ii = 1,nel
        ll = abs(ix(ii))
        if(ll.gt.0) then
          xsji     = xsj*shp(3,ii)*wg(l)
          dt(ll)   = dt(ll)   + xsji
          st(ll,1) = st(ll,1) + q1*xsji
          st(ll,3) = st(ll,3) + q2*xsji
          st(ll,4) = st(ll,4) + qm*xsji
          st(ll,5) = st(ll,5) + ul(1,ii)*xsji
        end if
      end do ! ii
   end do ! l
end
