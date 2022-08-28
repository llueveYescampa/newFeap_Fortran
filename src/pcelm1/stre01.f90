subroutine stre01(d,xl,ul,tl,shp,eps,sig,xx,yy,ndm,ndf,nel,ityp)
implicit none
integer                                        ndm,ndf,nel,ityp
double precision  d(*),xl(ndm,*),ul(ndf,*),tl(*),shp(3,*), &
                  eps(4),sig(4),xx,yy

!  Purpose: Compute stress and strain for point

!  Inputs:
!     d(*)     - Material set parameters
!     xl(ndm,*)- Nodal coordinates for element
!     ul(ndf,*) - Nodal solution parameters for element
!     tl(*)     - Nodal temperatures for element
!     shp(3,*)  - SHape functions and derivatives for point
!     ndm       - Spatial dimension of FEM mesh
!     ndf       - Number dof/node
!     nel       - Number of nodes on element
!     ityp      - Problem type

!  Outputs:
!     eps(*)    - Strains at point
!     sig(*)    - Stresses at point
!     xx        - x-coordinate at point
!     yy        - y-coordinate at point


   integer  j
   double precision   ta
   
!  Compute strains and coordinates

   call pconsd(eps,4,0.0d0)
   xx = 0.0
   yy = 0.0
   ta = -d(9)
   do j = 1,nel
     xx = xx + shp(3,j)*xl(1,j)
     yy = yy + shp(3,j)*xl(2,j)
     ta = ta + shp(3,j)*tl(j)
     eps(1) = eps(1) + shp(1,j)*ul(1,j)
     eps(2) = eps(2) + shp(1,j)*ul(2,j) + shp(2,j)*ul(1,j)
     eps(3) = eps(3) + shp(2,j)*ul(2,j)
     if(ityp.eq.3) then 
       eps(4) = eps(4) + shp(3,j)*ul(1,j)
     end if  
   end do
   ta = ta*d(10)

!  Compute stresses

   if(ityp.gt.2) then
     if(xx.ne.0.0d0) then
       eps(4) = eps(4)/xx
     else
       eps(4) = eps(1)
     end if
     sig(4) = d(1)*eps(4) + d(2)*(eps(1) + eps(3)) - ta
   else
     sig(4) = d(13)*(eps(1) + eps(3)) - ta
   end if
   sig(1) = d(1)*eps(1) + d(2)*(eps(3) + eps(4)) - ta
   sig(2) = d(3)*eps(2)
   sig(3) = d(1)*eps(3) + d(2)*(eps(1) + eps(4)) - ta
   if(ityp.eq.1) then 
     eps(4) = d(18)*(sig(1) + sig(3))
   end if  

end
