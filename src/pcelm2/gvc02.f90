subroutine gvc02(shp,shp3,xsj,wg,xl,type,ndm)
implicit none
integer  ndm
double precision  shp(3,4),shp3(4),xsj,wg,xl(ndm,4),type

!  Purpose: Compute volumetric integrals for b-bar

!  Inputs:
!     shp(3,*)    - Shape fucntions and derivatives
!     shp3(*)     - Shape function for axisymmetric, otherwise 1.0
!     xsj         - Jacobian
!     xl(ndm,*)   - Nodal coordinates for element
!     type        - Type of problem
!     ndm         - Spatial dimension of mesh

!  Outputs:
!     none        - stored in common block


   integer i
   double precision  rr   

   include 'elcom2.h'


   if(type.ne.0.0d0) then
     rr = 0.0d0
     do i = 1,4
       rr = rr + shp(3,i)*xl(1,i)
     end do ! i
     xsj = xsj*rr
   end if
   xsj = xsj*wg
   do i = 1,4
     shp3(i) = 0.0d0
     if(type.ne.0.0d0) then
       shp3(i) = shp(3,i)/rr
     end if  
     g(1,i) = g(1,i) + (shp(1,i) + shp3(i))*xsj
     g(2,i) = g(2,i) +  shp(2,i)*xsj
   end do ! i
end
