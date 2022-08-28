subroutine strn02(shp,xl,ul,type,xr0,xz0,ndm,ndf,eps)
implicit none
integer   ndm,ndf
double precision shp(3,1),xl(ndm,1),ul(ndf,1),type,xr0,xz0,eps(4)

!  Purpose: Compute strain for point

!  Inputs:
!     shp(3,*)  - Shape functions and derivatives for point
!     xl(ndm,*) - Nodal coordinates for element
!     ul(ndf,*) - Nodal solution parameters for element
!     type      - Type of problem
!     ndm       - Spatial dimension of FEM mesh
!     ndf       - Number dof/node

!  Outputs:
!     xr0       - x-(r)-coordinate at point
!     xz0       - y-(z)-coordinate at point
!     eps(*)    - Strains at point


   integer k

!  Compute strain and incremental tensors for constitutive equations

   xr0   = 0.0d0
   xz0   = 0.0d0
   call pconsd(eps,4,0.0d0)
   do k = 1,4
     xr0   = xr0   + shp(3,k)*xl(1,k)
     xz0   = xz0   + shp(3,k)*xl(2,k)
     eps(1) = eps(1) + ul(1,k)*shp(1,k)
     eps(2) = eps(2) + ul(2,k)*shp(2,k)
     eps(3) = eps(3) + shp(3,k)*ul(1,k)
     eps(4) = eps(4) + (ul(2,k)*shp(1,k) + ul(1,k)*shp(2,k))/2.
   end do ! k
   eps(3) = type*eps(3)/xr0

end
