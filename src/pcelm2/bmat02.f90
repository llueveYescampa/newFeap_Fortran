subroutine bmat02(sh3,shp,g,bbar,ib)
implicit none
  integer          :: ib
  double precision :: sh3,shp(3),g(*),bbar(4,2)

!  Purpose: Form B-bar matrix for plane and axisymmetric problems

!  Inputs:
!     sh3       - Shape functioin for axisymmetry, otherwise 1.0
!     shp(30    - Shape functions and derivatives for node
!     g(*)      - b-bar for volumetric terms
!     ib        - Switch:  =0 for B-bar; not 0 for displacement

!  Outputs:
!     bbar(4,2) - B-bar matrix for node


   integer    i
   double precision     bb1,bb2

   bbar(1,1) = shp(1)
   bbar(2,1) = 0.0d0
   bbar(3,1) = sh3
   bbar(4,1) = shp(2)
   bbar(1,2) = 0.0d0
   bbar(2,2) = shp(2)
   bbar(3,2) = 0.0d0
   bbar(4,2) = shp(1)

!  Correct for B-bar effects

   if(ib.eq.0) then
     bb1 = (g(1) - shp(1) - sh3)/3.0d0
     bb2 = (g(2) - shp(2))/3.0d0
     do i = 1,3
       bbar(i,1) = bbar(i,1) + bb1
       bbar(i,2) = bbar(i,2) + bb2
     end do ! i
   end if
end
