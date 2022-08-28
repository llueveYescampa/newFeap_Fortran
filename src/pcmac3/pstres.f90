subroutine pstres(sig,p1,p2,p3)
 implicit none
 double precision   sig(3),p1,p2,p3

!  Purpose: Principal stresses (2 dimensions): sig = sig-xx,sig-xy,sig-yy

!  Inputs:
!     sig(*)   - Stress state (as above)

!  Outputs:
!     p1,p2    - Principal stresses
!     p3       - Angle p1 makes with x-1 axis


   double precision   xi1,xi2,rho

   xi1 = (sig(1) + sig(3)) * 0.5
   xi2 = (sig(1) - sig(3)) * 0.5 
   rho = sqrt(xi2*xi2 + sig(2)*sig(2))
   p1 = xi1 + rho
   p2 = xi1 - rho
   p3 = 45.0
   if(xi2.ne.0.0d0) then
     p3 = 22.5*atan2(sig(2),xi2)/atan(1.0d0)
   end if  

end
