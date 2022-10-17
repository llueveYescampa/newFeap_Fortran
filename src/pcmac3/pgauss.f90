subroutine pgauss(l,lint,r,z,w)
implicit  none
integer   l,lint
double precision    r(*),z(*),w(*)
      

!  Purpose: Gauss points and weights for two dimensions

!  Inputs:
!     l         - Order of quadrature

!  Outputs:
!     lint      - Number of points
!     r(*),z(*) - Points
!     w(*)      - Weights


   integer   lr(9),lz(9),lw(9),i
   double precision    g,h

   include 'eldata.h'

   data lr/-1,1,1,-1,0,1,0,-1,0/
   data lz/-1,-1,1,1,-1,0,1,0,0/
   data lw/4*25,4*40,64/

   lint = l*l

!  1x1 integration

   if(l.eq.1) then
     r(1) = 0.
     z(1) = 0.
     w(1) = 4.

!  2x2 integration

   else if(l.eq.2) then
     g = 1.0/sqrt(3.d0)
     do i = 1,4
       r(i) = g*lr(i)
       z(i) = g*lz(i)
       w(i) = 1.
     end do  

!  3x3 integration

   else if(l.eq.3) then
     g = sqrt(0.60d0)
     h = 1.0/81.0d0
     do i = 1,9
       r(i) = g*lr(i)
       z(i) = g*lz(i)
       w(i) = h*lw(i)
     end do  
   end if

end
