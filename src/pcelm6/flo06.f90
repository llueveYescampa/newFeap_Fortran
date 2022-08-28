subroutine flo06(d,shp,ul,q1,q2,qm,uu,ndf,nel)
implicit none
integer ndf,nel
double precision d,shp(3,*),ul(ndf,*),q1,q2,qm,uu

   integer i
!  Compute flows at current point

   q1 = 0.0d0
   q2 = 0.0d0
   uu = 0.0d0
   do i = 1,nel
     q1 = q1 - d*shp(1,i)*ul(1,i)
     q2 = q2 - d*shp(2,i)*ul(1,i)
     uu = uu +   shp(3,i)*ul(1,i)
   end do ! i
   qm =    sqrt(q1*q1 + q2*q2)
end
