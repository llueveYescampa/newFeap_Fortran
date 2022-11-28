subroutine bmat05(bm,shp,sn,cs,recr,wks)
implicit none
  double precision :: bm(5,3),shp,sn,cs,recr,wks

!  Nonlinear B-matrix for beams and axisymmetric shells

   bm(1,1) = shp
   bm(1,2) = shp*wks
   bm(2,1) = 0.5d0*cs*recr
   bm(2,2) =-0.5d0*sn*recr
   bm(3,3) = shp
   bm(4,1) = bm(2,1)*sn*recr
   bm(4,2) = bm(2,2)*sn*recr
   bm(4,3) = bm(2,1)
   bm(5,2) = shp
   bm(5,3) =-0.5d0

end
