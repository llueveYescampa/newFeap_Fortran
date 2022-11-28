subroutine modl05(sig,vl,dd,d,sn,cs,sl,recr,wks)
implicit none
  double precision :: sig(5),vl(3,2),dd(5,5),d(*),sn,cs,sl,recr,wks

   double precision    v1,eps(5)
   
!  Elasticity matrix for beam/axisymmetric shells

   dd(1,1) = d(7)
   dd(2,2) = d(7)
   dd(1,2) = d(8)
   dd(2,1) = d(8)
   dd(3,3) = d(9)
   dd(4,4) = d(9)
   dd(3,4) = d(10)
   dd(4,3) = d(10)
   dd(5,5) = 5.*d(7)*(1.0 - d(2))/12.0

!  Strains

   v1     = 0.5*((vl(1,1) + vl(1,2))*cs - (vl(2,1) + vl(2,2))*sn)
   eps(1) = (vl(1,2) - vl(1,1))/sl + 0.5*wks*wks
   eps(2) =  v1*recr
   eps(3) = (vl(3,2) - vl(3,1))/sl
   eps(4) = (v1*sn*recr +0.5*(vl(3,1) + vl(3,2))*cs)*recr
   eps(5) = (vl(2,2) - vl(2,1))/sl - 0.5*(vl(3,1) + vl(3,2))

!  Stresses

   sig(1) = dd(1,1)*eps(1) + dd(1,2)*eps(2)
   sig(2) = dd(2,1)*eps(1) + dd(2,2)*eps(2)
   sig(3) = dd(3,3)*eps(3) + dd(3,4)*eps(4)
   sig(4) = dd(4,3)*eps(3) + dd(4,4)*eps(4)
   sig(5) = dd(5,5)*eps(5)

end
