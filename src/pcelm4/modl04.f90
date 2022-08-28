subroutine modl04(d,eps, sig,ad)
implicit none
double precision d(10), eps, sig, ad
      
   double precision yld, mySum, gam, h
   real phi
      
   integer         nh1,nh2
   common /hdata/  nh1,nh2
      
   integer maxa
include 'maxa.h'      
   common h(maxa)
      
!  trial stress
   ad  = d(7)
   sig = d(7)*(eps - h(nh1))
   if(d(8) .gt. 0.0d0) then
!     compute plastic corrections
      yld = d(8) + d(9)*h(nh1+2)
      phi = abs(sig - h(nh1+1))
      if(phi.gt.yld) then
         mySum      = d(7) + d(9) + d(10)
         gam      = (phi-yld-1.d-08*d(8))/mySum
         sig      = sig - d(7)*gam
         ad       = ad -  d(7)**2/mySum
!        update the history terms
         mySum      = gam*(sig - h(nh1+1))/phi
         h(nh1  ) = h(nh1  ) + gam*mySum
         h(nh1+1) = h(nh1+1) + gam*mySum*d(10)
         h(nh1+2) = h(nh1+2) + gam
      end if
   end if
end
