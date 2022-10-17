subroutine modl04(d,eps, sig,ad)
implicit none
double precision d(10), eps, sig, ad
      
   double precision yld, mySum, gam, phi
      
   include 'maxa.h'      

   include 'ddata.h'
   include 'hdata.h'
   
      
!  trial stress
   ad  = d(7)
   sig = d(7)*(eps - dm(nh1))
   if(d(8) .gt. 0.0d0) then
!     compute plastic corrections
      yld = d(8) + d(9)*dm(nh1+2)
      phi = abs(sig - dm(nh1+1))
      if(phi.gt.yld) then
         mySum      = d(7) + d(9) + d(10)
         gam      = (phi-yld-1.d-08*d(8))/mySum
         sig      = sig - d(7)*gam
         ad       = ad -  d(7)**2/mySum
!        update the history terms
         mySum      = gam*(sig - dm(nh1+1))/phi
         dm(nh1  ) = dm(nh1  ) + gam*mySum
         dm(nh1+1) = dm(nh1+1) + gam*mySum*d(10)
         dm(nh1+2) = dm(nh1+2) + gam
      end if
   end if
end
