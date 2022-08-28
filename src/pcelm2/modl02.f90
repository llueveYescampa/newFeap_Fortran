subroutine modl02(d,ul,eps,sig,xsj,ndf,ib)
implicit none
integer  ndf,ib
double precision     d(*),ul(ndf,*),eps(4),sig(5),xsj

!  Purpose:  Form material tangent and stress state at point

!  Inputs:
!       d(*)    -  material parameters
!       eps(4)  -  current strains at point
!       h(*)    -  history terms at point
!       nh      -  number of history terms

!  Ouputs:
!       ad(4,4) -  current material tangent moduli
!       sig(4)  -  stresses at point.
!       sig(5)  -  yield state at point.


   integer i,j
   double precision temp

   double precision g     ,ad
   common /elcom2/  g(2,4),ad(4,4)

   integer         nh1,nh2
   common /hdata/  nh1,nh2

   integer         maxa
include 'maxa.h'      
   double precision h
   common          h(maxa)

!  Compute material moduli and stresses

   call pconsd(ad,16,0.0d0)
   call elpl02(d,eps,h(nh2),h(nh2+4),h(nh2+8),sig,ul,ndf,ib)
   nh2 = nh2 + 9

!  Multiply by jacobian

   do i = 1,4
     sig(i) = sig(i)*xsj
     do j = 1,4
        ad(i,j) = ad(i,j)*xsj
     end do ! j
   end do ! i

!  Reorder stresses for stress divergence calculations and prints

   temp   = sig(4)
   sig(4) = sig(3)
   sig(3) = sig(2)
   sig(2) = temp

end
