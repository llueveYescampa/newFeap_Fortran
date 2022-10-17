subroutine stcn03(ix,d,ssa,tta,dt,st,nnp)
implicit none
integer    ix(*),              nnp
double precision d(*),ssa,tta,dt(nnp),st(nnp,*)

   integer   jj,ll
   double precision     sig1,sig2,sig3, xsj,ssg,ttg
   
   double precision    ss(4),tt(4)
   
   include 'elcom3.h'
   
   data ss/-1.0,1.0,1.0,-1.0/
   data tt/-1.0,-1.0,1.0,1.0/
   
!  Compute stress projections
   
   do jj = 1,4
     ll = abs(ix(jj))
     if(ll.gt.0) then
   
!      Compute weighted stresses at nodes
       xsj = xj0 + ss(jj)*xj1 + tt(jj)*xj2
       ssg = (ss(jj) - ssa)*beta(5)
       ttg = (tt(jj) - tta)*beta(4)
       dt(ll)   = dt(ll)   + xsj
       sig1     = (beta(1) + a1(1)*ttg + a2(1)*ssg)*xsj
       sig2     = (beta(2) + a1(2)*ttg + a2(2)*ssg)*xsj
       sig3     = (beta(3) + a1(3)*ttg + a2(3)*ssg)*xsj
       st(ll,1) = st(ll,1) + sig1
       st(ll,2) = st(ll,2) + sig3
       st(ll,3) = st(ll,3) + sig2
       st(ll,4) = st(ll,4) + d(13)*(sig1+sig2)
     end if
   end do ! jj
end
